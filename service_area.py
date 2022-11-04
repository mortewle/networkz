import geopandas as gpd
import pandas as pd
import numpy as np
from networkz.tilpass_graf import bestem_ids, naermeste_noder, oppdater_graf
from networkz.stottefunksjoner import gdf_concat


def service_area(graf: tuple, # grafen som lages i funksjonen 'graf' eller 'graf_fra_geometri'
    startpunkter: gpd.GeoDataFrame, 
    kostnad, # én eller flere antall meter/minutter
    kostnadskolonne: str = "kostnad", #navnet på kostnadskolonnen som returneres
    id_kolonne = None, # hvis ikke id-kolonne oppgis, brukes startpunktenes geometri som id
    search_tolerance = 5000, # maks avstand til nærmeste node
    forsok = 1,
    bufferdist_prosent = 10):

    if search_tolerance is None:
        search_tolerance = 100000000
        
    g, v, n = graf
    G, veger, noder = g.copy(), v.copy(), n.copy()

    startpunkter_kopi = startpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne, startpunkter_kopi)
    if id_kolonne is None:
        id_kolonne = "fra"

    #bruker geometrien som id
    startpunkter_kopi["wkt"] = [punkt.wkt for punkt in startpunkter_kopi.geometry]

    #lagre dict med id-ene for å kunne få tilbake opprinnelige id-er 
    if not "geom_wkt" in id_kolonner:
        id_dict_start = {wkt: idd  for idd, wkt in zip(startpunkter_kopi[id_kolonner[0]], startpunkter_kopi["wkt"])}
        
    #finn avstand til nærmeste noder
    startpunkter_kopi, for_langt_unna = naermeste_noder(noder, startpunkter_kopi, search_tolerance=search_tolerance)

    # loop for hver kostnad og hvert startpunkt
    alle_service_areas = []
    if isinstance(kostnad, int) or isinstance(kostnad, str):
        kostnad = [kostnad]
    for kost in kostnad:
        for i in startpunkter_kopi["wkt"]:

            startpunkt = startpunkter_kopi[startpunkter_kopi["wkt"]==i]

            # legg til lenker mellom startpunkt og nodene innen ønsket avstand
            G2, edgelist, len_naa = oppdater_graf(G, noder, 
                startpunkter = startpunkt, 
                sluttpunkter = None,
                edgelist = None,
                search_tolerance=search_tolerance,
                forsok = forsok,
                bufferdist_prosent = bufferdist_prosent)

            # hvis ingen punkter er innen search_tolerance, gå videre til neste startpunkt
            if G2 is None:
                continue
            
            # beregn alle kostnader fra startpunktet
            resultat = G2.distances(weights='weight', source = startpunkt["wkt"].iloc[0])

            # lag tabell av resultatene og fjern alt over ønsket kostnad
            df = pd.DataFrame(data={"name": np.array(G2.vs["name"]), kostnadskolonne: resultat[0]})
            df = df[df[kostnadskolonne] < kost]

            if len(df) == 0:
                continue
            
            # velg ut vegene som er i dataframen vi nettopp lagde. Og dissolve til én rad.
            sa = (veger
                .loc[veger.target.isin(df.name), ["geometry"]]
                .dissolve()
                .reset_index(drop=True)
            )
            # lag kolonner for id, kostnad og evt. node-info
            sa[id_kolonne] = i
            sa[kostnadskolonne] = kost
            sa["naermeste_node_fra"] = startpunkt["node_id"].iloc[0]
            sa["bufferdist"] = startpunkt["dist"].iloc[0] * (1 + bufferdist_prosent/100) ** forsok + 5*forsok
            alle_service_areas.append(sa)

    alle_service_areas = gdf_concat(alle_service_areas)

    if for_langt_unna is not None:
        alle_service_areas = gdf_concat([alle_service_areas, for_langt_unna])
        
    if not "geom_wkt" in id_kolonner:
        alle_service_areas[id_kolonne] = alle_service_areas[id_kolonne].map(id_dict_start)

    return alle_service_areas


