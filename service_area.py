import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import LineString
from networkz.stottefunksjoner import gdf_concat, bestem_ids, map_ids
from networkz.lag_graf import lag_graf


def service_area(G,
                 startpunkter: gpd.GeoDataFrame,
                 kostnad,
                 id_kolonne = None, # hvis ikke id-kolonne oppgis, brukes startpunktenes geometri som id
                 ):
        
    veger = G.nettverk
    startpunkter = startpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne, startpunkter)
    if id_kolonne is None:
        id_kolonne = "fra"
    
    if isinstance(kostnad, int) or isinstance(kostnad, str):
        kostnad = [kostnad]
        
    G2, startpunkter = lag_graf(G, G.kostnad, startpunkter)

    # loop for hver kostnad og hvert startpunkt
    alle_service_areas = []   
    for i in startpunkter["nz_idx"]:
        for kost in kostnad:

            startpunkt = startpunkter[startpunkter["nz_idx"]==i]
            
            # beregn alle kostnader fra startpunktet
            resultat = G2.distances(weights='weight', source = startpunkt["nz_idx"].iloc[0])

            # lag tabell av resultatene og fjern alt over ønsket kostnad
            df = pd.DataFrame(data={"name": np.array(G2.vs["name"]), G.kostnad: resultat[0]})
            df = df[df[G.kostnad] < kost]

            if len(df) == 0:
                alle_service_areas.append(gpd.GeoDataFrame(pd.DataFrame({id_kolonne: [i], G.kostnad: kost, "geometry": LineString()}), geometry="geometry", crs=25833))
                continue
            
            # velg ut vegene som er i dataframen vi nettopp lagde. Og dissolve til én rad.
            sa = (veger
                .loc[veger.target.isin(df.name), ["geometry"]]
                .dissolve()
                .reset_index(drop=True)
            )
            # lag kolonner for id, kostnad og evt. node-info
            sa[id_kolonne] = i
            sa[G.kostnad] = kost
            alle_service_areas.append(sa)

    alle_service_areas = gdf_concat(alle_service_areas)

    # få tilbake opprinnelige id-er
    alle_service_areas = map_ids(alle_service_areas, id_kolonner, startpunkter)

    alle_service_areas = alle_service_areas[[id_kolonne, G.kostnad, "geometry"]]
    
    return alle_service_areas


