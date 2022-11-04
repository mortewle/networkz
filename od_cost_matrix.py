import geopandas as gpd
import pandas as pd
import numpy as np
import pygeos
import time
from networkz.tilpass_graf import bestem_ids, naermeste_noder, oppdater_graf


def od_cost_matrix(graf: tuple,
                    startpunkter: gpd.GeoDataFrame, 
                    sluttpunkter: gpd.GeoDataFrame,
                    kostnad: str = "kostnad",
                    id_kolonne = None,
                    search_tolerance: int = 5000, 
                    returner_linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                    radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                    maks_distanse: int = None, # TODO?
                    maks_antall_naboer: int = None, # TODO?
                    forsok: int = 1, # hvis mer enn ett forsøk, lages nye lenker mellom de manglende punktene og noder innen en radius som øker for hvert forsøk
                    bufferdist_prosent = 10 # hvor mange prosent lenger enn nærmeste node man vil koble noder. Øker med x prosent hver gang, pluss 5 meter hver gang.
                    ):

    if search_tolerance is None:
        search_tolerance = 100000000
        
    g, v, n = graf
    G, noder = g.copy(), n.copy()
    
    startpunkter_kopi = startpunkter.copy().to_crs(25833)
    sluttpunkter_kopi = sluttpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne, startpunkter_kopi, sluttpunkter_kopi)

    startpunkter_kopi = startpunkter_kopi.loc[:, ["geometry", id_kolonner[0]]]
    sluttpunkter_kopi = sluttpunkter_kopi.loc[:, ["geometry", id_kolonner[1]]]
        
    #TODO
    if maks_distanse is not None:
        joinet = startpunkter_kopi.sjoin_nearest(sluttpunkter_kopi, distance_col="distanse")
    
    #TODO
    if maks_antall_naboer is not None:
        "putt noe her"
        
    # bruker geometrien som id. Vurdere å endre for å få opp farten?
    startpunkter_kopi["wkt"] = [punkt.wkt for punkt in startpunkter_kopi.geometry]
    sluttpunkter_kopi["wkt"] = [punkt.wkt for punkt in sluttpunkter_kopi.geometry]

    # lag dict som lagrer opprinnelige id-er (hvis ikke geometrien brukes som id-er)
    if not "geom_wkt" in id_kolonner:
        id_dict_start = {wkt: idd  for idd, wkt in zip(startpunkter_kopi[id_kolonner[0]], startpunkter_kopi["wkt"])}
        id_dict_slutt = {wkt: idd  for idd, wkt in zip(sluttpunkter_kopi[id_kolonner[1]], sluttpunkter_kopi["wkt"])}
    
    #finn avstand til nærmeste node
    antall_forsok = 1
    startpunkter_kopi, sluttpunkter_kopi, for_langt_unna = naermeste_noder(noder, startpunkter_kopi, sluttpunkter_kopi, search_tolerance)
    
    # legg til lenker mellom punktene og noder innen gitt avstand
    G2, edgelist, len_naa = oppdater_graf(
            graf = G, 
            noder = noder, 
            startpunkter = startpunkter_kopi, 
            sluttpunkter = sluttpunkter_kopi, 
            edgelist = None,
            search_tolerance=search_tolerance,
            forsok = antall_forsok,
            bufferdist_prosent = bufferdist_prosent)
    
    # lager en funksjon av selve beregningen siden denne gjentas for OD-parene som mangler
    def kjor_od_cost(G2, startpunkter_kopi, sluttpunkter_kopi, radvis, antall_forsok):
        
        if radvis: 
            fra_wkt, til_wkt, kostnader = [], [], []
            for f_wkt, t_wkt in zip(startpunkter_kopi["wkt"], sluttpunkter_kopi["wkt"]):
                resultat = G2.distances(weights='weight', 
                                        source = f_wkt, 
                                        target = t_wkt)
                fra_wkt.append(f_wkt)
                til_wkt.append(t_wkt)
                kostnader.append(resultat[0][0])    
        else:
            
            resultat = G2.distances(weights='weight', 
                                    source = startpunkter_kopi["wkt"], 
                                    target = sluttpunkter_kopi["wkt"])

            fra_wkt, til_wkt, kostnader = [], [], []
            for i, f_wkt in enumerate(startpunkter_kopi["wkt"]):
                for ii, t_wkt in enumerate(sluttpunkter_kopi["wkt"]):
                    fra_wkt.append(f_wkt)
                    til_wkt.append(t_wkt)
                    kostnader.append(resultat[i][ii])
              
        return pd.DataFrame(data = {"fra": fra_wkt, "til": til_wkt, kostnad: kostnader, "forsok": antall_forsok}).replace([np.inf, -np.inf], np.nan)
    
    df = kjor_od_cost(G2, startpunkter_kopi, sluttpunkter_kopi, radvis, antall_forsok)

    print("prosent som mangler nå: " + str(len(df[df.kostnad.isna()]) / len(df)*100))

    df["fra_til"] = df.fra.astype(str) + df.til.astype(str)
        
    #kjør igjen for de som mangler fram til ønsket antall forsøk er nådd
    mangler = df[df[kostnad].isna()]
    mangler = mangler.loc[mangler["fra"] != mangler["til"]]

    df_mangler = "tekst_som_ikke_har_lengde_paa_0"
    while len(mangler) > 0 and antall_forsok < forsok:
        
        antall_forsok += 1          

        tid = time.perf_counter()

        # hent punktene som mangler minst en rute
        mangler_start = startpunkter_kopi[startpunkter_kopi["wkt"].isin(mangler["fra"])]
        mangler_slutt = sluttpunkter_kopi[sluttpunkter_kopi["wkt"].isin(mangler["til"])]            
            
        #finn punktene som skaper problemer ved å telle opp hvor mange ganger hvert punkt er representert i 'mangler'
        mangler["antall_mangler_fra"] = mangler["fra"].map(mangler["fra"].value_counts())
        mangler["antall_mangler_til"] = mangler["til"].map(mangler["til"].value_counts())

        # velg ut punktene som er over gjennomsnittlig problematiske eller alle med mer enn 1 manglende rute hvis siste forsøk
        if antall_forsok == forsok:
            problempunktene_start = mangler.loc[mangler["antall_mangler_fra"] > 1, "fra"]
            problempunktene_slutt = mangler.loc[mangler["antall_mangler_til"] > 1, "til"]
        else:
            problempunktene_start = mangler.loc[mangler["antall_mangler_fra"] >= np.mean(mangler["antall_mangler_fra"]), "fra"]
            problempunktene_slutt = mangler.loc[mangler["antall_mangler_til"] >= np.mean(mangler["antall_mangler_til"]), "til"]
        
        # bruk kun problempunktene i oppdater_graf() for å gjøre det kjappere
        mangler_start_problempunktene = mangler_start[mangler_start["wkt"].isin(problempunktene_start)]
        mangler_slutt_problempunktene = mangler_slutt[mangler_slutt["wkt"].isin(problempunktene_slutt)]
        
        G2, edgelist, len_naa = oppdater_graf(G2, noder, 
                startpunkter = mangler_start_problempunktene, 
                sluttpunkter = mangler_slutt_problempunktene, 
                edgelist = edgelist,
                search_tolerance = search_tolerance,
                forsok = antall_forsok, 
                bufferdist_prosent = bufferdist_prosent, 
                len_forrige_edgelist = len_naa)
        
        # avslutt letingen hvis ingen nye lenker er funnet tre ganger på rad
        if G2 is None:
            df.loc[df.kostnad.isna(), "forsok"] = antall_forsok
            break
        
        df_mangler = kjor_od_cost(G2, mangler_start, mangler_slutt, radvis=radvis, antall_forsok=antall_forsok)
        
        df_mangler = df_mangler.loc[df_mangler["fra"] != df_mangler["til"]]

        # velg først ut de som manglet forrige gang
        df_mangler["fra_til"] = df_mangler.fra.astype(str) + df_mangler.til.astype(str)
        df_mangler = df_mangler[(df_mangler["fra_til"].isin(df.loc[df[kostnad].isna(), "fra_til"]))]        

        # så velg ut de som mangler nå også
        mangler = df_mangler[df_mangler[kostnad].isna()]

        # legg til de nye resultatene
        if len(df_mangler)>0 and len(mangler)>0:
            df_mangler = df_mangler[~df_mangler[kostnad].isna()]
            df = df[(~df["fra_til"].isin(df_mangler["fra_til"]))]
            df = pd.concat([df, df_mangler], axis=0, ignore_index=True)

#        print(f"forsok {antall_forsok}. antall nye: {len(df_mangler)}. antall mangler: {len(mangler)}. tid: {time.perf_counter()-tid}")

    # koble på info om nærmeste node-id og avstand
    df = (df
          .drop("fra_til", axis=1)
          .merge(startpunkter_kopi[["wkt", "node_id", "dist"]], left_on = "fra", right_on = "wkt")
          .rename(columns={"node_id": "naermeste_node_fra"})
    )
    
    df = df.loc[df["fra"] != df["til"]]
    
    # beregn så bufferdistansen som ble brukt
    df["bufferdist_fra"] = df["dist"] * (1 + bufferdist_prosent/100) **  df["forsok"] + 5*df["forsok"]

    # samme for sluttpunktene
    df = (df
          .drop(["wkt", "dist"], axis=1)
          .merge(sluttpunkter_kopi[["wkt", "node_id", "dist"]], left_on = "til", right_on = "wkt")
          .rename(columns={"node_id": "naermeste_node_til"})
    )
    df["bufferdist_til"] = df["dist"] * (1 + bufferdist_prosent/100) **  df["forsok"] + 5*df["forsok"]
    
    # omkod bufferdistansen til laveste verdi for hvert start- og sluttpunkt
    # for å få et bedre inntrykk av hvem som egentlig er problempunktene
    df["bufferdist_fra"] = df["fra"].map(df.groupby("fra").agg(np.min)["bufferdist_fra"])
    df["bufferdist_til"] = df["til"].map(df.groupby("til").agg(np.min)["bufferdist_til"])
            
    df = df.drop(["wkt", "dist"], axis=1) 
    
    # legg til radene som var for langt unna fra starten av        
    if for_langt_unna is not None:
        df = pd.concat([df, for_langt_unna], axis=0, ignore_index=True)

    # lag linjer mellom start og slutt
    if returner_linjer:
        fra = pygeos.from_shapely(gpd.GeoSeries.from_wkt(df["fra"], crs=25833))
        til =  pygeos.from_shapely(gpd.GeoSeries.from_wkt(df["til"], crs=25833))
        df["geometry"] = pygeos.shortest_line(fra, til)
        df = gpd.GeoDataFrame(df, geometry="geometry", crs=25833)
        
    # gi tilbake de opprinnelige id-ene
    if not "geom_wkt" in id_kolonner:
        df["fra"] = df["fra"].map(id_dict_start)
        df["til"] = df["til"].map(id_dict_slutt)

    return df
