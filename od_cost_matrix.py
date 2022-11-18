import geopandas as gpd
import pandas as pd
import numpy as np
import pygeos
from networkz.stottefunksjoner import bestem_ids, lag_midlr_id, map_ids
from networkz.lag_igraph import lag_graf, m_til_min, m_til_treige_min



def od_cost_matrix(G,
                   nettverk,
                    startpunkter: gpd.GeoDataFrame, 
                    sluttpunkter: gpd.GeoDataFrame,
                    id_kolonne = None,
                    linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                    radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                    cutoff: int = None,
                    destination_count: int = None,
                    ):
    
    startpunkter = startpunkter.copy().to_crs(25833)
    sluttpunkter = sluttpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne, startpunkter, sluttpunkter)

    startpunkter, sluttpunkter = lag_midlr_id(G.noder, startpunkter, sluttpunkter)
     
    def kjor_od_cost(G, nettverk, kostnad, startpunkter, sluttpunkter, radvis):
        
        G2, startpunkter, sluttpunkter = lag_graf(G,
                                                  nettverk,
                                                  kostnad,
                                                  startpunkter, 
                                                  sluttpunkter)
            
        if not radvis: 
            
            # selve avstandsberegningen her:
            resultat = G2.distances(weights='weight', 
                                    source = startpunkter["nz_idx"], 
                                    target = sluttpunkter["nz_idx"])
            
            # lag lister med avstander og id-er
            fra_wkt, til_wkt, kostnader = [], [], []
            for i, f_wkt in enumerate(startpunkter["nz_idx"]):
                for ii, t_wkt in enumerate(sluttpunkter["nz_idx"]):
                    fra_wkt.append(f_wkt)
                    til_wkt.append(t_wkt)
                    kostnader.append(resultat[i][ii])

        else:                
            fra_wkt, til_wkt, kostnader = [], [], []
            for f_wkt, t_wkt in zip(startpunkter["nz_idx"], sluttpunkter["nz_idx"]):
                resultat = G2.distances(weights='weight', 
                                        source = f_wkt, 
                                        target = t_wkt)
                fra_wkt.append(f_wkt)
                til_wkt.append(t_wkt)
                kostnader.append(resultat[0][0])
        
        df = pd.DataFrame(data = {"fra": fra_wkt, "til": til_wkt, kostnad: kostnader})

        # litt opprydning og koble på info om avstand til noder
        df = (df
            .replace([np.inf, -np.inf], np.nan)
            .loc[(df[kostnad] > 0) | (df[kostnad].isna())]
            .reset_index(drop=True)
            .merge(startpunkter[["dist_node_start", "nz_idx"]], left_on="fra", right_on="nz_idx", how="left")
            .drop("nz_idx", axis=1)
            .merge(sluttpunkter[["dist_node_slutt", "nz_idx"]], left_on="til", right_on="nz_idx", how="left")
            .drop("nz_idx", axis=1)
        )
        
        # gi tilbake angitt fart til nodene
        if "minutter" in kostnad and G.kost_til_nodene:
            df[kostnad] = df[kostnad] - m_til_treige_min(df["dist_node_start"], G.kjoretoy) - m_til_treige_min(df["dist_node_slutt"], G.kjoretoy)
            df[kostnad] = df[kostnad] + m_til_min(df["dist_node_start"], G.kjoretoy) + m_til_min(df["dist_node_slutt"], G.kjoretoy)        
        
        return df
        
    # kjør funksjonen én eller to ganger avhengig av hvor mange kostnader man har
    if isinstance(G.kostnad, str): 
        out = kjor_od_cost(G, nettverk, G.kostnad, startpunkter, sluttpunkter, radvis)
        
    else:
        if len(G.kostnad)==2:
            df = kjor_od_cost(G, nettverk, G.kostnad[0], startpunkter, sluttpunkter, radvis)
            df2 = kjor_od_cost(G, nettverk, G.kostnad[1], startpunkter, sluttpunkter, radvis)
            out = df.merge(df2[["fra", "til", G.kostnad[1]]], on = ["fra", "til"], how="outer")
        else:
            raise ValueError("Kan bare beregne 'meter' og 'minutter'")

    # fjern rutene som går fra og til samme punkt
    wkt_dict_start = {idd: geom.wkt for idd, geom in zip(startpunkter["nz_idx"], startpunkter.geometry)}
    wkt_dict_slutt = {idd: geom.wkt for idd, geom in zip(sluttpunkter["nz_idx"], sluttpunkter.geometry)}
    out["wkt_fra"] = out["fra"].map(wkt_dict_start)
    out["wkt_til"] = out["til"].map(wkt_dict_slutt)
    
    out = out[out["wkt_fra"] != out["wkt_til"]]
    
    if isinstance(G.kostnad, str):
        kost = G.kostnad
    else:
        kost = G.kostnad[0]
              
    if cutoff is not None:
        out = out[out[kost] < cutoff].reset_index(drop=True)
        
    if destination_count is not None:
        out = out.loc[~out[kost].isna()]
        out = out.loc[out.groupby('fra')[kost].idxmin()].reset_index(drop=True)

    # lag linjer mellom start og slutt
    if linjer:        
        fra = pygeos.from_shapely(gpd.GeoSeries.from_wkt(out["wkt_fra"], crs=25833))
        til =  pygeos.from_shapely(gpd.GeoSeries.from_wkt(out["wkt_til"], crs=25833))
        out["geometry"] = pygeos.shortest_line(fra, til)
        out = gpd.GeoDataFrame(out, geometry="geometry", crs=25833)
    
    out = out.drop(["wkt_fra", "wkt_til"], axis=1, errors="ignore")
        
    # få tilbake opprinnelige id-er
    out = map_ids(out, id_kolonner, startpunkter, sluttpunkter)

    return out

