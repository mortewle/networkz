import geopandas as gpd
import pandas as pd
import numpy as np
from networkz.stottefunksjoner import gdf_concat, bestem_ids, lag_midlr_id, map_ids
from networkz.lag_igraph import lag_graf, m_til_min, m_til_treige_min
from shapely.geometry import LineString



def shortest_path(G,
                  nettverk,
                    startpunkter: gpd.GeoDataFrame, 
                    sluttpunkter: gpd.GeoDataFrame,
                    id_kolonne = None,
                    cutoff: int = None,
                    destination_count: int = None,
                    ):

    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) # ignorer tåpelig advarsel
            
    startpunkter = startpunkter.copy().to_crs(25833)
    sluttpunkter = sluttpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne, startpunkter, sluttpunkter)

    startpunkter, sluttpunkter = lag_midlr_id(G.noder, startpunkter, sluttpunkter)
    
    G2, startpunkter, sluttpunkter = lag_graf(G,
                                              nettverk,
                                              G.kostnad,
                                              startpunkter, 
                                              sluttpunkter)
    
    # funksjonen get_shortest_paths() må loopes for hver fra-til-kombinasjon
    linjer = []
    for fra_id in startpunkter["nz_idx"]:
        for til_id in sluttpunkter["nz_idx"]:
            
            start = startpunkter[startpunkter["nz_idx"]==fra_id]
            slutt = sluttpunkter[sluttpunkter["nz_idx"]==til_id]
            
            if len(start)==0 or len(slutt)==0:
                continue
            
            res = G2.get_shortest_paths(weights='weight', 
                                        v = fra_id, 
                                        to = til_id)
            
            if len(res[0])==0:
                linjer.append(gpd.GeoDataFrame(pd.DataFrame({"fra": [fra_id], "til": [til_id], G.kostnad: [np.nan], "geometry": LineString()}), geometry="geometry", crs=25833))
                continue
            
            path = G2.vs[res[0]]["name"]
            
            linje = nettverk.loc[(nettverk.source.isin(path)) & (nettverk.target.isin(path)), ["geometry"]]
            
            linje = linje.dissolve()
            
            linje["fra"] = fra_id
            linje["til"] = til_id
            
            # for å få kostnaden også:
            kost = G2.distances(weights='weight', source = start["nz_idx"], target = slutt["nz_idx"])
            linje[G.kostnad] = kost[0][0]
            
            linjer.append(linje)
    
    linjer = gdf_concat(linjer)
    
    linjer = (linjer
          .replace([np.inf, -np.inf], np.nan)
          .loc[(linjer[G.kostnad] > 0) | (linjer[G.kostnad].isna())]
          .reset_index(drop=True)
          .merge(startpunkter[["dist_node_start", "nz_idx"]], left_on="fra", right_on="nz_idx", how="left")
          .drop("nz_idx", axis=1)
          .merge(sluttpunkter[["dist_node_slutt", "nz_idx"]], left_on="til", right_on="nz_idx", how="left")
          .drop("nz_idx", axis=1)
    )

    linjer = map_ids(linjer, id_kolonner, startpunkter, sluttpunkter)

    if cutoff is not None:
        linjer = linjer[linjer[G.kostnad] < cutoff]
        
    if destination_count is not None:
        linjer = linjer.loc[linjer.groupby('fra')[G.kostnad].idxmin()].reset_index(drop=True)
        
    if G.kostnad=="minutter" and G.kost_til_nodene is True:
        linjer[G.kostnad] = linjer[G.kostnad] - m_til_treige_min(linjer["dist_node_start"], G.kjoretoy) - m_til_treige_min(linjer["dist_node_slutt"], G.kjoretoy)
        linjer[G.kostnad] = linjer[G.kostnad] + m_til_min(linjer["dist_node_start"], G.kjoretoy) + m_til_min(linjer["dist_node_slutt"], G.kjoretoy)
    
    linjer = linjer[["fra", "til", G.kostnad, "geometry"]]
    
    return linjer

