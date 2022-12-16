import numpy as np
import pandas as pd
import geopandas as gpd
import pygeos
from networkz.id_greier import bestem_ids, lag_midlr_id, map_ids
from networkz.lag_igraph import lag_graf


def od_cost_matrix(G,
                    startpunkter: gpd.GeoDataFrame, 
                    sluttpunkter: gpd.GeoDataFrame,
                    id_kolonne = None,
                    linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                    radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. radvis=True går mye treigere.
                    cutoff: int = None,
                    destination_count: int = None,
                    ):
        
    startpunkter = startpunkter.copy().to_crs(25833)
    sluttpunkter = sluttpunkter.copy().to_crs(25833)

    id_kolonner = bestem_ids(id_kolonne,
                             startpunkter,
                             sluttpunkter)

    # må lage midlertidig id som går fra høyeste node_id+1
    startpunkter["nz_idx"], sluttpunkter["nz_idx"] = lag_midlr_id(G.noder,
                                                                  startpunkter,
                                                                  sluttpunkter)
        
    # så loop nettverksberegningen for hver kostnad (hvis flere)
    kostnadene = G.kostnad
    if isinstance(kostnadene, str):
         kostnadene = [kostnadene]

    out = []
    for kostnad in kostnadene:
                
        G2, startpunkter, sluttpunkter = lag_graf(G,
                                                  kostnad,
                                                  startpunkter, 
                                                  sluttpunkter)
        
        if not radvis:
            
            # selve avstandsberegningen her:
            resultat = G2.distances(weights='weight',
                                    source = startpunkter["nz_idx"],
                                    target = sluttpunkter["nz_idx"])
            
            fra_idx, til_idx, kostnader = [], [], []
            for i, f_idx in enumerate(startpunkter["nz_idx"]):
                for ii, t_idx in enumerate(sluttpunkter["nz_idx"]):
                    fra_idx.append(f_idx)
                    til_idx.append(t_idx)
                    kostnader.append(resultat[i][ii])

        else:                
            fra_idx, til_idx, kostnader = [], [], []
            for f_idx, t_idx in zip(startpunkter["nz_idx"], sluttpunkter["nz_idx"]):
                resultat = G2.distances(weights='weight', 
                                        source = f_idx, 
                                        target = t_idx)
                fra_idx.append(f_idx)
                til_idx.append(t_idx)
                kostnader.append(resultat[0][0])
        
        df = pd.DataFrame(data = {"fra": fra_idx, "til": til_idx, kostnad: kostnader})

        # litt opprydning
        df = (df
            .replace([np.inf, -np.inf], np.nan)
            .loc[(df[kostnad] > 0) | (df[kostnad].isna())]
            .reset_index(drop=True)
        )

        out.append(df)
            
    # samle resultatene for ulike kostnader, med fra og til som index
    out = pd.concat([df.set_index(["fra", "til"]) for df in out], ignore_index=False, axis=1)
    out = out.reset_index()

    # gi dataene kolonner med avstand til nærmeste node
    startpunkter = startpunkter.loc[:, ~startpunkter.columns.duplicated()]
    sluttpunkter = sluttpunkter.loc[:, ~sluttpunkter.columns.duplicated()]
    out = (out
           .merge(startpunkter[["dist_node_start", "nz_idx"]], left_on="fra", right_on="nz_idx", how="left")
           .drop("nz_idx", axis=1)
           .merge(sluttpunkter[["dist_node_slutt", "nz_idx"]], left_on="til", right_on="nz_idx", how="left")
           .drop("nz_idx", axis=1)
    )

    # hvis cutoff og/eller destination_count er True, brukes første kostnad i filtreringen
    if isinstance(G.kostnad, str):
        kostnad1 = G.kostnad
    else:
        kostnad1 = G.kostnad[0]
              
    if cutoff:
        out = out[out[kostnad1] < cutoff].reset_index(drop=True)
        
    if destination_count:
        out = out.loc[~out[kostnad1].isna()]
        out = out.loc[out.groupby('fra')[kostnad1].idxmin()].reset_index(drop=True)

    # 
    wkt_dict_start = {idd: geom.wkt for idd, geom in zip(startpunkter["nz_idx"], startpunkter.geometry)}
    wkt_dict_slutt = {idd: geom.wkt for idd, geom in zip(sluttpunkter["nz_idx"], sluttpunkter.geometry)}
    out["wkt_fra"] = out["fra"].map(wkt_dict_start)
    out["wkt_til"] = out["til"].map(wkt_dict_slutt)
    for kostnad in kostnadene:
        out[kostnad] = [0 if fra==til else out[kostnad].iloc[i] for i, (fra, til) in enumerate(zip(out.wkt_fra, out.wkt_til))]

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

