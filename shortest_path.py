import geopandas as gpd
import pandas as pd
import numpy as np
from networkz.stottefunksjoner import gdf_concat
from networkz.tilpass_graf import bestem_ids, naermeste_noder, oppdater_graf


def shortest_path(graf: tuple,
                    startpunkter: gpd.GeoDataFrame, 
                    sluttpunkter: gpd.GeoDataFrame,
                    kostnad: str = "kostnad", #navn på kolonnen
                    id_kolonne = None, 
                    search_tolerance: int = 5000,
                    forsok = 1,
                    bufferdist_prosent = 10):

    if search_tolerance is None:
        search_tolerance = 100000000
        
    g, v, n = graf
    G, veger, noder = g.copy(), v.copy(), n.copy()

    from shapely.geometry import LineString
    
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) # ignorer en tåpelig advarsel
        
    startpunkter_kopi = startpunkter.copy().to_crs(25833)
    sluttpunkter_kopi = sluttpunkter.copy().to_crs(25833)
    
    id_kolonner = bestem_ids(id_kolonne, startpunkter_kopi, sluttpunkter_kopi)

    startpunkter_kopi = startpunkter_kopi.loc[:, ["geometry", id_kolonner[0]]]
    sluttpunkter_kopi = sluttpunkter_kopi.loc[:, ["geometry", id_kolonner[1]]]

    startpunkter_kopi["wkt"] = [punkt.wkt for punkt in startpunkter_kopi.geometry]
    sluttpunkter_kopi["wkt"] = [punkt.wkt for punkt in sluttpunkter_kopi.geometry]

    #lag dict med id-er hvis ikke geometrien brukes
    if not "geom_wkt" in id_kolonner:
        id_dict_start = {wkt: idd  for idd, wkt in zip(startpunkter_kopi[id_kolonner[0]], startpunkter_kopi["wkt"])}
        id_dict_slutt = {wkt: idd  for idd, wkt in zip(sluttpunkter_kopi[id_kolonner[1]], sluttpunkter_kopi["wkt"])}
        
    #finn avstand til nærmeste noder (fjerner også de over search_tolerance)
    startpunkter_kopi, sluttpunkter_kopi, for_langt_unna = naermeste_noder(noder, startpunkter_kopi, sluttpunkter_kopi, search_tolerance)
    
    # funksjonen get_shortest_paths() må loopes for hver fra-til-kombinasjon
    linjer = []
    for fra_id in startpunkter_kopi["wkt"]:
        for til_id in sluttpunkter_kopi["wkt"]:
            antall_forsok = 1
            start = startpunkter_kopi[startpunkter_kopi["wkt"]==fra_id]
            slutt = sluttpunkter_kopi[sluttpunkter_kopi["wkt"]==til_id]
            G2, edgelist, len_naa = oppdater_graf(G, noder, 
                startpunkter = start, 
                sluttpunkter = slutt, 
                edgelist = None,
                search_tolerance=search_tolerance,
                forsok = antall_forsok,
                bufferdist_prosent = bufferdist_prosent)
            if G2 is None:
                continue
            res = G2.get_shortest_paths(weights='weight', v = fra_id, to = til_id)
            if len(res[0])==0:
                while len(res[0]) == 0 and antall_forsok < forsok:
                    antall_forsok += 1
                    G2, edgelist, len_naa = oppdater_graf(G2, noder, 
                            startpunkter = start,
                            sluttpunkter = slutt, 
                            edgelist = edgelist,
                            search_tolerance=search_tolerance,
                            forsok = antall_forsok,
                            bufferdist_prosent = bufferdist_prosent, 
                            len_forrige_edgelist = len_naa)       
                    if G2 is None:
                        break
                    res = G2.get_shortest_paths(weights='weight', v = fra_id, to = til_id)
            if len(res[0])==0:
                linjer.append(gpd.GeoDataFrame(pd.DataFrame({"fra": [fra_id], "til": [til_id], kostnad: [np.nan], "geometry": [LineString()]}), geometry="geometry", crs=25833))
                continue
            path = G2.vs[res[0]]["name"]
            linje = veger.loc[(veger.source.isin(path)) & (veger.target.isin(path)), ["geometry"]]
            linje = linje.dissolve()
            linje["fra"] = fra_id
            linje["til"] = til_id
            kost = G2.distances(weights='weight', source = start["wkt"], target = slutt["wkt"])
            linje[kostnad] = kost[0][0]
            start = start[["wkt", "node_id", "dist"]].rename(columns={"node_id": "naermeste_node_fra", "dist": "fra_node_avstand"})
            linje = linje.merge(start, left_on = "fra", right_on = "wkt", how="left")
            slutt = slutt[["wkt", "node_id", "dist"]].rename(columns={"node_id": "naermeste_node_til", "dist": "til_node_avstand"})
            linje = linje.merge(slutt, left_on = "til", right_on = "wkt", how="left")
            linje["forsok"] = antall_forsok
            linje["bufferdist_fra"] = linje["fra_node_avstand"] * (1 + bufferdist_prosent/100) **  linje["forsok"] + 5*linje["forsok"]
            linje["bufferdist_til"] = linje["til_node_avstand"] * (1 + bufferdist_prosent/100) **  linje["forsok"] + 5*linje["forsok"]
            linje = linje.loc[:, ~linje.columns.str.contains("wkt|_node_avstand")]
            linjer.append(linje)
    
    linjer = gdf_concat(linjer)

    if for_langt_unna is not None:
        linjer = gdf_concat([linjer, for_langt_unna])
        
    if not "geom_wkt" in id_kolonner:
        linjer["fra"] = linjer["fra"].map(id_dict_start)
        linjer["til"] = linjer["til"].map(id_dict_slutt)
    
    return linjer

