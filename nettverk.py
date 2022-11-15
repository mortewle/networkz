import numpy as np
import pandas as pd
import geopandas as gpd
import pygeos
from networkz.stottefunksjoner import fjern_tomme_geometrier, gdf_concat, les_geoparquet, kutt_linjer


def network_from_geometry(veger, 
                          minute_col = None):
    
    veger_kopi = veger.copy()
    
    veger_kopi["idx"] = veger_kopi.index
    
    veger_kopi = veger_kopi.to_crs(25833)
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = pygeos.line_merge(pygeos.force_2d(pygeos.from_shapely(veger_kopi.geometry)))   
    n = len(veger_kopi)
    veger_kopi = veger_kopi.explode(ignore_index=True)
    if len(veger_kopi)<n and minute_col is not None:
        print(f"Advarsel: {len(veger_kopi)-n} multigeometrier ble splittet. Minutt-kostnader blir feil for disse.")
        #TODO: lag ny minutt-kolonne manuelt her
    
    # nye node-id-er fra endepunkter som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder())   
    ytterpunkter = veger_kopi.copy()
    ytterpunkter["geometry"] = ytterpunkter.geometry.boundary 
    sirkler = ytterpunkter.loc[ytterpunkter.is_empty, "idx"] #sirkler har tom boundary
    veger_kopi = veger_kopi[~veger_kopi.idx.isin(sirkler)]

    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True) # boundary gir multipunkt med to ytterpunkter for hver linje, og explode() (til singlepart) gir en dobbelt så lang tabell med enkeltpunkter
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2] # gjør annenhvert ytterpunkt til source_wkt og target_wkt
    
    
    veger_kopi = make_node_ids(veger_kopi)
        
    if minute_col in veger_kopi.columns:
        veger_kopi["minutter"] = veger_kopi[minute_col]
        veger_kopi = veger_kopi[veger_kopi["minutter"] >= 0, 
                                ["source_wkt", "target_wkt", "minutter", "meter", "geometry"]]    
    else:
        veger_kopi = veger_kopi[veger_kopi["meter"] >= 0,
                                ["source_wkt", "target_wkt", "meter", "geometry"]]

    for col in ["source", "target", "source_wkt", "target_wkt"]:
        veger_kopi[col] = veger_kopi[col].astype(str).astype("category")
    
    veger_kopi["turn_restriction"] = False
    
    return veger_kopi


# nye node-id-er som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder())
def make_node_ids(veger_kopi):

#    while np.max(veger_kopi.length)>50:
 #       veger_kopi = kutt_linjer(veger_kopi, 50)
  #      print(np.max(veger_kopi.length))
    
    sources = veger_kopi[["source_wkt"]].rename(columns={"source_wkt":"wkt"})  
    targets = veger_kopi[["target_wkt"]].rename(columns={"target_wkt":"wkt"})
    noder = (pd.concat([sources, targets], axis=0, ignore_index=True)
                .drop_duplicates(subset=["wkt"])
                .reset_index(drop=True) # viktig at node_id følger index for å kunne indekse på numpy arrays i graf()
    )
    noder["node_id"] = noder.index
    noder["node_id"] = noder["node_id"].astype(str) # funker ikke med numeriske node-navn i igraph, pussig nok...
    
    #koble på de nye node-id-ene
    veger_kopi = (veger_kopi
            .drop(["source", "target", "nz_idx"], axis=1, errors="ignore")
            .merge(noder, left_on = "source_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"source"})
            .drop("wkt",axis=1)
            .merge(noder, left_on = "target_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"target"})
            .drop("wkt",axis=1)
    )
    
    veger_kopi["meter"] = veger_kopi.length
    
    # fjern duplikatlenkene med høyest kostnad (siden disse aldri vil bli brukt)
    if "minutter" in veger_kopi.columns:
        veger_kopi = veger_kopi.sort_values("minutter", ascending=True)
    else:     
        veger_kopi = veger_kopi.sort_values("meter", ascending=True)
             
    veger_kopi.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)

    return veger_kopi


# funksjon som tilpasser vegnettet til funksjonen graf(), som bygger grafen.
# outputen her kan lagres i parquet så man slipper å lage nettverket hver gang.
def lag_nettverk(veger,
                 source: str = "fromnodeid",
                 target: str = "tonodeid",
                 linkid: str = "linkid",
                 sperring = False,
                 turn_restrictions = None):
    
    veger_kopi = veger.copy()
    
    veger_kopi["idx"] = veger_kopi.index

    # endre til små bokstaver
    veger_kopi.columns = [col.lower() for col in veger_kopi.columns]

    # hvis ikke angitte kolonner finnes i vegdataene, sjekk om andre kolonner matcher. 
    # lager ny kolonne hvis ingen matcher. Gir feilmelding hvis flere enn en matcher.
    if not source in veger_kopi.columns:
        n = 0
        for col in veger_kopi.columns:
            if "from" in col and "node" in col or "source" in col:
                source = col
                n += 1
        if n == 1:
            print(f"Bruker '{source}' som source-kolonne")
        elif n == 0:
            veger_kopi[source] = np.nan
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde source-id-er")
    if not target in veger_kopi.columns:
        n = 0
        for col in veger_kopi.columns:
            if "to" in col and "node" in col or "target" in col:
                target = col
                n += 1
        if n == 1:
            print(f"Bruker '{target}' som target-kolonne")
        elif n == 0:
            veger_kopi[target] = np.nan
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde target-id-er")
    if not linkid in veger_kopi.columns:
        n = 0
        for col in veger_kopi.columns:
            if "link" in col and "id" in col:
                linkid = col
                n += 1
        if n == 1:
            print(f"Bruker '{linkid}' som linkid-kolonne")
        elif n == 0:
            veger_kopi[linkid] = np.nan #godta dette eller raise ValueError("Finner ikke linkid-kolonne") ?
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde linkid-id-er")

    # bestem om minuttkolonnen er 2022 eller tidligere
    if "drivetime_fw" in veger_kopi.columns and "drivetime_bw" in veger_kopi.columns:
        drivetime_fw, drivetime_bw = "drivetime_fw", "drivetime_bw"
    elif "ft_minutes" in veger_kopi.columns and "tf_minutes" in veger_kopi.columns:
        drivetime_fw, drivetime_bw = "ft_minutes", "tf_minutes"
    else:
        raise ValueError("Finner ikke kolonner med minutter")

    # finn eventuell kommunekolonne
    n = 0
    for col in veger_kopi.columns:
        if "komm" in col or "muni" in col:
            komm_col = col
            n += 1
    if n == 1:
        veger_kopi["KOMMUNENR"] = veger_kopi[komm_col].map(lambda x: str(int(x)).zfill(4)).astype("category")
    else:
        veger_kopi["KOMMUNENR"] = 0
    
    if sperring:
        if "sperring" in veger_kopi.columns:
            veger_kopi = veger_kopi[veger_kopi.sperring != 1] 
    
    
    # litt opprydning (til utm33-koordinater, endre navn på kolonner, beholde kun relevante kolonner, resette index)
    veger_kopi = (veger_kopi
                    .to_crs(25833)
                    .rename(columns={source: "source", target: "target", linkid: "linkid"}) 
                    [["idx", "source", "target", "linkid", drivetime_fw, drivetime_bw, "oneway", "KOMMUNENR", "geometry"]]
                    .reset_index(drop=True) 
                    )
    
    # fra multilinestring til linestring. Og fjerne z-koordinat fordi de ikke trengs
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = pygeos.line_merge(pygeos.force_2d(pygeos.from_shapely(veger_kopi.geometry)))   
      
    #hvis noen lenker fortsatt er multilinestrings, må de splittes for å ikke ha flere enn to ytterpunkter 
    n = len(veger_kopi)
    veger_kopi = veger_kopi.explode(ignore_index=True)
    if len(veger_kopi)<n:
        print(f"Advarsel: {len(veger_kopi)-n} multigeometrier ble splittet. Minutt-kostnader blir feil for disse.")
        #TODO: lag ny minutt-kolonne manuelt
    
    # fjern små veger som er kutta av fra resten av vegnettet.
    # gjerne privatveger bak bom.
    print(len(veger_kopi))
    import time
    print(time.perf_counter())
    import dask_geopandas as dg
    dask_gdf = dg.from_geopandas(veger_kopi[["geometry"]], npartitions=8)
    dask_gdf["geometry"] = dask_gdf.buffer(0.001, resolution = 1).compute() # lavest mulig resolution for å få det fort
    dissolvet = dask_gdf.dissolve().compute()
    singlepart = dissolvet.explode(ignore_index=True)
    store_nettverksdeler = singlepart[singlepart.area>5]
    veger_kopi = veger_kopi.sjoin(store_nettverksdeler, how="inner")
    print(len(veger_kopi))
    print(time.perf_counter())
    print(sum(veger_kopi.area))
    
    if sperring:
        while np.max(veger_kopi.length) > 1001:
            veger_kopi = kutt_linjer(veger_kopi, 1000)
        while np.max(veger_kopi.length) > 301:
            veger_kopi = kutt_linjer(veger_kopi, 300)
        while np.max(veger_kopi.length) > 101:
            veger_kopi = kutt_linjer(veger_kopi, 100)
        while np.max(veger_kopi.length) > 51:
            veger_kopi = kutt_linjer(veger_kopi, 50)
    
    # hent ut linjenes ytterpunkter. 
    # men først: sirkler har ingen ytterpunkter og heller ingen funksjon. Disse må fjernes
    ytterpunkter = veger_kopi.copy()
    ytterpunkter["geometry"] = ytterpunkter.geometry.boundary 
    sirkler = ytterpunkter.loc[ytterpunkter.is_empty, "idx"] #sirkler har tom boundary
    veger_kopi = veger_kopi[~veger_kopi.idx.isin(sirkler)]
    
    # lag kolonner med geometritekst (wkt) for source og target
    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True) # boundary gir multipunkt med to ytterpunkter for hver linje, og explode() (til singlepart) gir en dobbelt så lang tabell med enkeltpunkter
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2] # gjør annenhvert ytterpunkt til source_wkt og target_wkt

    #velg ut de enveiskjørte og snu source og target for lenkene som går "feil" vei
    ft = veger_kopi[(veger_kopi.oneway=="FT") | (veger_kopi.oneway=="F")] 
    tf = veger_kopi[(veger_kopi.oneway=="TF") | (veger_kopi.oneway=="T")]
    tf = tf.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})

    #dupliser lenkene som går begge veier og snu source og target i den ene
    begge_retninger1 = veger_kopi[veger_kopi.oneway=="B"]
    begge_retninger2 = begge_retninger1.copy()
    begge_retninger2 = begge_retninger2.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
    
    # lag minutt-kolonne
    begge_retninger1["minutter_bil"] = begge_retninger1[drivetime_fw]
    begge_retninger2["minutter_bil"] = begge_retninger2[drivetime_bw]
    ft["minutter_bil"] = ft[drivetime_fw]
    tf["minutter_bil"] = tf[drivetime_bw]

    veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf])

    # lag meter-kolonne
    veger_edges["meter"] = veger_edges.length
    
    # TODO: fullfør denne
    if turn_restrictions is not None:
        veger_edges = turn_restr(veger_edges, turn_restrictions)
    else:
        veger_edges["turn_restriction"] = False

    # rydd opp (fjern eventuelle 0-verdier, velg ut kolonner, fjern duplikat-lenkene med høyest kostnad, reset index)
    veger_edges = veger_edges.loc[veger_edges["minutter_bil"] >= 0, 
                                  ["source", "target", "source_wkt", "target_wkt", "minutter_bil", "meter", "turn_restriction", "KOMMUNENR", "geometry"]]
    
    #nye node-id-er som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder())
    veger_edges = make_node_ids(veger_edges)

    for col in ["source", "target", "source_wkt", "target_wkt", "turn_restriction", "KOMMUNENR"]:
        veger_edges[col] = veger_edges[col].astype(str).astype("category")
                    
    return veger_edges



def turn_restr(veger_edges, turn_restrictions):
    
    turn_restrictions.columns = [col.lower() for col in turn_restrictions.columns]

    # hvis 2021-data
    if "edge1fid" in turn_restrictions.columns:
        turn_restrictions1 = turn_restrictions.loc[turn_restrictions.edge1end=="Y", ["edge1fid", "edge2fid"]]
        turn_restrictions2 = turn_restrictions.loc[turn_restrictions.edge1end=="N", ["edge1fid", "edge2fid"]].rename(columns={"edge1fid":"edge2fid", "edge2fid":"edge1fid"})
        turn_restrictions = pd.concat([turn_restrictions1, turn_restrictions2], axis=0, ignore_index=True)
        lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = "idx", right_on = "edge1fid", how = "inner")
#            lenker_med_restriction = veger_edges[veger_edges.idx.isin(turn_restrictions.edge1fid)]
    # hvis 2022
    else:    
        lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = "linkid", right_on = "fromlinkid", how = "inner")
#      lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = ["source", "target"], right_on = ["fromfromnode", "fromtonode"], how = "inner")
    #    lenker_med_restriction2 = veger_edges.merge(turn_restrictions, left_on = ["source", "target", "linkid"], right_on = ["fromfromnode", "fromtonode", "fromlinkid"], how = "inner")
        
    # gjør lenkene med restrictions til første del av nye dobbellenker som skal lages
    lenker_med_restriction = (lenker_med_restriction
            .rename(columns={"target": "middlenode", "minutter_bil": "minutter1", "meter": "meter1", "geometry": "geom1", "idx": "edge1fid"}) 
            .loc[:, ["source", "source_wkt", "middlenode", "minutter1", "meter1", "geom1", "edge1fid"]] )
    # klargjør tabell som skal bli andre del av dobbellenkene
    restrictions = (veger_edges
            .copy()
            .rename(columns={"source": "middlenode", "minutter_bil": "minutter2", "meter": "meter2", "geometry": "geom2", "idx": "edge2fid"})
            .loc[:, ["middlenode","target", "target_wkt", "minutter2", "meter2", "geom2", "edge2fid"]] )

    # koble basert på den nye kolonnen 'middlenode', som blir midterste node i dobbellenkene
    fra_noder_med_restriction = lenker_med_restriction.merge(restrictions, on = "middlenode", how = "inner")

    # vi har nå alle dobbellenker som starter der et svingforbud starter. 
    # velg nå ut kun dobbellenkene det faktisk er svingforbud ()
    if "edge1fid" in turn_restrictions.columns:
        dobbellenker = fra_noder_med_restriction[~((fra_noder_med_restriction["edge1fid"].isin(turn_restrictions["edge1fid"])) & 
                                (fra_noder_med_restriction["edge2fid"].isin(turn_restrictions["edge2fid"])))]
    else:
        dobbellenker = fra_noder_med_restriction[~((fra_noder_med_restriction["source"].isin(turn_restrictions["fromfromnode"])) & 
#                                   (fra_noder_med_restriction["middlenode"].isin(turn_restrictions["fromtonode"])) &
                                (fra_noder_med_restriction["target"].isin(turn_restrictions["totonode"])))]

    # smelt lenkeparene sammen
    dobbellenker["minutter_bil"] = dobbellenker["minutter1"] + dobbellenker["minutter2"]
    dobbellenker["meter"] = dobbellenker["meter1"] + dobbellenker["meter2"]
    dobbellenker["geometry"] = pygeos.line_merge(pygeos.union(pygeos.from_shapely(dobbellenker.geom1), pygeos.from_shapely(dobbellenker.geom2)))
    dobbellenker = gpd.GeoDataFrame(dobbellenker, geometry = "geometry", crs = 25833)
    dobbellenker["turn_restriction"] = True

    if "edge1fid" in turn_restrictions.columns:
        veger_edges.loc[(veger_edges["idx"].isin(turn_restrictions["edge1fid"])), "turn_restriction"] = False
    else:
        veger_edges.loc[(veger_edges["linkid"].isin(turn_restrictions["fromlinkid"])), "turn_restriction"] = False

    return gdf_concat([veger_edges, dobbellenker])
       

    