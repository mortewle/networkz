import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from sklearn.neighbors import NearestNeighbors
from networkz.hjelpefunksjoner import fjern_tomme_geometrier, gdf_concat


# funksjon som tilpasser vegnettet til classen Graf().
def lag_nettverk(veger,
                 directed = True,
                 source: str = "fromnodeid",
                 target: str = "tonodeid",
                 linkid: str = "linkid",
                 minutter = ("drivetime_fw", "drivetime_bw"),
                 vegkategori: str = "category",
                 kommunekolonne = "municipality",
                 finn_isolerte = True,
                 utvid: int = None, # meter
                 turn_restrictions = None,
                 behold_kolonner: list = None
                 ):
    
    if behold_kolonner:
        if isinstance(behold_kolonner, str):
            behold_kolonner = [behold_kolonner]
        assert isinstance(behold_kolonner, (list, tuple)), "'behold_kolonner' må være en liste/tuple med kolonnenavn."
    else:
        behold_kolonner = []
    
    veger_kopi = veger.copy()
    
    veger_kopi["idx"] = veger_kopi.index

    veger_kopi.columns = [col.lower() for col in veger_kopi.columns]
    
    # hvis ikke angitte kolonner finnes i vegdataene, sjekk om andre kolonner matcher. 
    # lager ny kolonne hvis ingen matcher. Gir feilmelding hvis flere enn én matcher.
    veger_kopi["source"] = finn_source(veger_kopi, source)
    veger_kopi["target"] = finn_target(veger_kopi, target)
    veger_kopi["linkid"] = finn_linkid(veger_kopi, linkid)
    veger_kopi["category"] = finn_vegkategori(veger_kopi, vegkategori)   
    veger_kopi["KOMMUNENR"] = finn_og_omkod_kommunekolonne(veger_kopi, kommunekolonne)
    veger_kopi["drivetime_fw"], veger_kopi["drivetime_bw"] = finn_minutter(veger_kopi, minutter)
    
    if not "oneway" in veger_kopi.columns:
        veger_kopi["oneway"] = np.nan
    if not "sperring" in veger_kopi.columns:
        veger_kopi["sperring"] = -1
    
    # litt opprydning (til utm33-koordinater, endre navn på kolonner, beholde kun relevante kolonner, resette index)
    veger_kopi = (veger_kopi
                  .to_crs(25833)
                  .loc[:, ["idx", "source", "target", "linkid", "drivetime_fw", "drivetime_bw", "oneway", "sperring", "category", "KOMMUNENR", "geometry"] + behold_kolonner]
                  .reset_index(drop=True)
                  )
    
    # fra multilinestring til linestring. Og fjerne z-koordinat fordi de ikke trengs
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = shapely.line_merge(shapely.force_2d(veger_kopi.geometry)) 
    
    assert len(veger_kopi)>0

    #hvis noen lenker fortsatt er multilinestrings, må de splittes for å ikke ha flere enn to ytterpunkter snart
    n = len(veger_kopi)
    veger_kopi = veger_kopi.explode(ignore_index=True)
    if len(veger_kopi)<n:
        print(f"Advarsel: {n-len(veger_kopi)} multigeometrier ble splittet. Minutt-kostnader blir feil for disse.")
        #TODO: lag ny minutt-kolonne manuelt
    
    if finn_isolerte:
        veger_kopi = finn_isolerte_nettverk(veger_kopi, 
                                            lengde=10000,
                                            ruteloop_m=2250)
    else:
        veger_kopi["isolert"] = np.nan
    
    # hent ut linjenes ytterpunkter. 
    # men først: sirkler har ingen ytterpunkter og heller ingen funksjon. Disse må fjernes
    ytterpunkter = veger_kopi.copy()
    ytterpunkter["geometry"] = ytterpunkter.geometry.boundary 
    sirkler = ytterpunkter.loc[ytterpunkter.is_empty, "idx"] #sirkler har tom boundary
    veger_kopi = veger_kopi[~veger_kopi.idx.isin(sirkler)]
    
    assert len(veger_kopi)>0

    """
    while np.max(veger_edges.length) > 1001:
        veger_edges = kutt_linjer(veger_edges, 1000)
    while np.max(veger_edges.length) > 301:
        veger_edges = kutt_linjer(veger_edges, 300)
    while np.max(veger_edges.length) > 101:
        veger_edges = kutt_linjer(veger_edges, 100)
    while np.max(veger_edges.length) > 51:
        veger_edges = kutt_linjer(veger_edges, 50)
    """
    
    # lag kolonner med geometritekst (wkt) for source og target
    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True) # boundary gir multipunkt med to ytterpunkter for hver linje, og explode() (til singlepart) gir en dobbelt så lang tabell med enkeltpunkter
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2] # gjør annenhvert ytterpunkt til source_wkt og target_wkt

    if not directed or all(veger_kopi["oneway"].isna()):
        veger_edges = veger_kopi.copy()
        veger_edges["minutter"] = np.where(veger_edges["drivetime_fw"].fillna(0) > 0,
                                           veger_edges["drivetime_fw"], 
                                           veger_edges["drivetime_bw"]) 
        assert len(veger_edges)>0

    else:
        #velg ut de enveiskjørte og snu source og target for lenkene som går "feil" vei
        ft = veger_kopi.loc[(veger_kopi.oneway=="FT") | (veger_kopi.oneway=="F")] 
        tf = veger_kopi.loc[(veger_kopi.oneway=="TF") | (veger_kopi.oneway=="T")]
        tf = tf.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})       
        
        #dupliser lenkene som går begge veier og snu source og target i den ene
        begge_retninger1 = veger_kopi[veger_kopi.oneway=="B"]
        begge_retninger2 = begge_retninger1.copy()
        begge_retninger2 = begge_retninger2.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
        
        # lag minutt-kolonne
        begge_retninger1 = begge_retninger1.rename(columns={"drivetime_fw": "minutter"})
        begge_retninger2 = begge_retninger2.rename(columns={"drivetime_bw": "minutter"})
        ft = ft.rename(columns={"drivetime_fw": "minutter"})
        tf = tf.rename(columns={"drivetime_bw": "minutter"})

        # oneway=="N" er sperringer fram til og med 2021
        n = veger_kopi[(veger_kopi.oneway=="N")]
        if len(n)>0:
            n["minutter"] = np.where((n["drivetime_fw"].isna()) | (n["drivetime_fw"]==0) | (n["drivetime_fw"]==""),
                                        n["drivetime_bw"],
                                        n["drivetime_fw"])
            n2 = n.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
            veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf, n, n2])
        else:
            veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf])
    
        assert len(veger_edges)>0

    if turn_restrictions:
        veger_edges = lag_turn_restrictions(veger_edges, turn_restrictions)
    else:
        veger_edges["turn_restriction"] = np.nan
    
    #nye node-id-er som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder())
    veger_edges, noder = lag_node_ids(veger_edges)

    if utvid:
        if not isinstance(utvid, (float, int)):
            raise ValueError("utvid må være et tall (antall meter man vil utvide linjer)")
        tettede_hull = tett_nettverkshull(noder, veger_edges, utvid)
        veger_edges = gdf_concat([veger_edges, tettede_hull])
    
    veger_edges["meter"] = veger_edges.length

    veger_edges = veger_edges.loc[(veger_edges.minutter > 0) | (veger_edges.minutter.isna()),
                                  ["source", "target", "minutter", "meter", "turn_restriction", "oneway", "sperring", vegkategori, "isolert", "KOMMUNENR", "source_wkt", "target_wkt", "geometry"] + behold_kolonner]
    
    # fjern kolonner som ikke ble brukt        
    for col in ["minutter", "turn_restriction", "isolert", vegkategori, "sperring", "KOMMUNENR", "oneway"]:
        if len(veger_edges[~((veger_edges[col].isna()) | (veger_edges[col]==0.02) | (veger_edges[col]==-1))])==0:
            veger_edges = veger_edges.drop(col, axis=1)
    
    return veger_edges


# looping er treigt, men buffer+dissolve for store områder er mye treigere. 
# "deler" derfor dataene i ruter ved å gi kolonner koordinat-kategorier
def koor_kat(gdf, 
             meter, # minst mulig ruter er ikke alltid raskest. Varierer med hvor tunge dataene er.
             x2 = False # x2=True gir en kolonne til med ruter 1/2 hakk nedover og bortover. Hvis grensetilfeller er viktig
             ):
    
    # rund ned koordinatene og sett sammen til kolonne
    gdf["koor_kat"] = round(gdf.geometry.bounds.minx/meter,1).astype(int).astype(str) + "_" + round(gdf.geometry.bounds.miny/meter,1).astype(int).astype(str)
    
    if x2:

        gdf["koor_kat_x"] = gdf.geometry.bounds.minx / meter
        
        unike_x = gdf["koor_kat_x"].astype(int).unique()
        unike_x.sort()
        
        for x in unike_x:
            gdf.loc[(gdf["koor_kat_x"] >= x-0.5) & (gdf["koor_kat_x"] < x+0.5), "koor_kat_x2"] = x+0.5

        # samme for y
        gdf["koor_kat_y"] = gdf.geometry.bounds.miny/meter
        unike_y = gdf["koor_kat_y"].astype(int).unique()
        unike_y.sort()
        for y in unike_y:
            gdf.loc[(gdf["koor_kat_y"] >= y-0.5) & (gdf["koor_kat_y"] < y+0.5), "koor_kat_y2"] = y+0.5

        gdf["koor_kat2"] = gdf["koor_kat_x2"].astype(str) + "_" + gdf["koor_kat_y2"].astype(str)

        gdf = gdf.drop(["koor_kat_x","koor_kat_y","koor_kat_x2","koor_kat_y2"], axis=1)
        
    return gdf


def finn_isolerte_nettverk(veger, lengde: int, ruteloop_m: int):
    
    """ Gir vegdataene kolonnen 'isolert', hvor 0 betyr ikke isolert, og andre tall betyr størrelsen på det isolerte nettverket i meter.
    Finner de isolerte ved å bufre 0.001 meter, dissolve og explode (til singlepart). 
    Dette gjøres i loop for en rute av gangen, hvor index-verdier (idx) under angitt lengde lagres i tuple (fordi det tar mindre plass enn lister).
    Så gjentas loopen pga grensetilfeller basert på ruter som er halvveis forskjøvet. idx-er som er med i begge loop-runder, anses som isolerte. """
    
    # gir vegdataene to kolonner med koordinatkategorier. Den andre kolonnen inneholder rutekategorier som er halvveis forskøvet
    veger = koor_kat(veger, meter = ruteloop_m, x2 = True)
    
    if "sperring" in veger.columns:
        ikke_sperringer = veger[veger.sperring.astype(int) != 1]
        sperringer = veger[veger.sperring.astype(int) == 1]
    else:
        ikke_sperringer = veger.copy()
        sperringer = veger.copy()
        
    def kat_loop(veger, sperringer, kolonne):
        kanskje_isolerte = ()
        storrelsene = ()
        for kat in veger[kolonne].unique():
            sperringene = sperringer.loc[sperringer[kolonne] == kat, ["geometry", "idx"]]
            vegene = veger.loc[veger[kolonne] == kat, ["geometry"]]
            vegene["geometry"] = vegene.buffer(0.001, resolution = 1) # lavest mulig resolution for å få det fort
            dissolvet = vegene.dissolve()
            sum_lengde = dissolvet.length.sum()
            singlepart = dissolvet.explode(ignore_index=True)
            singlepart["utstrekning"] = singlepart.convex_hull.area
            lite_nettverk = singlepart[(singlepart.length < lengde*2) & 
                                       (singlepart.length < sum_lengde*0.5) &
                                       (singlepart["utstrekning"] < sum_lengde)
                                       ]
            for geom in lite_nettverk.geometry:
                nye_kanskje_isolerte =  tuple(veger.loc[veger.within(geom), "idx"]) + tuple(sperringene.loc[sperringene.intersects(geom), "idx"])
                nye_kanskje_isolerte = tuple(x for x in nye_kanskje_isolerte if x not in kanskje_isolerte)
                kanskje_isolerte = kanskje_isolerte + nye_kanskje_isolerte
                storrelsene = storrelsene + tuple(geom.length for _ in range(len(nye_kanskje_isolerte)))
        return kanskje_isolerte, storrelsene
       
    kanskje_isolerte, storrelsene = kat_loop(ikke_sperringer, sperringer, "koor_kat")
    kanskje_isolerte2, storrelsene2 = kat_loop(ikke_sperringer, sperringer, "koor_kat2")
    
    # dict med id-er og lengder
    kanskje_isolerte = {x: lengde for x, lengde in zip(kanskje_isolerte, storrelsene) }

    # behold bare id-ene som er små nettverk i begge rutelooper. Og velg maks. lengde for hver id
    isolerte_nettverk = {i: max([lengde, kanskje_isolerte[i]]) 
                        for i, lengde in  zip(kanskje_isolerte2, storrelsene2)
                        if i in kanskje_isolerte}

    # lag kolonne med lengden for id-ene som er isolerte. Ellers blir isolert 0
    veger.loc[veger.idx.isin(isolerte_nettverk), "isolert"] = [isolerte_nettverk[x] for x in veger.loc[veger.idx.isin(isolerte_nettverk), "idx"] ]

    veger.loc[~veger.idx.isin(isolerte_nettverk), "isolert"] = 0
    
    return veger
    
    
def tett_nettverkshull(noder, veger, avstand, crs=25833):
    
    """ Lager rette linjer der det er små hull i nettverket."""
        
    blindveger = noder[noder["n"] <= 1]
    blindveger = blindveger.reset_index(drop=True)

    if len(blindveger) <= 1:
        blindveger["minutter"] = -1
        return blindveger
    
    # koordinater til numpy array
    blindveger_array = np.array([(x, y) for x, y in zip(blindveger.geometry.x, blindveger.geometry.y)])
    
    # finn nærmeste to naboer og velg ut nest nærmeste (nærmeste er fra og til samme punkt)
    nbr = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(blindveger_array)
    avstander, idxs = nbr.kneighbors(blindveger_array)
    avstander = avstander[:,1]
    idxs = idxs[:,1]
    
    # start- og sluttgeometrier for de nye lenkene, hvis under ønsket avstand
    fra = np.array([geom.wkt for dist, geom in zip(avstander, blindveger.geometry)
                            if dist < avstand and dist > 0])
    
    til = np.array([blindveger.loc[blindveger.index==idx, "wkt"].iloc[0]
                            for dist, idx in zip(avstander, idxs)
                            if dist < avstand and dist > 0])

    # lag GeoDataFrame med rette linjer
    fra =  gpd.GeoSeries.from_wkt(fra, crs=crs)
    til =  gpd.GeoSeries.from_wkt(til, crs=crs)
    nye_lenker = pd.DataFrame()
    nye_lenker["geometry"] = shapely.shortest_line(fra, til)
    nye_lenker = gpd.GeoDataFrame(nye_lenker, geometry="geometry", crs=crs)
    
    # lag alle de andre kolonnene
    nye_lenker["sperring"] = 1
    nye_lenker["minutter"] = 0.02
    nye_lenker["turn_restriction"] = np.nan

    ytterpunkter = nye_lenker.geometry.boundary.explode(ignore_index=True) 
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    nye_lenker["source_wkt"], nye_lenker["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2]
    
    wkt_id_dict = {wkt: id for wkt, id in zip(blindveger["wkt"], blindveger["node_id"])}
    nye_lenker["source"] = nye_lenker["source_wkt"].map(wkt_id_dict)
    nye_lenker["target"] = nye_lenker["target_wkt"].map(wkt_id_dict)
    
    wkt_kat_dict = {wkt: kat for wkt, kat in zip(veger["source_wkt"], veger["category"])}
    nye_lenker["category"] = nye_lenker["source_wkt"].map(wkt_kat_dict)
    
    return nye_lenker


def lag_turn_restrictions(veger, turn_restrictions):
    
    veger_edges = veger.copy()
    
    # FID starter på  1
    veger_edges["idx"] = veger_edges["idx"] + 1

    turn_restrictions.columns = [col.lower() for col in turn_restrictions.columns]

    for col in turn_restrictions.columns:
        turn_restrictions[col] = turn_restrictions[col].astype(str)
                
    # hvis 2021-data
    if "edge1fid" in turn_restrictions.columns:
        veger_edges["idx"] = veger_edges["idx"].astype(str)
        turn_restrictions1 = turn_restrictions.loc[turn_restrictions.edge1end=="Y", ["edge1fid", "edge2fid"]].rename(columns={"edge1fid":"edge2fid", "edge2fid":"edge1fid"})
        turn_restrictions2 = turn_restrictions.loc[turn_restrictions.edge1end=="N", ["edge1fid", "edge2fid"]]
        turn_restrictions = pd.concat([turn_restrictions1, turn_restrictions2], axis=0, ignore_index=True)
        lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = "idx", right_on = "edge1fid", how = "inner")
        
    # hvis 2022
    else:    
        veger_edges["linkid"] = veger_edges["linkid"].astype(str)
        lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = "linkid", right_on = "fromlinkid", how = "inner")
#      lenker_med_restriction = veger_edges.merge(turn_restrictions, left_on = ["source", "target"], right_on = ["fromfromnode", "fromtonode"], how = "inner")
    #    lenker_med_restriction2 = veger_edges.merge(turn_restrictions, left_on = ["source", "target", "linkid"], right_on = ["fromfromnode", "fromtonode", "fromlinkid"], how = "inner")
        
    # gjør lenkene med restrictions til første del av nye dobbellenker som skal lages
    lenker_med_restriction = (lenker_med_restriction
                              .drop("edge1fid", axis=1, errors="ignore")
                              .rename(columns={"target": "middlenode", "minutter": "minutter1", "meter": "meter1", "geometry": "geom1", "idx": "edge1fid"})
                              .loc[:, ["source", "source_wkt", "middlenode", "minutter1", "meter1", "geom1", "edge1fid"]] )
    # klargjør tabell som skal bli andre del av dobbellenkene
    restrictions = (veger_edges
            .copy()
            .rename(columns={"source": "middlenode", "minutter": "minutter2", "meter": "meter2", "geometry": "geom2", "idx": "edge2fid"})
            .loc[:, ["middlenode","target", "target_wkt", "minutter2", "meter2", "geom2", "edge2fid"]] )

    # koble basert på den nye kolonnen 'middlenode', som blir midterste node i dobbellenkene
    fra_noder_med_restriction = lenker_med_restriction.merge(restrictions, 
                                                             on = "middlenode", 
                                                             how = "inner")
    
    # vi har nå alle dobbellenker som starter der et svingforbud starter. 
    # fjern nå dobbellenkene det faktisk er svingforbud
    if "edge1fid" in turn_restrictions.columns:
        dobbellenker = fra_noder_med_restriction[
                            ~((fra_noder_med_restriction["edge1fid"].isin(turn_restrictions["edge1fid"])) &
                                (fra_noder_med_restriction["edge2fid"].isin(turn_restrictions["edge2fid"])))]
    else:
        dobbellenker = fra_noder_med_restriction[
                            ~((fra_noder_med_restriction["source"].isin(turn_restrictions["fromfromnode"])) &
#                                   (fra_noder_med_restriction["middlenode"].isin(turn_restrictions["fromtonode"])) &
                                (fra_noder_med_restriction["target"].isin(turn_restrictions["totonode"])))]
    
    # smelt lenkeparene sammen
    dobbellenker["minutter"] = dobbellenker["minutter1"] + dobbellenker["minutter2"]
    dobbellenker["meter"] = dobbellenker["meter1"] + dobbellenker["meter2"]
    dobbellenker["geometry"] = shapely.line_merge(shapely.union(dobbellenker.geom1, dobbellenker.geom2))
    dobbellenker = gpd.GeoDataFrame(dobbellenker, geometry = "geometry", crs = 25833)
    dobbellenker["turn_restriction"] = True
    
    if "edge1fid" in turn_restrictions.columns:
        veger_edges.loc[(veger_edges["idx"].isin(turn_restrictions["edge1fid"])), "turn_restriction"] = False
    else:
        veger_edges.loc[(veger_edges["linkid"].isin(turn_restrictions["fromlinkid"])), "turn_restriction"] = False

    return gdf_concat([veger_edges, dobbellenker])


def lag_node_ids(veger, crs=25833):
    """ Nye node-id-er som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder()) """
    
    if not "oneway" in veger.columns:
        veger["oneway"] = 0
        fjern = True
    else:
        fjern = False
    
    sources = veger[["source_wkt", "oneway"]].rename(columns={"source_wkt":"wkt"})
    targets = veger[["target_wkt", "oneway"]].rename(columns={"target_wkt":"wkt"})
        
    noder = pd.concat([sources, targets], axis=0, ignore_index=True)
    
    # tell opp antall forekomster av hver node (del toveiskjørte på to)
    noder["n"] = np.where((noder["oneway"]=="B") | (noder["oneway"]=="N"), 
                          0.5, 1)
    noder["n"] = noder["wkt"].map(noder[["wkt", "n"]].groupby("wkt").sum()["n"])
    
    noder = noder.drop_duplicates(subset=["wkt"]).reset_index(drop=True) # viktig at node_id følger index for å kunne indekse på numpy arrays i graf()
    
    noder["node_id"] = noder.index
    noder["node_id"] = noder["node_id"].astype(str) # funker ikke med numeriske node-navn i igraph, pussig nok...

    noder = noder[["node_id", "wkt", "n"]]
    
    #koble på de nye node-id-ene
    veger = (veger
            .drop(["source", "target", "nz_idx"], axis=1, errors="ignore")
            .merge(noder, left_on = "source_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"source"})
            .drop(["wkt", "n"], axis=1)
            .merge(noder, left_on = "target_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"target"})
            .drop(["wkt", "n"], axis=1)
    )
    
    veger["meter"] = veger.length
    
    # fjern duplikatlenkene med høyest kostnad (siden disse aldri vil bli brukt)
    if "minutter" in veger.columns:
        veger = veger.sort_values("minutter", ascending=True)
    else:     
        veger = veger.sort_values("meter", ascending=True)
             
    veger = veger.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)
    
    if fjern:
        veger = veger.drop("oneway", axis=1)
        
    noder["geometry"] = gpd.GeoSeries.from_wkt(noder.wkt, crs=crs)
    noder = gpd.GeoDataFrame(noder, geometry="geometry", crs=crs)
    noder = noder.reset_index(drop=True)
    
    return veger, noder


def finn_source(veger, source):
    if not source in veger.columns:
        mulige = []
        for col in veger.columns:
            if "from" in col and "node" in col or "source" in col:
                source = col
                mulige.append(col)
            elif "start" in col and "node" in col or "fra" in col and "node" in col:
                source = col
                mulige.append(col)
        if len(mulige) == 1:
            print(f"Bruker '{source}' som source-kolonne")
        elif len(mulige) == 0:
            veger[source] = np.nan
        elif len(mulige) > 1:
            raise ValueError(f"Flere kolonner kan inneholde source-id-er: {', '.join(mulige)}")
    return veger[source]


def finn_target(veger, target):
    if not target in veger.columns:
        mulige = []
        for col in veger.columns:
            if "to" in col and "node" in col or "target" in col:
                target = col
                mulige.append(col)
            elif "end" in col and "node" in col or "slutt" in col and "node" in col:
                target = col
                mulige.append(col)
        if len(mulige) == 1:
            print(f"Bruker '{target}' som target-kolonne")
        elif len(mulige) == 0:
            veger[target] = np.nan
        elif len(mulige) > 1:
            raise ValueError(f"Flere kolonner kan inneholde target-id-er: {', '.join(mulige)}")
    return veger[target]


def finn_linkid(veger, linkid):
    if not linkid in veger.columns:
        n = 0
        for col in veger.columns:
            if "link" in col and "id" in col:
                linkid = col
                n += 1
        if n == 1:
            print(f"Bruker '{linkid}' som linkid-kolonne")
        elif n == 0:
            veger[linkid] = np.nan
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde linkid-id-er")
    return veger[linkid]


def finn_minutter(veger, minutter) -> tuple:
    if isinstance(minutter, str) and minutter in veger.columns:
        return veger[minutter], veger[minutter]
    elif minutter[0] in veger.columns and minutter[1] in veger.columns:
        return veger[minutter[0]], veger[minutter[1]]
    elif "drivetime_fw" in veger.columns and "drivetime_bw" in veger.columns:
        return veger["drivetime_fw"], veger["drivetime_bw"]
    elif "ft_minutes" in veger.columns and "tf_minutes" in veger.columns:
        return veger["ft_minutes"], veger["tf_minutes"]
    else:
        return np.nan, np.nan


def finn_vegkategori(veger, vegkategori):
    if vegkategori in veger.columns:
        return veger[vegkategori]
    elif "category" in veger.columns:
        return veger["category"]
    elif "vegtype" in veger.columns:
        return veger["vegtype"]
    elif "roadid" in veger.columns:
        veger["category"] = veger["roadid"].map(lambda x: x.replace('{','').replace('}','')[0])
        return veger["category"] 
    else:
        return np.nan


def finn_og_omkod_kommunekolonne(veger, kommunekolonne):
    if kommunekolonne in veger.columns:
        return veger[kommunekolonne].map(lambda x: str(int(x)).zfill(4)).astype("category")
    n = 0
    for col in veger.columns:
        if "komm" in col or "muni" in col:
            komm_col = col
            n += 1
    if n == 1:
        return veger[komm_col].map(lambda x: str(int(x)).zfill(4)).astype("category")
    else:
        return np.nan