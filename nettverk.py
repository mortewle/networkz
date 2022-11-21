import numpy as np
import pandas as pd
import geopandas as gpd
import pygeos
from networkz.stottefunksjoner import fjern_tomme_geometrier, gdf_concat, kutt_linjer
       

# funksjon som tilpasser vegnettet til funksjonen graf(), som bygger grafen.
# outputen her kan lagres i parquet så man slipper å lage nettverket hver gang.
def lag_nettverk(veger,
                 source: str = "fromnodeid",
                 target: str = "tonodeid",
                 linkid: str = "linkid",
                 minutter = ("drivetime_fw", "drivetime_bw"),
                 vegkategori: str = "category",
                 kommunekolonne = "municipality",
                 isolerte_nettverk = True,
                 turn_restrictions = None):
    
    veger_kopi = veger.copy()
    
    veger_kopi["idx"] = veger_kopi.index

    # endre til små bokstaver
    veger_kopi.columns = [col.lower() for col in veger_kopi.columns]

    # finn kolonner
    veger_kopi["source"] = finn_source(veger_kopi, source)
    veger_kopi["target"] = finn_target(veger_kopi, target)
    veger_kopi["linkid"] = finn_linkid(veger_kopi, linkid)
    veger_kopi["category"] = finn_vegkategori(veger_kopi, vegkategori)   
    veger_kopi["KOMMUNENR"] = finn_og_omkod_kommunekolonne(veger_kopi, kommunekolonne)
    veger_kopi["drivetime_fw"], veger_kopi["drivetime_bw"] = finn_minutter(veger_kopi, minutter)
    
    if not "sperring" in veger_kopi.columns:
        veger_kopi["sperring"] = -1
        
    # litt opprydning (til utm33-koordinater, endre navn på kolonner, beholde kun relevante kolonner, resette index)
    veger_kopi = (veger_kopi
                  .to_crs(25833)
                  [["idx", "source", "target", "linkid", "drivetime_fw", "drivetime_bw", "oneway", "sperring", "category", "KOMMUNENR", "geometry"]]
                  .reset_index(drop=True) 
                  )
    
    # fra multilinestring til linestring. Og fjerne z-koordinat fordi de ikke trengs
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = pygeos.line_merge(pygeos.force_2d(pygeos.from_shapely(veger_kopi.geometry)))   
    
    #hvis noen lenker fortsatt er multilinestrings, må de splittes for å ikke ha flere enn to ytterpunkter snart
    n = len(veger_kopi)
    veger_kopi = veger_kopi.explode(ignore_index=True)
    if len(veger_kopi)<n:
        print(f"Advarsel: {len(veger_kopi)-n} multigeometrier ble splittet. Minutt-kostnader blir feil for disse.")
        #TODO: lag ny minutt-kolonne manuelt
    
    if isolerte_nettverk:
        veger_kopi = finn_isolerte_nettverk(veger_kopi, 
                                            storrelse=10000, 
                                            ruteloop=2250)
    else:
        veger_kopi["isolert"] = np.nan
        
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
    begge_retninger1["minutter"] = begge_retninger1["drivetime_fw"]
    begge_retninger2["minutter"] = begge_retninger2["drivetime_bw"]
    ft["minutter"] = ft["drivetime_fw"]
    tf["minutter"] = tf["drivetime_bw"]

    n = veger_kopi[(veger_kopi.oneway=="N")]
    if len(n)>0:
        n["minutter"] = np.where((n["drivetime_fw"].isna()) | (n["drivetime_fw"]==0) | (n["drivetime_fw"]==""),
                                     n["drivetime_bw"],
                                     n["drivetime_fw"])
        n2 = n.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
        veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf, n, n2])
    else:
        veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf])

    # lag meter-kolonne
    veger_edges["meter"] = veger_edges.length
    
    # TODO: fullfør denne
    if turn_restrictions is not None:
        veger_edges = turn_restr(veger_edges, turn_restrictions)
    else:
        veger_edges["turn_restriction"] = False

    # rydd opp (fjern eventuelle 0-verdier, velg ut kolonner, fjern duplikat-lenkene med høyest kostnad, reset index)
    veger_edges = veger_edges.loc[veger_edges["minutter"] >= 0, 
                                  ["idx", "source", "target", "source_wkt", "target_wkt", "minutter", "meter", "turn_restriction", "sperring", vegkategori, "isolert", "KOMMUNENR", "geometry"]]
    
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
         
    #nye node-id-er som følger index (fordi jeg indexer med numpy arrays i avstand_til_noder())
    veger_edges = make_node_ids(veger_edges)

    # category-type for å spare plass
    for col in ["source", "target", "source_wkt", "target_wkt", "turn_restriction", vegkategori, "KOMMUNENR"]:
        veger_edges[col] = veger_edges[col].astype(str).astype("category")
            
    return veger_edges


# looping er treigt, men buffer+dissolve for store områder er mye treigere. 
# "deler" derfor dataene i ruter ved å gi kolonner koordinat-kategorier
def koor_kat(gdf, 
             meter, # minst mulig ruter er ikke alltid raskest. Varierer med hvor tunge dataene er.
             x2 = False # x2=True gir en kolonne til med ruter 1/2 hakk nedover og bortover. Hvis grensetilfeller er viktig
             ):
    
    # rund ned
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


def finn_isolerte_nettverk(veger, storrelse, ruteloop):
    
    veger = koor_kat(veger, meter = ruteloop, x2 = True)

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
            lite_nettverk = singlepart[(singlepart.length < storrelse*2) & 
                                       (singlepart.length < sum_lengde*0.5) &
                                       (singlepart["utstrekning"] < sum_lengde)
                                       ] #.unary_union
            for geom in lite_nettverk.geometry:
                nye_kanskje_isolerte =  tuple(veger.loc[veger.within(geom), "idx"]) + tuple(sperringene.loc[sperringene.intersects(geom), "idx"])
                nye_kanskje_isolerte = tuple(x for x in nye_kanskje_isolerte if x not in kanskje_isolerte)
                kanskje_isolerte = kanskje_isolerte + nye_kanskje_isolerte
                storrelsene = storrelsene + tuple(geom.length for _ in range(len(nye_kanskje_isolerte)))
        return kanskje_isolerte, storrelsene
       
    kanskje_isolerte, storrelsene = kat_loop(ikke_sperringer, sperringer, "koor_kat")
    kanskje_isolerte2, storrelsene2 = kat_loop(ikke_sperringer, sperringer, "koor_kat2")
    
#    isolerte_nettverk = [x for x in kanskje_isolerte if x in kanskje_isolerte2]

    # dict med id-er og lengder
    kanskje_isolerte = {x: lengde for x, lengde in zip(kanskje_isolerte, storrelsene) }

    # behold bare id-ene som er små nettverk i begge rutelooper. Og velg maks. lengde for hver id
    isolerte_nettverk = {i: max([lengde, kanskje_isolerte[i]]) 
                        for i, lengde in  zip(kanskje_isolerte2, storrelsene2)
                        if i in kanskje_isolerte}

                
#    veger.loc[veger.idx.isin(isolerte_nettverk), "isolert"] = 1

    # lag kolonne med lengden for id-ene som er isolerte. Ellers blir isolert 0
    veger.loc[veger.idx.isin(isolerte_nettverk), "isolert"] = [isolerte_nettverk[x] for x in veger.loc[veger.idx.isin(isolerte_nettverk), "idx"] ]

    veger.loc[~veger.idx.isin(isolerte_nettverk), "isolert"] = 0
    
    return veger
    

def turn_restr(veger, turn_restrictions):
    
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
    dobbellenker["geometry"] = pygeos.line_merge(pygeos.union(pygeos.from_shapely(dobbellenker.geom1), pygeos.from_shapely(dobbellenker.geom2)))
    dobbellenker = gpd.GeoDataFrame(dobbellenker, geometry = "geometry", crs = 25833)
    dobbellenker["turn_restriction"] = True
    
    if "edge1fid" in turn_restrictions.columns:
        veger_edges.loc[(veger_edges["idx"].isin(turn_restrictions["edge1fid"])), "turn_restriction"] = False
    else:
        veger_edges.loc[(veger_edges["linkid"].isin(turn_restrictions["fromlinkid"])), "turn_restriction"] = False

    return gdf_concat([veger_edges, dobbellenker])
 

# lag retningsløst nettverk
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


def finn_source(veger, source):
    # hvis ikke angitte kolonner finnes i vegdataene, sjekk om andre kolonner matcher. 
    # lager ny kolonne hvis ingen matcher. Gir feilmelding hvis flere enn en matcher.
    if not source in veger.columns:
        n = 0
        for col in veger.columns:
            if "from" in col and "node" in col or "source" in col:
                source = col
                n += 1
        if n == 1:
            print(f"Bruker '{source}' som source-kolonne")
        elif n == 0:
            veger[source] = np.nan
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde source-id-er")
    return veger[source]


def finn_target(veger, target):
    if not target in veger.columns:
        n = 0
        for col in veger.columns:
            if "to" in col and "node" in col or "target" in col:
                target = col
                n += 1
        if n == 1:
            print(f"Bruker '{target}' som target-kolonne")
        elif n == 0:
            veger[target] = np.nan
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde target-id-er")
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
            veger[linkid] = np.nan #godta dette eller raise ValueError("Finner ikke linkid-kolonne") ?
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
        raise ValueError("Finner ikke kolonner med minutter")


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
        raise ValueError("Finner ikke vegkategori-kolonne")


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
        return 0