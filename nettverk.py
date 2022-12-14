import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import line_merge, union, shortest_line, force_2d
from sklearn.neighbors import NearestNeighbors
from networkz.hjelpefunksjoner import fjern_tomme_geometrier, gdf_concat


# funksjon som tilpasser vegnettet til classen Graf().
# mange parametre for at den skal kunne brukes for alle nettverk, men burde nok vært splitta litt opp.
def lag_nettverk(veger,
                 directed = True,
                 
                 source: str = "fromnode",
                 target: str = "tonode",
                 linkid: str = "linkid",
                 minutter = ("drivetime_fw", "drivetime_bw"),
                 vegkategori: str = "category",
                 kommunekolonne = "municipality",
                 
                 finn_isolerte = True,
                 utvid: int = None, # meter
                 turn_restrictions = None,
                 stigningsprosent = False,
                 behold_kolonner: list = None,
                 copy = True
                 ):
    
    if behold_kolonner:
        if isinstance(behold_kolonner, str):
            behold_kolonner = [behold_kolonner]
        assert isinstance(behold_kolonner, (list, tuple)), "'behold_kolonner' må være string/liste/tuple med kolonnenavn."
    else:
        behold_kolonner = []
    
    if copy:
        veger_kopi = veger.copy()
    else:
        veger_kopi = veger
        
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
    
    # litt opprydning
    veger_kopi = (veger_kopi
                  .to_crs(25833)
                  .loc[:, ["idx", "source", "target", "linkid", "drivetime_fw", "drivetime_bw", "oneway", "sperring", "category", "KOMMUNENR", "geometry"] + behold_kolonner]
                  .reset_index(drop=True)
                  )
    
    # fra multilinestring til linestring. Og fjerne z-koordinat fordi de ikke trengs
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = line_merge(veger_kopi.geometry)
    
    assert len(veger_kopi), "vegene har 0 rader"
    
    #hvis noen lenker fortsatt er multilinestrings, må de splittes for å ikke ha flere enn to ytterpunkter
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
    
    del ytterpunkter, sirkler
    
    assert len(veger_kopi), "vegene har 0 rader"

    """
    del opp lange linjer?
    while np.max(veglenker.length) > 1001:
        veglenker = kutt_linjer(veglenker, 1000)
    while np.max(veglenker.length) > 301:
        veglenker = kutt_linjer(veglenker, 300)
    while np.max(veglenker.length) > 101:
        veglenker = kutt_linjer(veglenker, 100)
    while np.max(veglenker.length) > 51:
        veglenker = kutt_linjer(veglenker, 50)
    """
    
    # lag kolonner med geometritekst (wkt) for source og target
    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True) # boundary gir multipunkt med to ytterpunkter for hver linje, og explode() (til singlepart) gir en dobbelt så lang tabell med enkeltpunkter
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2] # gjør annenhvert ytterpunkt til source_wkt og target_wkt
    
    if stigningsprosent:
        veger_kopi = beregn_stigningsprosent(veger_kopi, ytterpunkter)       
        behold_kolonner = behold_kolonner + ["stigningsprosent"]

    veger_kopi["geometry"] = force_2d(veger_kopi.geometry)
    
    if not directed or all(veger_kopi["oneway"].isna()):
        veglenker = veger_kopi
        veglenker["minutter"] = np.where(veglenker["drivetime_fw"].fillna(0) > 0,
                                         veglenker["drivetime_fw"], 
                                         veglenker["drivetime_bw"])        
    else:
        veglenker = snu_og_samle_veglenker(veger_kopi, stigningsprosent)
    
    assert len(veglenker), "vegene har 0 rader"
    
    if turn_restrictions:
        veglenker = lag_turn_restrictions(veglenker, turn_restrictions)
    else:
        veglenker["turn_restriction"] = np.nan
    
    #nye node-id-er som følger index (fordi det indexes med numpy arrays i avstand_til_noder())
    veglenker, noder = lag_node_ids(veglenker)
    
    if utvid:
        if not isinstance(utvid, (float, int)):
            raise ValueError("utvid må være et tall (antall meter man vil utvide linjer)")
        tettede_hull = tett_nettverkshull(noder, veglenker, utvid)
        veglenker = gdf_concat([veglenker, tettede_hull])
    
    veglenker["meter"] = veglenker.length

    veglenker = veglenker.loc[(veglenker.minutter > 0) | (veglenker.minutter.isna()),
                                  ["source", "target", "minutter", "meter", "turn_restriction", "oneway", "sperring", vegkategori, "isolert", "KOMMUNENR", 
                                   "source_wkt", "target_wkt", "geometry"] + behold_kolonner]
    
    # fjern kolonner som ikke ble brukt        
    for col in ["minutter", "turn_restriction", "isolert", vegkategori, "sperring", "KOMMUNENR", "oneway"]:
        if len(veglenker[~((veglenker[col].isna()) | (veglenker[col]==0.02) | (veglenker[col]==-1))])==0:
            veglenker = veglenker.drop(col, axis=1)
    
    return veglenker


def beregn_stigningsprosent(veger, ytterpunkter):
    
    assert all(ytterpunkter.has_z), "Vegdataene må ha z-koordinater for å kunne beregne stigning."
    
    hoyde = [z for z in ytterpunkter.geometry.z]
    veger["hoyde_source"], veger["hoyde_target"] = hoyde[0::2], hoyde[1::2]
    
    veger["stigningsprosent"] = (veger.hoyde_target - veger.hoyde_source) / veger.length * 100
    
    veger.loc[(veger.stigningsprosent>100) | (veger.stigningsprosent<-100), "stigningsprosent"] = 0
    
    return veger


def snu_og_samle_veglenker(veger, stigningsprosent):
    #velg ut de enveiskjørte og snu source og target for lenkene som går "feil" vei
    ft = veger.loc[(veger.oneway=="FT") | (veger.oneway=="F")] 
    tf = veger.loc[(veger.oneway=="TF") | (veger.oneway=="T")]
    tf = tf.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})           
    
    #dupliser lenkene som går begge veier og snu source og target i den ene
    begge_retninger1 = veger[veger.oneway=="B"]
    begge_retninger2 = begge_retninger1.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
    
    # lag minutt-kolonne
    begge_retninger1 = begge_retninger1.rename(columns={"drivetime_fw": "minutter"})
    begge_retninger2 = begge_retninger2.rename(columns={"drivetime_bw": "minutter"})
    ft = ft.rename(columns={"drivetime_fw": "minutter"})
    tf = tf.rename(columns={"drivetime_bw": "minutter"})
    
    if stigningsprosent:
        tf["stigningsprosent"] = tf["stigningsprosent"] * -1
        begge_retninger2["stigningsprosent"] = begge_retninger2["stigningsprosent"] * -1
        
    # oneway=="N" er sperringer fram til og med 2021
    n = veger[(veger.oneway=="N")]
    if len(n)>0:
        n["minutter"] = np.where((n["drivetime_fw"].isna()) | (n["drivetime_fw"]==0) | (n["drivetime_fw"]==""),
                                    n["drivetime_bw"],
                                    n["drivetime_fw"])
        n2 = n.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
        veglenker = gdf_concat([begge_retninger1, begge_retninger2, ft, tf, n, n2])
    else:
        veglenker = gdf_concat([begge_retninger1, begge_retninger2, ft, tf])
    
    return veglenker


def gridish(gdf, meter, x2 = False):
    """
    Enkel rutedeling av dataene. For å fleksibelt kunne loope for små områder sånn at ting blir håndterbart. 
    Gir dataene kolonne med avrundede xy-koordinater. Rundes av til valgfritt antall meter.
    x2=True gir en kolonne til med ruter 1/2 hakk nedover og bortover. Hvis grensetilfeller er viktig, kan man loope en gang per rutekategorikolonne. 
    """
    
    # rund ned koordinatene og sett sammen til kolonne
    gdf["gridish"] = round(gdf.geometry.bounds.minx/meter,1).astype(int).astype(str) + "_" + round(gdf.geometry.bounds.miny/meter,1).astype(int).astype(str)
    
    if x2:

        gdf["gridish_x"] = gdf.geometry.bounds.minx / meter
        
        unike_x = gdf["gridish_x"].astype(int).unique()
        unike_x.sort()
        
        for x in unike_x:
            gdf.loc[(gdf["gridish_x"] >= x-0.5) & (gdf["gridish_x"] < x+0.5), "gridish_x2"] = x+0.5

        # samme for y
        gdf["gridish_y"] = gdf.geometry.bounds.miny/meter
        unike_y = gdf["gridish_y"].astype(int).unique()
        unike_y.sort()
        for y in unike_y:
            gdf.loc[(gdf["gridish_y"] >= y-0.5) & (gdf["gridish_y"] < y+0.5), "gridish_y2"] = y+0.5

        gdf["gridish2"] = gdf["gridish_x2"].astype(str) + "_" + gdf["gridish_y2"].astype(str)

        gdf = gdf.drop(["gridish_x","gridish_y","gridish_x2","gridish_y2"], axis=1)
        
    return gdf


def finn_isolerte_nettverk(veger, lengde: int, ruteloop_m: int):
    """ 
    Gir vegdataene kolonnen 'isolert', som indikerer hvilke veger som er adskilt fra hovedvegnettet med bom/sperring e.l..
    Dette gjøres ved å samle veger som nesten overlapper (innen 0.001 meter), så velge ut ansamlingene som er under en viss størrelse og utstrekning.
    Dette er fryktelig tungt, så det gjøres i loop for et lite område av gangen. 
    Så gjentas loopen for områder som er halvveis forskjøvet på grunn av grensetilfeller. 
    Vegene som er med i begge loops, anses som isolerte. 
    """
    
    # gir vegdataene to kolonner med koordinatkategorier. Den andre kolonnen inneholder rutekategorier som er halvveis forskøvet
    veger = gridish(veger, meter = ruteloop_m, x2 = True)
    
    # fjerner sperringer før beregningen
    if "sperring" in veger.columns:
        ikke_sperringer = veger.loc[veger.sperring.astype(int) != 1]
        sperringer = veger.loc[veger.sperring.astype(int) == 1]
    else:
        ikke_sperringer = veger
        sperringer = None
    
    # samle nesten overlappende med buffer, dissolve og explode. Loopes for hver rute.
    def buffdissexp_gridish_loop(veger, sperringer, lengde, kolonne):
        
        # lagrer veg-indeksene (idx) i tuple fordi det tar mindre plass enn lister
        kanskje_isolerte = ()
        for i in veger[kolonne].unique():
            
            vegene = veger.loc[veger[kolonne] == i, ["geometry"]]
            if sperringer is not None:
                sperringene = sperringer.loc[sperringer[kolonne] == i, ["geometry", "idx"]]
            
            vegene["geometry"] = vegene.buffer(0.001, resolution = 1) # lavest mulig oppløsning for å få det fort
            dissolvet = vegene.dissolve()
            singlepart = dissolvet.explode(ignore_index=True)
            
            # velger nettverk under gitt lengde - hvis de er under halvparten av total lengde og mindre utstrekning enn total lengde
            sum_lengde = dissolvet.length.sum()
            singlepart["utstrekning"] = singlepart.convex_hull.area # fordi lange, isolerte veger, gjerne skogsbilveger, kan være brukbare for turer innad i skogen.
            lite_nettverk = singlepart[(singlepart.length < lengde*2) & 
                                       (singlepart.length < sum_lengde*0.5) &
                                       (singlepart["utstrekning"] < sum_lengde) ]
            
            # legg til nye idx-er i tuple-en med kanskje isolerte nettverk
            for geom in lite_nettverk.geometry:
                
                nye_kanskje_isolerte = tuple(veger.loc[veger.within(geom), "idx"])
                
                if sperringer is not None:
                    nye_kanskje_isolerte = nye_kanskje_isolerte + tuple(sperringene.loc[sperringene.intersects(geom), "idx"])
                                
                nye_kanskje_isolerte = tuple(x for x in nye_kanskje_isolerte if x not in kanskje_isolerte)
                
                kanskje_isolerte = kanskje_isolerte + nye_kanskje_isolerte
                                
        return kanskje_isolerte
       
    kanskje_isolerte = buffdissexp_gridish_loop(ikke_sperringer, sperringer, lengde, "gridish")
    kanskje_isolerte2 = buffdissexp_gridish_loop(ikke_sperringer, sperringer, lengde, "gridish2")

    isolerte_nettverk = [idx for idx in kanskje_isolerte 
                         if idx in kanskje_isolerte2]

    veger["isolert"] = np.where(veger.idx.isin(isolerte_nettverk), 
                                1, 0)

    return veger.loc[:, ~veger.columns.str.contains("gridish")]
    
    
def tett_nettverkshull(noder, veger, avstand):
    """ 
    Lager rette linjer der det er små hull i nettverket. 
    Bruker NearestNeighbors fra scikit-learn, fordi det er utrolig raskt. Man trenger ikke loope for områder en gang. 
    scikit-learn bruker numpy arrays, som må konverteres tilbake til geopandas via index-ene.
    """
    
    crs = noder.crs
    
    # velger ut nodene som kun finnes i én lenke. Altså blindveier i en nettverksanalyse.
    blindveger = noder[noder["n"] <= 1]
    
    # viktig å nullstille index siden sklearn kneighbors gir oss en numpy.array med indekser
    blindveger = blindveger.reset_index(drop=True)
    
    if len(blindveger) <= 1:
        blindveger["minutter"] = -1
        return blindveger
    
    # koordinater i tuple som numpy array
    blindveger_array = np.array([(x, y) for x, y in zip(blindveger.geometry.x, blindveger.geometry.y)])
    
    # finn nærmeste to naboer
    nbr = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(blindveger_array)
    avstander, idxs = nbr.kneighbors(blindveger_array)
    
    # velg ut nest nærmeste (nærmeste er fra og til samme punkt)
    avstander = avstander[:,1]
    idxs = idxs[:,1]
    
    """
    Nå har vi 1d-numpy arrays av lik lengde som blindvegene. 
    'idxs' inneholder numpy-indeksen for vegen som er nærmest, altså endepunktene for de nye lenkene.
    For å konvertere dette fra numpy til geopandas, trengs geometri og node-id. 
    """
    
    # fra-geometrien kan hentes direkte siden avstandene og blindvegene har samme rekkefølge
    fra = np.array([geom.wkt for dist, geom in zip(avstander, blindveger.geometry)
                            if dist < avstand and dist > 0])
    
    # til-geometrien må hentes via index-en
    til = np.array([blindveger.loc[blindveger.index==idx, "wkt"].iloc[0]
                            for dist, idx in zip(avstander, idxs)
                            if dist < avstand and dist > 0])


    # lag GeoDataFrame med rette linjer
    fra = gpd.GeoSeries.from_wkt(fra, crs=crs)
    til = gpd.GeoSeries.from_wkt(til, crs=crs)
    nye_lenker = shortest_line(fra, til)
    nye_lenker = gpd.GeoDataFrame({"geometry": nye_lenker}, geometry="geometry", crs=crs)
    
    # så lage resten av kolonnene
    
    # lager wkt på samme måte som for hele nettverket. Prøvde å bruke array-ene over, men tror floats med .0 der blir til int.
    ytterpunkter = nye_lenker.geometry.boundary.explode(ignore_index=True)
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    nye_lenker["source_wkt"], nye_lenker["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2]

    wkt_id_dict = {wkt: id for wkt, id in zip(blindveger["wkt"], blindveger["node_id"])}
    nye_lenker["source"] = nye_lenker["source_wkt"].map(wkt_id_dict)
    nye_lenker["target"] = nye_lenker["target_wkt"].map(wkt_id_dict)

    assert not any(nye_lenker.source.isna()) and not any(nye_lenker.target.isna()), "Manglende source/target-id-er."
    
    wkt_kat_dict = {wkt: kat for wkt, kat in zip(veger["source_wkt"], veger["category"])}
    nye_lenker["category"] = nye_lenker["source_wkt"].map(wkt_kat_dict)
    
    nye_lenker["sperring"] = 1
    nye_lenker["minutter"] = 0.02
    
    return nye_lenker


def tilpass_veger_sykkelfot(veger):
    veger = veger.loc[veger.sykkelforbud != "Ja"]
    veger = veger.rename(columns={"trafikkretning": "oneway"})
    veger["oneway"] = veger.oneway.map({"MED": "FT", "MOT": "TF", "BEGGE": "B"})
    return veger


def lag_turn_restrictions(veger, turn_restrictions):
        
    # FID starter på  1
    veger["fid"] = veger["idx"] + 1

    turn_restrictions.columns = [col.lower() for col in turn_restrictions.columns]

    for col in turn_restrictions.columns:
        turn_restrictions[col] = turn_restrictions[col].astype(str)
                
    # hvis 2021-data
    if "edge1fid" in turn_restrictions.columns:
        veger["fid"] = veger["fid"].astype(str)
        turn_restrictions1 = turn_restrictions.loc[turn_restrictions.edge1end=="Y", ["edge1fid", "edge2fid"]].rename(columns={"edge1fid":"edge2fid", "edge2fid":"edge1fid"})
        turn_restrictions2 = turn_restrictions.loc[turn_restrictions.edge1end=="N", ["edge1fid", "edge2fid"]]
        turn_restrictions = pd.concat([turn_restrictions1, turn_restrictions2], axis=0, ignore_index=True)
        lenker_med_restriction = veger.merge(turn_restrictions, left_on = "fid", right_on = "edge1fid", how = "inner")
        
    # hvis 2022
    else:    
        veger["linkid"] = veger["linkid"].astype(str)
        lenker_med_restriction = veger.merge(turn_restrictions, left_on = "linkid", right_on = "fromlinkid", how = "inner")
#      lenker_med_restriction = veger.merge(turn_restrictions, left_on = ["source", "target"], right_on = ["fromfromnode", "fromtonode"], how = "inner")
    #    lenker_med_restriction2 = veger.merge(turn_restrictions, left_on = ["source", "target", "linkid"], right_on = ["fromfromnode", "fromtonode", "fromlinkid"], how = "inner")
        
    # gjør lenkene med restrictions til første del av nye dobbellenker som skal lages
    lenker_med_restriction = (lenker_med_restriction
                              .drop("edge1fid", axis=1, errors="ignore")
                              .rename(columns={"target": "middlenode", "minutter": "minutter1", "meter": "meter1", "geometry": "geom1", "fid": "edge1fid"})
                              .loc[:, ["source", "source_wkt", "middlenode", "minutter1", "meter1", "geom1", "edge1fid"]] )
    # klargjør tabell som skal bli andre del av dobbellenkene
    restrictions = (veger
                    .rename(columns={"source": "middlenode", "minutter": "minutter2", "meter": "meter2", "geometry": "geom2", "fid": "edge2fid"})
                    .loc[:, ["middlenode","target", "target_wkt", "minutter2", "meter2", "geom2", "edge2fid"]] 
                    )

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
    dobbellenker["geometry"] = line_merge(union(dobbellenker.geom1, dobbellenker.geom2))
    dobbellenker = gpd.GeoDataFrame(dobbellenker, geometry = "geometry", crs = 25833)
    dobbellenker["turn_restriction"] = True
    
    if "edge1fid" in turn_restrictions.columns:
        veger.loc[(veger["fid"].isin(turn_restrictions["edge1fid"])), "turn_restriction"] = False
    else:
        veger.loc[(veger["linkid"].isin(turn_restrictions["fromlinkid"])), "turn_restriction"] = False

    return gdf_concat([veger, dobbellenker])


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
    
    #koble på de nye node-id-ene
    nodeordbok = {wkt: node_id for wkt, node_id in zip(noder["wkt"], noder["node_id"])}
    veger["source"] = veger["source_wkt"].map(nodeordbok)
    veger["target"] = veger["target_wkt"].map(nodeordbok)
    
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