import geopandas as gpd
import pandas as pd
import numpy as np
import pygeos
import igraph
from networkz.stottefunksjoner import fjern_tomme_geometrier, gdf_concat


# funksjon som tilpasser vegnettet til funksjonen graf(), som bygger grafen.
# outputen her kan lagres i parquet så man slipper å lage nettverket hver gang.
def lag_nettverk(veger, 
                      source: str = "fromnodeid",
                      target: str = "tonodeid",
                      linkid: str = "linkid",
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
            if "from" in col and "node" in col:
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
            if "to" in col and "node" in col:
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
            veger_kopi[linkid] = np.nan #midlr.?
            #raise ValueError("Finner ikke linkid-kolonne")
        elif n > 1:
            raise ValueError("Flere kolonner kan inneholde linkid-id-er")

    # bestem om minuttkolonnen er 2022 eller tidligere
    if "drivetime_fw" in veger_kopi.columns and "drivetime_bw" in veger_kopi.columns:
        drivetime_fw, drivetime_bw = "drivetime_fw", "drivetime_bw"
    elif "ft_minutes" in veger_kopi.columns and "tf_minutes" in veger_kopi.columns:
        drivetime_fw, drivetime_bw = "ft_minutes", "tf_minutes"
    else:
        raise ValueError("Finner ikke kolonner med minutter")

    # litt opprydning (til utm33-koordinater, endre navn på kolonner, beholde kun relevante kolonner, resette index)
    veger_kopi = (veger_kopi
                    .to_crs(25833)
                    .rename(columns={source: "source", target: "target", linkid: "linkid"}) 
                    [["idx", "source", "target", "linkid", drivetime_fw, drivetime_bw, "oneway", "geometry"]]
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
    
    # hent ut linjenes ytterpunkter. 
    # men først: sirkler har ingen ytterpunkter og heller ingen funksjon (ellers hadde de vært delt i flere lenker). Disse må fjernes
    # finner sirklene ved å isolere linjene som har EMTPY boundary (siden sirkler ikke har yttergrenser)
    ytterpunkter = veger_kopi.copy()
    ytterpunkter["geometry"] = ytterpunkter.geometry.boundary 
    sirkler = ytterpunkter.loc[ytterpunkter.is_empty, "idx"]
    veger_kopi = veger_kopi[~veger_kopi.idx.isin(sirkler)]
    
    #boundary gir multipunkt med to ytterpunkter for hver linje, og explode() (til singlepart) gir en dobbelt så lang tabell med enkeltpunkter
    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True)
    #lag liste med wkt-geometri (well-known text). TODO: vurder å endre til tuple med koordinater for å spare plass
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(ytterpunkter.x, ytterpunkter.y)]
    #gjør annenhvert ytterpunkt til source_wkt og target_wkt
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2]

    #hvis noen lenker mangler node-id, lag nye node-id-er
    mangler_node_id = veger_kopi[(veger_kopi["source"].isna()) | (veger_kopi["source"]=='') | (veger_kopi["target"].isna()) | (veger_kopi["target"]=='')]
    if len(mangler_node_id)>0:
        #lag node-tabell
        sources = veger_kopi[["source_wkt"]].rename(columns={"source_wkt":"wkt"})  
        targets = veger_kopi[["target_wkt"]].rename(columns={"target_wkt":"wkt"})
        noder = (pd.concat([sources, targets], axis=0, ignore_index=True)
                 .drop_duplicates(subset=["wkt"])
                 .reset_index(drop=True)
        )

        #nye node-id-er
        noder["node_id"] = noder.index
        noder["node_id"] = noder["node_id"].astype(str) #funker ikke med numeriske node-navn i igraph, irriterende nok...
        
        #koble node-id-ene med source og target i vegene
        veger_kopi = (veger_kopi
            .drop(["source", "target"], axis=1)
            .merge(noder, left_on = "source_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"source"})
            .merge(noder, left_on = "target_wkt", right_on = "wkt", how="inner")
            .rename(columns={"node_id":"target"})
        )

    #velg ut de enveiskjørte
    ft = veger_kopi[(veger_kopi.oneway=="FT") | (veger_kopi.oneway=="F")] 
    tf = veger_kopi[(veger_kopi.oneway=="TF") | (veger_kopi.oneway=="T")]
    #snu source og target for lenkene som går "feil" vei
    tf = tf.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})

    #dupliser lenkene som går begge veier og snu source og target i den ene
    begge_retninger1 = veger_kopi[veger_kopi.oneway=="B"]
    begge_retninger2 = begge_retninger1.copy()
    begge_retninger2 = begge_retninger2.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
    
    # lag minutt-kolonne
    begge_retninger1["minutter"] = begge_retninger1[drivetime_fw]
    begge_retninger2["minutter"] = begge_retninger2[drivetime_bw]
    ft["minutter"] = ft[drivetime_fw]
    tf["minutter"] = tf[drivetime_bw]

    #hva er oneway=='N'? Antar at det er toveiskjørt.
    n = veger_kopi[veger_kopi.oneway=="N"]
    if len(n)>0:
        n2 = n.copy()
        n2 = n2.rename(columns={"source": "target", "target": "source", "source_wkt": "target_wkt", "target_wkt": "source_wkt"})
        n["minutter"] = n[drivetime_fw]
        n2["minutter"] = n2[drivetime_bw]
        n = n[~((n.minutter.isna()) | (n.minutter==""))]
        n2 = n2[~((n2.minutter.isna()) | (n2.minutter==""))]
        veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf, n, n2])
    else:
        veger_edges = gdf_concat([begge_retninger1, begge_retninger2, ft, tf])

    # lag meter-kolonne
    veger_edges["meter"] = veger_edges.length

    # rydd opp (fjern eventuelle 0-verdier, velg ut kolonner, fjern duplikat-lenkene med høyest kostnad, reset index)
    veger_edges = (veger_edges
        .loc[veger_edges["minutter"] >= 0, ["idx", "source", "target", "source_wkt", "target_wkt", "linkid", "minutter", "meter", "geometry"]]
        .sort_values("minutter", ascending=True)
        .drop_duplicates(subset=["source", "target"])
        .reset_index(drop=True)
    )
    
    # TODO: fullfør denne
    if turn_restrictions is not None:
        
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
                .rename(columns={"target": "middlenode", "minutter": "minutter1", "meter": "meter1", "geometry": "geom1", "idx": "edge1fid"}) 
                .loc[:, ["source", "source_wkt", "middlenode", "minutter1", "meter1", "geom1", "edge1fid"]] )
        # klargjør tabell som skal bli andre del av dobbellenkene
        restrictions = (veger_edges
                .copy()
                .rename(columns={"source": "middlenode", "minutter": "minutter2", "meter": "meter2", "geometry": "geom2", "idx": "edge2fid"})
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

        dobbellenker["minutter"] = dobbellenker["minutter1"] + dobbellenker["minutter2"]
        dobbellenker["meter"] = dobbellenker["meter1"] + dobbellenker["meter2"]
        dobbellenker["geometry"] = pygeos.line_merge(pygeos.union(pygeos.from_shapely(dobbellenker.geom1), pygeos.from_shapely(dobbellenker.geom2)))
        dobbellenker = gpd.GeoDataFrame(dobbellenker, geometry = "geometry", crs = 25833)
        dobbellenker["turn_restriction"] = True

        if "edge1fid" in turn_restrictions.columns:
            veger_edges.loc[(veger_edges["idx"].isin(turn_restrictions["edge1fid"])), "turn_restriction"] = False
        else:
            veger_edges.loc[(veger_edges["linkid"].isin(turn_restrictions["fromlinkid"])), "turn_restriction"] = False

        veger_edges = gdf_concat([veger_edges, dobbellenker])

        veger_edges = veger_edges[["source", "target", "source_wkt", "target_wkt", "linkid", "minutter", "meter", "turn_restriction", "geometry"]]
    else:
        veger_edges["turn_restriction"] = False

    return veger_edges


def graf(nettverk = r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\python\nettverk_landet.parquet", 
         directed=True, 
         turn_restrictions=False, #midlr satt til false fordi det ikke funker
         weight: str = "minutter", 
         kommune = None, # TODO: legg til kommunekolonne i vegdataene så man kan gjøre mindre analyser
         aar = None):

    # hvis man har gitt filsti som nettverk, prøv å les filen, eventuelt prøv i dapla, linux, windows-prod
    if isinstance(nettverk, str):
        try:
            nettverk = gpd.read_parquet(nettverk)
        except Exception:
            try:
                nettverk = gpd.read_file(nettverk)
            except Exception:
                try:
                    nettverk = gpd.read_parquet(r"dapla_sti\nettverk_landet.parquet")
                except Exception:
                    try:
                        nettverk = gpd.read_parquet(r"linux_sti\nettverk_landet.parquet")
                    except Exception:
                        raise ValueError(f"Finner ikke {nettverk}")

    # fjern enten lenker med svingforbud eller lenker som hindrer svinger hvis man ikke vil ha turn_restrictions på
    if directed and turn_restrictions:
        nettverk = nettverk[nettverk["turn_restriction"] != False]
    else:
        nettverk# = nettverk[nettverk["turn_restriction"] != True] midlr. fordi turn_restrictions ikke funker

    # bestem om kostnad er minutter, meter eller evt km
    if "min" in weight:
        weight = "minutter"
    elif "dist" in weight or "meter" in weight:
        weight = "meter"
    elif "kilome" in weight or "km" in weight:
        weight = "km"    
        nettverk = nettverk.to_crs(25833)
        nettverk[weight] = nettverk.length / 1000
    else:
        raise ValueError("weight må være en string som inneholder 'min', 'dist', 'meter', 'km' eller 'kilome'")

    # lag liste med tuples med lenker og legg dem til i grafen 
    edges = [(str(fra), str(til)) for fra, til in zip(nettverk["source"], nettverk["target"])]
    G = igraph.Graph.TupleList(edges, directed=directed)
    G.es['weight'] = nettverk[weight]

    #lag noder som inneholder både sources og targets (fordi startpunkter skal kobles til nærmeste source, sluttpunkter til target)
    sources = (nettverk[["source", "source_wkt"]]
        .drop_duplicates(subset=["source"])
        .rename(columns={"source": "node_id", "source_wkt": "wkt"})  
    )
    targets = (nettverk[["target", "target_wkt"]]
        .drop_duplicates(subset=["target"])
        .rename(columns={"target": "node_id", "target_wkt": "wkt"})
    )
    sources["hva"] = "source"  
    targets["hva"] = "target"  
    noder = pd.concat([sources, targets], axis=0, ignore_index=True)
    noder["geometry"] = gpd.GeoSeries.from_wkt(noder.wkt, crs=25833)
    noder = gpd.GeoDataFrame(noder, geometry="geometry", crs=25833)

    return G, nettverk.copy(), noder


# funksjon som lager graf kun basert på linjenes endepunkter.
# meter som kostnad
# kun undirected fordi man sannsynligvis ikke vet hva som er start og hva som er slutt
def graf_fra_geometri(veger,
                      samle_innen_dist = 1,
                      utm_desimaler = None):

    # til utm33, fjerne zm-koordinater og fra multilinestring til linestring
    veger_kopi = veger.copy().to_crs(25833).loc[:, ["geometry"]]
    veger_kopi = fjern_tomme_geometrier(veger_kopi)
    veger_kopi["geometry"] = pygeos.line_merge(pygeos.force_2d(pygeos.from_shapely(veger_kopi.geometry)))    
    veger_kopi = veger_kopi.explode(ignore_index=True)

    # hent ut linjenes ytterpunkter
    # først: sirkler har ingen ytterpunkter. Disse må fjernes.
    veger_kopi["idx"] = veger_kopi.index
    ytterpunkter = veger_kopi.copy()
    ytterpunkter["geometry"] = ytterpunkter.geometry.boundary
    ytterpunkter = ytterpunkter[~ytterpunkter.is_empty]
    veger_kopi = veger_kopi[veger_kopi.index.isin(ytterpunkter["idx"])].reset_index()

    ytterpunkter = veger_kopi.geometry.boundary.explode(ignore_index=True)
    if utm_desimaler is not None:
        x, y = round(ytterpunkter.x, utm_desimaler), round(ytterpunkter.y, utm_desimaler)
    else:
        x, y = ytterpunkter.x, ytterpunkter.y
    
    wkt_geom = [f"POINT ({x} {y})" for x, y in zip(x, y)]
    veger_kopi["source_wkt"], veger_kopi["target_wkt"] = wkt_geom[0::2], wkt_geom[1::2]

    veger_kopi = omkod_wkt(veger_kopi, "target_wkt", "source_wkt", samle_innen_dist)
    veger_kopi = omkod_wkt(veger_kopi, "source_wkt", "target_wkt", samle_innen_dist)
           
    veger_kopi['source'], veger_kopi['target'] = veger_kopi.source_wkt, veger_kopi.target_wkt    

    # lag graf fra liste med tuples (edges) og bruk meter som kostnad
    edges = [(fra, til) for fra, til in zip(veger_kopi.source, veger_kopi.target)]
    G = igraph.Graph.TupleList(edges, directed=False)
    G.es['weight'] = round(veger_kopi.length, 2)

    #lag node-gdf
    node_ids = [x for x in np.unique([veger_kopi.source_wkt, veger_kopi.target_wkt])]
    noder = gpd.GeoDataFrame({"geometry": gpd.GeoSeries.from_wkt(node_ids)}, geometry="geometry", crs=25833)
    noder["node_id"] = node_ids

    return G, veger_kopi, noder


# forsøk på å "flytte" noder som kun er koblet til én lenke. 
def omkod_wkt(veger, wkt_kolonne1, wkt_kolonne2, maks_dist):
    veger["idx2"] = veger.index
    veger["n"] = veger[wkt_kolonne1].map(veger[wkt_kolonne1].value_counts()).fillna(0).astype(int)
    ensom_wkt = veger[veger["n"]==1].copy()
    ensom_wkt["geometry"] = gpd.GeoSeries.from_wkt(ensom_wkt[wkt_kolonne1], crs=25833)
    ensom_wkt = gpd.GeoDataFrame(ensom_wkt, geometry="geometry", crs=25833).drop(wkt_kolonne1, axis=1)
    noder = veger[[wkt_kolonne2]].drop_duplicates(subset=[wkt_kolonne2]).rename(columns={wkt_kolonne2:"node_id"})
    noder["geometry"] = gpd.GeoSeries.from_wkt(noder["node_id"], crs=25833)
    noder = gpd.GeoDataFrame(noder, geometry="geometry", crs=25833)
    ensom_wkt = ensom_wkt.sjoin_nearest(noder[["node_id","geometry"]], distance_col="dist")
    ensom_wkt = ensom_wkt.rename(columns={"node_id": wkt_kolonne1})
    ensom_wkt = ensom_wkt.loc[ensom_wkt.dist<maks_dist, ~ensom_wkt.columns.str.contains("index|level|dist")]
    return veger[(veger["idx2"].isin(ensom_wkt["idx2"])) | (veger["n"]>1)].reset_index(drop=True)
