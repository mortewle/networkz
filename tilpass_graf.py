import pandas as pd


# funksjon som gir start- og sluttpunktene en avstand til nærmeste node navnet på node-id-en som er nærmest
def naermeste_noder(noder, startpunkter, sluttpunkter = None, search_tolerance = None):
        
    # hvis directed graf, er nodene delt i sources og targets
    if "hva" in noder.columns:
        noder_start = noder.loc[noder["hva"]=="source", ["node_id", "geometry"]].copy()
        noder_slutt = noder.loc[noder["hva"]=="target", ["node_id", "geometry"]].copy()
    else:
        noder_start = noder[["node_id", "geometry"]]
        noder_slutt = noder[["node_id", "geometry"]]

    # spatial join nearest mellom nodene og punktene
    startpunkter_joinet = (startpunkter[["wkt", "geometry"]]
                           .sjoin_nearest(noder_start, 
                                          distance_col="dist")
                           .drop_duplicates("wkt") #hvis flere er like nærme, blir det duplikatrader
                           .reset_index(drop=True)
    )
    # fjern rader som er for langt unna nærmeste node
    startpunkter_joinet2 = startpunkter_joinet.loc[startpunkter_joinet.dist <= search_tolerance,  
                                                  ["wkt", "geometry", "dist", "node_id"]]

    # og lag en tabell som kun inneholder de som er for langt unna
    if len(startpunkter_joinet2) < len(startpunkter):
        for_langt_unna = startpunkter[~startpunkter.wkt.isin(startpunkter_joinet2.wkt)]
    else:
        for_langt_unna = None
        
    # for service_area(), avsluttes funksjonen her
    if sluttpunkter is None:
        return startpunkter_joinet2, for_langt_unna

    #gjør det samme for sluttpunktene
    sluttpunkter_joinet = (sluttpunkter
                           [["wkt", "geometry"]]
                           .sjoin_nearest(noder_slutt, distance_col="dist")
                           .drop_duplicates("wkt")
                           .reset_index(drop=True)
    )
    sluttpunkter_joinet2 = sluttpunkter_joinet.loc[sluttpunkter_joinet.dist <= search_tolerance,  
                                                  ["wkt", "geometry", "dist", "node_id"]]

    # og lag en lang dataframe med alle fra-til-kombinasjoner av punktene som er for langt unna nærmeste node
    if for_langt_unna is not None:
        for_langt_liste = [(x, y) for y in sluttpunkter["wkt"] for x in for_langt_unna["wkt"]]
        for_langt_unna = pd.DataFrame(data=for_langt_liste, columns=['fra','til'])
        for_langt_unna = for_langt_unna.merge(startpunkter_joinet[["wkt", "geometry", "dist", "node_id"]], 
                                              left_on = "fra", right_on = "wkt", how = "left")
        
    if len(sluttpunkter_joinet2) < len(sluttpunkter):
        for_langt_unna_slutt = sluttpunkter[~sluttpunkter.wkt.isin(sluttpunkter_joinet2.wkt)]
        for_langt_liste = [(x, y) for y in startpunkter["wkt"] for x in for_langt_unna_slutt["wkt"]]
        for_langt_unna_slutt = pd.DataFrame(data=for_langt_liste, columns=['fra','til'])
        for_langt_unna_slutt = for_langt_unna_slutt.merge(sluttpunkter_joinet[["wkt", "geometry", "dist", "node_id"]], 
                                              left_on = "til", right_on = "wkt", how = "left")
        for_langt_unna = pd.concat([for_langt_unna, for_langt_unna_slutt], axis=0, ignore_index=True)
    
    return startpunkter_joinet2, sluttpunkter_joinet2, for_langt_unna


# legger til lenker mellom start- og sluttpunktene og noder innen en viss avstand 
def oppdater_graf(graf, noder, 
                  startpunkter, sluttpunkter, 
                  edgelist, # liste over alle lenker som er lagt til i grafen hittil
                  search_tolerance, 
                  forsok, # hvilket forsøk vi er på nå. 
                  bufferdist_prosent, # hvor mye lenger enn nærmeste punkt man vil koble noder.
                  len_forrige_edgelist=0 # lengden på edgelisten forrige gang. Hvis ingen flere noder blir koblet til tre ganger på rad, avsluttes letingen
                  ):
    
    G = graf.copy()

    # hvis directed graf, er nodene delt i sources og targets
    if "hva" in noder.columns:
        noder_start = noder.loc[noder["hva"]=="source", ["node_id", "geometry"]].copy()
        noder_slutt = noder.loc[noder["hva"]=="target", ["node_id", "geometry"]].copy()
    else:
        noder_start = noder.copy()
        noder_slutt = noder.copy()

    # første forsøk vil det ikke være noen lenker lagt til. Lag da en tom liste.
    if edgelist is None:
        edgelist = []
    
    # legg til startpunktene som noder
    G.add_vertices([wkt for wkt in startpunkter["wkt"]])
    
    startpunkter_kopi = startpunkter[["geometry", "wkt", "dist"]].copy()
    # bufre med ønsket antall prosent mer enn nærmeste node. Bufrer mer og mer for hvert forsøk. 
    # Opphøyer det i antall forsøk heller enn å gange det for å få det eksponentielt
    # legger også til 5 meter for hvert forsøk siden korte avstander nesten ikke vil øke hvis man kun baserer seg på prosent
    
    startpunkter_kopi["bufferdist"] = startpunkter_kopi["dist"] * (1 + bufferdist_prosent/100) ** forsok + 5*forsok
    # hvis distansen er over search_tolerance: gjør den lik search_tolerance
    startpunkter_kopi["bufferdist"] = [dist if dist < search_tolerance else search_tolerance for dist in startpunkter_kopi["bufferdist"]]
    # bufre 
    startpunkter_kopi["geometry"] = startpunkter_kopi.buffer(startpunkter_kopi["bufferdist"])
    # koble med nodene igjen. Startpunkt-bufferne dupliseres da med én rad per overlappende node
    startpunkter_joinet = startpunkter_kopi.sjoin(noder_start, how="inner")
    # lag lenker mellom startpunkt-id og node_id hvis de ikke er i edgelist allerede
    edges = [(fra, til) for fra, til in zip(startpunkter_joinet["wkt"], startpunkter_joinet["node_id"]) if (til, fra) not in edgelist]
    # 0 i kostnad for lenkene fram til startnoden
    kostnader = G.es["weight"] + [0 for x in range(len(startpunkter_joinet))]

    # samme opplegg for sluttpunktene
    if sluttpunkter is not None:
        G.add_vertices([wkt for wkt in sluttpunkter["wkt"]])
        sluttpunkter_kopi = sluttpunkter[["geometry", "wkt", "dist"]].copy()
        sluttpunkter_kopi["bufferdist"] = sluttpunkter_kopi["dist"] * (1 + bufferdist_prosent/100) ** forsok + 5*forsok#*np.log(bufferdist_prosent+3)
        sluttpunkter_kopi["bufferdist"] = [dist if dist< search_tolerance else search_tolerance for dist in sluttpunkter_kopi["bufferdist"]]
        sluttpunkter_kopi["geometry"] = sluttpunkter_kopi.buffer(sluttpunkter_kopi["bufferdist"])
        sluttpunkter_joinet = sluttpunkter_kopi.sjoin(noder_slutt, how="inner")
        #motsatt rekkefølge på edgene her
        edges = edges + [(til, fra) for fra, til in zip(sluttpunkter_joinet["wkt"], sluttpunkter_joinet["node_id"]) if (til, fra) not in edgelist] 
        kostnader = kostnader + [0 for x in range(len(sluttpunkter_joinet))]

    # oppdater grafen
    G.add_edges(edges)
    G.es["weight"] = kostnader

    # hvis man ikke finner noen nye punkter tre ganger på rad, gi opp letingen
    # TODO: vurdere å droppe dette, men da kan man jo fanges i en langvarig loop...
    antall_edges_naa = len(edgelist)
    edgelist = edgelist + edges
    if len(edgelist)==antall_edges_naa==len_forrige_edgelist:
        return None, None, None
    
    # returnerer oppdatert graf og liste over edges som er lagt til nå og antallet før de nye ble lagt til gang
    return G, edgelist, antall_edges_naa
    

# funksjon som sjekker om id-kolonnene finnes, eller om geometri (wkt) skal brukes
# returnerer tuple med kolonnenavn for start- og sluttpunktene
# unødvendig komplisert kanskje, men funker
def bestem_ids(id_kolonne, startpunkter, sluttpunkter=None):
    if id_kolonne is None:
        startpunkter["geom_wkt"] = [x.wkt for x in startpunkter.geometry]
        if sluttpunkter is None:
            print("Bruker startpunktenes geometri som id")
        else:
            sluttpunkter["geom_wkt"] = [x.wkt for x in sluttpunkter.geometry]
            print("Bruker start- og sluttpunktenes geometri som id")
        return ("geom_wkt", "geom_wkt")
    elif isinstance(id_kolonne, str):
        if sluttpunkter is not None:
            if id_kolonne in startpunkter.columns and id_kolonne in sluttpunkter.columns:
                return (id_kolonne, id_kolonne)
            elif "geom" in id_kolonne:
                startpunkter["geom_wkt"] = [x.wkt for x in startpunkter.geometry]
                sluttpunkter["geom_wkt"] = [x.wkt for x in sluttpunkter.geometry]
                return ("geom_wkt", "geom_wkt")
            else:
                raise ValueError("id_kolonne finnes ikke i start- og/eller sluttpunkt-dataene")
        if id_kolonne in startpunkter.columns:
            return (id_kolonne, id_kolonne)
        elif "geom" in id_kolonne:
            startpunkter["geom_wkt"] = [x.wkt for x in startpunkter.geometry]
            return ("geom_wkt", "geom_wkt")
        else:
            raise ValueError("id_kolonne finnes ikke i startpunkt-dataene")
    elif isinstance(id_kolonne, list) or isinstance(id_kolonne, tuple) and len(id_kolonne)==2:
        if id_kolonne[0] in startpunkter.columns and id_kolonne[1] in sluttpunkter.columns:
            return (id_kolonne[0], id_kolonne[1])
        else:
            raise ValueError("id_kolonne finnes ikke i start- eller sluttpunkt-dataene")
    else:
        raise ValueError("id_kolonne er verken None, string, liste eller tuple.")
