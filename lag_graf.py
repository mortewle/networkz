import geopandas as gpd
import pandas as pd
import numpy as np
import pygeos
import igraph
from sklearn.neighbors import NearestNeighbors
from networkz.stottefunksjoner import fjern_tomme_geometrier


def lag_graf(G,
             kostnad,
             startpunkter,
             sluttpunkter=None
             ):
    
    nettverk = G.nettverk.copy()
    noder = G.noder
                                    
    # fjern enten lenker med svingforbud eller lenker som hindrer svinger hvis man ikke vil ha turn_restrictions på
    if G.directed and G.turn_restrictions:
        nettverk = nettverk[nettverk["turn_restriction"] != False]
    else:
        nettverk# = nettverk[nettverk["turn_restriction"] != True] midlr. fordi turn_restrictions ikke funker
    
    startpunkter["nz_idx"] = range(len(startpunkter))
    startpunkter["nz_idx"] = startpunkter["nz_idx"] + np.max(noder.node_id.astype(int)) + 1
    startpunkter["nz_idx"] = startpunkter["nz_idx"].astype(str)
    
    if sluttpunkter is not None:
        sluttpunkter["nz_idx"] = range(len(sluttpunkter))
        sluttpunkter["nz_idx"] = sluttpunkter["nz_idx"] + np.max(startpunkter.nz_idx.astype(int)) + 1
        sluttpunkter["nz_idx"] = sluttpunkter["nz_idx"].astype(str)

    # alle lenkene og kostnadene i nettverket 
    edges = [(str(fra), str(til)) for fra, til in zip(nettverk["source"], nettverk["target"])]
    kostnader = list(nettverk[kostnad])

    # lenker mellom startpunktene og nærmeste noder
    edges_start, avstander_start, startpunkter = avstand_til_noder(startpunkter, 
                                                                   G, 
                                                                   hva = "start")
    
    if len(startpunkter)==0:
        raise ValueError("Ingen startpunkter innen search_tolerance")

    # omkod meter til minutter
    if kostnad=="minutter" and G.kost_til_nodene:
        if sluttpunkter is not None:
            avstander_start = [m_til_treige_min(x, G.kjoretoy) for x in avstander_start]
        else:
            avstander_start = [m_til_min(x, G.kjoretoy) for x in avstander_start]
    elif G.kost_til_nodene is False:
        avstander_start = [0 for _ in avstander_start]
      
    edges = edges + edges_start
    kostnader = kostnader + avstander_start
    
    if sluttpunkter is not None:
        edges_slutt, avstander_slutt, sluttpunkter = avstand_til_noder(sluttpunkter, G, hva = "slutt")  

        if len(sluttpunkter)==0:
            raise ValueError("Ingen sluttpunkter innen search_tolerance")
        
        if kostnad=="minutter" and G.kost_til_nodene:
            avstander_slutt = [m_til_treige_min(x, G.kjoretoy) for x in avstander_slutt]
        elif G.kost_til_nodene is False:
            avstander_slutt = [0 for _ in avstander_slutt]
        
        edges = edges + edges_slutt
        kostnader = kostnader + avstander_slutt
    
    # lag liste med tuples med lenker og legg dem til i grafen 
    G2 = igraph.Graph.TupleList(edges, directed=G.directed)
    G2.es['weight'] = kostnader
    
    if sluttpunkter is None:
        return G2, startpunkter
    return G2, startpunkter, sluttpunkter


# funksjon som gjør om meter til minutter for lenkene mellom punktene og nabonodene
# ganger luftlinjeavstanden med 1.5 siden det alltid er svinger i Norge
def m_til_min(x, kjoretoy):
    if kjoretoy=="sykkel":
        return (x*1.5) / 166.6666667 # = 10 km/t
    if kjoretoy=="fot":
        return (x*1.5) / 50 # = 3 km/t
    if kjoretoy=="bil":
        return (x*1.5) / 500 # = 30 km/t
    

# funksjon som gjør om meter til minutter med lav hastighet
def m_til_treige_min(x, kjoretoy):
    if kjoretoy=="sykkel":
        return (x*1.5) / 100 # = 6 km/t
    if kjoretoy=="fot":
        return (x*1.5) / 25 # = 1.5 km/t
    if kjoretoy=="bil":
        return (x*1.5) / 150 # = 9 km/t
      
    
def avstand_til_noder(punkter, G, hva):
        
    """
    Her finner man avstanden til de n nærmeste nodene for hvert start-/sluttpunkt.
    Gjør om punktene og nodene til 1d numpy arrays bestående av koordinat-tuples
    sklearn kneighbors returnerer 2d numpy arrays med avstander og tilhørende indexer fra node-arrayen
    Derfor må node_id-kolonnen være identisk med index, altså gå fra 0 og oppover uten mellomrom
    """
    
    noder = G.noder
    
    punkter = punkter.reset_index(drop=True)
    
    # sørg for at nodene er i riktig, numerisk rekkefølge
    noder.node_id = noder.node_id.astype(int)
    noder = noder.sort_values("node_id")
    noder.node_id = noder.node_id.astype(str)
    
    if (len(noder)-1) != np.max(noder.node_id.astype(int)):
        raise ValueError("Nodenes node_id er ikke identisk med index.")
    
    # arrays med koordinat-tupler
    punkter_array = np.array([(x, y) for x, y in zip(punkter.geometry.x, punkter.geometry.y)])
    noder_array = np.array([(x, y) for x, y in zip(noder.geometry.x, noder.geometry.y)])
    
    # avstand fra punktene til 100 nærmeste noder
    # de langt unna vil ikke være attraktive pga lav hastighet fram til nodene 
    nbr = NearestNeighbors(n_neighbors=G.naboer, algorithm='ball_tree').fit(noder_array)
    avstander, idxs = nbr.kneighbors(punkter_array)
    
    punkter["dist_node"] = np.min(avstander, axis=1)
    
    # så en stygg list comprehension fram til jeg finner noe enklere
    # når punktene er sluttpunkter, må lenkene snus, altså fra node_id til sluttpunkt
    if hva=="start":
        edges = np.array([[(nz_idx, node_id)
                for node_id in idxs[i]]
                for i, nz_idx in zip(punkter.index, punkter.nz_idx)])
    else:
        edges = np.array([[(node_id, nz_idx) 
                for node_id in idxs[i]]
                for i, nz_idx in zip(punkter.index, punkter.nz_idx)])

    avstander = np.array([[dist
                          if dist <= G.search_tolerance and dist <= (dist_min*(1+G.dist_konstant/100)+G.dist_konstant) else 0
                          for dist in avstander[i]] 
                          for i, dist_min in zip(punkter.index, punkter.dist_node)])
    
    # velg ut alt som er under search_tolerance
    edges = edges[avstander!=0]
    avstander = avstander[avstander!=0]
    punkter = punkter[punkter.dist_node <= G.search_tolerance]
    
    # flat ut listene
    edges = list(edges)
    avstander = list(avstander)

    if hva=="start":
        punkter = punkter.rename(columns={"dist_node": "dist_node_start"})
    else:
        punkter = punkter.rename(columns={"dist_node": "dist_node_slutt"})
        
    return edges, avstander, punkter

