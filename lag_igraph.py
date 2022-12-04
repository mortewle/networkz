import numpy as np
import igraph
from sklearn.neighbors import NearestNeighbors


# lager igraph-graf som inkluderer lenker til/fra start-/sluttpunktene
def lag_graf(G,
             kostnad,
             startpunkter,
             sluttpunkter=None
             ):
                 
    # alle lenkene og kostnadene i nettverket
    edges = [(str(fra), str(til)) for fra, til in zip(G.nettverk["source"], G.nettverk["target"])]
    kostnader = list(G.nettverk[kostnad])

    # lenker mellom startpunktene og nærmeste noder
    edges_start, avstander_start, startpunkter = avstand_til_noder(startpunkter, 
                                                                   G, 
                                                                   hva = "start")
    
    if len(startpunkter)==0:
        raise ValueError("Ingen startpunkter innen search_tolerance")

    # omkod meter til minutter
    if G.kost_til_nodene==0:
        avstander_start = [0 for _ in avstander_start]
    elif kostnad=="minutter":
        avstander_start = [m_til_min(x, G) for x in avstander_start]
    
    edges = edges + edges_start
    kostnader = kostnader + avstander_start
    
    # samme for sluttpunktene
    if sluttpunkter is not None:
        edges_slutt, avstander_slutt, sluttpunkter = avstand_til_noder(sluttpunkter, G, hva = "slutt")

        if len(sluttpunkter)==0:
            raise ValueError("Ingen sluttpunkter innen search_tolerance")
        
        if G.kost_til_nodene==0:
            avstander_slutt = [0 for _ in avstander_slutt]
        elif kostnad=="minutter":
            avstander_slutt = [m_til_min(x, G) for x in avstander_slutt]

        edges = edges + edges_slutt
        kostnader = kostnader + avstander_slutt
    
    # lag liste med tuples med lenker og legg dem til i grafen 
    G2 = igraph.Graph.TupleList(edges, directed=G.directed)
    G2.es['weight'] = kostnader

    if sluttpunkter is None:
        return G2, startpunkter

    return G2, startpunkter, sluttpunkter


def m_til_min(x, G):
    """ 
    gjør om meter til minutter for lenkene mellom punktene og nabonodene.
    ganger luftlinjeavstanden med 1.5 siden det alltid er svinger i Norge. """
    
    return (x * 1.5) / (16.666667 * G.kost_til_nodene)


def avstand_til_noder(punkter, G, hva):
        
    """
    Her finner man avstanden til de n nærmeste nodene for hvert start-/sluttpunkt.
    Gjør om punktene og nodene til 1d numpy arrays bestående av koordinat-tuples
    sklearn kneighbors returnerer 2d numpy arrays med avstander og tilhørende indexer fra node-arrayen
    Derfor må node_id-kolonnen være identisk med index, altså gå fra 0 og oppover uten mellomrom
    """
    
    noder = G.noder
    
    punkter = punkter.reset_index(drop=True)
        
    # arrays med koordinat-tupler
    punkter_array = np.array([(x, y) for x, y in zip(punkter.geometry.x, punkter.geometry.y)])
    noder_array = np.array([(x, y) for x, y in zip(noder.geometry.x, noder.geometry.y)])
    
    # avstand fra punktene til 50 nærmeste noder
    # de langt unna vil ikke være attraktive pga lav hastighet fram til nodene 
    nbr = NearestNeighbors(n_neighbors=50, algorithm='ball_tree').fit(noder_array)
    avstander, idxs = nbr.kneighbors(punkter_array)
    
    punkter["dist_node"] = np.min(avstander, axis=1)
    
    # lag lenker fra punktene til nodene
    if hva=="start":
        edges = np.array([[(nz_idx, node_id) for node_id in idxs[i]] 
                          for i, nz_idx in zip(punkter.index, punkter.nz_idx)])
    # motsatt retning for sluttpunktene
    else:
        edges = np.array([[(node_id, nz_idx) for node_id in idxs[i]]
                          for i, nz_idx in zip(punkter.index, punkter.nz_idx)])

    # lag array med avstandene. 0 hvis mer enn search_tolerance eller dist_faktor-en
    avstander = np.array([[dist
                          if dist <= G.search_tolerance and dist <= (dist_min*(1+G.dist_faktor/100)+G.dist_faktor) else 0
                          for dist in avstander[i]] 
                          for i, dist_min in zip(punkter.index, punkter.dist_node)])
         
    # velg ut alt som er under search_tolerance og innenfor dist_faktor-en
    edges = edges[avstander != 0]
    avstander = avstander[avstander != 0]
    punkter = punkter[punkter.dist_node <= G.search_tolerance]
    
    edges = list(edges)
    avstander = list(avstander)

    if hva=="start":
        punkter = punkter.rename(columns={"dist_node": "dist_node_start"})
    else:
        punkter = punkter.rename(columns={"dist_node": "dist_node_slutt"})
        
    return edges, avstander, punkter

