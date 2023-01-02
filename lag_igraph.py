import numpy as np
import igraph
from sklearn.neighbors import NearestNeighbors


def lag_graf(G, kostnad, # kostnad er eget parameter siden od_cost_matrix looper gjennom kostnadene.
             startpunkter, sluttpunkter=None):
    """ 
    Lager igraph.Graph som inkluderer lenker til/fra start-/sluttpunktene. """
    
    # alle lenkene og kostnadene i nettverket
    lenker = [(str(fra), str(til)) for fra, til in zip(G.nettverk["source"], G.nettverk["target"])]
    kostnader = list(G.nettverk[kostnad])
        
    # lenker mellom startpunktene og nærmeste noder
    lenker_start, avstander_start, startpunkter = avstand_til_noder(startpunkter,
                                                                    G, 
                                                                    hva = "start")
        
    if len(startpunkter)==0:
        raise ValueError("Ingen startpunkter innen search_tolerance")
    
    # omkod meter til minutter
    avstander_start = beregn_kostnad(avstander_start, kostnad, G.kost_til_nodene)
    
    lenker = lenker + lenker_start
    kostnader = kostnader + avstander_start
    
    # samme for sluttpunktene
    if sluttpunkter is not None:
        lenker_slutt, avstander_slutt, sluttpunkter = avstand_til_noder(sluttpunkter, G, hva = "slutt")
    
        if len(sluttpunkter)==0:
            raise ValueError("Ingen sluttpunkter innen search_tolerance")
        
        avstander_slutt = beregn_kostnad(avstander_slutt, kostnad, G.kost_til_nodene)

        lenker = lenker + lenker_slutt
        kostnader = kostnader + avstander_slutt

    # lag liste med tuples med lenker og legg dem til i grafen
    G2 = igraph.Graph.TupleList(lenker, directed=G.directed)
    G2.es['weight'] = kostnader
        
    if sluttpunkter is None:
        return G2, startpunkter

    return G2, startpunkter, sluttpunkter


def avstand_til_noder(punkter, G, hva):
    """ 
    Her finner man avstanden til de n nærmeste nodene for hvert start-/sluttpunkt.
    Gjør om punktene og nodene til 1d numpy arrays bestående av koordinat-tuples
    sklearn kneighbors returnerer 2d numpy arrays med avstander og tilhørende indexer fra node-arrayen
    Derfor må node_id-kolonnen være identisk med index, altså gå fra 0 og oppover uten mellomrom. 
    """
    
    noder = G.noder
    
    punkter = punkter.reset_index(drop=True)
    
    # arrays med koordinat-tupler
    punkter_array = np.array([(x, y) for x, y in zip(punkter.geometry.x, punkter.geometry.y)])
    noder_array = np.array([(x, y) for x, y in zip(noder.geometry.x, noder.geometry.y)])
    
    # avstand fra punktene til 50 nærmeste noder (gjerne vil bare de nærmeste være attraktive pga lav hastighet fram til nodene) 
    nbr = NearestNeighbors(n_neighbors=50, algorithm='ball_tree').fit(noder_array)
    avstander, idxs = nbr.kneighbors(punkter_array)
    
    punkter["dist_node"] = np.min(avstander, axis=1)
    
    # lag lenker fra punktene til nodene
    if hva=="start":
        lenker = np.array([[(nz_idx, node_id) for node_id in idxs[i]] 
                          for i, nz_idx in zip(punkter.index, punkter.nz_idx)])
    # motsatt retning for sluttpunktene
    else:
        lenker = np.array([[(node_id, nz_idx) for node_id in idxs[i]]
                          for i, nz_idx in zip(punkter.index, punkter.nz_idx)])

    # lag array med avstandene. 0 hvis mer enn search_tolerance eller dist_faktor-en
    avstander = np.array([[dist
                          if dist <= G.search_tolerance and dist <= dist_faktor_avstand(dist_min, G.dist_faktor) else 0
                          for dist in avstander[i]] 
                          for i, dist_min in zip(punkter.index, punkter.dist_node)])
    
    # velg ut alt som er under search_tolerance og innenfor dist_faktor-en
    lenker = lenker[avstander != 0]
    avstander = avstander[avstander != 0]
    punkter = punkter[punkter.dist_node <= G.search_tolerance]
    
    lenker = [tuple(arr) for arr in lenker]
    avstander = [arr for arr in avstander]
    
    if hva=="start":
        punkter = punkter.rename(columns={"dist_node": "dist_node_start"})
    else:
        punkter = punkter.rename(columns={"dist_node": "dist_node_slutt"})
        
    return lenker, avstander, punkter


def beregn_kostnad(avstander, kostnad, kost_til_nodene):
    """ 
    Gjør om meter til minutter for lenkene mellom punktene og nabonodene.
    og ganger luftlinjeavstanden med 1.5 siden det alltid er svinger i Norge. 
    Gjør ellers ingenting.
    """
        
    if kost_til_nodene==0:
        return [0 for _ in avstander]

    elif kostnad=="meter":
        return [x*1.5 for x in avstander]
        
    elif kostnad=="minutter":
        return [(x * 1.5) / (16.666667 * kost_til_nodene) for x in avstander]
    
    
def dist_faktor_avstand(dist_min: int, dist_faktor: int) -> int: 
    """ 
    Finner terskelavstanden for lagingen av lenker mellom start- og sluttpunktene og nodene. Alle noder innenfor denne avstanden kobles med punktene.
    Terskelen er avstanden fra hvert punkt til nærmeste node pluss x prosent pluss x meter, hvor x==dist_faktor.
    Så hvis dist_faktor=10 og avstanden til node er 100, blir terskelen 120 meter (100*1.10 + 10). """
    
    return (dist_min * (1 + dist_faktor/100) + dist_faktor)