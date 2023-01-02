import numpy as np


def lag_midlr_id(noder, startpunkter, sluttpunkter=None):
    """ 
    Lager id-kolonne som brukes som node-id-er i igraph.Graph(). 
    Fordi start- og sluttpunktene må ha node-id-er som ikke finnes i nettverket.
    """
    
    startpunkter["nz_idx"] = range(len(startpunkter))
    startpunkter["nz_idx"] = startpunkter["nz_idx"] + np.max(noder.node_id.astype(int)) + 1
    startpunkter["nz_idx"] = startpunkter["nz_idx"].astype(str)
        
    if sluttpunkter is None:
        return startpunkter["nz_idx"]

    sluttpunkter["nz_idx"] = range(len(sluttpunkter))
    sluttpunkter["nz_idx"] = sluttpunkter["nz_idx"] + np.max(startpunkter.nz_idx.astype(int)) + 1
    sluttpunkter["nz_idx"] = sluttpunkter["nz_idx"].astype(str)
    
    return startpunkter["nz_idx"], sluttpunkter["nz_idx"]


def map_ids(df, id_kolonner, startpunkter, sluttpunkter=None):    
    """ Fra midlertidig til opprinnelig id. """
    
    if "wkt" in id_kolonner:
        id_dict_start = {nz_idx: idd.wkt  for idd, nz_idx in zip(startpunkter.geometry, startpunkter["nz_idx"])}
        df["fra"] = df["fra"].map(id_dict_start)
        if sluttpunkter is not None:
            id_dict_slutt = {nz_idx: idd.wkt  for idd, nz_idx in zip(sluttpunkter.geometry, sluttpunkter["nz_idx"])}
            df["til"] = df["til"].map(id_dict_slutt)
        return df
        
    if sluttpunkter is None:
        id_dict_start = {nz_idx: idd  for idd, nz_idx in zip(startpunkter[id_kolonner], startpunkter["nz_idx"])}
        df[id_kolonner] = df[id_kolonner].map(id_dict_start)
        return df
    
    if isinstance(id_kolonner, str):
        id_kolonner = (id_kolonner, id_kolonner)

    id_dict_start = {nz_idx: idd  for idd, nz_idx in zip(startpunkter[id_kolonner[0]], startpunkter["nz_idx"])}
    df["fra"] = df["fra"].map(id_dict_start)
    
    id_dict_slutt = {nz_idx: idd  for idd, nz_idx in zip(sluttpunkter[id_kolonner[1]], sluttpunkter["nz_idx"])}
    df["til"] = df["til"].map(id_dict_slutt)
        
    return df
 

def bestem_ids(id_kolonne, startpunkter, sluttpunkter=None):
    """
    Sjekker om id-kolonnene fins.
    Returnerer start- og sluttpunktene med kun geometri og id-kolonne.
    Returnerer også navnet på id-kolonnen(e) for å kunne mappe senere. 
    """
    
    if not id_kolonne or id_kolonne=="geometry":
        id_kolonne = "wkt"
        
    if sluttpunkter is None:
        assert isinstance(id_kolonne, str), "Kan bare ha én id-kolonne som string i service_area."
        if id_kolonne=="wkt":
            return startpunkter[["geometry"]], id_kolonne
        return startpunkter[[id_kolonne, "geometry"]], id_kolonne
    
    if isinstance(id_kolonne, str):
        if id_kolonne=="wkt":
            return startpunkter[["geometry"]], sluttpunkter[["geometry"]], id_kolonne
        return startpunkter[[id_kolonne, "geometry"]], sluttpunkter[[id_kolonne, "geometry"]], id_kolonne
    
    elif isinstance(id_kolonne, (list, tuple)) and len(id_kolonne)==2:
        if not id_kolonne[0] in startpunkter and id_kolonne[1] in sluttpunkter:
            raise ValueError("id_kolonne finnes ikke i start- eller sluttpunkt-dataene")
        return startpunkter[[id_kolonne[0], "geometry"]], sluttpunkter[[id_kolonne[0], "geometry"]], id_kolonne
        
    else:
        raise ValueError("id_kolonne må være string, liste, tuple eller None/False/0.")
