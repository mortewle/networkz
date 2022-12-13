import numpy as np


def lag_midlr_id(noder, startpunkter, sluttpunkter=None):
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

    # hvis id_kolonne er oppgitt, map/koble tilbake denne id-en    
    if not "geom_wkt" in id_kolonner:
        id_dict_start = {nz_idx: idd  for idd, nz_idx in zip(startpunkter[id_kolonner[0]], startpunkter["nz_idx"])}
        if sluttpunkter is None:
            df[id_kolonner[0]] = df[id_kolonner[0]].map(id_dict_start)
        else:
            df["fra"] = df["fra"].map(id_dict_start)
            id_dict_slutt = {nz_idx: idd  for idd, nz_idx in zip(sluttpunkter[id_kolonner[1]], sluttpunkter["nz_idx"])}
            df["til"] = df["til"].map(id_dict_slutt)
    
    # hvis ingen id_kolonne er oppgitt, brukes geometrien i wkt-format
    else:
        id_dict_start = {nz_idx: idd.wkt  for idd, nz_idx in zip(startpunkter.geometry, startpunkter["nz_idx"])}
        df["fra"] = df["fra"].map(id_dict_start)
        if sluttpunkter is not None:
            id_dict_slutt = {nz_idx: idd.wkt  for idd, nz_idx in zip(sluttpunkter.geometry, sluttpunkter["nz_idx"])}
            df["til"] = df["til"].map(id_dict_slutt)
            
    return df.reset_index(drop=True)


# funksjon som sjekker om id-kolonnene finnes, eller om geometri (wkt) skal brukes
# returnerer tuple med kolonnenavn for start- og sluttpunktene
def bestem_ids(id_kolonne, startpunkter, sluttpunkter=None) -> tuple:
    if id_kolonne is None:
        return ("geom_wkt", "geom_wkt")
    
    elif isinstance(id_kolonne, str):
        if id_kolonne=="geometry":
            return ("geom_wkt", "geom_wkt")
        if sluttpunkter is not None:
            if id_kolonne in startpunkter.columns and id_kolonne in sluttpunkter.columns:
                return (id_kolonne, id_kolonne)
            else:
                raise ValueError("id_kolonne finnes ikke i start- og/eller sluttpunkt-dataene")
        if id_kolonne in startpunkter.columns:
            return (id_kolonne, id_kolonne)
        else:
            raise ValueError("id_kolonne finnes ikke i startpunkt-dataene")
    
    elif isinstance(id_kolonne, list) or isinstance(id_kolonne, tuple) and len(id_kolonne)==2:
        if id_kolonne[0] in startpunkter.columns and id_kolonne[1] in sluttpunkter.columns:
            return (id_kolonne[0], id_kolonne[1])
        else:
            raise ValueError("id_kolonne finnes ikke i start- eller sluttpunkt-dataene")
    
    else:
        raise ValueError("id_kolonne er verken None, string, liste eller tuple.")