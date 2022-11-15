import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.wkt import loads
import pygeos


def les_geoparquet(fil, crs=25833):
    try:
        import dapla as dp
        df = dp.read_pandas(sti)
        df["geometry"] = gpd.GeoSeries.from_wkb(df.geometry)
        return gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
    except Exception:
        return gpd.read_parquet(fil).to_crs(25833)
    
    
# konverterer til geodataframe fra geoseries, shapely-objekt, wkt, liste med shapely-objekter eller shapely-sekvenser 
# OBS: når man har shapely-objekter eller wkt, bør man sette crs. 
def til_gdf(geom, set_crs=None, **qwargs) -> gpd.GeoDataFrame:

    if isinstance(geom, str):
        from shapely.wkt import loads
        geom = loads(geom)
        gdf = gpd.GeoDataFrame({"geometry": gpd.GeoSeries(geom)}, **qwargs)
    elif isinstance(geom, tuple):
        gdf, geom_kolonne = geom
        print(gdf, geom_kolonne)
        gdf["geometry"] = gpd.GeoSeries.from_wkt(gdf[geom_kolonne], crs=set_crs, **qwargs)
        gdf = gpd.GeoDataFrame(gdf, geometry="geometry", **qwargs)
    else:
        gdf = gpd.GeoDataFrame({"geometry": gpd.GeoSeries(geom)}, **qwargs)

    if set_crs:
        gdf = gdf.set_crs(set_crs)
    
    return gdf


# fjerner tomme geometrier og NaN-geometrier
def fjern_tomme_geometrier(gdf):
    if isinstance(gdf, gpd.GeoDataFrame):
        kopi = gdf.copy(deep=True)
        kopi = kopi[~kopi.geometry.is_empty]
        kopi = kopi.dropna(subset = ["geometry"])
    elif isinstance(gdf, gpd.GeoSeries):
        kopi = gdf.copy(deep=True)
        kopi = kopi[~kopi.is_empty]
        kopi = kopi.dropna()
    else:
        raise ValueError("Input må være GeoDataFrame eller GeoSeries")
    return kopi
    

# samler liste med geodataframes til en lang geodataframe
def gdf_concat(gdf_liste: list, crs=None, axis=0, ignore_index=True, geometry="geometry", **concat_qwargs) -> gpd.GeoDataFrame:

    if crs is not None:        
        #prøv å transformere alle gdf-ene til ønsket crs. Hvis det ikke funker, er det nok fordi man har naive crs. Gir da advarsel, men setter likevel crs.
        try:
            gdf_liste = [gdf.to_crs(crs) for gdf in gdf_liste]
        except ValueError:
            print("OBS: ikke alle gdf-ene dine har crs. Hvis du nå samler latlon og utm, må du først bestemme crs med set_crs(), så gi dem samme crs med to_crs()")

        gdf = gpd.GeoDataFrame(pd.concat(gdf_liste, axis=axis, ignore_index=ignore_index, **concat_qwargs), geometry=geometry, crs=crs)
            
    else:
        #hvis mer enn ett unikt crs, prøv å endre til samme som første i lista. Hvis valueerror har man nok naive crs. Gir da advarsel og lager gdf uten å spesifisere crs.
        if len(set([str(x.crs) for x in gdf_liste])) > 1:
            try:
                gdf_liste = [x.to_crs(gdf_liste[0].crs) for x in gdf_liste if len(x)>0]
            except ValueError:
                print("OBS: ikke alle gdf-ene dine har crs. Sorg for at alle dataene hadde samme crs i utgangspunktet")

        gdf = gpd.GeoDataFrame(pd.concat(gdf_liste, axis=axis, ignore_index=ignore_index, **concat_qwargs), geometry=geometry)
            
    return(gdf)


# lager n tilfeldige punkter innenfor et gitt område (mask)
def tilfeldige_punkter(n, mask=None):
    import random
    if mask is None:
        x = np.array([random.random()*10**7 for _ in range(n*1000)])
        y = np.array([random.random()*10**8 for _ in range(n*1000)])
        punkter = til_gdf([loads(f"POINT ({x} {y})") for x, y in zip(x, y)], crs=25833)
        return punkter
    mask_kopi = mask.copy()
    mask_kopi = mask_kopi.to_crs(25833)
    out = gpd.GeoDataFrame({"geometry":[]}, geometry="geometry", crs=25833)
    while len(out) < n:
        x = np.array([random.random()*10**7 for _ in range(n*1000)])
        x = x[(x > mask_kopi.bounds.minx.iloc[0]) & (x < mask_kopi.bounds.maxx.iloc[0])]
        
        y = np.array([random.random()*10**8 for _ in range(n*1000)])
        y = y[(y > mask_kopi.bounds.miny.iloc[0]) & (y < mask_kopi.bounds.maxy.iloc[0])]
        
        punkter = til_gdf([loads(f"POINT ({x} {y})") for x, y in zip(x, y)], crs=25833)
        overlapper = punkter.clip(mask_kopi)
        out = gdf_concat([out, overlapper])
    out = out.sample(n).reset_index(drop=True).to_crs(mask.crs)
    out["idx"] = out.index
    return out


# funksjon som sjekker om id-kolonnene finnes, eller om geometri (wkt) skal brukes
# returnerer tuple med kolonnenavn for start- og sluttpunktene
# unødvendig komplisert kanskje, men funker
def bestem_ids(id_kolonne, startpunkter, sluttpunkter=None) -> tuple:
    if id_kolonne is None:
        if sluttpunkter is None:
            print("Bruker startpunktenes geometri som id")
        else:
            print("Bruker start- og sluttpunktenes geometri som id")
        return ("geom_wkt", "geom_wkt")
    elif isinstance(id_kolonne, str):
        if sluttpunkter is not None:
            if id_kolonne in startpunkter.columns and id_kolonne in sluttpunkter.columns:
                return (id_kolonne, id_kolonne)
            elif "geom" in id_kolonne:
                return ("geom_wkt", "geom_wkt")
            else:
                raise ValueError("id_kolonne finnes ikke i start- og/eller sluttpunkt-dataene")
        if id_kolonne in startpunkter.columns:
            return (id_kolonne, id_kolonne)
        elif "geom" in id_kolonne:
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


def map_ids(df, id_kolonner, startpunkter, sluttpunkter=None):
    
    if not "geom_wkt" in id_kolonner:
        id_dict_start = {nz_idx: idd  for idd, nz_idx in zip(startpunkter[id_kolonner[0]], startpunkter["nz_idx"])}
        if sluttpunkter is None:
            df[id_kolonner[0]] = df[id_kolonner[0]].map(id_dict_start)
        else:
            df["fra"] = df["fra"].map(id_dict_start)
            id_dict_slutt = {nz_idx: idd  for idd, nz_idx in zip(sluttpunkter[id_kolonner[1]], sluttpunkter["nz_idx"])}
            df["til"] = df["til"].map(id_dict_slutt)
    
    else:
        id_dict_start = {nz_idx: idd.wkt  for idd, nz_idx in zip(startpunkter.geometry, startpunkter["nz_idx"])}
        df["fra"] = df["fra"].map(id_dict_start)
        if sluttpunkter is not None:
            id_dict_slutt = {nz_idx: idd.wkt  for idd, nz_idx in zip(sluttpunkter.geometry, sluttpunkter["nz_idx"])}
            df["til"] = df["til"].map(id_dict_slutt)
            
    return df


#kutter linjer x meter fra startpunktet.
#må gjentas hvis man har linjer som er mer enn dobbelt så lange som spesifisert lengde
def kutt_linjer(gdf, lengde):
    from shapely.geometry import Point, LineString
    from shapely.ops import unary_union
    
    #cut med shapely av sgilles
    def cut(line, distance):
        # Cuts a line in two at a distance from its starting point
        if distance <= 0.0 or distance >= line.length:
            return LineString(line)
        coords = list(line.coords)
        for i, p in enumerate(coords):
            pd = line.project(Point(p))
            if pd == distance:
                return unary_union([
                    LineString(coords[:i+1]),
                    LineString(coords[i:])])
            if pd > distance:
                cp = line.interpolate(distance)
                return unary_union([
                    LineString(coords[:i] + [(cp.x, cp.y)]),
                    LineString([(cp.x, cp.y)] + coords[i:])])
            
    #singlepart, kutt for hver rad, singlepart
    singlepart = gpd.GeoDataFrame(gdf, geometry="geometry", crs=gdf.crs).explode(ignore_index=True)
    singlepart["geometry"] = singlepart.geometry.apply(lambda x: cut(x, lengde))
    singlepart = gpd.GeoDataFrame(singlepart, geometry="geometry", crs=gdf.crs)
    singlepart = singlepart.explode(ignore_index=True)
    singlepart = gpd.GeoDataFrame(singlepart, geometry="geometry", crs=gdf.crs)
    singlepart = singlepart[singlepart.length>0.01]
    return singlepart
