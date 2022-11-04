import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.wkt import loads
import pygeos


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

    if crs:        
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
                print("OBS: gdf-ene dine har ulik crs. Gir alle samme crs som forste gdf i lista")
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
