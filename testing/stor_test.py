# %%
import geopandas as gpd
import pandas as pd
import numpy as np
import time, sys
from shapely.wkt import loads

sys.path.append(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\Dokumenter\GitHub")
import networkz as nz
from networkz.stottefunksjoner import *

# %%
resultater = []
for aar in [2021, 2022]:
    
    print("\n\n\n" + str(aar) + "\n")
    
    if aar==2021:
        veger = gpd.read_parquet(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\data\veger_oslo_2021.parquet")
        #veger = nz.les_fil(r"C:\data\vegnett2021.gdb\ERFKPS")
    if aar==2022:
        veger = gpd.read_parquet(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\data\veger_oslo_og_naboer2.parquet")
    
    nettverk = nz.lag_nettverk(veger, turn_restrictions = None)
    
    G = nz.graf(nettverk, weight = "minutter")

    for n in [1000]:#, 1000]:#, 10000]:
        
        punkter = gpd.read_parquet(f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/tilfeldige_adresser_{n}.parquet")
        
        hovedoya = nz.til_gdf("POINT (261319.30000000013 6647824.800000001)", 
                              crs=punkter.crs)
        hovedoya["geometry"] = hovedoya.buffer(10)
        hovedoya["hovedoya"] = 1
        
        punkter = punkter.sjoin(hovedoya, how="left")
        punkter = punkter[punkter.hovedoya != 1]

        print("\n", len(punkter), "\n")
        
        for bufferdist_prosent in [100,50,25,10,0]:
            print(bufferdist_prosent)
            
            tid = time.perf_counter()
            
            od = nz.od_cost_matrix(G, punkter, punkter, 
                             #      id_kolonne="idx",
                                forsok = 30, search_tolerance=None,
                                bufferdist_prosent = bufferdist_prosent,
                                returner_linjer=False)

            df = pd.DataFrame({
                "aar": aar,
                "n": n,
                "bufferdist_prosent": bufferdist_prosent,
                "tid": time.perf_counter()-tid,
                "forsok_maks": np.max(od.forsok),
                "mangler_prosent": len(od[od.kostnad.isna()]) / len(od)*100,
                "kostnad_mean": np.mean(od.loc[~od.kostnad.isna(), "kostnad"]),
                "kostnad_median": np.median(od.loc[~od.kostnad.isna(), "kostnad"]),
                               }, index=[0])
            resultater.append(df)
            
resultater = pd.concat(resultater, axis=0, ignore_index=True)
resultater.to_csv("C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/resultater_adresser.csv", sep=";", index=False)
resultater
# %%
aar = 2021
n = 1000
bufferdist = 20
# %%
if aar==2021:
    veger = gpd.read_parquet(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\data\veger_oslo_2021.parquet")
#    veger = nz.les_fil(r"C:\data\vegnett2021.gdb\ERFKPS")
if aar==2022:
    veger = gpd.read_parquet(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\data\veger_oslo_og_naboer.parquet")

nettverk = nz.lag_nettverk(veger, turn_restrictions = None)

G = nz.graf(nettverk, weight = "minutter")
    
punkter = gpd.read_parquet(f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/tilfeldige_adresser_{n}.parquet")

tid = time.perf_counter()

od = nz.od_cost_matrix(G, punkter, punkter, 
                       id_kolonne="idx",
                    forsok = 10, search_tolerance=None,
                    bufferdist_prosent = bufferdist_prosent,
                    returner_linjer=False)

df = pd.DataFrame({
    "aar": aar,
    "n": n,
    "bufferdist_prosent": bufferdist_prosent,
    "tid": time.perf_counter()-tid,
    "forsok_maks": np.max(od.forsok),
    "mangler_prosent": len(od[od.kostnad.isna()]) / len(od)*100 
                    }, index=[0])

df       
#%%
def problempunkter(od):
    od["fra_mangler"] = od.fra.map(od[od.kostnad.isna()].fra.value_counts())
    problempunkter = od[(od.fra_mangler>np.mean(od.fra_mangler))].drop_duplicates("fra")
    problempunkter["idx"] = problempunkter["fra"]
    
    od["til_mangler"] = od.til.map(od[od.kostnad.isna()].til.value_counts())
    problempunkter2 = od[(od.til_mangler>np.mean(od.til_mangler))].drop_duplicates("til")
    problempunkter2["idx"] = problempunkter2["til"]
    
    problempunkter = pd.concat([problempunkter, problempunkter2],axis=0, ignore_index=True)
    return problempunkter[["idx"]].drop_duplicates("idx")
    
problempunkter(od)



# %%
od["fra_mangler"] = od.fra.map(od[od.kostnad.isna()].fra.value_counts())
problempunkter = od[(od.fra_mangler>900)]
print(problempunkter)
problempunkter["geometry"] = [loads(x) for x in problempunkter.fra]
problempunkter = gpd.GeoDataFrame(problempunkter, geometry="geometry", crs=25833)
problempunkter.explore()
# %%
# %%
od["til_mangler"] = od.til.map(od[od.kostnad.isna()].til.value_counts())
problempunkter = od[(od.til_mangler>900)]
problempunkter["geometry"] = [loads(x) for x in problempunkter.til]
problempunkter = gpd.GeoDataFrame(problempunkter, geometry="geometry", crs=25833)
problempunkter.explore()
# %%
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

od["til_mangler"] = od.til.map(od[od.kostnad.isna()].til.value_counts())
problempunkter = od[(od.til_mangler>900)]
problempunkter = nz.til_gdf(problempunkter, "fra", crs=25833)#.plot()
problempunkter.explore()

# %%
import geopandas as gpd
gpd.__version__
