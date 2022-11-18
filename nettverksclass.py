import pandas as pd
import geopandas as gpd
from networkz.nettverk import make_node_ids
from networkz.stottefunksjoner import les_geoparquet


# årene det ligger tilrettelagte vegnettverk på Dapla
# hvis du vil bruke et annet nettverk, kan du kjøre det gjennom lag_nettverk()
NYESTE_AAR = 2022
ELDSTE_AAR = 2019


#NETTVERK = f"ssb-prod-dapla-felles-data-delt/GIS/Vegnett/{NYESTE_AAR}/vegnett_{NYESTE_AAR}.parquet"
NETTVERKSSTI = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"


class Nettverk:
           
    # TODO: forenkle dette
    def hent_nettverk(self):
        
        if int(self.aar)>NYESTE_AAR or int(self.aar)<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        
        if isinstance(self.nettverk, gpd.GeoDataFrame):
            return self.nettverk.to_crs(25833), None
        
        if self.nettverk is None:
            if self.aar is None:
                raise ValueError("Både 'aar' og 'nettverk' kan ikke være None")
            
        if isinstance(self.nettverk, str):
            if self.aar is not None:
                self._nettverk = self._nettverk.replace(str(NYESTE_AAR), str(self.aar))
            try:
                if "parquet" in self._nettverk:
                    return les_geoparquet(self._nettverk).to_crs(25833), self._nettverk
                return gpd.read_file(self._nettverk).to_crs(25833), self._nettverk
            except Exception:
                raise ValueError(f"Finner ikke {self._nettverk}")
        else:
            raise ValueError("'nettverk' må enten være filsti, None eller GeoDataFrame")


    def lag_noder(self, nettverk):
        
        sources = (nettverk
                [["source", "source_wkt"]]
                .drop_duplicates(subset=["source"])
                .rename(columns={"source": "node_id", "source_wkt": "wkt"})
                .reset_index(drop=True)
        )
        
        targets = (nettverk
                [["target", "target_wkt"]]
                .loc[~nettverk.target.isin(sources.node_id)]
                .drop_duplicates(subset=["target"])
                .rename(columns={"target": "node_id", "target_wkt": "wkt"})
                .reset_index(drop=True)
        )
        
        noder = pd.concat([sources, 
                        targets], axis=0, ignore_index=True)
        
        noder["geometry"] = gpd.GeoSeries.from_wkt(noder.wkt, crs=25833)
        noder = gpd.GeoDataFrame(noder, geometry="geometry", crs=25833)
        
        noder.node_id = noder.node_id.astype(int)
        noder = noder.sort_values("node_id").reset_index(drop=True)
        noder.node_id = noder.node_id.astype(str)
        
        return noder.drop("wkt", axis=1).drop_duplicates(subset=["node_id"])
    
    
    def velg_kommuner_eller_ikke(self):
        if self._kommuner is not None:
            self.nettverk.KOMMUNENR = self.nettverk.KOMMUNENR.map(lambda x: str(int(x)).zfill(4))
            if isinstance(self._kommuner, str) or isinstance(self._kommuner, int) or isinstance(self._kommuner, float):
                self._kommuner = str(int(self._kommuner)).zfill(4)       
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR==self._kommuner]
            else:
                self._kommuner = [str(int(k)).zfill(4) for k in self._kommuner]   
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR.isin(list(self._kommuner))]
            if len(self.nettverk)==0:
                raise ValueError("Ingen rader matcher med kommunenumrene dine.")
        
        return self.nettverk, self._kommuner
    
        
    def filtrer_nettverket(self):
        
        if self.directed and self.turn_restrictions:
            self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != False]
        else:
            self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != True]

        if self.sperring is not None:
            self.nettverk = self.fjern_sperringer()

        if len(self.nettverk)==0:
            raise ValueError("Nettverket har 0 rader")      
        
        if (len(self.noder)-1) != self.noder.node_id.astype(int).max():
            self.nettverk = make_node_ids(self.nettverk)
                        
        return self.nettverk
         
         
    def omkod_minutter(self):
        if self.kostnad=="minutter" or self.kostnad[0]=="minutter":
            try:
                if "bil" in self.kjoretoy or "car" in self.kjoretoy or "auto" in self.kjoretoy:
                    self._kjoretoy = "bil"
                    self.nettverk["minutter"] = self.nettverk.minutter_bil
                elif "sykkel" in self.kjoretoy or "bike" in self.kjoretoy or "bicyc" in self.kjoretoy:
                    self._kjoretoy = "sykkel"
                    self.nettverk["minutter"] = self.nettverk.length / (self.fart_sykkel * 16.6666666)
                elif "fot" in self.kjoretoy or "foot" in self.kjoretoy:
                    self._kjoretoy = "fot"
                    self.nettverk["minutter"] = self.nettverk.length / (self.fart_fot * 16.6666666)
                else:
                    raise ValueError("kjoretoy må være bil, sykkel eller fot")
            except AttributeError:
                raise AttributeError("Finner ikke minuttkolonne.")
            
        return self.nettverk
        
    
    def fjern_sperringer(self):
        
        if self.sperring is not None:
            
            for kat in [kat.lower() for kat in self.sperring]:
                self.nettverk = self.nettverk[~((self.nettverk["sperring"].fillna(0).astype(int) == 1) & (self.nettverk["category"].str.lower()==kat))]
                if self.fjern_isolerte:
                    self.nettverk = self.nettverk[~((self.nettverk["isolert"].astype(int) == 1) & (self.nettverk["category"].str.lower()==kat))]
            
        return self.nettverk
              

    # TODO: forenkle dette mye
    def bestem_kostnad(self, kostnad):
        
        if isinstance(kostnad, str):
            kostnad = [kostnad]
        if not isinstance(kostnad, list) and not isinstance(kostnad, tuple):
            raise ValueError("kostnad må være string, liste eller tuple")

        kostnader = []
        for kost in kostnad:
            if kost in self.nettverk:
                kostnader.append(kost)
        
        # hjelp til med å finne kostnadskolonnen     
        if len(kostnader) != len(kostnad):
            for col in self.nettverk.columns:
                if "minutter" in col:
                    kostnader.append("minutter")
                elif "minutes" in col:
                    kostnader.append("minutes")
                elif "meter" in kostnad:
                    self.nettverk["meter"] = self.nettverk.length
                    kostnader.append("meter")
                elif "dist" in kostnad:
                    self.nettverk["meter"] = self.nettverk.length
                    kostnader.append("meter")
                    
        if len(kostnader) == 0:
            raise ValueError("Finner ikke kostnadskolonne")

        kostnader = list(set(kostnader))
        
        if len(kostnader) > len(kostnad):
            raise ValueError(f"Flere enn {len(kostnad)} kolonner kan inneholde kostnaden{'e' if len(kostnad)>1 else ''} {', '.join(kostnad)}")
        
        if len(kostnader)==1:
            kostnader = kostnader[0]
        
        return kostnader  