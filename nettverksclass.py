import geopandas as gpd
from networkz.hjelpefunksjoner import les_geoparquet


# årene det ligger tilrettelagte vegnettverk på Dapla
# hvis man vil bruke et annet nettverk, kan man kjøre det gjennom lag_nettverk(). OBS: tar tid å finne isolerte nettverksøyer, så det bør man unngå å gjøre mer enn én gang.
NYESTE_AAR = 2022
ELDSTE_AAR = 2019


#NETTVERK = f"ssb-prod-dapla-felles-data-delt/GIS/Vegnett/{NYESTE_AAR}/vegnett_{NYESTE_AAR}.parquet"
NETTVERKSSTI_BIL = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"
NETTVERKSSTI_SYKKELFOT = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"


class Nettverk:
    """ Class med metoder for henting, filtrering og omkoding av vegnettverket."""  
    
    def hent_nettverk(self) -> gpd.GeoDataFrame:
        
        if int(self.aar)>NYESTE_AAR or int(self.aar)<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        
        self.nettverkssti = self.nettverkssti.replace(str(NYESTE_AAR), str(self.aar))
        nettverk = les_geoparquet(self.nettverkssti).to_crs(25833)
        
        return nettverk
           
    
    def velg_kommuner_eller_ikke(self) -> tuple[gpd.GeoDataFrame, list]:
        if self._kommuner:
            
            self.nettverk.KOMMUNENR = self.nettverk.KOMMUNENR.map(lambda x: str(int(x)).zfill(4))
            
            if isinstance(self._kommuner, (str, int, float)):
                self._kommuner = str(int(self._kommuner)).zfill(4)       
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR==self._kommuner]
            else:
                self._kommuner = [str(int(k)).zfill(4) for k in self._kommuner]   
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR.isin(list(self._kommuner))]
                
            if len(self.nettverk)==0:
                raise ValueError("Ingen rader matcher med kommunenumrene dine.")
        
        return self.nettverk, self._kommuner
    
    
    def filtrer_nettverk(self) -> gpd.GeoDataFrame:
        
        if self.kjoretoy=="bil":
            if self.directed and self.turn_restrictions:
                self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != False]
            else:
                if "turn_restriction" in self.nettverk.columns:
                    self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != True]

            if self.sperring:
                self.nettverk = self.sett_opp_sperringer()
            
        if self.fjern_isolerte:
            self.nettverk = self.nettverk[self.nettverk["isolert"].fillna(0) == 0]
                
        if len(self.nettverk)==0:
            raise ValueError("Nettverket har 0 rader")      
        
        return self.nettverk
    
    
    def sett_opp_sperringer(self) -> gpd.GeoDataFrame:
        
        if self.sperring is True or not "category" in self.nettverk.columns:
            return self.nettverk[self.nettverk["sperring"].astype(int) != 1]

        for vegkat in [vegkat.lower() for vegkat in self.sperring]:
            self.nettverk = self.nettverk[~((self.nettverk["sperring"].astype(int) == 1) & (self.nettverk["category"].str.lower() == vegkat))]
            
        return self.nettverk
        