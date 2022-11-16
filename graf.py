if __name__=="__main__":
    import sys
    sys.path.append(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\Dokumenter\GitHub")

import geopandas as gpd
import pandas as pd
from dataclasses import dataclass
from networkz.od_cost_matrix import od_cost_matrix
from networkz.service_area import service_area
from networkz.shortest_path import shortest_path
from networkz.nettverk import make_node_ids
from networkz.stottefunksjoner import les_geoparquet


# årene det ligger tilrettelagte vegnettverk på Dapla
# hvis du vil bruke et annet nettverk, kan du kjøre det gjennom lagnettverk()
NYESTE_AAR = 2022
ELDSTE_AAR = 2019


#NETTVERK = f"ssb-prod-dapla-felles-data-delt/GIS/Vegnett/{NYESTE_AAR}/vegnett_{NYESTE_AAR}.parquet"
NETTVERKSSTI = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"


@dataclass
class Nettverk:
    
    aar: int = NYESTE_AAR
    kjoretoy: str = "bil"
    nettverk: str = NETTVERKSSTI
    kommuner: list[str] = None
    
    # TODO: forenkle dette
    def hent_nettverk(self, NYESTE_AAR):
        
        if self.aar>NYESTE_AAR or self.aar<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        
        if isinstance(self.nettverk, gpd.GeoDataFrame):
            return self.nettverk.to_crs(25833), None
        
        if self.nettverk is None:
            if self.aar is None:
                raise ValueError("Både 'aar' og 'nettverk' kan ikke være None")
            
        if isinstance(self.nettverk, str):
            if self.aar is not None:
                self.nettverk = self.nettverk.replace(str(NYESTE_AAR), str(self.aar))
            try:
                if "parquet" in self.nettverk:
                    return les_geoparquet(self.nettverk).to_crs(25833), self.nettverk
                return gpd.read_file(self.nettverk).to_crs(25833), self.nettverk
            except Exception:
                raise ValueError(f"Finner ikke {self.nettverk}")
        else:
            raise ValueError("'nettverk' må enten være filsti, None eller GeoDataFrame")

    def lag_noder(self):
        
        sources = (self.nettverk
                [["source", "source_wkt"]]
                .drop_duplicates(subset=["source"])
                .rename(columns={"source": "node_id", "source_wkt": "wkt"})
                .reset_index(drop=True)
        )
        
        targets = (self.nettverk
                [["target", "target_wkt"]]
                .loc[~self.nettverk.target.isin(sources.node_id)]
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
        if self.kommuner is not None:
            self.nettverk.KOMMUNENR = self.nettverk.KOMMUNENR.map(lambda x: str(int(x)).zfill(4))
            if isinstance(self.kommuner, str) or isinstance(self.kommuner, int) or isinstance(self.kommuner, float):
                self.kommuner = str(int(self.kommuner)).zfill(4)       
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR==self.kommuner]
            else:
                self.kommuner = [str(int(k)).zfill(4) for k in self.kommuner]   
                self.nettverk = self.nettverk[self.nettverk.KOMMUNENR.isin(list(self.kommuner))]
            if len(self.nettverk)==0:
                raise ValueError("Ingen rader matcher med kommunenumrene dine.")
        
        return self.nettverk, self.kommuner
    
    
@dataclass
class ReglerNettverk:
    kostnad: str = "minutter"
    directed: bool = True 
    turn_restrictions: bool = False #midlr false
    sperring: str = "ERFKPS" # hvilke vegkategorier hvor vegbommer skal være til hinder. Så hvis sperring=='ERFK_S', er det lov å kjøre gjennom private bommer. Hvis sperring er None, er alle bommer lov.
    fjern_isolerte: bool = True
    
    
    def filtrer_nettverk(self):
                
        if self.directed and self.turn_restrictions:
            self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != False]
        else:
            self.nettverk = self.nettverk[self.nettverk["turn_restriction"] != True]
            
        if self.sperring is not None:
            for kat in [kat for kat in self.sperring]:
                self.nettverk = self.nettverk[~((self.nettverk["sperring"].fillna(0).astype(int) == 1) & (self.nettverk["category"]==kat))]
                if self.fjern_isolerte:
                    self.nettverk = self.nettverk[~((self.nettverk["isolert"].astype(int) == 1) & (self.nettverk["category"]==kat))]
                    
        return self.nettverk
    
    
    # TODO: forenkle dette mye
    def bestem_kostnad(self):
        
        if isinstance(self.kostnad, str):
            kostnad = [self.kostnad]
        elif isinstance(self.kostnad, list) and not isinstance(self.kostnad, tuple):
            kostnad = self.kostnad
        else:
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
                elif "dist" in kostnad or "meter" in kostnad:
                    self.nettverk["meter"] = self.nettverk.length
                    kostnader.append("meter")

        if len(kostnader) == 0:
            raise ValueError("Finner ikke kostnadskolonne")

        if len(kostnader) > len(kostnad):
            raise ValueError(f"Flere enn {len(kostnad)} kolonner kan inneholde kostnaden{'e' if len(kostnad)>1 else ''} {', '.join(kostnad)}")
        
        if len(kostnader)==1:
            kostnader = kostnader[0]
        
        return kostnader  


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
    
    
@dataclass
class ReglerNodeKobling:
    search_tolerance: int = 5000
    dist_konstant: int = 25
    kost_til_nodene: bool = True
    naboer: int = 100


@dataclass
class ReglerSykkelFot:
    fart_sykkel: int = 20
    fart_fot: int = 5
    
    # TODO
    def omkod_fart(self):
        return self.nettverk
 
    
@dataclass(order=True)
class Graf(Nettverk, ReglerNettverk, ReglerNodeKobling, ReglerSykkelFot):
        
        
    def __post_init__(self):

        self.nettverk, self.nettverkssti = self.hent_nettverk(NYESTE_AAR)

        if self.kjoretoy!="bil":
            self.nettverk = self.omkod_fart()
        
        self.nettverk, self.kommuner = self.velg_kommuner_eller_ikke()

        self.noder = self.lag_noder()

        self.kostnad = self.bestem_kostnad()
        
        self.nettverk = self.omkod_minutter()


    def finn_feil(self):
        
        if len(self.nettverk)==0:
            raise ValueError("Nettverket har 0 rader")        
        
        if (len(self.noder)-1) != self.noder.node_id.astype(int).max():
            self.nettverk = make_node_ids(self.nettverk)
        self.noder = self.lag_noder()
       
        
    def od_cost_matrix(self, 
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        returner_linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                        radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                        cutoff: int = None,
                        destination_count: int = None,
                        ):
        
        self.nettverk = self.filtrer_nettverk()
        
        self.finn_feil()
                 
        return od_cost_matrix(self, startpunkter, sluttpunkter, id_kolonne, returner_linjer, radvis, cutoff, destination_count)


    def service_area(self,
                     startpunkter: gpd.GeoDataFrame, 
                     kostnad,
                     id_kolonne = None
                     ):
        
        self.finn_feil()

        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad i service_area")
        
        return service_area(self, startpunkter, kostnad, id_kolonne)


    def shortest_path(self,
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        cutoff: int = None,
                        destination_count: int = None,
                        ):

        self.finn_feil()
        
        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad i shortest_path")
        
        return shortest_path(self,startpunkter, sluttpunkter, id_kolonne, cutoff, destination_count)


    # for å på gjeldende attributter når man printer class-objektet
    def __repr__(self):
        for attr, val in self.__dict__.items():
            if attr=="nettverkssti" or attr=="_noder":
                continue
            elif attr=="nettverk" and isinstance(val, gpd.GeoDataFrame):
                try:
                    print(f"{attr} = GeoDataFrame hentet fra: {self.nettverkssti},")    
                except AttributeError:
                    print(f"{attr} = GeoDataFrame,")
            elif isinstance(val, str):
                print(f"{attr.strip('_')} = '{val}',")
            elif isinstance(val, gpd.GeoDataFrame):
                print(f"{attr.strip('_')} = GeoDataFrame,")                
            else:
                print(f"{attr.strip('_')} = {val},")
        return "\n"
 
 
    # mer info om attributtene
    def info(self):
        print("aar: ")
        print("nettverk: ")
        print("kostnad: ")
        print("kjoretoy: ")
        print("directed: ")
        print("turn_restrictions: ")
        print("search_tolerance: ")
        print("dist_konstant: ")
        print("km_t_til_nodene: ")
 
 

if __name__=="__main__":
    
    G = Graf(kjoretoy="sykkel", aar = 2022)
    
    print(G)
    print("")
    print(type(G.nettverk))
    print((G.fart_sykkel))
    print("")
    
    punkter = gpd.read_parquet(f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/tilfeldige_adresser_1000.parquet")

    print(G.od_cost_matrix(punkter.sample(1), punkter))
    
    G.aar = 2021
    
    print(G.od_cost_matrix(punkter.sample(1), punkter))
    
    print("ferdig")
