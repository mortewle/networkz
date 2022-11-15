import geopandas as gpd
import pandas as pd
from networkz.od_cost_matrix import od_cost_matrix
from networkz.service_area import service_area
from networkz.shortest_path import shortest_path
from networkz.nettverk import make_node_ids


# årene det ligger tilrettelagte vegnettverk på Dapla
# hvis du vil bruke et annet nettverk, kan du kjøre det gjennom lag_nettverk()
NYESTE_AAR = 2022
ELDSTE_AAR = 2019


#NETTVERK = f"ssb-prod-dapla-felles-data-delt/GIS/Vegnett/{NYESTE_AAR}/vegnett_{NYESTE_AAR}.parquet"
NETTVERKSSTI = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"


# class som inneholder vegnettet, nodene og regler for hvordan nettverksanalysen skal gjennomføres
class Graf:
        
    def __init__(self,
                 
                 # disse handler om hvilket nettverk som skal leses inn fra fellesbøtta, og hvilket område som skal beholdes
                 aar=NYESTE_AAR,
                 nettverk = NETTVERKSSTI,
                 kommuner = None,
                 kjoretoy="bil", 
                 
                 # regler for selve ruteberegningene
                 kostnad="minutter", 
                 directed=True, 
                 turn_restrictions=False, #midlr false
                 
                 # regler for hvordan man vil koble punkter til noder. Om man vil ha mest mulig fullstendige eller mest mulig riktige resultater.
                 search_tolerance = 5000, 
                 dist_konstant = 25,
                 kost_til_nodene = True,
                 naboer = 100,
                 
                 # konstant fart i km/t når sykkel og fot er kjoretoy og kostnad er minutter
                 fart_sykkel = 20,
                 fart_fot = 5,
                 ):

        if aar>NYESTE_AAR or aar<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        
        self._aar = aar

        self._nettverk, self.nettverkssti = self.hent_nettverk(nettverk, aar, NYESTE_AAR)
        
        self._kommuner = self.velg_kommuner_eller_ikke(kommuner)
               
        self.naboer = naboer
                
        self._noder = self.lag_noder()
        
        self._kjoretoy = kjoretoy
        
        self.fart_sykkel = fart_sykkel
        self.fart_fot = fart_fot
        
        self._kostnad = self.bestem_kostnad(kostnad)
           
        self.directed = directed
        
        self.turn_restrictions = turn_restrictions
        
        self.search_tolerance = search_tolerance if search_tolerance is not None else 100000000
        
        self.dist_konstant = dist_konstant
        
        self.kost_til_nodene = kost_til_nodene
    

    # for å printe relevante attributter når man skriver navnet på class-objektet
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
 
        
    # sørg for at nodene oppdaterer seg når nettverket endres
    @property
    def nettverk(self):
        return self._nettverk
    
    @nettverk.setter
    def nettverk(self, endret_nettverk):  
        if len(self._nettverk)!=len(endret_nettverk):    
            self._nettverk = make_node_ids(endret_nettverk)
        self._noder = self.lag_noder()
        return self._nettverk

    
    # sørg for at kostnad er enten "minutter", "meter" eller begge
    @property
    def kostnad(self):
        return self._kostnad

    @kostnad.setter       
    def kostnad(self, kostnad):
        self._kostnad = self.bestem_kostnad(kostnad)
        return self._kostnad
    
    
    # disse skal være immutable for å gjøre ting enklere
    @property
    def aar(self):
        return self._aar
    @property
    def noder(self):
        return self._noder    
    @property
    def kjoretoy(self):
        return self._kjoretoy
    @property
    def kommuner(self):
        return self._kommuner        
     
    
    def od_cost_matrix(self, 
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        returner_linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                        radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                        cutoff: int = None,
                        destination_count: int = None,
                        ):
        
        return od_cost_matrix(self, startpunkter, sluttpunkter, id_kolonne, returner_linjer, radvis, cutoff, destination_count)


    def service_area(self,
                     startpunkter: gpd.GeoDataFrame, 
                     kostnad,
                     id_kolonne = None
                     ):
        
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
        
        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad i shortest_path")
        
        return shortest_path(self,startpunkter, sluttpunkter, id_kolonne, cutoff, destination_count)


    def velg_kommuner_eller_ikke(self, kommuner):
        if kommuner is not None:
            self._nettverk.KOMMUNENR = self._nettverk.KOMMUNENR.map(lambda x: str(int(x)).zfill(4))
            if isinstance(kommuner, str) or isinstance(kommuner, int) or isinstance(kommuner, float):
                kommuner = str(int(kommuner)).zfill(4)       
                self._nettverk = self._nettverk[self._nettverk.KOMMUNENR==kommuner]
            else:
                kommuner = [str(int(k)).zfill(4) for k in kommuner]   
                self._nettverk = self._nettverk[self._nettverk.KOMMUNENR.isin(list(kommuner))]
            if len(self._nettverk)==0:
                raise ValueError("Ingen rader matcher med kommunenumrene dine.")


    def hent_nettverk(self, nettverk, aar, NYESTE_AAR):
        
        if isinstance(nettverk, gpd.GeoDataFrame):
            return nettverk.to_crs(25833), None
        
        if nettverk is None:
            if aar is None:
                raise ValueError("Både 'aar' og 'nettverk' kan ikke være None")
            
        if isinstance(nettverk, str):
            if aar is not None:
                nettverk = nettverk.replace(str(NYESTE_AAR), str(aar))
            try:
                if "parquet" in nettverk:
                    return les_geoparquet(nettverk).to_crs(25833), nettverk
                return gpd.read_file(nettverk).to_crs(25833), nettverk
            except Exception:
                raise ValueError(f"Finner ikke {nettverk}")
        else:
            raise ValueError("'nettverk' må enten være filsti, None eller GeoDataFrame")


    def bestem_kostnad(self, kostnad):
        
        if isinstance(kostnad, str):
            kostnad = [kostnad]
        
        if not isinstance(kostnad, list) and not isinstance(kostnad, tuple):
            raise ValueError("kostnad må være string, liste eller tuple")

        kostnader = []
        for kost in kostnad:
            if kost in self._nettverk:
                kostnader.append(kost)
        
        # hjelp til med å finne kostnadskolonnen     
        if len(kostnader) != len(kostnad):
            if "minutter" in self._nettverk:
                kostnader.append("minutter")
            elif "minutes" in self._nettverk:
                kostnader.append("minutes")
            elif "dist" in kostnad or "meter" in kostnad:
                kostnader.append("meter")
            else:
                raise ValueError("Finner ikke kostnadskolonne")

        if len(kostnader) > len(kostnad):
            raise ValueError(f"Flere enn {len(kostnad)} kolonner kan inneholde kostnaden{'e' if len(kostnad)>1 else ''} {', '.join(kostnad)}")
        
        if len(kostnader)==1:
            kostnader = kostnader[0]
            
        # omkod minutter
        if kostnad=="minutter" or kostnad[0]=="minutter":
            try:
                if "bil" in self.kjoretoy or "car" in self.kjoretoy or "auto" in self.kjoretoy:
                    self._kjoretoy = "bil"
                    self._nettverk["minutter"] = self._nettverk.minutter_bil
                elif "sykkel" in self.kjoretoy or "bike" in self.kjoretoy or "bicyc" in self.kjoretoy:
                    self._kjoretoy = "sykkel"
                    self._nettverk["minutter"] = self._nettverk.length / (self.fart_sykkel * 16.6666666)
                elif "fot" in self.kjoretoy or "foot" in self.kjoretoy:
                    self._kjoretoy = "fot"
                    self._nettverk["minutter"] = self._nettverk.length / (self.fart_fot * 16.6666666)
                else:
                    raise ValueError("kjoretoy må være bil, sykkel eller fot")
            except AttributeError:
                raise AttributeError("Finner ikke minuttkolonne.")
            
        return kostnader
 
 
    def lag_noder(self):
        
        sources = (self._nettverk
                [["source", "source_wkt"]]
                .drop_duplicates(subset=["source"])
                .rename(columns={"source": "node_id", "source_wkt": "wkt"})
                .reset_index(drop=True)
        )
        
        targets = (self._nettverk
                [["target", "target_wkt"]]
                .loc[~self._nettverk.target.isin(sources.node_id)]
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
       
        
class Nettverk:
    def __init__(self,
                aar = NYESTE_AAR,
                nettverk = NETTVERKSSTI,
                kommuner = None,
                kjoretoy = "bil"
    ):
        
        if aar>NYESTE_AAR or aar<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        
        self._aar = aar

        self._nettverk, self.nettverkssti = self.hent_nettverk(nettverk, 
                                                          aar, 
                                                          NYESTE_AAR)
        
        self._kommuner = self.velg_kommuner_eller_ikke(kommuner)       
    
        self._kjoretoy = kjoretoy

        self._noder = self.lag_noder()


class ReglerRuteBeregning:
    def __init__(self,
                 kostnad="minutter", 
                 directed=True, 
                 turn_restrictions=False, #midlr false
    ):
        
        self._kostnad = bestem_kostnad(self, kostnad)
           
        self.directed = directed
        
        self.turn_restrictions = turn_restrictions
        
from dataclasses import dataclass

@dataclass
class ReglerNodeKobling:
    search_tolerance = 5000
    dist_konstant = 25
    kost_til_nodene = True
    naboer = 100
    
    def __post_init__(self):
        self.search_tolerance = self.search_tolerance if self.search_tolerance is not None else 100000000

        
@dataclass
class ReglerSykkelFot:
    fart_sykkel = 20
    fart_fot = 5
      
  
@dataclass
class Graf:
        
    nettverk = Nettverk
    regler_rute = ReglerRuteBeregning
    regler_noder = ReglerNodeKobling
    regler_sykkel_fot = ReglerSykkelFot

    # for å printe relevante attributter når man skriver navnet på class-objektet
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
 
 
    # mer info
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
 
        
    # gjør sånn at nodene oppdaterer seg når nettverket endres
    @property
    def nettverk(self):        
        return self._nettverk
    
    @nettverk.setter
    def nettverk(self, val):  
        if len(self._nettverk)!=len(val):    
            self._nettverk = make_node_ids(val)
        self._noder = lag_noder(self._nettverk)
        return self._nettverk

    
    # sørg for at kostnad er enten "minutter", "meter" eller begge
    @property
    def kostnad(self):
        return self._kostnad

    @kostnad.setter       
    def kostnad(self, kostnad):
        self._kostnad = bestem_kostnad(self, kostnad)
        return self._kostnad
    
    
    # disse skal være immutable for å gjøre ting enklere
    @property
    def aar(self):
        return self._aar
    @property
    def noder(self):
        return self._noder    
    @property
    def kjoretoy(self):
        return self._kjoretoy
    @property
    def kommuner(self):
        return self._kommuner        
     
    
    def od_cost_matrix(self, 
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        returner_linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                        radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                        cutoff: int = None,
                        destination_count: int = None,
                        ):
        
        return od_cost_matrix(self, startpunkter, sluttpunkter, id_kolonne, returner_linjer, radvis, cutoff, destination_count)


    def service_area(self,
                     startpunkter: gpd.GeoDataFrame, 
                     kostnad,
                     id_kolonne = None
                     ):
        
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
        
        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad i shortest_path")
        
        return shortest_path(self,startpunkter, sluttpunkter, id_kolonne, cutoff, destination_count)


    
"""

    Tenker å droppe disse fordi de gjør ting unødvendig komplisert
    
    @aar.setter
    def aar(self, val):
        
        if self.nettverkssti is not None:
            self.nettverkssti = self.nettverkssti.replace(str(self._aar), str(val))
            self._nettverk, self.nettverkssti = hent_nettverk(self.nettverkssti, NETTVERKSSTI, NYESTE_AAR, val)

        self._noder = lag_noder(self._nettverk)
                
        return self._aar


    @kjoretoy.setter
    def kjoretoy(self, val):
        self._kjoretoy = val
        if self._kjoretoy=="sykkel":
            self.nettverk["minutter"] = self.nettverk.length / (self.fart_sykkel * 16.6666666)
        elif self._kjoretoy=="fot":
            self.nettverk["minutter"] = self.nettverk.length / (self.fart_fot * 16.6666666)
        elif self._kjoretoy=="bil":
            self.nettverk["minutter"] = self.nettverk.minutter_bil
        self.nettverk = self.nettverk[["source", "target", "minutter", "minutter_bil", "meter", "turn_restriction", "geometry"]]
        return self._kjoretoy
"""