import geopandas as gpd
from dataclasses import dataclass
from copy import copy, deepcopy

from networkz.nettverksclass import Nettverk, NYESTE_AAR, NETTVERKSSTI_BIL, NETTVERKSSTI_SYKKELFOT
from networkz.od_cost_matrix import od_cost_matrix
from networkz.shortest_path import shortest_path
from networkz.service_area import service_area
from networkz.nettverk import lag_node_ids


#TODO: sykkelfot, 
# beregn min. manuelt med ny metode


"""
Først regler for kjøretøyene. I hver sin class for å gjøre det litt mer oversiktlig. 
Class-ene har ingen effekt utover å definere standardreglene for kjøretøyet som velges i Graf-class-en.
Reglene/parametrene kan settes i Graf(). """

@dataclass
class ReglerFot:
    directed: bool = False
    fart: int = 5
    kost_til_nodene: int = 3
    nettverkssti: str = NETTVERKSSTI_SYKKELFOT
    max_aadt: int = None
    max_fartsgrense: int = None
    
@dataclass
class ReglerSykkel:
    directed: bool = True
    fart: int = 20
    kost_til_nodene: int = 5
    nettverkssti: str = NETTVERKSSTI_SYKKELFOT
    max_aadt: int = None
    max_fartsgrense: int = None

@dataclass
class ReglerBil:
    directed: bool = True
    turn_restrictions: bool = False #svingforbud. midlr false
    sperring: str = "ERFKPS" # hvilke vegkategorier hvor vegbommer skal være til hinder. Hvis sperring er None, er alle bommer lov. Hvis sperring=='ERFKS', er det lov å kjøre gjennom private bommer.    
    kost_til_nodene: int = 10
    nettverkssti: str = NETTVERKSSTI_BIL
    

class Graf(Nettverk):
    """ 
    Class som inneholder vegnettet, nodene og generelle regler for hvordan nettverksanalysen skal gjennomføres.
    Regler knyttet til kjøretøyet går via ReglerSykkel, ReglerFot og ReglerBil, men parametrene godtas her.
    Super-class-en Nettverk inneholder metoder for å hente, filtrere og omkode nettverket. """
    
    def __init__(self,
                 
                 # disse handler om hvilket nettverk som skal leses inn, og hvilket område som skal beholdes
                 aar=NYESTE_AAR,
                 *,
                 kjoretoy="bil",
                 kommuner = None,

                 # hvis man vil bruke et annet nettverk
                 nettverk: gpd.GeoDataFrame = None, 
                 
                 # generelle regler for nettverksanalysen
                 kostnad = "minutter",
                 fjern_isolerte = True, 
                 dist_faktor = 15,
                 search_tolerance = 1000,
                 
                 # regler knyttet til kjøretøyet (parametrene i kjøretøy-class-ene)
                 **qwargs 
                 ):
        
        self._aar = aar
        self.nettverk = nettverk
        self._kommuner = kommuner
        
        self._kjoretoy = self.bestem_kjoretoy(kjoretoy)

        self._kostnad = kostnad
                    
        if self._kjoretoy=="bil":
            self.regler = ReglerBil(**qwargs)
        elif self._kjoretoy=="sykkel":
            self.regler = ReglerSykkel(**qwargs)
        elif self._kjoretoy=="fot":
            self.regler = ReglerFot(**qwargs)
        else:
            raise ValueError("kjortetoy må være bil, sykkel eller fot.")
        
        self.search_tolerance = search_tolerance if search_tolerance else 100000000
        self._fjern_isolerte = fjern_isolerte
        self.dist_faktor = dist_faktor
        self.kost_til_nodene = self.kost_til_nodene if self.kost_til_nodene else 0
        
        # hent og klargjør nettverket for året og kommunene som er angitt
        if nettverk is None:
            self.nettverk = self.hent_nettverk()
        else:
            if not "source" in nettverk.columns or not "target" in nettverk.columns:
                raise ValueError("Finner ikke kolonnene 'source' og/eller 'target'. Kjør nettverket gjennom lag_nettverk() før Graf()")
            self.nettverk = nettverk
        
        self.nettverk, self._kommuner = self.velg_kommuner_eller_ikke()
        
        self.nettverk, self._noder = lag_node_ids(self.nettverk)
        
        self._kostnad = self.sjekk_kostnad(kostnad)
        
        if not "isolert" in self.nettverk.columns:
            self._fjern_isolerte = False
            
        if self._kjoretoy != "bil":
            self.nettverk["minutter"] = self.nettverk.length / (self.fart * 16.666667)
                
        if self._kjoretoy=="bil":
            if not "sperring" in self.nettverk.columns:
                self.sperring = None

            if not "turn_restriction" in self.nettverk.columns:
                self.turn_restrictions = False
                   
    
    def od_cost_matrix(self, 
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                        radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                        cutoff: int = None,
                        destination_count: int = None,
                        ):
        
        self.nettverk = self.filtrer_nettverk()
        self.nettverk, self._noder = lag_node_ids(self.nettverk)
 
        return od_cost_matrix(self, startpunkter, sluttpunkter, id_kolonne, linjer, radvis, cutoff, destination_count)


    def service_area(self,
                     startpunkter: gpd.GeoDataFrame, 
                     impedance,
                     id_kolonne = None
                     ):
        
        self.nettverk = self.filtrer_nettverk()
        self.nettverk, self._noder = lag_node_ids(self.nettverk)

        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad (str) i shortest_path")
        
        return service_area(self, startpunkter, impedance, id_kolonne)


    def shortest_path(self,
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        cutoff: int = None,
                        destination_count: int = None,
                        ):

        self.nettverk = self.filtrer_nettverk()
        self.nettverk, self._noder = lag_node_ids(self.nettverk)
        
        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad (str) i shortest_path")
        
        return shortest_path(self, startpunkter, sluttpunkter, id_kolonne, cutoff, destination_count)
    
 
    def info(self) -> None:
        print("aar: ")
        print("nettverk: ")
        print("kostnad: ")
        print("kjoretoy: ")
        print("directed: ")
        print("turn_restrictions: ")
        print("search_tolerance: ")
        print("dist_faktor: ")
        print("kost_til_nodene: km/t ")
 
    
    def bestem_kjoretoy(self, kjoretoy) -> str:
        
        kjoretoy = kjoretoy.lower()
        
        if "bil" in kjoretoy or "car" in kjoretoy or "auto" in kjoretoy:
            self._kjoretoy = "bil"
        elif "sykkel" in kjoretoy or "bike" in kjoretoy or "bicyc" in kjoretoy:
            self._kjoretoy = "sykkel"
        elif "fot" in kjoretoy or "foot" in kjoretoy:
            self._kjoretoy = "fot"
        else:
            self._kjoretoy = kjoretoy
            
        return self._kjoretoy
    
    
    def sjekk_kostnad(self, kostnad):
        """ sjekk om kostnadskolonnen finnes i dataene """
        
        if isinstance(kostnad, str):
            kostnad = [kostnad]
        if not isinstance(kostnad, (list, tuple)):
            raise ValueError("kostnad må være string, liste eller tuple")
        
        kostnader = []
        for kost in kostnad:
            for kolonne in self.nettverk.columns:
                if kost in kolonne or kolonne in kost:
                    kostnader.append(kolonne)
                elif "min" in kost and "min" in kolonne:
                    kostnader.append(kolonne)
                elif "meter" in kost or "dist" in kost:
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
    
    
    # sørg for at kostnaden er brukbar
    @property
    def kostnad(self):
        return self._kostnad

    @kostnad.setter       
    def kostnad(self, ny_kostnad):
        self._kostnad = self.sjekk_kostnad(ny_kostnad)
        return self._kostnad


    # hindre at fjern_isolerte settes til True når ikke isolert-kolonnen finnes 
    @property
    def fjern_isolerte(self):
        return self._fjern_isolerte

    @fjern_isolerte.setter       
    def fjern_isolerte(self, ny_verdi):
        if ny_verdi and not "isolert" in self.nettverk.columns:
            raise ValueError("Kan ikke sette fjern_isolerte til True når kolonnen 'isolert' ikke finnes.")
        return self._fjern_isolerte    

    
    # disse skal det ikke være lov å endre
    @property
    def aar(self):
        return self._aar
    @property
    def kjoretoy(self):
        return self._kjoretoy
    @property
    def kommuner(self):
        return self._kommuner
    @property
    def noder(self):
        return self._noder                


    def copy(self):
        return copy(self)
    
    def deepcopy(self):
        return deepcopy(self)
    
    
    # for å printe attributtene
    def __repr__(self) -> None:
        attrs = []
        for attr, val in self.__dict__.items():
            if attr in attrs:
                continue
            elif attr=="nettverkssti" or attr=="_noder":
                continue
            elif isinstance(val, (ReglerSykkel, ReglerFot, ReglerBil)):
                for attr, val in val.__dict__.items():
                    if attr=="nettverkssti":
                        pass
                    elif isinstance(val, str):
                        print(f"{attr.strip('_')} = '{val}',")
                    else:
                        print(f"{attr.strip('_')} = {val},")
                    attrs.append(attr)
                continue
            elif attr=="nettverk" and self.nettverkssti:
                print(f"{attr} = GeoDataFrame hentet fra: {self.nettverkssti},")    
            elif isinstance(val, gpd.GeoDataFrame):
                print(f"{attr.strip('_')} = GeoDataFrame,")
            elif isinstance(val, str):
                print(f"{attr.strip('_')} = '{val}',")
            else:
                print(f"{attr.strip('_')} = {val},")
            attrs.append(attr)
        del attrs
        return ""
 
 
    # for å gjøre regel-atributtene tilgjengelig direkte i Graf-objektet.
    def __getattr__(self, navn):
        return self.regler.__getattribute__(navn)        