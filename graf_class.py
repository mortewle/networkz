if __name__=="__main__":
    import sys
    sys.path.append(r"C:\Users\ort\OneDrive - Statistisk sentralbyrå\Dokumenter\GitHub")

import geopandas as gpd
import pandas as pd
from networkz.od_cost_matrix import od_cost_matrix
from networkz.service_area import service_area
from networkz.shortest_path import shortest_path
from networkz.nettverk import make_node_ids
from networkz.nettverksclass import Nettverk, NETTVERKSSTI, NYESTE_AAR
from copy import copy, deepcopy


class Graf(Nettverk):
    """ class som inneholder vegnettet, nodene og regler for hvordan nettverksanalysen skal gjennomføres """
    
    def __init__(self,
                 
                 # disse handler om hvilket nettverk som skal leses inn, og hvilket område som skal beholdes
                 aar=NYESTE_AAR,
                 nettverk = NETTVERKSSTI,
                 kommuner = None,
                 kjoretoy="bil", 
                 
                 # regler for selve ruteberegningene
                 kostnad="minutter", 
                 directed=True, 
                 turn_restrictions=False, #midlr false
                 sperring = "ERFKPS", # hvilke vegkategorier hvor vegbommer skal være til hinder. Så hvis sperring=='ERFK_S', er det lov å kjøre gjennom private bommer. Hvis sperring er None, er alle bommer lov.
                 fjern_isolerte = True, 
                 
                 # regler for hvordan man vil koble punkter til noder. Om man vil ha mest mulig fullstendige eller mest mulig riktige resultater.
                 search_tolerance = 5000, 
                 dist_konstant = 25,
                 kost_til_nodene = True,
                 naboer = 100,
                 
                 # konstant fart i km/t når sykkel og fot er kjoretoy og kostnad er minutter
                 fart_sykkel = 20,
                 fart_fot = 5,
                 ):
        
        self._aar = aar
        self.nettverk = nettverk               
        self._kjoretoy = kjoretoy
        self._kommuner = kommuner
        self.directed = directed
        self.turn_restrictions = turn_restrictions
        self.sperring = sperring
        self.fjern_isolerte = fjern_isolerte
        self.search_tolerance = search_tolerance if search_tolerance is not None else 100000000
        self.dist_konstant = dist_konstant
        self.kost_til_nodene = kost_til_nodene
        self.naboer = naboer          
        self.fart_sykkel = fart_sykkel
        self.fart_fot = fart_fot
        
        # hent og klargjør nettverket for året og kommunene som er angitt
        self.nettverk, self.nettverkssti = self.hent_nettverk()
        self.nettverk, self._kommuner = self.velg_kommuner_eller_ikke()
        self._kostnad = self.bestem_kostnad(kostnad)
        self.nettverk = self.omkod_minutter()
        self._noder = self.lag_noder()


    def od_cost_matrix(self, 
                        startpunkter: gpd.GeoDataFrame, 
                        sluttpunkter: gpd.GeoDataFrame,
                        id_kolonne = None,
                        linjer = False, # om man vil at rette linjer mellom start- og sluttpunktene returnert
                        radvis = False, # hvis False beregnes kostnaden fra alle startpunkter til alle sluttpunkter. 
                        cutoff: int = None,
                        destination_count: int = None,
                        ):
        
        # lager kopi av det filtrerte nettverket for å ikke endre på det opprinnelige
        self.nettverk = self.filtrer_nettverket()
        self._noder = self.lag_noder()
                         
        return od_cost_matrix(self, startpunkter, sluttpunkter, id_kolonne, linjer, radvis, cutoff, destination_count)


    def service_area(self,
                     startpunkter: gpd.GeoDataFrame, 
                     kostnad,
                     id_kolonne = None
                     ):
        
        self.nettverk = self.filtrer_nettverket()
        self._noder = self.lag_noder()

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

        self.nettverk = self.filtrer_nettverket()
        self._noder = self.lag_noder()
        
        if not isinstance(self.kostnad, str):
            raise ValueError("Kan bare ha én kostnad i shortest_path")
        
        return shortest_path(self, startpunkter, sluttpunkter, id_kolonne, cutoff, destination_count)
    
    
    # for å få gjeldende attributter når man printer class-objektet
    def __repr__(self):
        for attr, val in self.__dict__.items():
            if attr=="nettverkssti" or attr=="_noder":
                continue
            elif attr=="nettverk" and self.nettverkssti:
                print(f"{attr} = GeoDataFrame hentet fra: {self.nettverkssti},")    
            elif isinstance(val, gpd.GeoDataFrame):
                print(f"{attr.strip('_')} = GeoDataFrame,")                
            elif isinstance(val, str):
                print(f"{attr.strip('_')} = '{val}',")
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
 
    
    # sørg for at kostnaden funker
    @property
    def kostnad(self):
        return self._kostnad

    @kostnad.setter       
    def kostnad(self, kostnad):
        self._kostnad = self.bestem_kostnad(kostnad)
        return self._kostnad
    
    
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
    



def main():
    """ testing """
    
    import numpy as np

    punkter = gpd.read_parquet(f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/tilfeldige_adresser_1000.parquet")

    G = Graf(aar = 2021,
                sperring=None,
                fjern_isolerte=False)
    
    punkter_sample = punkter.sample(1)
    
    resultater = []
    for sperring in [None, "P", "ERFKS", "ERFKPS"]:
        for fjern_isolerte in [True, False]:
        
            G2 = G.copy()
            
            G2.fjern_isolerte = fjern_isolerte
            G2.sperring = sperring
            
            od = G.od_cost_matrix(punkter_sample,
                                    punkter, 
                                    id_kolonne="idx")

            df = pd.DataFrame({
                "fjern_isolerte":fjern_isolerte,
                "sperring": sperring,
                "mangler_prosent": len(od[od[G.kostnad].isna()]) / len(od)*100,
                "kostnad_median": np.median(od.loc[~od[G.kostnad].isna(), G.kostnad]),
                                }, index=[0])
            resultater.append(df)

    resultater = pd.concat(resultater, axis=0, ignore_index=True)
    print(resultater)

    G = Graf()    
    G = Graf(kjoretoy="sykkel", kommuner="0301", aar = 2022, fjern_isolerte=False, sperring=None)
    G.kostnad = ["minutter", "meter"]
    print(G)
    
    print(G.od_cost_matrix(punkter.sample(1), punkter, id_kolonne="idx", linjer=True, cutoff=60, destination_count=10).sample(1))
    
    print(G.nettverk.isolert.value_counts())
    print(G.nettverk.sperring.value_counts())
    
    G2 = G.copy()

    G.fjern_isolerte = True
    G.sperring = "ERFKSP"
    
    G.nettverk = G.nettverk[G.nettverk.category != "S"]
    
    G.kostnad = "minutter"

    print(G.shortest_path(punkter.sample(1), punkter.sample(5)).head(1))
    print(G.service_area(punkter.sample(1), 10).head(1))
    
    print(G2.nettverk.isolert.value_counts())
    print(G2.nettverk.sperring.value_counts())
    
    print("ferdig")


if __name__=="__main__":
    main()
