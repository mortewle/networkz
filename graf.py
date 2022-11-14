import geopandas as gpd
from networkz.od_cost_matrix import od_cost_matrix
from networkz.service_area import service_area
from networkz.shortest_path import shortest_path
from networkz.nettverk import hent_nettverk, lag_noder, make_node_ids
from networkz.stottefunksjoner import bestem_kostnad


# årene det ligger tilrettelagte vegnettverk på Dapla
# hvis du vil bruke et annet nettverk, kan du kjøre det gjennom lag_nettverk()
NYESTE_AAR = 2022
ELDSTE_AAR = 2019


#NETTVERK = f"ssb-prod-dapla-felles-data-delt/GIS/Vegnett/{NYESTE_AAR}/vegnett_{NYESTE_AAR}.parquet"
NETTVERKSSTI = f"C:/Users/ort/OneDrive - Statistisk sentralbyrå/data/nettverk_{NYESTE_AAR}.parquet"


# class som inneholder vegnettet, nodene og regler for hvordan nettverksanalysen skal gjennomføres
# kjør nz.Graf().attributter() for å printe alle attributtene og detault-verdier
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
                 
                 # konstant fart i km/t når sykkel og fot er kjoretoy og kostnad er minutter
                 fart_sykkel = 20,
                 fart_fot = 5,
                 ):

        if aar>NYESTE_AAR or aar<ELDSTE_AAR:
            raise ValueError(f"aar må være mellom {ELDSTE_AAR} og {NYESTE_AAR}")
        self._aar = aar

        self._nettverk, self.nettverkssti = hent_nettverk(nettverk, 
                                                          aar, 
                                                          NYESTE_AAR)
        
        # filtrer på KOMMUNENR hvis kommuner ikke er None
        self._kommuner = kommuner
        if kommuner is not None:
            if isinstance(kommuner, str):
                self._nettverk = self._nettverk[self._nettverk.KOMMUNENR==kommuner]
            else:
                self._nettverk = self._nettverk[self._nettverk.KOMMUNENR.isin(list(kommuner))]
        
        self.fart_sykkel = fart_sykkel
        self.fart_fot = fart_fot
                
        self._noder = lag_noder(self._nettverk)
        
        self._kjoretoy = kjoretoy
        #if self.kjoretoy=="sykkel":
        #   self.nettverk.sykkelforbud=="Nei"
        
        self._kostnad = bestem_kostnad(self, kostnad)
           
        self.directed = directed
        
        self.turn_restrictions = turn_restrictions
        
        self.search_tolerance = search_tolerance if search_tolerance is not None else 100000000
        
        if dist_konstant<0:
            raise ValueError("'naboer' må være minst 0")
        self.dist_konstant = dist_konstant
        
        if isinstance(True, bool):
            self.kost_til_nodene = kost_til_nodene
        else:
            raise ValueError("kost_til_nodene må være True eller False")
    
        self.fart_sykkel = fart_sykkel
        self.fart_fot = fart_fot

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