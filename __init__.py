
__version__ = "0.2.0"

# sjekk om manglende dependencies (kopiert fra pandas' __init__)
_hard_dependencies = ("numpy", "pandas", "geopandas", "igraph", "sklearn")
_missing_dependencies = []

for _dependency in _hard_dependencies:
    try:
        __import__(_dependency)
    except ImportError as _e:
        _missing_dependencies.append(f"{_dependency}: {_e}")

if _missing_dependencies:
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(_missing_dependencies)
    )
del _hard_dependencies, _dependency, _missing_dependencies

from networkz.graf_class import Graf
from networkz.nettverksclass import Nettverk

from networkz.nettverk import lag_nettverk, network_from_geometry

from networkz.service_area import service_area

from networkz.shortest_path import shortest_path

from networkz.od_cost_matrix import od_cost_matrix

from networkz.stottefunksjoner import (
    til_gdf, 
    fjern_tomme_geometrier, 
    gdf_concat,
    tilfeldige_punkter,
    les_geoparquet,
    kutt_linjer
    )
