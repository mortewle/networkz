
# (kopiert fra pandas' __init__)
# Let users know if they're missing any of our hard dependencies
_hard_dependencies = ("numpy", "pandas", "geopandas", "igraph")
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


from networkz.od_cost_matrix import od_cost_matrix

from networkz.service_area import service_area

from networkz.shortest_path import shortest_path

from networkz.lag_graf import (
    lag_nettverk,
    graf,
    graf_fra_geometri,
    omkod_wkt
    )

from networkz.tilpass_graf import (
    naermeste_noder,
    oppdater_graf,
    bestem_ids
    )

from networkz.stottefunksjoner import (
    til_gdf, 
    fjern_tomme_geometrier, 
    gdf_concat,
    tilfeldige_punkter
    )
