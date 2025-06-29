from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.patches import Polygon
import os

def plot_mpas(ax, grid, field, **kwargs):

    # Read grid info
    latCell = np.array(grid.variables['latCell'][:]) * 180 / np.pi
    lonCell = np.array(grid.variables['lonCell'][:]) * 180 / np.pi
    latVertex = np.array(grid.variables['latVertex'][:]) * 180 / np.pi
    lonVertex = np.array(grid.variables['lonVertex'][:]) * 180 / np.pi
    verticesOnCell = np.array(grid.variables['verticesOnCell'][:])

    maxLatVerticesOnCell = np.array([latVertex[np.trim_zeros(np.array(vertices)) - 1].max() for vertices in verticesOnCell])
    minLatVerticesOnCell = np.array([latVertex[np.trim_zeros(np.array(vertices)) - 1].min() for vertices in verticesOnCell])
    maxLonVerticesOnCell = np.array([lonVertex[np.trim_zeros(np.array(vertices)) - 1].max() for vertices in verticesOnCell])
    minLonVerticesOnCell = np.array([lonVertex[np.trim_zeros(np.array(vertices)) - 1].min() for vertices in verticesOnCell])


    # Get cells internal to extent of axes projection
    extent = ax.get_extent()
    if extent:
        [lon_lower, lon_upper, lat_lower, lat_upper] = extent
        lat_cells = np.where((lat_lower <= maxLatVerticesOnCell) & (minLatVerticesOnCell <= lat_upper))[0]
        lon_cells = np.where((lon_lower <= maxLonVerticesOnCell) & (minLonVerticesOnCell <= lon_upper))[0]
        internal_cells = np.array([cell for cell in lat_cells if cell in lon_cells])
    else:
        internal_cells = lon_cells

    # separate callable kwargs
    callable_kwargs = {}
    callable_keys = []
    for key in kwargs:
        if callable(kwargs[key]):
            callable_kwargs[key] = kwargs[key]
            callable_keys.append(key)
    for key in callable_keys:
        del kwargs[key]

    # Get array of values index at internal_cell
    values = field[internal_cells]

    # add a null value indicator (padding) at end of latVertex and lonVertex
    paddedLatVertex = np.append(latVertex, [-999999])
    paddedLonVertex = np.append(lonVertex, [-999999])

    vertices = verticesOnCell[internal_cells] - 1
    lats = np.array([paddedLatVertex[cell] for cell in vertices])
    lons = np.array([paddedLonVertex[cell] for cell in vertices])

    # Create polygon for each internal cell
    for i in range(len(vertices)):
        la= lats[i][lats[i] != -999999]
        lo= lons[i][lons[i] != -999999]

        called_kwargs = {}
        for key in callable_kwargs:
            called_kwargs[key] = callable_kwargs[key](values[i])
        if abs(lo.max() - lo.min()) < 100:
            coords = np.array([lo, la]).T
            polygon = Polygon(
                coords,
                **kwargs,
                **called_kwargs
            )
            ax.add_patch(polygon)

    return ax
