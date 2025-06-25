from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.patches import Polygon

def plot_mpas(ax, grid, field, **kwargs):

    # Read grid info
    latCell = np.array(grid.variables['latCell'][:]) * 180 / np.pi
    lonCell = np.array(grid.variables['lonCell'][:]) * 180 / np.pi
    latVertex = np.array(grid.variables['latVertex'][:]) * 180 / np.pi
    lonVertex = np.array(grid.variables['lonVertex'][:]) * 180 / np.pi
    verticesOnCell = np.array(grid.variables['verticesOnCell'][:])

    # Get cells internal to extent of axes projection
    extent = ax.get_extent()
    if extent:
        [lon_lower, lon_upper, lat_lower, lat_upper] = extent
        lat_cells = np.where((lat_lower <= latCell) & (latCell <= lat_upper))[0]
        lon_cells = np.where((lon_lower <= lonCell) & (lonCell <= lon_upper))[0]
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

    for cell in internal_cells:
        vertices = np.trim_zeros(verticesOnCell[cell]) - 1
        lats = latVertex[vertices]
        lons = lonVertex[vertices]
        called_kwargs = {}
        for key in callable_kwargs:
            called_kwargs[key] = callable_kwargs[key](field[cell])
        if abs(lons.max() - lons.min()) < 100:
            coords = np.array([lons, lats]).T
            polygon = Polygon(
                coords,
                **kwargs,
                **called_kwargs
            )
            ax.add_patch(polygon)

    return ax
