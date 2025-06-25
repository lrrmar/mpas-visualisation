import cartopy
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from netCDF4 import Dataset
import numpy as np

from plot_mpas import plot_mpas

grid = Dataset('x1.10242.init.nc')
data = Dataset('diag.2014-09-10_00.00.00.nc')
var = np.array(data.variables['t2m'][:][0])


vmin, vmax = var.min(), var.max()
norm = Normalize(vmin=vmin, vmax=vmax)
colormap = cm.jet
def get_color(val):
    return colormap(norm(val))


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.PlateCarree())
extent = [90, 130, -25, 20] # x0, x1, y0, y1
ax.set_extent(extent, crs=cartopy.crs.PlateCarree())

ax = plot_mpas(ax, grid, var, alpha=1, color=get_color)

fig.colorbar(cm.ScalarMappable(norm=norm, cmap=colormap), ax=ax)


ax.coastlines()
plt.show()
