from plot_mpas import PlottingSession
from matplotlib.colors import Normalize,ListedColormap
from matplotlib import cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy
data = Dataset('diag.2014-09-10_00.00.00.nc')
var = np.array(data.variables['t2m'][:][0])

session = PlottingSession()
vmin, vmax = var.min(), var.max()
norm = Normalize(vmin=vmin, vmax=vmax)
colormap = cm.viridis
def get_color(val):
    return colormap(norm(val))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.PlateCarree())
extent = [-180, 180, -90, 90] # x0, x1, y0, y1
ax.set_extent(extent, crs=cartopy.crs.PlateCarree())

ax = session.plot_to_ax(ax, var, color=get_color)
ax.coastlines()
plt.show()
