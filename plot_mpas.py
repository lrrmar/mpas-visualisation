from matplotlib import cm
from netCDF4 import Dataset
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import os
import sys
import time
import matplotlib
from multiprocessing import shared_memory

class PlottingSession:
    def __init__(
        self,
        latCell = False,
        lonCell = False,
        latVertex = False,
        lonVertex = False,
        verticesOnCell = False,
        mesh_file = False,
        num_cells = False,
    ):
        if __name__ == "__main__":
            self.mesh_info = MPASMeshInfo(
                latCell,
                lonCell,
                latVertex,
                lonVertex,
                verticesOnCell,
                sys.argv[1]
            )
            self._allocate_memory()
            try:
                print('Session is ready for plotting')
                while True:
                    time.sleep(5)
            except KeyboardInterrupt:
                self.shm.unlink()
                exit()
        else:
            self.mesh_info = MPASMeshInfo(
                latCell,
                lonCell,
                latVertex,
                lonVertex,
                verticesOnCell,
                mesh_file,
            )


    def plot_to_ax(self, ax, field, **kwargs):
        if self.mesh_info.waiting_for_server:
            self._read_memory(field)
        self.mpas_axes = MPASAxes(self.mesh_info, ax)
        self.mpas_plot = MPASPlot(self.mesh_info, self.mpas_axes, field, **kwargs)
        return self.mpas_plot.ax.ax

    def _allocate_memory(self):

        float_data = np.array(
            [
                self.mesh_info.paddedLatVertex,
                self.mesh_info.paddedLonVertex
            ]
        )

        int_data_a = np.array(
            [
                self.mesh_info.maxLatVerticesOnCellID,
                self.mesh_info.minLatVerticesOnCellID,
                self.mesh_info.maxLonVerticesOnCellID,
                self.mesh_info.minLonVerticesOnCellID
            ]
        )
        int_data_b = self.mesh_info.verticesOnCell

        shm_flt = shared_memory.SharedMemory(create=True, size=float_data.nbytes, name='float_data')
        self.shm_flt = shm_flt
        float_buffer = np.ndarray(float_data.shape, dtype=float_data.dtype,buffer=shm_flt.buf)
        float_buffer[:] = float_data[:]

        shm_int_a = shared_memory.SharedMemory(create=True, size=int_data_a.nbytes, name='int_data_a')
        self.shm_int_a = shm_int_a
        int_buffer_a = np.ndarray(int_data_a.shape, dtype=int_data_a.dtype,buffer=shm_int_a.buf)
        int_buffer_a[:] = int_data_a[:]

        shm_int_b = shared_memory.SharedMemory(create=True, size=int_data_b.nbytes, name='int_data_b')
        self.shm_int_b = shm_int_b
        int_buffer_b = np.ndarray(int_data_b.shape, dtype=int_data_b.dtype,buffer=shm_int_b.buf)
        int_buffer_b[:] = int_data_b[:]


    def _read_memory(self, field):
        num_cells = len(field)
        self.shm_flt = shared_memory.SharedMemory(name='float_data')
        float_data_shape = (2, num_cells,)
        float_buffer = np.ndarray(float_data_shape, dtype=np.float64, buffer=self.shm_flt.buf)[:]
        paddedLatVertex = float_buffer[0]
        paddedLonVertex = float_buffer[1]

        self.shm_int_a = shared_memory.SharedMemory(name='int_data_a')
        int_data_a_shape = (4, num_cells,)
        int_buffer_a = np.ndarray(int_data_a_shape, dtype=np.int32, buffer=self.shm_int_a.buf)[:]
        maxLatVerticesOnCell = int_buffer_a[0]
        minLatVerticesOnCell = int_buffer_a[1]
        maxLonVerticesOnCell = int_buffer_a[2]
        minLonVerticesOnCell = int_buffer_a[3]

        self.shm_int_b = shared_memory.SharedMemory(name='int_data_b')
        int_data_b_shape = (num_cells, 10, )
        int_buffer_b = np.ndarray(int_data_b_shape, dtype=np.int32, buffer=self.shm_int_b.buf)[:]
        verticesOnCell = int_buffer_b
        self.mesh_info = MPASMeshRead(
            paddedLonVertex,
            paddedLatVertex,
            verticesOnCell,
            maxLatVerticesOnCell,
            minLatVerticesOnCell,
            maxLonVerticesOnCell,
            minLonVerticesOnCell,
        )

class MPASMeshRead:

    def __init__(
        self,
        paddedLonVertex,
        paddedLatVertex,
        verticesOnCell,
        maxLatVerticesOnCellID,
        minLatVerticesOnCellID,
        maxLonVerticesOnCellID,
        minLonVerticesOnCellID,
    ):
        self.paddedLonVertex = paddedLonVertex
        self.paddedLatVertex = paddedLatVertex
        self.verticesOnCell = verticesOnCell
        self.maxLatVerticesOnCellID = maxLatVerticesOnCellID
        self.minLatVerticesOnCellID = minLatVerticesOnCellID
        self.maxLonVerticesOnCellID = maxLonVerticesOnCellID
        self.minLonVerticesOnCellID = minLonVerticesOnCellID

class MPASMeshInfo:

    def __init__(
        self,
        latCell = False,
        lonCell = False,
        latVertex = False,
        lonVertex = False,
        verticesOnCell = False,
        mesh_file = False
    ):
        if not latCell and not lonCell and not latVertex and not lonVertex and not mesh_file:
            self.waiting_for_server = True
            return
        mesh_info = {
            'latCell': latCell,
            'lonCell': lonCell,
            'latVertex': latVertex,
            'lonVertex': lonVertex,
            'verticesOnCell': verticesOnCell,
        }

        # Handle missing mesh info arays
        missing_mesh_info = [key for key in mesh_info if not mesh_info[key]]
        if len(missing_mesh_info) > 0:
            if not mesh_file:
                print(f"missing {missing_mesh_info}, please supply a mesh file or wait for connection to server")
            print(f"missing {missing_mesh_info}, loading {mesh_file}")
            try:
                self.mesh_data = Dataset(mesh_file)
            except Exception as e:
                print('Could not load mesh file')
                print(e)
                return None

            for var_name in missing_mesh_info:
                try:
                    mesh_info[var_name] = self.mesh_data.variables[var_name][:]
                except Exception as e:
                    print(f"Could not get {var_name} from {mesh_file}")
                    print(e)
                    return None

        # Add mesh info
        self.latCell = mesh_info['latCell'] * 180 / np.pi
        self.lonCell = mesh_info['lonCell'] * 180 / np.pi
        self.latVertex = mesh_info['latVertex'] * 180 / np.pi
        self.lonVertex = mesh_info['lonVertex'] * 180 / np.pi
        self.verticesOnCell = mesh_info['verticesOnCell']

        # Add max coordinate info for fitting within extent
        self.paddedLatVertex = np.append(self.latVertex, [-999999])
        self.paddedLonVertex = np.append(self.lonVertex, [-999999])
        vertices = self.verticesOnCell - 1 # fortran to python indexing
        self.maxLatVerticesOnCell = np.array([self.paddedLatVertex[vertices][self.paddedLatVertex[vertices] > -999999].max() for vertices in self.verticesOnCell])
        self.minLatVerticesOnCell = np.array([self.paddedLatVertex[vertices][self.paddedLatVertex[vertices] > -999999].min() for vertices in self.verticesOnCell])
        self.maxLonVerticesOnCell = np.array([self.paddedLonVertex[vertices][self.paddedLonVertex[vertices] > -999999].max() for vertices in self.verticesOnCell])
        self.minLonVerticesOnCell = np.array([self.paddedLonVertex[vertices][self.paddedLonVertex[vertices] > -999999].min() for vertices in self.verticesOnCell])

        # Find the ID of the maximum and minimum latitude and longitude vertex per cell
        self.maxLatVerticesOnCellID = np.array([
            vertices[np.argmax(self.paddedLatVertex[vertices][self.paddedLatVertex[vertices] > -999999]
                               )] for vertices in self.verticesOnCell])
        self.minLatVerticesOnCellID = np.array([
            vertices[np.argmin(self.paddedLatVertex[vertices][self.paddedLatVertex[vertices] > -999999]
                               )] for vertices in self.verticesOnCell])
        self.maxLonVerticesOnCellID = np.array([
            vertices[np.argmax(self.paddedLonVertex[vertices][self.paddedLonVertex[vertices] > -999999]
                               )] for vertices in self.verticesOnCell])
        self.minLonVerticesOnCellID = np.array([
            vertices[np.argmin(self.paddedLonVertex[vertices][self.paddedLonVertex[vertices] > -999999]
                               )] for vertices in self.verticesOnCell])


class MPASAxes:

    def __init__(self, mesh, ax):
        self.ax = ax
        self.mesh_info = mesh
        self.internal_cells = self._cells_in_extent()

    extent_lookup = {}

    def _cells_in_extent(self):
        new_extent = self.ax.get_extent()
        if not new_extent:
            print('axes do not have an extent')
            return None

        lookup_key = self._generate_extent_lookup_key(new_extent)
        if lookup_key not in MPASAxes.extent_lookup:
            [lon_lower, lon_upper, lat_lower, lat_upper] = new_extent
            lat_cells = np.where((lat_lower <= self.mesh_info.paddedLatVertex[self.mesh_info.maxLatVerticesOnCellID]) & (self.mesh_info.paddedLatVertex[self.mesh_info.minLatVerticesOnCellID] <= lat_upper))[0]
            lon_cells = np.where((lon_lower <= self.mesh_info.paddedLonVertex[self.mesh_info.maxLonVerticesOnCellID]) & (self.mesh_info.paddedLonVertex[self.mesh_info.minLonVerticesOnCellID] <= lon_upper))[0]
            cells_in_extent = np.array([cell for cell in lat_cells if cell in lon_cells])
            MPASAxes.extent_lookup[lookup_key] = cells_in_extent
        return MPASAxes.extent_lookup[lookup_key]

    def _generate_extent_lookup_key(self, extent):
        return str(extent)

    def ax(self):
        return self.ax

class MPASPlot:

    def __init__(self, mesh, ax, field, **kwargs):
        self.mesh = mesh
        self.ax = ax
        self.field = field
        self.kwargs = kwargs
        self.callable_kwargs = None
        self._separate_callable_cells()
        self._create_polygons()
        self._render_on_ax()

    def _separate_callable_cells(self):
        kwargs = self.kwargs.copy()
        callable_kwargs = {}
        callable_keys = []
        for key in kwargs:
            if callable(kwargs[key]):
                callable_kwargs[key] = kwargs[key]
                callable_keys.append(key)
        for key in callable_keys:
            del kwargs[key]
        self.kwargs = kwargs
        self.callable_kwargs = callable_kwargs


    def _create_polygons(self):
        # Get array of values index at internal_cell
        values = self.field[self.ax.internal_cells]

        vertices = self.mesh.verticesOnCell[self.ax.internal_cells] - 1
        # add a null value indicator (padding) at end of latVertex and lonVertex
        lats = np.array([self.mesh.paddedLatVertex[cell] for cell in vertices])
        lons = np.array([self.mesh.paddedLonVertex[cell] for cell in vertices])
        # Create polygon for each internal cell
        polygons = []
        for i in range(len(vertices)):
            init_time = time.time()
            la= lats[i][lats[i] != -999999]
            lo= lons[i][lons[i] != -999999]

            called_kwargs = {}
            for key in self.callable_kwargs:
                called_kwargs[key] = self.callable_kwargs[key](values[i])
            if abs(lo.max() - lo.min()) < 100:
                coords = np.array([lo, la]).T
                polygon = Polygon(
                    coords,
                    rasterized = True,
                    **self.kwargs,
                    **called_kwargs
                )
                polygons.append(polygon)
        self.polygons = polygons

    def _render_on_ax(self):
        if len(self.polygons) > 0:
            self.ax.ax.add_collection(PatchCollection(self.polygons, match_original=True))

def plot_mpas(ax, mesh, field, **kwargs):

    init_time = time.time()

    # Read mesh info
    latCell = np.array(mesh.variables['latCell'][:]) * 180 / np.pi
    lonCell = np.array(mesh.variables['lonCell'][:]) * 180 / np.pi
    latVertex = np.array(mesh.variables['latVertex'][:]) * 180 / np.pi
    lonVertex = np.array(mesh.variables['lonVertex'][:]) * 180 / np.pi
    verticesOnCell = np.array(mesh.variables['verticesOnCell'][:])
    print(f"Reading cell info: {time.time() - init_time}")
    init_time = time.time()

    paddedLatVertex = np.append(latVertex, [-999999])
    paddedLonVertex = np.append(lonVertex, [-999999])
    vertices = verticesOnCell - 1
    maxLatVerticesOnCell = np.array([paddedLatVertex[vertices][paddedLatVertex[vertices] > -999999].max() for vertices in verticesOnCell])
    minLatVerticesOnCell = np.array([paddedLatVertex[vertices][paddedLatVertex[vertices] > -999999].min() for vertices in verticesOnCell])
    maxLonVerticesOnCell = np.array([paddedLonVertex[vertices][paddedLonVertex[vertices] > -999999].max() for vertices in verticesOnCell])
    minLonVerticesOnCell = np.array([paddedLonVertex[vertices][paddedLonVertex[vertices] > -999999].min() for vertices in verticesOnCell])
    print(f"Finding max coords: {time.time() - init_time}")
    init_time = time.time()

    # Get cells internal to extent of axes projection
    extent = ax.get_extent()
    if extent:
        [lon_lower, lon_upper, lat_lower, lat_upper] = extent
        lat_cells = np.where((lat_lower <= maxLatVerticesOnCell) & (minLatVerticesOnCell <= lat_upper))[0]
        lon_cells = np.where((lon_lower <= maxLonVerticesOnCell) & (minLonVerticesOnCell <= lon_upper))[0]
        internal_cells = np.array([cell for cell in lat_cells if cell in lon_cells])
    else:
        internal_cells = lon_cells
    print(f"Finding internal cells: {time.time() - init_time}")
    init_time = time.time()

    # separate callable kwargs
    callable_kwargs = {}
    callable_keys = []
    for key in kwargs:
        if callable(kwargs[key]):
            callable_kwargs[key] = kwargs[key]
            callable_keys.append(key)
    for key in callable_keys:
        del kwargs[key]
    print(f"Separating callables: {time.time() - init_time}")
    init_time = time.time()

    # Get array of values index at internal_cell
    values = field[internal_cells]

    vertices = verticesOnCell[internal_cells] - 1
    # add a null value indicator (padding) at end of latVertex and lonVertex
    lats = np.array([paddedLatVertex[cell] for cell in vertices])
    lons = np.array([paddedLonVertex[cell] for cell in vertices])
    print(f"padding coords: {time.time() - init_time}")
    init_time = time.time()

    # Create polygon for each internal cell
    polygons = []
    for i in range(len(vertices)):
        init_time = time.time()
        la= lats[i][lats[i] != -999999]
        lo= lons[i][lons[i] != -999999]

        called_kwargs = {}
        for key in callable_kwargs:
            called_kwargs[key] = callable_kwargs[key](values[i])
        if abs(lo.max() - lo.min()) < 100:
            coords = np.array([lo, la]).T
            polygon = Polygon(
                coords,
                rasterized = True,
                **self.kwargs,
                **called_kwargs
            )
            polygons.append(polygon)
    print(f"Creating polygons: {time.time() - init_time}")

    init_time = time.time()
    ax.add_collection(PatchCollection(polygons), match_original=True)
    print(f"Rendering polygons: {time.time() - init_time}")

    return ax

if __name__ == '__main__':
    PlottingSession()
