"""xorca_lonlat2ij

This module finds the i and j indeces on an ORCA grid that correspond to
pairs of longitude and latitude.

Example
-------
>>> import numpy as np
>>> import xarray as xr
>>> llat_r = np.arange(-90, 91, 45)
>>> llon_r = np.arange(-180, 181, 60)
>>> llon_rr, llat_rr = np.meshgrid(llon_r, llat_r * (-1.))
>>> data = xr.Dataset({'time_counter': (['t'], [1]),
...                    'glamf': (['t', 'y', 'x'], llon_rr[None,...]),
...                    'gphif': (['t', 'y', 'x'], llat_rr[None,...])})
>>> lon_points = np.array([170, 10, -80])
>>> lat_points = np.array([85, -60, 7])
>>> points = np.vstack((lat_points, lon_points)).T
>>> points
array([[ 85, 170],
       [-60,  10],
       [  7, -80]])
>>> grid = 'F'
>>> xgcm = False
>>> xarray_out = False
>>> index = get_ij(data, points, grid=grid, xgcm=xgcm, xarray_out=xarray_out)
>>> index
array([[0, 6],
       [3, 3],
       [2, 2]])

"""

import numpy as np
import xarray as xr
from scipy import spatial


def prepare_coords(coords):
    """Convert coords to an n*2 array."""
    coords = np.asarray(coords).astype(np.float)
    if coords.ndim == 1:
        coords = np.array([coords])
    return coords


def transform_coordinates(coords):
    """Transform coordinates to cartesian

    Parameters
    ----------
    coords : array
        Numpy array of shape=(n, 2). First column is latitudes, second column
        is longitudes.
        `array([[lat_1, lon_1], [lat_2, lon_2], ..., [lat_n, lon_n]])`

    Returns
    -------
    X : array
        Numpy array of shape=(n, 2) containing the cartesian coordinates
        corresponding to `coords`.

    """
    # WGS 84 reference coordinate system parameters
    A = 6378.137  # major axis [km]
    E2 = 6.69437999014e-3  # eccentricity squared

    coords = prepare_coords(coords)

    # convert to radiants
    lat_rad = np.radians(coords[:, 0])
    lon_rad = np.radians(coords[:, 1])

    # convert to cartesian coordinates
    r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
    x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
    y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
    z = r_n * (1 - E2) * np.sin(lat_rad)

    return np.column_stack((x, y, z))


def construct_tree(data, lat_name, lon_name):
    """ Construct a KD-tree for point lookup using `scipy.spatial.cKDTree`

    Parameters
    ----------
    data : xarray.DataArray
        Containing any variable with the latitude and longitude coordinates to
        construct the KD-tree from.
    lat_name : str
        Name of the latitude variable.
    lon_name : str
        Name of the longitude variable.

    Returns
    -------
    tree : scipy.spatial.ckdtree.cKDTree
        The KD-tree constructed based on the latitude and longitude coordinates
        in `data`.

    """
    # reshape and stack coordinates
    coords = np.column_stack((data[lat_name].values.ravel(),
                              data[lon_name].values.ravel()))

    # construct KD-tree
    tree = spatial.cKDTree(transform_coordinates(coords))

    return tree


def query_tree(tree, points, data_shape):
    """ Query the KD-tree for nearest neighbors

    Parameters
    ----------
    tree : scipy.spatial.ckdtree.cKDTree
        Output of construct_tree().
    points : array
        Numpy array of shape=(n, 2). First column is latitudes, second column
        is longitudes.
        `array([[lat_1, lon_1], [lat_2, lon_2], ..., [lat_n, lon_n]])`
    data_shape : tuple
        Shape of the original dataset (needed for regridding the index).

    Returns
    -------
    X1, X2 : xarray.DataArray
        `X1` is the j indeces and `X2` the i indeces. `X1` and `X2` have the
        dimension `location` of len=n.

    """
    _, index = tree.query(transform_coordinates(points))

    # regrid to 2D grid
    index = np.unravel_index(index, data_shape)

    # return DataArray indexers
    return xr.DataArray(index[0], dims='location'), \
        xr.DataArray(index[1], dims='location')


def specify_grid(grid):
    """Specify the grid to search on (`F` for use of xorca_brokenline).

    Parameters
    ----------
    grid : str
        Name of the grid (T, U, V, F).

    Returns
    -------
    S : str
        Name of the dummy variable living on the specified grid.

    """
    # make sure the letter describing the grid is lowercase
    lgrid = str(grid).lower()

    # define name of dummy variable (e1?)
    if (lgrid == 't') or (lgrid == 'u') \
            or (lgrid == 'v') or (lgrid == 'f'):
        var = 'e1' + lgrid
    else:
        raise ValueError(str(grid) + ' is not a valid grid-name')

    return var


def define_latlon(grid, xgcm=False):
    """Defining the names of the latitude and longitude variables.

    Parameters
    ----------
    grid : str
        Name of the grid (T, U, V, F)
    xgcm : bool, optional
        If `True`, dataset has to be compatible with xgcm. Default is `False`.

    Returns
    -------
    S1, S2 : str
        `S1` is the name of the latitude variable and `S2` the name of the
        longitude variable.

    """
    # make sure the letter describing the grid is lowercase
    lgrid = str(grid).lower()

    # put together the names
    if xgcm:
        if (lgrid == 't'):
            lat_name = 'llat_cc'
            lon_name = 'llon_cc'
        elif (lgrid == 'u'):
            lat_name = 'llat_cr'
            lon_name = 'llon_cr'
        elif (lgrid == 'v'):
            lat_name = 'llat_rc'
            lon_name = 'llon_rc'
        elif (lgrid == 'f'):
            lat_name = 'llat_rr'
            lon_name = 'llon_rr'
        else:
            raise ValueError(str(grid + ' is not a valid grid-name'))
    else:
        if lgrid not in ['t', 'u', 'v', 'f']:
            raise ValueError(str(grid + ' is not a valid grid-name'))
        else:
            lat_name = 'gphi' + lgrid
            lon_name = 'glam' + lgrid

    return lat_name, lon_name


def get_ij(data, points, grid='F', xgcm=False, xarray_out=False):
    """Get i and j indeces of ORCA grid of specified lon and lat points.

    Parameters
    ----------
    data : xarray.DataSet
        Dataset of any variable(s) living on the specified grid if `xgcm=True`
        (i.e. DataArray has been loaded with load_orca_dataset) or DataSet of
        original `mesh_hgr.nc/meshmask.nc` if `xgcm=False` (default).
    points : array
        Numpy array of shape=(n, 2). First column is latitudes, second column
        is longitudes.
        `array([[lat_1, lon_1], [lat_2, lon_2], ..., [lat_n, lon_n]])`
    grid : str, optional
        Name of the grid (T, U, V, F). Default is `'F'` (which is the grid
        used for `xorca_brokenline`)
    xgcm : bool, optional
        If `True`, dataset has to be compatible with xgcm. Default is `False`.
    xarray_out : bool, optional
        If `True` output will be of type `xarray.DataArray`. Default is
        `False`.

    Returns
    -------
    X : array or tuple
        `X` is a numpy array of shape=(n, 2), if `xarray_out=False`. First
        column is the j indeces corresponding to the latitudes of the first
        column of `points`, second column is the i indeces corresponding to the
        longitudes of the second column of `points`.
        `array([[j_1, i_1], [j_2, i_2], ..., [j_n, i_n]])`
        `X` is a tuple of two `xarray.DataArray`, if `xarray_out=True`. The
        first entry is the j indeces and the second entry, the i indeces. Each
        `xarray.DataArray` has the dimension `location` of len=n.
        `(<xarray.DataArray (location: n)>
        array([j_1, j_2, ..., j_n])
        Dimensions without coordinates: location,
        <xarray.DataArray (location: n)>
        array([i_1, i_2, ..., i_n])
        Dimensions without coordinates: location)`

    """

    # define the names of latitude and longitude variables
    lat_name, lon_name = define_latlon(grid, xgcm=xgcm)

    # get data from dataset (different for xgcm-compatible or not)
    if xgcm:
        var = specify_grid(grid)
        data_in = data[var]
        data_shape = data[var].shape
    else:
        data_in = data
        data_shape = data[lat_name].squeeze().shape

    # construct the KD-tree
    tree = construct_tree(data_in, lat_name, lon_name)

    # query the indeces from the KD-tree
    index = query_tree(tree, points, data_shape)

    # return the indeces as array or xarray.DataArray
    if xarray_out:
        return index
    else:
        return np.vstack((index[0].values, index[1].values)).T
