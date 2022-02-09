from itertools import chain, islice
import numpy as np
from functools import reduce
import xarray as xr
import xarray.ufuncs as uf


def heading(ll_pairs):
    """
    Compute the heading of a section leg defined by a two latlon pairs

    Parameter
    ---------
    ll_pairs : list of tuples
        List of tuples of length 2 that defines section leg.
        `[(j0,i0),(j1,i1)]`

    Returns
    --------
    heading: Float

    """

    conv = np.pi/180  # for degree to radian conversion
    xa = ll_pairs[0][1]*conv
    xb = ll_pairs[1][1]*conv
    ya = -uf.log(np.tan(np.pi/4.-conv*ll_pairs[0][0]/2))
    yb = -uf.log(np.tan(np.pi/4.-conv*ll_pairs[1][0]/2))
    xb_xa = (xb-xa) % (2*np.pi)

    if (xb_xa >= np.pi):
        xb_xa -= 2 * np.pi
    # if (xb_xa <= -np.pi):  # Not needed because xb_xa >= 0
    #    xb_xa += 2 * np.pi

    heading = uf.arctan2(xb_xa, (yb-ya))
    if (heading < 0):
        heading += 2*np.pi
    return heading


def coordinate_rotation(lat, lon, dim):
    """
    Compute the local inclination of the coordinate system, i.e. the heading
    of every u or v interface. Which interface is processed is controlled by
    paramter dim. Note that for the y case the result is offset by -pi/2.

    Parameter
    ---------
    lat : xarray.DataArray
        Array of shape (1,y,x) giving latitude values
    lon : xarray.DataArray
        Array of shape (1,y,x) giving longitude values
    dim: string
        String to specify which dimension ('x' or 'y')
        is processed, i.e. inclination of v or u interface.

    Returns
    --------
    heading: xarray.DataArray
        Array of shape (1,y,x) giving local inclination in radian.

    """

    conv = np.pi/180  # for degree to radian conversion
    xa = lon*conv
    ya = -uf.log(uf.tan(np.pi/4.-conv*lat/2))
    if dim == 'y':
        xb = lon.roll(y=-1, roll_coords=True)*conv
        yb = -uf.log(uf.tan(np.pi/4.-conv*lat.roll(y=-1, roll_coords=True)/2))
    if dim == 'x':
        xb = lon.roll(x=-1, roll_coords=True)*conv
        yb = -uf.log(uf.tan(np.pi/4.-conv*lat.roll(x=-1, roll_coords=True)/2))
    xb_xa = (xb-xa) % (2*np.pi)
    xb_xa = (xb_xa - 2*np.pi).where(xb_xa >= np.pi, xb_xa)
    xb_xa = (xb_xa + 2*np.pi).where(xb_xa <= -np.pi, xb_xa)
    heading = uf.arctan2(xb_xa, (yb-ya))
    heading = uf.arctan2(xb_xa, (yb-ya))
    if dim == 'y':
        heading[0, -1, :] = heading[0, -2, :]
    if dim == 'x':
        heading -= np.pi/2
        heading[0, :, -1] = heading[0, :, -2]
    return heading  # *(-1)


def rotate_velocities(u, v, theta):
    """
    Rotate velocities by angel theta in clockwise direction.

        u_rot = u * xr.ufuncs.cos(theta) - v * np.sin(theta)
        v_rot = u * xr.ufuncs.sin(theta) + v * np.cos(theta)

    Parameter
    ---------
    v: xarray.Dataset
        Array of shape N providing meridional velocity component
    u: xarray.Dataset
        Array of shape N providing zonal velocity component
    theta: xarray.Dataset or scalar
        Angle [radian] used to rotate velocities. At least one dimension
        as to match the dimension of u and v or theta has to be a scalar
        value

    Returns
    -------
    v_rot: xarray.Dataset
        Array of shape N providing rotated meridional velocity component
    u_rot: xarray.Dataset
        Array of shape N providing rotated zonal velocity component
    -------
    """

    u_rot = u * xr.ufuncs.cos(theta) - v * np.sin(theta)
    v_rot = u * xr.ufuncs.sin(theta) + v * np.cos(theta)
    return u_rot, v_rot


def section_indices(ji_pairs=None, mesh_hgr=None):
    """Forms a staircase of consectuive segments from
    ij_pairs[0] to ij_pairs[-1]. Number of segments will be one
    less than number of input points.

    Parameters
    ----------
    ji_pairs : list of tuples
        List of tuples of lenth n that defines section nodes.
        `[(j0,i0),(j1,i1), ..., (jn,in)]`

    Return
    ------
    ji_pairs :  Iterator
        Iterator of tuples of length n

    """
    if mesh_hgr:
        ji_pairs_for_sub_sections = map(section_leg_indices,
                                        ji_pairs[:-1], ji_pairs[1:],
                                        [mesh_hgr]*(len(ji_pairs)-1))
    else:
        ji_pairs_for_sub_sections = map(section_leg_indices,
                                        ji_pairs[:-1], ji_pairs[1:])

    ji_pairs_stitched_together = reduce(
        lambda x, y: chain(x, islice(y, 1, None)),
        ji_pairs_for_sub_sections)

    return ji_pairs_stitched_together


def section_leg_indices(ji0=None, ji1=None, mesh_hgr=None):
    """Return sequence of (j,i) forming a staircase of from ji0 to ji1.

    Parameters
    ----------
    ji0 : tuple
        Tuple (j0, i0) defining the start point of the section leg.
    ji1 : tuple
        Tuple (j1, i1) defining the end point of the section leg.

    Return
    ------
    ji :  Iterator
        Iterator of tuples

    """
    #  Using np.linspace() for ii and jj can give diagonal steps.
    # The function gives "large" staircases instead as as may small steps as
    # possible when number of i-points is equal to number of j-points,
    # In this case the starting point is shifted to ii(i+1) and
    # the missing point is appended at the end.
    j0 = ji0[0]
    j1 = ji1[0]
    i0 = ji0[1]
    i1 = ji1[1]
    n_of_ipoints = abs(i1 - i0)
    n_of_jpoints = abs(j1 - j0)
    n_of_points = n_of_jpoints + n_of_ipoints
    assert n_of_points != 0, ("Leg includes zero points. ",
                              "Check start and endpoints.")

    if n_of_ipoints != n_of_jpoints:
        n_of_points += 1
        ii_leg = np.round(np.linspace(i0, i1, n_of_points))
        jj_leg = np.append([j0], np.cumsum(np.diff(ii_leg) == 0)
                           * np.sign(j1-j0)+j0)
    else:
        dx = np.sign(i1-i0)
        i0 += dx
        ii_leg = np.round(np.linspace(i0, i1, n_of_points))
        jj_leg = np.append([j0], np.cumsum(np.diff(ii_leg) == 0)
                           * np.sign(j1-j0)+j0)
        i0 -= dx
        ii_leg = np.append(i0, ii_leg)
        jj_leg = np.append(j0, jj_leg)
    # cast to list
    jj_leg, ii_leg = list(jj_leg.astype("int")), list(ii_leg.astype("int"))
    # Heading of section legs
    if mesh_hgr:
        latlon = [(mesh_hgr.gphif.isel(y=j0, x=i0).values.item(),
                  mesh_hgr.glamf.isel(y=j0, x=i0).values.item()),
                  (mesh_hgr.gphif.isel(y=j1, x=i1).values.item(),
                  mesh_hgr.glamf.isel(y=j1, x=i1).values.item())]
        a = heading(latlon)
        a = [a]*len(jj_leg)
        return zip(jj_leg, ii_leg, a)
    else:
        return zip(jj_leg, ii_leg)


def reduce_dataset(data, vars_to_keep, coords_to_keep):
    """ Removes all variables and coordinates from input dataset that are not
    listed in paramters.

    Parameters
    ----------
    data : xarray.Dataset
    vars_to_keep : sequence of strings
    dims_to_keep : sequence of strings

    Returns
    -------
    xarray.Dataset
    """
    for v in data.data_vars:
        if v not in vars_to_keep:
            data = data.drop(v)
    for c in data.coords:
        if c not in coords_to_keep:
            data = data.drop(c)
    return data


def apply_mask(data, mask):
    """ Mask all variables in a dataset with the provided mask. Routine uses
    mask(z=0) for variables without a vertical dimension.

    Parameters
    ----------
    data : xarray.Dataset
            Dataset containing 4D ('time_counter','x','y','z') and
            3D ('time_counter','x','y','z') variables.
    mask : xarray.Dataset
            Mask should not contain any dimensions that are not present in data
            to avoid broadcasting

    Returns
    -------
    xarray.Dataset
    """
    for v in data.data_vars:
        if 'z' in data[v].dims:
            data[v] = data[v].where(mask == 1)
        else:
            data[v] = data[v].where(mask.isel(z=0) == 1)
    return data


def interpolate_grids(data2, interpolation_scheme):
    """ Interpolate all variables in data to a different grid.
    Possible interpolations are:
    T -> U; T -> V; U -> V; V -> U.

    Parameter
    ---------
    data : xarray.Dataset
            Dataset containing 4D ('time_counter','x','y','z') and
            3D ('time_counter','x','y','z') variables.
    interpolation_scheme: strings
            Valied values are: t2u, t2v, u2v, v2u

    Returns
    -------
    xarray.Dataset
    """
    data = data2.copy()
    if interpolation_scheme not in ["t2u", "t2v", "u2v", "v2u"]:
        raise ValueError('Interpolation scheme is unknown. Valid values are:\
                          t2u, t2v, u2v, v2u')
    if interpolation_scheme == "t2v":
        for v in data.data_vars:
            data[v] = (data[v] + data[v].shift(y=-1)) * 0.5
    if interpolation_scheme == "t2u":
        for v in data.data_vars:
            data[v] = (data[v] + data[v].shift(x=-1)) * 0.5
    if interpolation_scheme == "u2v":
        for v in data.data_vars:
            data[v] = (data[v] + data[v].shift(x=1) +
                       data[v].shift(y=-1) +
                       data[v].shift(x=1, y=-1)) * 0.25
    if interpolation_scheme == "v2u":
        for v in data.data_vars:
            data[v] = (data[v] + data[v].shift(x=-1) +
                       data[v].shift(y=1) +
                       data[v].shift(x=-1, y=1)) * 0.25
    return data


def shift_grids(gridU, gridV, mesh_hgr, mesh_zgr, mask, gridT=None,
                coords_to_keep=('time_counter', ),
                vars_to_keep=('vozocrtx', 'sozotaux', 'vomecrty', 'sometauy',
                              'votemper', 'vosaline', 'vosigma0', 'sossheig'),
                var_e1v='e1v', var_e2v='e2v', var_e1u='e1u', var_e2u='e2u',
                var_e1t='e1t', var_e2t='e2t', var_e3u='e3u', var_e3v='e3v',
                var_e3t='e3t'):
    """Interpolate all input grids onto U-points and V-points and delete
    variables and coordinates not needed for selection of section.

    Parameters
    ----------
    gridU : xarray.Dataset
        Dataset of shape (t,z,y,x) containing variables living on the U-point.
    gridV : xarray.Dataset
        Dataset of shape (t,z,y,x) containing variables living on the V-point.
    mesh_hgr : xarray.Dataset
        Dataset of shape (1,y,x) with horizontal grid scale factors.
        Variable names can be given below.
    mesh_zgr : xarray.Dataset
        Dataset of shape (1,z,y,x) with vertical grid scale factors.
        Variable names can be given below.
    mask : xarray.Dataset
        Dataset of shape (1,z,y,x) containing umask, vmask. tmask is only
        required when gridT passed
    gridT : xarray.Dataset, optional
        Dataset of shape (t,z,y,x) containing variables living on the T-point.
    coords_to_keep : sequence of strings, optional
        All dimension not stated here will be dropped
        Default is: ('time_counter', )
    vars_to_keep : sequence of strings, optional
        All variables not stated here will be dropped.
        Default is: ('vozocrtx', 'vomecrty', 'votemper', 'vosaline')
    var_e1v: string, optional
        Variable name of the zonal scale factor of the V-grid.
        Default is e1v.
    var_e2v: string, optional
        Variable name of the meridional scale factor of the V-grid.
        Default is e2v.
    var_e1u: string, optional
        Variable name of the zonal scale factor of the U-grid.
        Default is e1u.
    var_e2u: string, optional
        Variable name of the meridional scale factor of the U-grid.
        Default is e2u.
    var_e1t: string, optional
        Variable name of the zonal scale factor of the T-grid.
        Default is e1u.
    var_e2t: string, optional
        Variable name of the meridional scale factor of the T-grid.
        Default is e2t.
    var_e3u: string, optional
        Variable name of the vertical scale factor of the U-grid.
        Default is e3u.
    var_e3v: string, optional
        Variable name of the vertical scale factor of the V-grid.
        Default is e3v.
    var_e3t: string, optional
        Variable name of the vertical scale factor of the T-grid.
        Default is e3t.

    Returns
    -------
    gridU : xarray.Dataset
        Dataset of shape (t,z,y,x) containing all variables
        interpolated to U-point
    gridV : xarray.Dataset
        Dataset of shape (t,z,y,x) containing all variables
        interpolated to V-point

    """

    gridU = reduce_dataset(gridU, vars_to_keep, coords_to_keep).\
        rename({'depthu': 'z'})
    gridV = reduce_dataset(gridV, vars_to_keep, coords_to_keep).\
        rename({'depthv': 'z'})
    gridU = apply_mask(gridU, mask['umask'].isel(t=0))
    gridV = apply_mask(gridV, mask['vmask'].isel(t=0))

    gridU_v = interpolate_grids(gridU, "u2v")
    gridV_u = interpolate_grids(gridV, "v2u")
    gridU = xr.merge([gridU, gridV_u])
    gridV = xr.merge([gridV, gridU_v])
    if gridT:
        gridT = reduce_dataset(gridT, vars_to_keep, coords_to_keep).\
            rename({'deptht': 'z'})
        gridT = apply_mask(gridT, mask['tmask'].isel(t=0))
        gridT_v = interpolate_grids(gridT, "t2v")
        gridT_u = interpolate_grids(gridT, "t2u")
        gridU = xr.merge([gridU, gridT_u])
        gridV = xr.merge([gridV, gridT_v])
    return gridU, gridV


def select_auxillary_variables(variable, jj, ii, selection_index, axis):
    # This function is meant to be used by select_section(). Passing the
    # ij iterator instead of jj,ii does not work beacuse at the time the
    # function is called by select_section() the iterator has already been
    # exhausted.
    """Select auxillar variables that have the shape (1,y,z) or (1,z,y,x)
    according to ji.
    A shifted value (variable.shift(axis=1)) is selected where
    selection_index=1

    Parameters
    ---------
    variable : xarray.DataArray
        Array of shape (1,y,x) or (1,z,y,x)
    jj,ii :xarray.DataArray
        Arrays of length n providing jj,ii coordinates of the section.
    selectin_index : xarray.DataArray, Boolean
        Array of length n
    axis : String
        Axis along which variable is shifted.

    Returns
    -------
        variable_selection : xarray.DataArray
            Array of legnth n or shape (n,z) depending on input.
    """
    shifted = xr.concat([variable.isel(t=0),
                        variable.isel(t=0).shift({axis: -1})], dim='s')
    shifted = shifted.isel(x=ii[:-1], y=jj[:-1])
    shifted = shifted.isel(s=1).\
        where(selection_index == 1, shifted.isel(s=0)).\
        where(selection_index, 0)
    return shifted


def select_section(ji, gridU, gridV, mesh_hgr=None, mesh_zgr=None, mask=None,
                   var_e1v='e1v', var_e2u='e2u', var_e3u='e3u', var_e3v='e3v'):
    """Form a section by selecting quantities at U or V, depending on the
        orientation of each segment.
        We define velocities and transports, that go to the right of the
        progression of the section, as positiv.
        This means that the sign of the velocities depends on the order of the
        input points.
        For each leg, velocities are converted as follows:
            i0_leg < i1_leg -> V*(-1)
            j0_leg > j1_leg -> U*(-1)

    Parameters
    ----------
    ji : Iterator of tuples.
        Iterator of n tuples of ji-coordinates defining the section.
        If tuples include angle of section leg as last lement this angle is
        used to rotate velocites.
    gridU : xarray.Dataset
        Dataset of shape (t,z,y,x) holding variables at U-point.
        Use dataset returned by broken_line.reduce_dataset()
    gridV : xarray.Dataset
        Dataset of shape (t,z,y,x) holding variables at V-point.
        Use dataset returned by broken_line.reduce_dataset()
    mesh_hgr : xarray.Dataset, optional
        Dataset of shape (1,y,x) with horizontal grid scale factors.
        Variable names can be given below.
        Only needed for calculation of volume transport.
    mesh_zgr : xarray.Dataset, optional
        Dataset of shape (1,z,y,x) with vertical grid scale factors.
        Variable names can be given below.
        Only needed for calculation of volume transport.
    mask : xarray.Dataset, optional
        Dataset must include mask files for U-points (umask), V-points (vmask)
        and T-points (tmask) when a T-grid file is processed and a vertical
        axis nav_lev.
        Masks need need shape (t,z,y,x) and nav_lev length z.
    var_e1v: string, optional
        Variable name of the zonal scale factor of the V-grid.
        Default is e1v.
    var_e2u: string, optional
        Variable name of the meridional scale factor of the U-grid.
        Default is e1v.
    var_e3u: string, optional
        Variable name of the vertical scale factor of the U-grid.
        Default is e3u.
    var_e3v: string, optional
        Variable name of the vertical scale factor of the V-grid.
        Default is e3v.

    Returns
    -------
    section: xarray.Dataset
        Dataset of shape (t,z,n)


    """

    ji = list(zip(*ji))
    jj = ji[0]
    ii = ji[1]
    ii = xr.DataArray(list(ii), dims='c')
    jj = xr.DataArray(list(jj), dims='c')
    if len(ji) == 3:
        a = ji[2]
        a = xr.DataArray(list(a[:-1]), dims='c')
    iidiff = ii.diff('c')
    jjdiff = jj.diff('c')
    gridU = xr.concat([gridU, gridU.shift(y=-1)], dim='s').\
        rename({'vozocrtx': 'u_normal'})
    gridV = xr.concat([gridV, gridV.shift(x=-1)], dim='s').\
        rename({'vomecrty': 'u_normal'})
    if (("vomecrty" in gridU.data_vars) & ("vozocrtx" in gridV.data_vars)):
        gridU = gridU.rename({'vomecrty': 'u_along'})
        gridV = gridV.rename({'vozocrtx': 'u_along'})
    if (("sozotaux" in gridU.data_vars) & ("sometauy" in gridV.data_vars)):
        gridU = gridU.rename({'sozotaux': 'tau_normal'})
        gridV = gridV.rename({'sometauy': 'tau_normal'})
    if (("sometauy" in gridU.data_vars) & ("sozotaux" in gridV.data_vars)):
        gridU = gridU.rename({'sometauy': 'tau_along'})
        gridV = gridV.rename({'sozotaux': 'tau_along'})
    # We have one point more than segments
    gridU = gridU.isel(x=ii[:-1], y=jj[:-1])
    gridV = gridV.isel(x=ii[:-1], y=jj[:-1])
    #    select_section.chose_shifted.set_to_zero_where_not_needed
    gridU = gridU.isel(s=1).where(jjdiff == 1, gridU.isel(s=0))\
        .where(jjdiff, 0)
    gridV = gridV.isel(s=1).where(iidiff == 1, gridV.isel(s=0))\
        .where(iidiff, 0)
    #
    if ((len(ji) == 3) and mesh_hgr):
        # Compute and select inclination of coordinate system
        a_coord_x =\
            coordinate_rotation(mesh_hgr.gphif, mesh_hgr.glamf, 'x')
        a_coord_y =\
            coordinate_rotation(mesh_hgr.gphif, mesh_hgr.glamf, 'y')
        a_coord_x =\
            select_auxillary_variables(a_coord_x, jj, ii, iidiff, 'x')
        a_coord_y =\
            select_auxillary_variables(a_coord_y, jj, ii, jjdiff, 'y')
        a_coord = a_coord_x + a_coord_y
        # correct angles of section legs by local
        # inclination of coordinate system
        a = a - a_coord
        # rotate velocties
        u_rot_gridu, v_rot_gridu =\
            rotate_velocities(gridU['u_normal'], gridU['u_along'], a)
        u_rot_gridv, v_rot_gridv =\
            rotate_velocities(gridV['u_along'], gridV['u_normal'], a)
        gridU = gridU.assign({'u_rot_normal': u_rot_gridu,
                             'u_rot_along': v_rot_gridu})
        gridV = gridV.assign({'u_rot_normal': u_rot_gridv,
                             'u_rot_along': v_rot_gridv})
    # Define direction of positiv transport
    gridU['u_normal'] = gridU['u_normal'] * jjdiff
    gridV['u_normal'] = gridV['u_normal'] * iidiff * (-1)
    if (("u_along" in gridU.data_vars) & ("u_along" in gridV.data_vars)):
        gridU['u_along'] = gridU['u_along'] * jjdiff
        gridV['u_along'] = gridV['u_along'] * iidiff
    if (("tau_normal" in gridU.data_vars) & ("tau_normal" in gridV.data_vars)):
        gridU['tau_normal'] = gridU['tau_normal'] * jjdiff
        gridV['tau_normal'] = gridV['tau_normal'] * iidiff * (-1)
    if (("tau_along" in gridU.data_vars) & ("tau_along" in gridV.data_vars)):
        gridU['tau_along'] = gridU['tau_along'] * jjdiff
        gridV['tau_along'] = gridV['tau_along'] * iidiff

    section = gridU + gridV
    section = section.assign({'ii': ii[:-1]})
    section = section.assign({'jj': jj[:-1]})
    if mesh_hgr:
        # Select dx
        e2u =\
            select_auxillary_variables(mesh_hgr[var_e2u], jj, ii, jjdiff, 'y')
        e1v =\
            select_auxillary_variables(mesh_hgr[var_e1v], jj, ii, iidiff, 'x')
        section = section.assign({'dx': e2u + e1v})
        section['c'] = section['dx'].cumsum()
        # Select lat
        gphiu =\
            select_auxillary_variables(mesh_hgr['gphiu'], jj, ii, jjdiff, 'y')
        gphiv =\
            select_auxillary_variables(mesh_hgr['gphiv'], jj, ii, iidiff, 'x')
        section = section.assign({'lat': gphiu + gphiv})
        # Select lon
        glamu =\
            select_auxillary_variables(mesh_hgr['glamu'], jj, ii, jjdiff, 'y')
        glamv =\
            select_auxillary_variables(mesh_hgr['glamv'], jj, ii, iidiff, 'x')
        section = section.assign({'lon': glamu + glamv})
    if mesh_zgr:
        # Select dz
        e3u =\
            select_auxillary_variables(mesh_zgr[var_e3u], jj, ii, jjdiff, 'y')
        e3v =\
            select_auxillary_variables(mesh_zgr[var_e3v], jj, ii, iidiff, 'x')
        section = section.assign({'dz': e3u + e3v})
        section['z'] = mesh_zgr['nav_lev'].values
    if mask:
        # Select mask
        umask = select_auxillary_variables(mask['umask'], jj, ii, jjdiff, 'y')
        vmask = select_auxillary_variables(mask['vmask'], jj, ii, iidiff, 'x')
        section = section.assign({'mask': umask + vmask})
    return section


def calculate_transport(section, S_ref=34.8):
    """ Calculate volume, freshwater transport and heat transport
    across a section (ignoring ssh layer).
    Due to the definition of positiv velocities to the right of the section,
    the barotropic volume transport should be zero on a closed broken line.

    Parameters
    ----------
    section: xarray.Dataset
        Dataset returned by broken_line.section()
    S_ref: float, optional
        Reference salinity used for calculation of freshwater transport.
        Default  is S_ref=34.8

    Returns
    -------
    trsp : xarray.Dataset
        Dataset of shape (t)
        - if Salt is given contains total Freshwater Transport [fwt_trsp]
                  across section refrenced to S_ref
        - if Temp is given contains total Heat Transport [ht_trsp]
                  across the section using CP_ref and RHO_ref

    """
    # Volume transport
    trsp = section['u_normal'] * section['dx'] * section['dz']
    trsp_sum = trsp.sum(('c', 'z')) / 1e6
    ds = xr.Dataset({'trsp': trsp_sum})
    # Freshwater transport
    if 'vosaline' in section.data_vars:
        fw_trsp = (1 - section['vosaline'] / S_ref) * trsp
        fw_trsp = fw_trsp.sum(('c', 'z')) / 1e6   # converted into Sverdrup
        ds.update({'fw_trsp': fw_trsp})
        ds['fw_trsp'].attrs['units'] = 'm³/s'
        ds['fw_trsp'].attrs['ref_salinity'] = S_ref
    # Heat transport
    # use constant specific heat capacity CP_ref and desnity RHO_ref
    CP_ref = 4000.0
    RHO_ref = 1000.0
    if 'votemper' in section.data_vars:
        ht_trsp = RHO_ref * CP_ref * section['votemper'] * trsp
        ht_trsp = ht_trsp.sum(('c', 'z'))
        ds.update({'ht_trsp': ht_trsp})
        ds['ht_trsp'].attrs['units'] = 'W'        # Watt or Joule/second
        ds['ht_trsp'].attrs['ref_density'] = RHO_ref
        ds['ht_trsp'].attrs['ref_spec_heat_capacity'] = CP_ref
    return ds
