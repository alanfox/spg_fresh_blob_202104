```python
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

import importlib
import sys
```


```python
from dask.distributed import Client
```


```python
client = Client(n_workers=1, threads_per_worker=8, memory_limit=10e9)
client
```


```python
# parameters
project_path = Path.cwd() / '..' / '..' 
project_path = project_path.resolve()

data_path = "/data/spg_fresh_blob_202104_data/raw/"

interim_data_path = Path('data/interim/endtracks/plusDist/')

sectionPath = Path('data/external/')
sectionFilename = 'osnap_pos_wp.txt'
sectionname = 'osnap'

year = 1990
```


```python
year_str = str(year)
```


```python
interim_data_filename = "endtracks_20211215_randomvel_mxl_osnap_backwards_"+year_str+"_subset_10percent.nc"
```


```python
# data_stores_subsets = list(sorted(Path(data_path).glob("*_????_subset.zarr/")))[:use_number_subset_years]
data_stores_subsets = list(sorted(Path(data_path).glob("*_"+year_str+"_subset_10percent.zarr/")))
```


```python
display(data_stores_subsets)
```


```python
ds_subsets = xr.concat(
    [xr.open_zarr(store) for store in data_stores_subsets],
    dim="traj",
)

display(ds_subsets)
print(ds_subsets.nbytes / 1e9, "GiB")
```

## add cumulative track lengths to file


```python
ds_subsets_lat_diff=ds_subsets.lat.diff(dim='obs',n=1)
ds_subsets_lon_diff=ds_subsets.lat.diff(dim='obs',n=1)

```


```python
# Calculate Great Circle distance between requeted point and 
# matched grid point. 
# Based on https://andrew.hedges.name/experiments/haversine/

# copied from internet

# approximate radius of earth in km
R = 6373000.0

lat1 = xr.ufuncs.radians(ds_subsets.lat)
lon1 = xr.ufuncs.radians(ds_subsets.lon)
lat2 = lat1.shift(obs=1)
lon2 = lon1.shift(obs=1)

dlat = lat2 - lat1
dlon = lon2 - lon1

a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

distance = R * c

```


```python
ds_subsets = ds_subsets.assign({'dist':distance.cumsum(dim='obs')})
```


```python
ds_subsets.dist.attrs = {'long_name':'alongtrack distance from OSNAP'
                        ,'units':'metres'}

```


```python
ds_subsets
```

## Extract data on osnap line (initialised positions)


```python
ds_subsets_osnap = ds_subsets.isel(obs=0)
```

## Update some attributes


```python
# ds_subsets_osnap['vol_trans_normal'] = 
ds_subsets_osnap.mxl.attrs = {'units':'m','long_name':'mixed layer depth'}
ds_subsets_osnap.salt.attrs = {'units':'PSU','long_name':'salinity'}
ds_subsets_osnap.temp.attrs = {'units':'degC','long_name':'temperature'}
ds_subsets_osnap.uvel.attrs = {'units':'degrees_east/second','long_name':'u velocity (raw)'}
ds_subsets_osnap.vvel.attrs = {'units':'degrees_north/second','long_name':'v velocity (raw)'}
```

## Flag tracks by source region and pathway


```python
def apply_left_of_line(ds, lon_1, lon_2, lat_1, lat_2):
    '''Apply an area crossing criterion.
    
    Larvae in ds selected while they are in a selected area.
    '''
    # particles are selected if they pass through given area.
    position =  ((lon_2 -lon_1) * (ds.lat - lat_1) - 
                     (ds.lon - lon_1) * (lat_2 - lat_1))
                        
    return position > 0.0, position < 0
```

#### from Labrador sea or from Gulf Stream


```python
# from labrador sea
ds_in1, ds_notin1 = apply_left_of_line(ds_subsets,-75,-40,40,65)
ds_in2, ds_notin2 = apply_left_of_line(ds_subsets,-100,-58.2,48,52)
ds_in3, ds_notin3 = apply_left_of_line(ds_subsets,-45,-45,60,70)
ds_lab_in = ds_in1*ds_in2*ds_in3
# from west of 60W, south of Flemish Cap (to test path from labrador sea)
ds_in1, ds_notin1 = apply_left_of_line(ds_subsets,-60,-60,33,63)
ds_in2, ds_notin2 = apply_left_of_line(ds_subsets,-58.2,-100,52,48)
ds_60w_in = ds_in1*ds_in2
# from gulf stream
ds_in1, ds_notin1 = apply_left_of_line(ds_subsets,-60,-100,33,33)
ds_in2, ds_notin2 = apply_left_of_line(ds_subsets,-44,-44,0,33)
ds_gst_in = ds_in1 * ds_in2

```


```python
# check trajectory routes
LabCu = ds_lab_in.max("obs")
LC60W = ds_60w_in.max("obs")
GulfS = ds_gst_in.max("obs")

# check when lef lab sea,crossed 60w or gulf stream. defaults to zero
LabCu_exit_index = ds_lab_in.argmax(axis=1)
LC60W_exit_index = ds_60w_in.argmax(axis=1)
GulfS_exit_index = ds_gst_in.argmax(axis=1)

LabCu_exit_index = LabCu_exit_index.where(LabCu_exit_index > 0,len(ds_subsets.obs)-1)
LC60W_exit_index = LC60W_exit_index.where(LC60W_exit_index > 0,len(ds_subsets.obs)-1)
GulfS_exit_index = GulfS_exit_index.where(GulfS_exit_index > 0,len(ds_subsets.obs)-1)
```


```python
# check if most recently left Lab Sea
LabCu_is_source = (LabCu_exit_index < GulfS_exit_index)
GulfS_is_source = (LabCu_exit_index > GulfS_exit_index)
LC60W_is_path = (LabCu_exit_index > LC60W_exit_index).where(LabCu_is_source,False)
```


```python
LCdir_is_path = LabCu_is_source.where(LC60W_is_path==False,False)
other_is_source = (LabCu_is_source==False).where(GulfS_is_source == False,False)
```

flag particles on osnap line by origin


```python
ds_subsets_osnap = ds_subsets_osnap.assign({'LabCu_is_source':LabCu_is_source})
ds_subsets_osnap = ds_subsets_osnap.assign({'LC60W_is_path':LC60W_is_path})
ds_subsets_osnap = ds_subsets_osnap.assign({'LCdir_is_path':LCdir_is_path})
ds_subsets_osnap = ds_subsets_osnap.assign({'GulfS_is_source':GulfS_is_source})
ds_subsets_osnap = ds_subsets_osnap.assign({'other_is_source':other_is_source})

```

flag by pathway for Labrador Current parcels


```python
ds_subsets_osnap.LabCu_is_source.attrs = {'long_name':'flag from Labrador Current'}
ds_subsets_osnap.LC60W_is_path.attrs   = {'long_name':'flag from LC via 60W'}
ds_subsets_osnap.LCdir_is_path.attrs   = {'long_name':'flag from LC direct'}
ds_subsets_osnap.GulfS_is_source.attrs = {'long_name':'flag from Gulf Stream'}
ds_subsets_osnap.other_is_source.attrs = {'long_name':'flag source not found'}
```

## Find the 'obs' index of the point where parcel leaves the source region

### known source regions


```python
# test individaul positions to see when the source was left
# defaults to zero if particular source was not on track
LabCu_exit_index = (ds_lab_in.where(ds_subsets_osnap.LabCu_is_source,False)).argmax(axis=1)
GulfS_exit_index = (ds_gst_in.where(ds_subsets_osnap.GulfS_is_source,False)).argmax(axis=1)

# combine in to one array. Nonzero numbers should not overlap for 
# Lab current and Gulf Stream
exit_index = LabCu_exit_index + GulfS_exit_index
# convert zeros to max dim obs
exit_index = exit_index.where(exit_index > 0,len(ds_subsets.obs)-1)
```

### index last non nan value for 'other' parcels


```python
a = ds_subsets.lat  # just a random selection of variable, nans the same for all variables
b = (~np.isnan(a)).cumsum(dim='obs').argmax(dim='obs') # finds last non-nan in dim 'obs'. nicked from stackoverflow search
exit_index = xr.ufuncs.minimum(exit_index,b).compute()
```

### extract source positions from full array


```python
ds_subsets_sourc = ds_subsets.isel(traj=xr.DataArray(range(len(ds_subsets.traj)),dims='traj'),obs=exit_index)
```

add the source and pathway flags to ds_subsets_sourc to match for xr.concat


```python
# flag particles on osnap line by origin and pathway
ds_subsets_sourc = ds_subsets_sourc.assign({'LabCu_is_source':LabCu_is_source})
ds_subsets_sourc = ds_subsets_sourc.assign({'LC60W_is_path':LC60W_is_path})
ds_subsets_sourc = ds_subsets_sourc.assign({'LCdir_is_path':LCdir_is_path})
ds_subsets_sourc = ds_subsets_sourc.assign({'GulfS_is_source':GulfS_is_source})
ds_subsets_sourc = ds_subsets_sourc.assign({'other_is_source':other_is_source})

```


```python
ds_subsets_sourc.LabCu_is_source.attrs = {'long_name':'flag from Labrador Current'}
ds_subsets_sourc.LC60W_is_path.attrs   = {'long_name':'flag from LC via 60W'}
ds_subsets_sourc.LCdir_is_path.attrs   = {'long_name':'flag from LC direct'}
ds_subsets_sourc.GulfS_is_source.attrs = {'long_name':'flag from Gulf Stream'}
ds_subsets_sourc.other_is_source.attrs = {'long_name':'flag source not found'}
```

### combine source and osnap positions and characteristics


```python
ds_subsets_paths = xr.concat([ds_subsets_osnap,
                             ds_subsets_sourc],
                             dim='ends')
```

### flag particles entering from north from Greenland Sea or Davis Strait


```python
b = b.compute()
```


```python
ds_subsets_domexi = ds_subsets.isel(traj=xr.DataArray(range(len(ds_subsets.traj)),dims='traj'),obs=b)
```


```python
# counts how many times enters or leaves lab sea
spgnoloop = (abs(ds_lab_in.astype(int).diff(dim='obs')).sum(dim='obs')<3).compute()
```


```python
# from Hudson Bay
ds_in1, ds_notin1 = apply_left_of_line(ds_subsets,-68,-68,33,63)
ds_in2, ds_notin2 = apply_left_of_line(ds_subsets,-95,-60,52,52)
ds_hud_in = ds_in1*ds_in2

HudBa = ds_hud_in.max("obs")

```


```python
HudBa.data.compute()
```


```python
Green_is_source = (ds_subsets_paths.isel(ends=0).LabCu_is_source &
                  (ds_subsets_domexi.lat > 65) & 
                  (ds_subsets_domexi.lon > -44) &
                  spgnoloop)
Davis_is_source = (ds_subsets_paths.isel(ends=0).LabCu_is_source &
                  (ds_subsets_domexi.lat > 65) & 
                  (ds_subsets_domexi.lon < -44) &
                  spgnoloop)
Hudba_is_source = (ds_subsets_paths.isel(ends=0).LabCu_is_source & 
                  HudBa &
                  spgnoloop)
```


```python
Green_is_source = xr.concat([(Green_is_source),(Green_is_source)],dim='ends')
ds_subsets_paths = ds_subsets_paths.assign({'Green_is_source':Green_is_source})
Davis_is_source = xr.concat([(Davis_is_source),(Davis_is_source)],dim='ends')
ds_subsets_paths = ds_subsets_paths.assign({'Davis_is_source':Davis_is_source})
Hudba_is_source = xr.concat([(Hudba_is_source),(Hudba_is_source)],dim='ends')
ds_subsets_paths = ds_subsets_paths.assign({'Hudba_is_source':Hudba_is_source})

```

add the source and pathway flags to ds_subsets_sourc to match for xr.concat


```python
ds_subsets_paths.Davis_is_source.attrs = {'long_name':'flag from Davis Strait'}
ds_subsets_paths.Green_is_source.attrs = {'long_name':'flag from Greenland Sea'}
ds_subsets_paths.Hudba_is_source.attrs = {'long_name':'flag from Hudson Bay'}

```

## We want to test for tracks which route north of osnap line between source and final times

This is because a common strategy is to remove these from analysis and only consider the 'direct' paths.

### Test particle positions

### section position data


```python
lonlat = xr.Dataset(pd.read_csv(project_path / sectionPath / sectionFilename,delim_whitespace=True))
```

#### south/north of osnap-e


```python
# do north and south separately because of missing values

south = xr.Dataset()
north = xr.Dataset()
epsilon = 0.05
for i in range(len(lonlat.lon)-1):
    south['subsect'+str(i)],north['subsect'+str(i)] = apply_left_of_line(ds_subsets,lonlat.lon[i+1],lonlat.lon[i],lonlat.lat[i+1]+epsilon,lonlat.lat[i]+epsilon)

# check in osnap east 
south_oe,north_oe = apply_left_of_line(ds_subsets,-44,-44,30,60)

# south_a = south.subsect0 + south.subsect1 + south.subsect2 
# south_b = south.subsect3 * south.subsect4 * south.subsect5
south_b = south.subsect4 * south.subsect5
south_c = south.subsect6 + south.subsect7 + south.subsect8 
south_d = south.subsect8 * south.subsect9 * south.subsect10 * south.subsect11 
# south_e = south.subsect12
# south_all = south_a * south_c * south_e * (south_b + south_d)
south_all = south_oe + (south_c * (south_b + south_d))

# north_a = north.subsect0 * north.subsect1 * north.subsect2 
# north_b = north.subsect3 + north.subsect4 + north.subsect5
north_b = north.subsect4 + north.subsect5
north_c = north.subsect6 * north.subsect7 * north.subsect8 
north_d = north.subsect8 + north.subsect9 + north.subsect10 + north.subsect11
# north_e = north.subsect12
# north_all = north_a + north_c + north_e + (north_b * north_d)
north_all = north_oe * (north_c  + (north_b * north_d))

```


```python
north_all = north_all.reset_coords(drop=True)
south_all = south_all.reset_coords(drop=True)
```


```python
# test individual positions to see when the parcel was first (in 'obs', last in time) north of osnap line
# defaults to zero if particular source was not on track
north_osnap_index = (north_all).argmax(axis=1)
# convert zeros to max dim obs
north_osnap_index = north_osnap_index.where(north_osnap_index > 0,len(ds_subsets.obs)-1)
```

Check if found north of osnap e between leaving source and arriving at osnap and flag to ds_subsets_paths


```python
north_osnap = xr.concat([(north_osnap_index < exit_index),(north_osnap_index < exit_index)],dim='ends')
ds_subsets_paths = ds_subsets_paths.assign({'north_of_osnap':north_osnap})
ds_subsets_paths.north_of_osnap.attrs = {'long_name':'flag path goes north of osnap-e'}
display(ds_subsets_paths)
```

## Output to netcdf


```python
ds_subsets_paths.to_netcdf(project_path / interim_data_path / interim_data_filename)
```


```python
conda list
```


```python

```


```python

```
