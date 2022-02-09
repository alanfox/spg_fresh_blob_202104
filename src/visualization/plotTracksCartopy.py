import cmocean as co
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime, timedelta
import seaborn as sns
import cmocean as co
import pandas as pd
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy

def plotTracksCartopy(ds1,cmap,title='',figname=''):
    
    central_lon, central_lat = -45, 47.5
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize = (12,8),subplot_kw={'projection': ccrs.Orthographic(central_lon, central_lat)})
    extent = [-80, -10, 25, 70]
    ax.set_extent(extent)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    
    dsmask = ds1

    pcm = ax.scatter(
        dsmask.lon.data.flatten(),
        dsmask.lat.data.flatten(),
        3,
#         (dsmask.time.data.flatten()-ds.time.isel(obs=0,traj=0).data).astype('timedelta64[D]')
        mdates.date2num(dsmask.time.data.flatten())
#         ,cmap= cmap
        ,zorder=2
        ,transform=ccrs.PlateCarree()
    #   ,alpha=0.3
    )
    cb = fig.colorbar(pcm,ax=ax,shrink=0.8,label = "date")
    loc = mdates.MonthLocator(bymonth=range(1,13,6))
    cb.ax.yaxis.set_major_locator(loc)
    cb.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
#     cb.ax.tick_params(labelsize=20)


    ax.scatter(
        dsmask.lon.isel(obs=0).data.flatten(),
        dsmask.lat.isel(obs=0).data.flatten(),6,zorder=5,
        transform=ccrs.PlateCarree()
    )
    ax.plot(np.linspace(-58.2,-40,num=20),np.linspace(52,65,num=20),
        zorder=5,color='C1',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    ax.plot(np.linspace(-80,-58.2,num=20),np.linspace(52,52,num=20),
        zorder=5,color='C1',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    ax.plot(np.linspace(-80,-60,num=20),np.linspace(51,51,num=20),
        zorder=5,color='C4',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    ax.plot(np.linspace(-60,-60,num=20),np.linspace(51,30,num=20),
        zorder=5,color='C4',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    
    ax.plot(np.linspace(-44,-100,num=20),np.linspace(33,33,num=20),
        zorder=5,color='C2',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    ax.plot(np.linspace(-44,-44,num=20),np.linspace(0,33,num=20),
        zorder=5,color='C2',
        linewidth=3,
        transform=ccrs.PlateCarree()
    )
    
    title=ax.set_title(title)
    
    #     plt.savefig(project_path / figure_path / 'plot_transports_by_source_37WtoScot_500M.eps', bbox_extra_artists=(lgd,))
    plt.savefig(figname + '.png', bbox_inches='tight')
    plt.savefig(figname + '.pdf', bbox_inches='tight')


    return
