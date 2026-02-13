import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
import os

path_ins = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-instant.nc"
ds = xr.open_dataset(path_ins)
ds_reg = ds.sel(latitude=slice(12, -55), longitude=slice(262, 328))

frames = []
for t in range(len(ds_reg.valid_time)):
    fig = plt.figure(figsize=(8, 10))
    ax = plt.axes([0.1, 0.1, 0.85, 0.85], projection=ccrs.PlateCarree())
    ax.set_extent([262, 328, -55, 12])
    ax.set_aspect('auto', adjustable='datalim')
    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    ax.gridlines(draw_labels=True, alpha=0.3)

    u, v = ds_reg['u10'].isel(valid_time=t), ds_reg['v10'].isel(valid_time=t)
    mag = np.sqrt(u**2 + v**2)
    
    cf = ax.contourf(ds_reg.longitude, ds_reg.latitude, mag, levels=np.linspace(0, 25, 50), cmap='YlGnBu')
    ax.quiver(ds_reg.longitude[::10], ds_reg.latitude[::10], u[::10,::10], v[::10,::10], 
              color='black', alpha=0.5, scale=400, transform=ccrs.PlateCarree())
    
    plt.colorbar(cf, label='Rapidez (m/s)', orientation='horizontal', pad=0.08, shrink=0.8)
    fecha = str(ds_reg.valid_time.values[t])[:16]
    plt.title(f'VIENTOS 10M - {fecha} UTC', fontweight='bold')
    
    tmp = f'uv_{t:02d}.png'
    plt.savefig(tmp, dpi=120, bbox_inches='tight', pad_inches=0.05)
    frames.append(tmp); plt.close()

imageio.mimsave('Sudamerica_Vientos_Final.gif', [imageio.imread(f) for f in frames], duration=400, loop=0)
for f in frames: os.remove(f)


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
import os

path_ins = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-instant.nc"
ds = xr.open_dataset(path_ins)
ds_reg = ds.sel(latitude=slice(12, -55), longitude=slice(262, 328))

frames = []
for t in range(len(ds_reg.valid_time)):
    fig = plt.figure(figsize=(8, 10))
    ax = plt.axes([0.1, 0.1, 0.85, 0.85], projection=ccrs.PlateCarree())
    ax.set_extent([262, 328, -55, 12])
    ax.set_aspect('auto', adjustable='datalim')
    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    ax.gridlines(draw_labels=True, alpha=0.3)

    data = ds_reg['ie'].isel(valid_time=t) * 1000
    cf = ax.contourf(ds_reg.longitude, ds_reg.latitude, data, levels=np.linspace(-0.2, 0.2, 51), cmap='RdYlBu', extend='both')
    
    plt.colorbar(cf, label='Flujo (mm)', orientation='horizontal', pad=0.08, shrink=0.8)
    fecha = str(ds_reg.valid_time.values[t])[:16]
    plt.title(f'FLUJO DE HUMEDAD - {fecha} UTC', fontweight='bold')
    
    tmp = f'ie_{t:02d}.png'
    plt.savefig(tmp, dpi=120, bbox_inches='tight', pad_inches=0.05)
    frames.append(tmp); plt.close()

imageio.mimsave('Sudamerica_Flujo_Humedad_Final.gif', [imageio.imread(f) for f in frames], duration=400, loop=0)
for f in frames: os.remove(f)



import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio.v2 as imageio
import os

path_acc = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-accum.nc"
ds = xr.open_dataset(path_acc)
ds_reg = ds.sel(latitude=slice(12, -55), longitude=slice(262, 328))

frames = []
for t in range(len(ds_reg.valid_time)):
    fig = plt.figure(figsize=(8, 10))
    ax = plt.axes([0.1, 0.1, 0.85, 0.85], projection=ccrs.PlateCarree())
    ax.set_extent([262, 328, -55, 12])
    ax.set_aspect('auto', adjustable='datalim')
    
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
    ax.gridlines(draw_labels=True, alpha=0.3)

    data = ds_reg['e'].isel(valid_time=t) * 1000
    cf = ax.contourf(ds_reg.longitude, ds_reg.latitude, data, levels=np.linspace(-1, 0.2, 50), cmap='RdYlGn', extend='both')
    
    plt.colorbar(cf, label='Evaporación (mm)', orientation='horizontal', pad=0.08, shrink=0.8)
    fecha = str(ds_reg.valid_time.values[t])[:16]
    plt.title(f'EVAPORACIÓN - {fecha} UTC', fontweight='bold')
    
    tmp = f'e_{t:02d}.png'
    plt.savefig(tmp, dpi=120, bbox_inches='tight', pad_inches=0.05)
    frames.append(tmp); plt.close()

imageio.mimsave('Sudamerica_Evaporacion_Final.gif', [imageio.imread(f) for f in frames], duration=400, loop=0)
for f in frames: os.remove(f)