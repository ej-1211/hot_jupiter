# Hot Jupiter

# MODES Package
### Vertical Structure Function
> - No input data for u v z yet
> - everthing is real (r8)
1. ```mod_numvsf.sigma.f90``` (Using ```/share/mod_adm.f90``` and ```/share/mod_const.f90```)
a. Read the surface stability $\texttt{stab(1)}$
b. Read and calcualte $d\sigma$ from $\texttt{vgrid}$


2. ```mod_vsf_driver.f90``` (Using ```/share/mod_adm.f90``` , ```/share/mod_const.f90```,```/share/mod_write_netcdf.f90``` and ```mod_numvsf.sigma.f90```)


Input : Stability, vgrid
Output : Equivalent Height ```equivalent_height.data```, Vertical Structure Function ```vsf.data.nc```

### Horizontal 