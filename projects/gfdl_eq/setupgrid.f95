SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------



  ! === Init local variables for the subroutine ===
  INTEGER                                    :: i ,j ,k ,kk
  CHARACTER (len=200)                        :: gridfile

! === Template for setting up grids. Move the code from readfile.f95
 ALLOCATE ( depth(imt,jmt) )  !LD: omits (:, eta_rho[-1]) and (xi_rho[-1],:) 
 ALLOCATE ( z_w(imt,jmt,km) ) !LD: made 3D; 
  !Order is   t  k  i  j 
  map2d    = [3, 4, 1, 2]
  map3d    = [2, 3, 4, 1]
 
  gridfile = "/nbhome/Liz.Drenkard/tracmass_stuff/test_grd.nc"
  ncTpos = 1
  print *, trim(gridfile)
  dxv = get2DfieldNC(trim(gridfile), 'dxt')
  dyu = get2DfieldNC(trim(gridfile), 'dyt')  
  dxdy = dyu*dxv
  depth = get2DfieldNC(trim(gridfile), 'deptho')
  mask = get2DfieldNC(trim(gridfile), 'wet')
  kmt = 1 

end SUBROUTINE setupgrid
