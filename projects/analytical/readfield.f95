SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  USE mod_seed, only: nff! LD ADDED, for nff  
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1,hr
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_r,Cs_r
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  ! ===   ===   ===

  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt), dzt0(imt,jmt) ) ! NOTE: omits (:, eta_rho[-1]) and (xi_rho[-1],:) 
     allocate ( sc_w(km), Cs_w(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
  initFieldcond: if(ints.eq.intstart) then
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
  end if initFieldcond

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_w = 0
  Cs_w = 0

  call updateclock
  call datasetswap
  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = 'test_xxxx_xx-xx.nc'
  write (dstamp(6:9),'(i4i2)') currYear
  write(dstamp(11:12),'(i2.2)')  currMon
  write(dstamp(14:15),'(i2.2)')  currDay
!  write(dstamp(35:36),'(i2.2)')  currHour
!  write(dstamp(38:39),'(i2.2)')  currMin

  dataprefix  = trim(inDataDir) // dstamp
!  print *, dataprefix
  tpos        = intpart1+1
  uvel        = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel        = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh         = get2dfieldNC(trim(dataprefix) ,'zeta')
#ifdef explicit_w
  wvel      = get3DfieldNC(trim(dataprefix) ,'omega')
#endif
  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where
  where (ssh > 1000)
     ssh = 0
  end where

  hs(:imt,:jmt,2) = ssh(:imt,:jmt)

#ifdef explicit_w
  wflux(:,:,:,2) = 0.
  do j=1,jmt
    do i=1,imt
      wflux(i,j,0:km-1,2) = wvel(i,j,1:km)*dxdy(i,j)
    end do
  end do
#endif

  Cs_w = get1DfieldNC (trim(dataprefix), 'Cs_w')
  sc_w = get1DfieldNC (trim(dataprefix), 's_w')
  hc   = getScalarNC (trim(dataprefix), 'hc')

  do k=1,km
     dzt0 = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
     z_w(:,:,k) = ssh + (ssh + depth) * dzt0
  end do
 
  dzt(:,:,1:km-1,2) = z_w(:,:,2:km) - z_w(:,:,1:km-1)
  dzt(:,:,km,2) = ssh - z_w(:,:,km)
  dzt(:,:,:,1)=dzt(:,:,:,2)

 ! NOTE: not filling last i and j positions of dzu and dzv, respectively... so = 0
  !       ergo, you do NOT want to be calculating particle trajectory with:
  !       1) uflux from points (imt,:) 
  !       2) vflux from points (:,jmt)
  ! SOLUTION: make kill zone for northern and eastern boundary 2 points wide in maphil.in

  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5

  do k=1,km
     uflux(:,:,k,2)   = uvel(:,:,k) * dzu(:,:,k) * dyu(:imt,:jmt)
     vflux(:,:,k,2)   = vvel(:,:,k) * dzv(:,:,k) * dxv(:imt,:jmt)
  end do

  if (nff .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

#ifdef tempsalt
  tem(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'temp')
  sal(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'salt')
  !rho(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'rho')
#ifdef larval_fish
  srflux(:,:,2)     = get2DfieldNC(trim(dataprefix) ,   'swrad')
! Note: this works as long as surface AKt is zero.
  ak2(:,:,:)        = get3DfieldNC(trim(dataprefix) ,   'AKt')
  akt(:,:,0:km-1,2) = ak2(:,:,:)
#endif
#endif
  return

end subroutine readfields
