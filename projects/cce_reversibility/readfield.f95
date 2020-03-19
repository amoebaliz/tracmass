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
  CHARACTER (len=200)                        :: dataprefix ,dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1, hr
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  
  ! ===   ===   ===

  alloCondUVW: if(.not. allocated (ssh)) then
    allocate ( ssh(imt,jmt), dzt0(imt,jmt) )
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

  dstamp      = CM4_his_daily_xxxx-xx-xx.nc

  write(dstamp(15:19),'(i4.4)')  currYear
  write(dstamp(21:22),'(i2.2)')  currMon
  write(dstamp(24:25),'(i2.2)')  currDay

  dataprefix  = trim(inDataDir) // dstamp

  tpos        = intpart1+1
  uvel        = get3DfieldNC(trim(dataprefix) , 'ssu')
  vvel        = get3DfieldNC(trim(dataprefix) , 'ssv')
  ssh         = get2dfieldNC(trim(dataprefix) , 'zos')

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

  hs(:,:,1) = ssh

  print *, 'MEEEEEP'

  do k=1,km
     uflux(:,:,k,2)   = uvel(:,:,k) * 1 * dyu(:imt,:jmt)
     vflux(:,:,k,2)   = vvel(:,:,k) * 1 * dxv(:imt,:jmt)
     dzt = 1
  end do

  if (nff .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if
  return

end subroutine readfields
