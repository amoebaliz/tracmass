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
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  !REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh
  !REAL,       ALLOCATABLE, DIMENSION(:,:,:)  :: dzt0
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

  dstamp      = 'xxxx/VIP-LD.HCo13T_avg_xxxx-xx-xxT00:00:00.nc'
!  dstamp      = 'xxxx/VIP-LD.HCo12T_avg_xxxx-xx-xxTxx:xx:00.nc'

  write (dstamp(1:4),'(i4.4)')   currYear
  write (dstamp(24:27),'(i4i2)') currYear
  write(dstamp(29:30),'(i2.2)')  currMon
  write(dstamp(32:33),'(i2.2)')  currDay
! write(dstamp(35:36),'(i2.2)')  currHour
! write(dstamp(38:39),'(i2.2)')  currMin
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  uvel        = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel        = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh         = get2dfieldNC(trim(dataprefix) ,'zeta')
  hs(:,:,2)   = ssh
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
!z_w(:,:,0) = depth
  do k=1,km
     dzt0 = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
     !dzt0(:,:,1) = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
     z_w(:,:,k) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
     !z_w(:,:,k,1) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:,1) 
 end do

 dzt(:,:,1:km-1,2) = z_w(:,:,2:km) - z_w(:,:,1:km-1)
 dzt(:,:,km,2) = ssh(:imt,:) - z_w(:,:,km)
 ! dzt(:,:,1:km-1,2) = z_w(:,:,2:km,1) - z_w(:,:,1:km-1,1)
 ! dzt(:,:,km,2) = ssh(:imt,:) - z_w(:,:,km,1)


  dzt(:,:,:,1)=dzt(:,:,:,2)
  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5
!  dzu(1:imt-1,:,:,1) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
!  dzv(:,1:jmt-1,:,1) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5
  do k=1,km
     !uflux(:,:,k,2)   = uvel(:,:,k) * dzu(:,:,k) * dyu
     !vflux(:,:,k,2)   = vvel(:,:,k) * dzv(:,:,k) * dxv
     uflux(:,:,k,2)   = uvel(:imt,:,k) * dzu(:,:,k) * dyu(:imt,:)
     vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k) * dxv(:imt,:)
 !    uflux(:,:,k,2)   = uvel(:imt,:,k) * dzu(:,:,k,1) * dyu(:imt,:)
  !   vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k,1) * dxv(:imt,:)

  end do

  if (nff .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

  return

end subroutine readfields
