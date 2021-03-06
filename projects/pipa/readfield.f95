SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_traj
  USE mod_getfile
  use mod_seed

#ifdef tempsalt
  USE mod_dens
#endif

  IMPLICIT none  
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER (len=200)                          :: dstamp,datauprefix,datavprefix
  CHARACTER (len=200)                          :: fieldFile
 
  ! = Loop variables
  INTEGER                                      :: i, j, k
  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:)     :: rhozvec, depthzvec, latvec
 REAL*4 ,ALLOCATABLE, DIMENSION(:)     :: tmpzvec, salzvec
#endif

 LOGICAL around

!---------------------------------------------------------------
 call datasetswap
 call updateclock

initFieldcond: if(ints.eq.intstart) then
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif

endif initFieldcond

! === EXTRACTING GCM FIELDS FOR GIVEN TIME ===
 dstamp      = 'GLORYS_xx_xx_'
 write(dstamp(8:9),'(i2.2)') currMon
 write(dstamp(11:12),'(i2.2)') currDay
 print *, (trim(inDataDir)//trim(dstamp)//'U.nc')
 uvel = get3DfieldNC(trim(inDataDir)//trim(dstamp)//'U.nc', 'vozocrtx')
 vvel = get3DfieldNC(trim(inDataDir)//trim(dstamp)//'V.nc', 'vomecrty')
!#ifdef tempsalt
! tem(:,:,:,2) = get3DfieldNC(trim(fieldFile)//'T.nc', 'votemper')
! sal(:,:,:,2) = get3DfieldNC(trim(fieldFile)//'S.nc', 'vosaline')
!#endif /*tempsalt*/

! Remove missing value place holders
where (uvel > 1000)
     uvel = 0
end where
where (vvel > 1000)
     vvel = 0
end where

do k=1,KM
  do j=1,JMT
     do i = 1,IMT
        uflux(i,j,km+1-k,nsp) = uvel(i,j,k) * dyu(i,j) * dzu(i,j,k)
        vflux(i,j,km+1-k,nsp) = vvel(i,j,k) * dxv(i,j) * dzv(i,j,k)

        ! Testing There are no velocity value assigned below sea floor at U and V points  
        if(k>kmv(i,j) .and. vvel(i,j,k)/=0. ) then
          print *,'vflux=',vflux(i,j,km+1-k,nsp), vvel(i,j,k),i,j,k,kmv(i,j),nsp
          stop 4966
        endif

        if(k>kmu(i,j) .and. uvel(i,j,k)/=0. ) then
          print *,'uflux=',uflux(i,j,k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
          stop 4966
        endif


     enddo
  enddo
enddo
return
end subroutine readfields
