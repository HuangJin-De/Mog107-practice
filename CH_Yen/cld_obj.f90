PROGRAM cloud_object
USE netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/03/24                                   !
! Thisi program aims to                                   !
! 1) identify cloud objects                               !
! 2) calculate the size and centroid of each cloud object.!
! Inputs:  satellite precipitation file                   !
! Outputs: cloud objects data (.txt)                      !
!          cloud ID and size (.nc)                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
INTEGER :: ncid, dimid, intvl=8, varid, varid1, varid2, varid3, varid4, t_dimid, i_dimid, j_dimid, accum, cnt, newcnt
INTEGER, DIMENSION(3) :: dimids, start, ccccc
INTEGER :: i, ii, j, jj, t, tt, ttt, y, yy, k, ni, nj, total = 0, ci2, cj2
                       ! now_index, centroid_index, cnt: counts of grids within
                       ! a cloud, idnum: idnumber of IDs of clouds
REAL :: radius, ci, cj, idnum
REAL, PARAMETER :: crit = 1.
REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lon
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PCP1, PCP2
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ID, ID2, SZE2
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: SZE ! accum, cnt, newcnt must be integer while calculating size
                                                    ! change to REAL beofre
                                                    ! nf90_put_var
! nc file dosen't accept integer type, neither dose it accept imcomplete array
! => duplicate ID and SZE to ID2 and SZE2 
INTEGER, DIMENSION(:,:), ALLOCATABLE :: CLD 

CHARACTER(len=50), PARAMETER :: path = '/data/dadm1/obs/TRMM/'
CHARACTER(len=50) :: filename, output

INTEGER :: t1, t2, tt1, tt2, clock_rate, clock_max
CHARACTER(len=9) :: dimnameii
CHARACTER(len=8) :: dimnamejj
CHARACTER(len=4) :: dimnamett

CHARACTER(len=4) :: year
INTEGER :: init_yr, finl_yr !1998!2015
PRINT*, 'enter the initial year and final year:'
READ(*,'(I4)') init_yr
READ(*,'(I4)') finl_yr
yy = finl_yr-init_yr+1

CALL system_clock ( t1, clock_rate, clock_max )
DO y = 1, yy
  CALL system_clock ( tt1, clock_rate, clock_max )
  ! read longitude and latitude first
  WRITE(year, '(I4)') init_yr+y-1
  PRINT*, 'year = ',year
  output = 'new_cld_obj/'//trim('cloud_')//trim(year)//trim('.txt')
  OPEN(10, FILE=output, FORM='FORMATTED', STATUS='UNKNOWN', ACCESS='SEQUENTIAL')
  filename = trim('TRMM3B42/3B42.')//trim(year)//trim('.3hr.nc')
  CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
  CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
  CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
  CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
  CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
  CALL check(nf90_inq_dimid(ncid, 'time', dimid))
  CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))
  
  ALLOCATE(Lon(ii))
  ALLOCATE(Lat(jj))
  ALLOCATE(PCP1(ii,jj,tt))
  ALLOCATE(PCP2(ii+2,jj+2,tt)); PCP2 = 0.
  ALLOCATE(ID(ii+2,jj+2,tt/intvl)); ID = 0.
  ALLOCATE(SZE(ii+2,jj+2,tt/intvl)); SZE = 0
  ALLOCATE(ID2(ii,jj,tt/intvl)); ID2 = 0.
  ALLOCATE(SZE2(ii,jj,tt/intvl)); SZE2 = 0.

  CALL check(nf90_inq_varid(ncid, 'pcp', varid))
  CALL check(nf90_get_var(ncid, varid, PCP1))
  CALL check(nf90_inq_varid(ncid, 'longitude', varid))
  CALL check(nf90_get_var(ncid, varid, Lon))
  CALL check(nf90_inq_varid(ncid, 'latitude', varid))
  CALL check(nf90_get_var(ncid, varid, Lat))


  PCP2(2:ii,2:jj,:) = PCP1
  DEALLOCATE(PCP1)
  
  PRINT*, 'START SEARCHING CLOUD!!'
  DO t = 1, tt, intvl
    ttt = (t-1)/intvl+1
    idnum = 0
    total = 0
    DO j = 2, jj
      DO i = 2, ii
        ! if finds a grid with pcp and yet to be counted
        ! => activate the counting of new cloud
        IF ((PCP2(i,j,t) > crit) .AND. ID(i,j,ttt) == 0) THEN
          ! store the information of the first grid of the new found cloud
          idnum = idnum+1
          !PRINT '(A7,I5,A10)', 'Cloud #', idnum, ' detected!'
          ni = i; nj = j
          ID(ni,nj,ttt) = idnum
          cnt = 1; accum = 1
          ALLOCATE(CLD(ii*jj-total,2)); CLD = 0
          CLD(cnt,1) = ni; CLD(cnt,2) = nj
          ! collecte the neighbor grids belonging to this new cloud and iterate the
          ! process until all grids collected
          DO WHILE (newcnt /= 0 .OR. accum == 1)
            ! discover more new grids step by step from already found grids and
            ! collect it
            !PRINT*, '1'
            newcnt = 0
            DO k = accum-cnt+1, accum
              !PRINT*, 'k=',k
              ni = CLD(k,1); nj = CLD(k,2)
              IF (PCP2(ni+1,nj,t) > crit .AND. ID(ni+1,nj,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni+1; CLD(accum+newcnt,2) = nj
                ID(ni+1,nj,ttt) = idnum
              END IF
              IF (PCP2(ni+1,nj+1,t) > crit .AND. ID(ni+1,nj+1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni+1; CLD(accum+newcnt,2) = nj+1
                ID(ni+1,nj+1,ttt) = idnum
              END IF
              IF (PCP2(ni,nj+1,t) > crit .AND. ID(ni,nj+1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni; CLD(accum+newcnt,2) = nj+1
                ID(ni,nj+1,ttt) = idnum
              END IF
              IF (PCP2(ni-1,nj+1,t) > crit .AND. ID(ni-1,nj+1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni-1; CLD(accum+newcnt,2) = nj+1
                ID(ni-1,nj+1,ttt) = idnum
              END IF
              IF (PCP2(ni-1,nj,t) > crit .AND. ID(ni-1,nj,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni-1; CLD(accum+newcnt,2) = nj
                ID(ni-1,nj,ttt) = idnum
              END IF
              IF (PCP2(ni-1,nj-1,t) > crit .AND. ID(ni-1,nj-1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni-1; CLD(accum+newcnt,2) = nj-1
                ID(ni-1,nj-1,ttt) = idnum
              END IF
              IF (PCP2(ni,nj-1,t) > crit .AND. ID(ni,nj-1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni; CLD(accum+newcnt,2) = nj-1
                ID(ni,nj-1,ttt) = idnum
              END IF
              IF (PCP2(ni+1,nj-1,t) > crit .AND. ID(ni+1,nj-1,ttt) == 0) THEN
                newcnt = newcnt+1
                CLD(accum+newcnt,1) = ni+1; CLD(accum+newcnt,2) = nj-1
                ID(ni+1,nj-1,ttt) = idnum
              END IF
              ! end checking the cross of a given old grid
            END DO
            ! end counting all old grids
            ! if there is no more grids to be counted, exit the while loop
            IF (newcnt == 0) THEN
              EXIT
            END IF
            cnt = newcnt ! all new grids are counted, total counts confirmed
            accum = accum+cnt ! add the new counts to the accumulation of cloud counts
          END DO
          ! end counting all grids
          DO k = 1, accum
            SZE(CLD(k,1),CLD(k,2),ttt) = accum
          END DO
          ci = float(SUM(CLD(:,1), CLD(:,1) > 0))/float(COUNT(CLD(:,1) > 0))-1.0 ! -1 because of the boundary added to data
          cj = float(SUM(CLD(:,2), CLD(:,2) > 0))/float(COUNT(CLD(:,2) > 0))-1.0
          ci2 = int(ci); cj2 = int(cj)
          ci = Lon(ci2)+(ci-ci2)*(Lon(ci2+1)-Lon(ci2))
          cj = Lat(cj2)+(cj-cj2)*(Lat(cj2+1)-Lat(cj2))
          DEALLOCATE(CLD)
          radius = SQRT(real(accum)*(111./4.)**2/3.14)
          WRITE(10, 100) t, int(idnum), radius, ci, cj 
          100 FORMAT(1X,I4,1X,I5,1X,F7.2,1X,F7.2,1X,F6.2)
          total = total+accum
        END IF
        ! end checking if there is a grid of a cloud yet discovered
      END DO
    END DO
  END DO
  DEALLOCATE(PCP2); 
  CALL check(nf90_close(ncid))
  CLOSE(10)
  SZE2=REAL(SZE) ! can NOT transfer type after the variable is defined
  ID2 = ID(2:ii+1,2:jj+1,:) ! nf90_put_var dosen't accept ID(2:ii+1,2:jj+1,:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  !! overwrite this file, if it already exists.
  print*,'here1'
  filename = 'cld_obj_'//trim(year)//trim('.nc') !'CWV_daily.nc'
  CALL check(nf90_create(filename, nf90_clobber, ncid))
  
  ! Define the dimensions. NetCDF will hand back an ID for each.
  CALL check(nf90_def_dim(ncid, 'latitude', jj, j_dimid)) ! not represent
  CALL check(nf90_def_dim(ncid, 'longitude', ii, i_dimid))
  CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid)) ! or tt/intvl (limited dimension)

  dimids = (/ i_dimid, j_dimid, t_dimid /)
  print*, 'here2'
  ! Define the variable.
  CALL check(nf90_def_var(ncid, "longitude", nf90_double,  i_dimid, varid1))
  CALL check(nf90_def_var(ncid, "latitude", nf90_double,  j_dimid, varid2))
  CALL check(nf90_def_var(ncid, "objid", nf90_double,  dimids, varid3)) ! netcdf dosen't accept integer
  CALL check(nf90_def_var(ncid, "objsize", nf90_double,  dimids, varid4))

  CALL check(nf90_enddef(ncid))
  print*, 'here3'
  CALL check(nf90_put_var(ncid, varid1, Lon))
  print*,'here4'
  CALL check(nf90_put_var(ncid, varid2, Lat))
  print*,'here5'
  !start = (/1,1,1/)
  !ccccc = (/ii,jj,5/)
  !DO t = 1, tt/intvl
   ! start(3) = t
  CALL check(nf90_put_var(ncid, varid3, ID2))!, start=start, count=ccccc)) !3/26 18:55 bug: max. 3 times
  !END DO
  CALL check(nf90_put_var(ncid, varid4, SZE2))
  print*,'here6'
  CALL check(nf90_close(ncid))
  PRINT *, 'SUCCESS writing ', filename
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL system_clock ( tt2, clock_rate, clock_max )
  PRINT*, 'Read and write the year ', trim(year), ' successfully! ',  'Elapsed CPU time = ',real(tt2-tt1)/real(clock_rate)
  DEALLOCATE(Lon); DEALLOCATE(Lat); DEALLOCATE(ID); DEALLOCATE(ID2); DEALLOCATE(SZE); DEALLOCATE(SZE2)
END DO
CALL system_clock ( t2, clock_rate, clock_max )
PRINT*, 'All years done! Elapsed CPU time = ',real(t2-t1)/real(clock_rate)


CONTAINS
        SUBROUTINE check(STATUS)
                INTEGER, INTENT( in) :: STATUS
                IF (STATUS /= nf90_noerr) THEN
                        PRINT*, trim(nf90_strerror(STATUS))
                        STOP 2
                END IF
        END SUBROUTINE

END PROGRAM cloud_object
