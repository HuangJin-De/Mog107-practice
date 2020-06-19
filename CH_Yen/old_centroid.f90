PROGRAM centroid
USE netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/03/24      
! Thisi program aims to calculate the centroids of clouds.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
INTEGER :: ncid, dimid, varid, t_dimid, i_dimid, j_dimid
INTEGER, DIMENSION(2) :: dimids
INTEGER :: i, ii, j, jj, tt, k, ni, nj, ci = 0, cj = 0, accum, cnt, newcnt, idnum = 0
                       ! now_index, centroid_index, cnt: counts of grids within
                       ! a cloud, idnum: idnumber of IDs of clouds
REAL, PARAMETER :: crit = 1.
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PCPP
REAL, DIMENSION(:,:), ALLOCATABLE :: PCP
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ID 
INTEGER, DIMENSION(:,:), ALLOCATABLE :: CLD 

CHARACTER(len=50), PARAMETER :: path = '/data/dadm1/obs/TRMM/'
CHARACTER(len=50) :: filename

INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=9) :: dimnameii
CHARACTER(len=8) :: dimnamejj
CHARACTER(len=4) :: dimnamett

CHARACTER(len=4) :: year
INTEGER, PARAMETER :: init_yr=1998, finl_yr=2015

! read longitude and latitude first
filename = 'TRMM3B42/3B42.1998.3hr.nc'
CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
CALL check(nf90_inq_dimid(ncid, 'time', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))

PRINT*, '-----------------------------------------'
PRINT*, 'Dimension: x = ',ii,', y = ',jj
PRINT*, '-----------------------------------------'

ALLOCATE(PCPP(ii,jj,tt))
ALLOCATE(ID(ii+2,jj+2))
ID = 0

CALL check(nf90_inq_varid(ncid, 'pcp', varid))
CALL check(nf90_get_var(ncid, varid, PCPP))

ALLOCATE(PCP(ii+2,jj+2)); PCP = 0
PCP(2:ii,2:jj) = PCPP(:,:,1)

PRINT*, 'START SEARCHING CLOUD!!'

DO j = 2, jj
  DO i = 2, ii
    ! if finds a grid with pcp and yet to be counted
    ! => activate the counting of new cloud
    IF ((PCP(i,j) > crit) .AND. ID(i,j) == 0) THEN
      ! store the information of the first grid of the new found cloud
      idnum = idnum+1
      PRINT*, 'Cloud #', idnum, ' detected! (', i, j, ')'
      ni = i; nj = j
      ID(ni,nj) = idnum
      cnt = 1; accum = 1
      ALLOCATE(CLD(ii*jj,2)); CLD = 0
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
          IF (PCP(ni+1,nj) > crit .AND. ID(ni+1,nj) == 0) THEN
            newcnt = newcnt+1
            CLD(accum+newcnt,1) = ni+1; CLD(accum+newcnt,2) = nj
            ID(ni+1,nj) = idnum
          END IF
          IF (PCP(ni,nj+1) > crit .AND. ID(ni,nj+1) == 0) THEN
            newcnt = newcnt+1
            CLD(accum+newcnt,1) = ni; CLD(accum+newcnt,2) = nj+1
            ID(ni,nj+1) = idnum
          END IF
          IF (PCP(ni-1,nj) > crit .AND. ID(ni-1,nj) == 0) THEN
            newcnt = newcnt+1
            CLD(accum+newcnt,1) = ni-1; CLD(accum+newcnt,2) = nj
            ID(ni-1,nj) = idnum
          END IF
          IF (PCP(ni,nj-1) > crit .AND. ID(ni,nj-1) == 0) THEN
            newcnt = newcnt+1
            CLD(accum+newcnt,1) = ni; CLD(accum+newcnt,2) = nj-1
            ID(ni,nj-1) = idnum
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
      PRINT*, '3'
      ! end counting all grids
      ci = float(SUM(CLD(:,1), CLD(:,1) > 0))/float(COUNT(CLD(:,1) > 0))
      cj = float(SUM(CLD(:,2), CLD(:,2) > 0))/float(COUNT(CLD(:,2) > 0))
      DEALLOCATE(CLD)
      PRINT*, 'Cloud ', idnum, ' finished!'
    END IF
    ! end checking if there is a grid of a cloud yet discovered
  END DO
END DO
CALL check(nf90_close(ncid))

! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
filename = 'test.nc' !'CWV_daily.nc'
CALL check(nf90_create(filename, nf90_clobber, ncid))
 
! Define the dimensions. NetCDF will hand back an ID for each.
!CALL check(nf90_def_dim(ncid, 'time', tt, t_dimid))
CALL check(nf90_def_dim(ncid, 'longitude', ii+2, i_dimid))
CALL check(nf90_def_dim(ncid, 'latitude', jj+2, j_dimid)) ! not represent z-direction

! The dimids array is used to pass the IDs of the dimensions of
! the variables. Note that in fortran arrays are stored in
! column-major format.
dimids =  (/ i_dimid, j_dimid /)

! Define the variable.
CALL check(nf90_def_var(ncid, "cloudID", nf90_double,  dimids, varid))

! End define mode. This tells netCDF we are done defining metadata.
CALL check(nf90_enddef(ncid))

! Write the pretend data to the file. Although netCDF supports
CALL check(nf90_put_var(ncid, varid, ID))

CALL check(nf90_close(ncid))
PRINT *, 'SUCCESS writing ', filename

CONTAINS
        SUBROUTINE check(STATUS)
                INTEGER, INTENT( in) :: STATUS
                IF (STATUS /= nf90_noerr) THEN
                        PRINT*, trim(nf90_strerror(STATUS))
                        STOP 2
                END IF
        END SUBROUTINE

END PROGRAM centroid
