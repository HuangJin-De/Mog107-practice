PROGRAM size_pdf
! Built by: C.-H. Yen, on 03/18/2020

USE netcdf

IMPLICIT NONE
INTEGER :: ncid, dimid, varid, varid1, varid2, varid3
CHARACTER(len=50), PARAMETER :: path = '/data/dadm1/obs/TRMM/'
CHARACTER(len=50) :: filename

INTEGER :: i, ii, j, jj, t, tt, y, yy, nth_y, sth_y, wst_x, est_x, y_dimid, d_dimid, b, bb, m
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=9) :: dimnameii
CHARACTER(len=8) :: dimnamejj
CHARACTER(len=4) :: dimnamett

CHARACTER(len=4) :: year
INTEGER, PARAMETER :: init_yr=1998, finl_yr=2015
REAL, PARAMETER :: nth = 30., sth = -10., wst = 100., est = 160.

REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lon
REAL, DIMENSION(:,:,:), ALLOCATABLE :: VAR, PCP, SZE
REAL, DIMENSION(:), ALLOCATABLE :: Bins ! bin scale of the cloud size distribution
REAL, DIMENSION(:,:), ALLOCATABLE :: DIST ! to record the counts of each bin of cloud size
INTEGER, DIMENSION(12) :: Hrs
INTEGER :: accum ! accumulation of 3 hours per month 


yy = finl_yr-init_yr+1

bb = 10 ! number of bins for the size distribution
ALLOCATE(Bins(bb))
ALLOCATE(DIST(12,bb))
DIST = 0

DO b = 1, bb
        Bins(b) = 10**(6./bb*b)
END DO

! read longitude and latitude first
filename = 'TRMM3B42/3B42.1998.3hr.nc'
CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))

PRINT*, '-----------------------------------------'
PRINT*, 'Dimension: x = ',ii,', y = ',jj
PRINT*, '-----------------------------------------'

ALLOCATE(Lat(jj))
ALLOCATE(Lon(ii))

! read latitude to get the interval of tropics
CALL check(nf90_inq_varid(ncid, 'latitude', varid))
CALL check(nf90_get_var(ncid, varid, Lat))
CALL check(nf90_inq_varid(ncid, 'longitude', varid))
CALL check(nf90_get_var(ncid, varid, Lon))

nth_y = MAXLOC(Lat, 1, Lat <= nth)
sth_y = MAXLOC(Lat, 1, Lat <= sth)
est_x = MAXLOC(Lon, 1, Lon <= est)
wst_x = MAXLOC(Lon, 1, Lon <= wst)
!print*, nth_y, sth_y, est_x, wst_x

! read precipitation rate and cloud size
DO y = 1, yy
        WRITE(year, '(I4)') init_yr+y-1
        filename = trim('TRMM3B42size/TRMMsize_3hrs_')//trim(year)//trim('.nc')
        CALL system_clock ( t1, clock_rate, clock_max )
        CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
        CALL check(nf90_inq_dimid(ncid, 'Time', dimid))
        CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))
        ALLOCATE(VAR(ii, jj, tt)) ! order of dimension must be reversed from the original nc file
        !ALLOCATE(PCP(est_x-wst_x+1, nth_y-sth_y+1, tt))
        ALLOCATE(SZE(est_x-wst_x+1, nth_y-sth_y+1, tt))
        IF (tt == 2920) THEN
                Hrs = 8*(/31,28,31,30,31,30,31,31,30,31,30,31/)
        ELSEIF (tt == 2928) THEN
                Hrs = 8*(/31,29,31,30,31,30,31,31,30,31,30,31/)
        END IF
        ! read values from ncid
        CALL check(nf90_inq_varid(ncid, 'objsize', varid))
        CALL check(nf90_get_var(ncid, varid, VAR))
        CALL system_clock ( t2, clock_rate, clock_max )
        PRINT*, 'Read ', trim(filename), ' successfully! ',  'Elapsed CPU time = ', real(t2-t1)/real(clock_rate)
        SZE = VAR(wst_x:est_x, sth_y:nth_y, :)
        CALL check(nf90_close(ncid))  
        ! 3/18 check the dimension of time of pcp (2928) vs size (2920)
        
        ! count the cloud size
        m = 1
        accum = Hrs(1)
        DO t = 1, tt
        print*, t
                IF (t > accum) THEN
                        m = m+1
                        accum = accum+Hrs(m)
                END IF
                DO b = 1, bb-1
                        DIST(m,b) = DIST(m,b)+COUNT(SZE >= Bins(b) .AND. SZE <= Bins(b+1))       
                END DO
        END DO
        print*, DIST            
        !DEALLOCATE(PCP)
        DEALLOCATE(SZE)
        DEALLOCATE(VAR)

END DO

CONTAINS
        SUBROUTINE check(STATUS)
                INTEGER, INTENT( in) :: STATUS
                IF (STATUS /= nf90_noerr) THEN
                        PRINT*, trim(nf90_strerror(STATUS))
                        STOP 2
                END IF
        END SUBROUTINE
END PROGRAM size_pdf
