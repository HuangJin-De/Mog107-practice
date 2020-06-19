PROGRAM cwv_progm
USE netcdf

IMPLICIT NONE

INTEGER :: ncid, dimid, varid
CHARACTER(len=50) :: path = '/data/dadm1/reanalysis/ECMWF/ITM/daily/'
CHARACTER(len=30) :: filename

INTEGER :: i, ii, j, jj, k, kk, d, dd, y, yy, nth_y, sth_y, y_dimid, d_dimid, z_dimid, i_dimid, j_dimid
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=9) :: dimnameii
CHARACTER(len=8) :: dimnamejj
CHARACTER(len=5) :: dimnamekk
CHARACTER(len=4) :: dimnamett

REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: VAR, Qv
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: CWV
REAL, DIMENSION(:,:,:), ALLOCATABLE :: DRY ! areal fraction of dry region with CWV lower than a criterion
REAL, DIMENSION(4), PARAMETER :: crtn = (/10., 20., 30., 40./) ! various criteria
REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lev
INTEGER, DIMENSION(3) :: dimids
REAL, PARAMETER :: g = 9.8, rho_w = 1000 ! density of liquid water (kg/m^3)
REAL, PARAMETER :: nth = 30., sth = -30.

CHARACTER(len=4) :: year
INTEGER, PARAMETER :: init_yr=1979, finl_yr=2017

yy = finl_yr-init_yr+1

! open ncid
filename = 'Q/daily_interim_Q_1979.nc'
CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))

! get the IDs of dimensions from ncid
! get the lengths of dimensions
CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
CALL check(nf90_inq_dimid(ncid, 'level', dimid))
CALL check(nf90_inquire_dimension(ncid, dimid, dimnamekk, kk))

PRINT*, '-----------------------------------------'
PRINT*, 'Dimension: x = ',ii,', y = ',jj,', z = ',kk
PRINT*, '-----------------------------------------'


ALLOCATE(Lat(jj))
ALLOCATE(Lev(kk))
ALLOCATE(DRY(4, 366, yy))

! read latitude to get the interval of tropics
CALL check(nf90_inq_varid(ncid, 'latitude', varid))
CALL check(nf90_get_var(ncid, varid, Lat))
CALL check(nf90_inq_varid(ncid, 'level', varid))
CALL check(nf90_get_var(ncid, varid, Lev))
CALL check(nf90_close(ncid))

nth_y = MAXLOC(Lat, 1, Lat <= nth)
sth_y = MAXLOC(Lat, 1, Lat <= sth)
ALLOCATE(CWV(ii,sth_y-nth_y+1,366, yy))

! import values of variables from ncid
DO y = 1, yy
        WRITE(year, '(I4)') init_yr+y-1
        filename = trim('Q/daily_interim_Q_')//trim(year)//trim('.nc')
        CALL system_clock ( t1, clock_rate, clock_max )
        CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
        CALL check(nf90_inq_dimid(ncid, 'time', dimid))
        CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, dd))

        ! allocate memory to the variables
        ALLOCATE(VAR(ii,jj,kk,dd))
        ALLOCATE(Qv(ii,sth_y-nth_y+1,kk,dd))

        CALL check(nf90_inq_varid(ncid, 'q', varid))
        CALL check(nf90_get_var(ncid, varid, VAR))
        CALL check(nf90_close(ncid))
        CALL system_clock ( t2, clock_rate, clock_max )
        PRINT*, 'Read ', trim(filename), ' successfully! ',  'Elapsed CPU time = ', real(t2-t1)/real(clock_rate)
        Qv = VAR(:,nth_y:sth_y,:,:)
        DEALLOCATE(VAR)
        ! calculate Column-Integrated_water_Vapor
        DO d = 1, dd
                DO k = 1, kk-1
                        CWV(:,:,d,y) = CWV(:,:,d,y)+Qv(:,:,k,d)*(Lev(k)-Lev(k+1))*100 ! hPa -> Pa
                END DO
                CWV(:,:,d,y) = CWV(:,:,d,y)/(g*rho_w)*1000
                DO i = 1, 4
                        DRY(i,d,y) = REAL(COUNT(CWV(:,:,d,y) < crtn(i)))/REAL(SIZE(CWV(:,:,d,y)))
                END DO
        END DO
        DEALLOCATE(Qv)
END DO


! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
filename = 'dry_daily.nc' !'CWV_daily.nc'
CALL check(nf90_create(filename, nf90_clobber, ncid))

! Define the dimensions. NetCDF will hand back an ID for each.
CALL check(nf90_def_dim(ncid, 'year', yy, y_dimid))
CALL check(nf90_def_dim(ncid, 'day', 366, d_dimid))
CALL check(nf90_def_dim(ncid, 'crtrn', 4, z_dimid)) ! not represent z-direction

! The dimids array is used to pass the IDs of the dimensions of
! the variables. Note that in fortran arrays are stored in
! column-major format.
dimids =  (/ z_dimid, d_dimid, y_dimid /)

! Define the variable.
CALL check(nf90_def_var(ncid, "dry", nf90_double,  dimids, varid))

! End define mode. This tells netCDF we are done defining metadata.
CALL check(nf90_enddef(ncid))

! Write the pretend data to the file. Although netCDF supports
CALL check(nf90_put_var(ncid, varid, DRY))
 
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
END PROGRAM cwv_progm
