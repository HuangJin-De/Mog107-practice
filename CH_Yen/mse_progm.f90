PROGRAM mse_progm
USE netcdf

IMPLICIT NONE

INTEGER :: ncid, dimid, varid, varid1, varid2, varid3
CHARACTER(len=50) :: path = '/data/dadm1/reanalysis/ECMWF/ITM/daily/'
CHARACTER(len=30) :: filename

INTEGER :: i, ii, j, jj, k, kk, t, tt, y, yy, nth_y, sth_y, y_dimid, d_dimid, z_dimid
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=9) :: dimnameii
CHARACTER(len=8) :: dimnamejj
CHARACTER(len=5) :: dimnamekk
CHARACTER(len=4) :: dimnamett

REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: VAR, TMP, GEO, Qv, TRASH
REAL, DIMENSION(:,:,:), ALLOCATABLE :: MSE, HGT
REAL, DIMENSION(:), ALLOCATABLE :: LAT, Lev
INTEGER, DIMENSION(3) :: dimids
REAL, PARAMETER :: g = 9.8, Cp = 1004., Lv = 2.5*10**6
REAL, PARAMETER :: nth = 30., sth = -30.

CHARACTER(len=4) :: year
INTEGER, PARAMETER :: init_yr=1979, finl_yr=2017
 
yy = finl_yr-init_yr+1
! open ncid
filename = 'T/daily_interim_T_1979.nc'
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

ALLOCATE(LAT(jj))
ALLOCATE(Lev(kk))
ALLOCATE(MSE(kk, 366, yy))
ALLOCATE(HGT(kk, 366, yy))

! read latitude to get the interval of tropics
CALL check(nf90_inq_varid(ncid, 'latitude', varid))
CALL check(nf90_get_var(ncid, varid, LAT))
CALL check(nf90_inq_varid(ncid, 'level', varid))
CALL check(nf90_get_var(ncid, varid, Lev))

nth_y = MAXLOC(LAT, 1, LAT <= nth)
sth_y = MAXLOC(LAT, 1, LAT <= sth)

! import values of variables from ncid
DO y = 1, finl_yr-init_yr+1
        WRITE(year, '(I4)') 1979+y-1
        filename = trim('T/daily_interim_T_')//trim(year)//trim('.nc')
        CALL system_clock ( t1, clock_rate, clock_max )
        CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
        CALL check(nf90_inq_dimid(ncid, 'time', dimid))
        CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))

        ! allocate memory to the variables
        ALLOCATE(VAR(ii,jj,kk,tt))
        ALLOCATE(TMP(ii,sth_y-nth_y+1,kk,tt))
        ALLOCATE(GEO(ii,sth_y-nth_y+1,kk,tt))
        ALLOCATE(Qv(ii,sth_y-nth_y+1,kk,tt))

        ! import values of temperature from ncid
        CALL check(nf90_inq_varid(ncid, 't', varid))
        CALL check(nf90_get_var(ncid, varid, VAR))
        CALL system_clock ( t2, clock_rate, clock_max )
        PRINT*, 'Read ', trim(filename), ' successfully! ',  'Elapsed CPU time = ', real(t2-t1)/real(clock_rate)
        TMP = VAR(:,nth_y:sth_y,:,:)
        CALL check(nf90_close(ncid))

        ! import values of geopotential from ncid
        filename = trim('Z/daily_interim_Z_')//trim(year)//trim('.nc')
        CALL system_clock ( t1, clock_rate, clock_max )
        CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
        CALL check(nf90_inq_varid(ncid, 'z', varid))
        CALL check(nf90_get_var(ncid, varid, VAR))
        CALL system_clock ( t2, clock_rate, clock_max )
        PRINT*, 'Read ', trim(filename), ' successfully! ',  'Elapsed CPU time = ', real(t2-t1)/real(clock_rate)
        GEO = VAR(:,nth_y:sth_y,:,:)
        CALL check(nf90_close(ncid))
        
        ! import values of water vapor from ncid
        filename = trim('Q/daily_interim_Q_')//trim(year)//trim('.nc')
        CALL system_clock ( t1, clock_rate, clock_max )
        CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
        CALL check(nf90_inq_varid(ncid, 'q', varid))
        CALL check(nf90_get_var(ncid, varid, VAR))
        CALL system_clock ( t2, clock_rate, clock_max )
        PRINT*, 'Read ', trim(filename), ' successfully! ',  'Elapsed CPU time = ', real(t2-t1)/real(clock_rate)
        Qv = VAR(:,nth_y:sth_y,:,:)
        CALL check(nf90_close(ncid))

        DO t = 1, tt
                DO k = 1, kk
                        MSE(k,t,y) = Cp*SUM(TMP(:,:,k,t))/SIZE(TMP(:,:,k,t)) &
                                   + SUM(GEO(:,:,k,t))/SIZE(GEO(:,:,k,t)) &
                                   + Lv*SUM(Qv(:,:,k,t))/SIZE(Qv(:,:,k,t))
                        HGT(k,t,y) = SUM(GEO(:,:,k,t))/SIZE(GEO(:,:,k,t))/g
                END DO
        END DO
        DEALLOCATE(VAR)
        DEALLOCATE(TMP)
        DEALLOCATE(GEO)
        DEALLOCATE(Qv)
END DO


! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
filename = 'mse_daily.nc'
CALL check(nf90_create(filename, nf90_clobber, ncid))

! Define the dimensions. NetCDF will hand back an ID for each.
CALL check(nf90_def_dim(ncid, 'year', finl_yr-init_yr+1, y_dimid))
CALL check(nf90_def_dim(ncid, 'day', 366, d_dimid))
CALL check(nf90_def_dim(ncid, 'level', kk, z_dimid))

! The dimids array is used to pass the IDs of the dimensions of
! the variables. Note that in fortran arrays are stored in
! column-major format.
dimids =  (/ z_dimid, d_dimid, y_dimid /)

! Define the variable.
CALL check(nf90_def_var(ncid, "mse", nf90_int,  dimids, varid1)) ! use tpye integer to store MSE to save memory
CALL check(nf90_def_var(ncid, 'height', nf90_float, dimids, varid2))
CALL check(nf90_def_var(ncid, "level", nf90_float, z_dimid, varid3))
! End define mode. This tells netCDF we are done defining metadata.
CALL check(nf90_enddef(ncid))

! Write the pretend data to the file. Although netCDF supports
CALL check(nf90_put_var(ncid, varid1, MSE))
CALL check(nf90_put_var(ncid, varid2, HGT))
CALL check(nf90_put_var(ncid, varid3, Lev))
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
END PROGRAM mse_progm
