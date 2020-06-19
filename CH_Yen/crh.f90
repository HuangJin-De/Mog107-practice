PROGRAM crh_olr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/04/23                                   !
! Thisi program aims to                                   !
! 1) calculate column-integrated RH                       !
! Inputs:  SPCAM outputs                                  !
! Outputs: CRH (nc)                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf
IMPLICIT NONE
INTEGER :: ncid, dimid, varid, varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8
INTEGER :: i_dimid, j_dimid, t_dimid, d_dimid, m_dimid
INTEGER, DIMENSION(3) :: dimids
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=50) :: path ='/data/dadm1/model_output/SPCAM/CPL64'
CHARACTER(len=50) :: filename
CHARACTER(len=15) :: dimnameii, dimnamejj, dimnamekk, dimnamett
REAL, PARAMETER ::eps=0.622,  es0=611., Lv=2.5*10.**6., Rv=461., T0=273.15, cwv_max = 100., g = 9.8, sth = -30., nth = 30.
CHARACTER(len=2) :: year, mth, day
INTEGER, DIMENSION(2) :: date
INTEGER, PARAMETER :: dd=365, mm=12,  nn=100
INTEGER, DIMENSION(12), PARAMETER :: days=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
INTEGER :: yy, init_yr, finl_yr
INTEGER :: sth_y, nth_y

INTEGER :: i, ii, iii, j, jj, jjj, k, kk, t, tt, d, y, n
REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lon, Lev

REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: TMP, Q, P, Qs, W, W_ddyy, W_dd, W_mm
REAL, DIMENSION(:,:,:), ALLOCATABLE :: W_clim

REAL, DIMENSION(:,:,:), ALLOCATABLE :: CRH, CRH_ddyy, CRH_dd, CRH_mm
REAL, DIMENSION(:,:), ALLOCATABLE :: CRH_clim

REAL, DIMENSION(nn) :: CWV

REAL, DIMENSION(nn,dd) :: CNT_ddyy
REAL, DIMENSION(nn,dd) :: CNT_dd
REAL, DIMENSION(nn,mm) :: CNT_mm
REAL, DIMENSION(nn) :: CNT_clim

REAL, DIMENSION(:,:,:), ALLOCATABLE :: MF_ddyy, MF_dd, MF_mm ! mass flux
REAL, DIMENSION(:,:), ALLOCATABLE :: MF_clim
 

CALL system_clock ( t1, clock_rate, clock_max )

PRINT*, 'init_yr:'
READ(*,'(I2)') init_yr
PRINT*, 'finl_yr:'
READ(*, '(I2)') finl_yr

yy = finl_yr-init_yr+1
DO y = 1, yy
  WRITE(year, '(I2.2)') init_yr+y-1
  PRINT*, 'year = ', year
  DO d = 1, dd
    CALL day2date(d,date)
    WRITE(mth, '(I2.2)') date(1)
    WRITE(day, '(I2.2)') date(2)
    !PRINT*, year,'/',mth,'/',day

    filename = '/CPL64.cam.h0.00'//year//'-'//mth//'-'//day//'-00000.nc'
    CALL check(nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid))
    CALL check(nf90_inq_dimid(ncid, 'lon', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
    CALL check(nf90_inq_dimid(ncid, 'lat', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
    CALL check(nf90_inq_dimid(ncid, 'lev', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamekk, kk))
    CALL check(nf90_inq_dimid(ncid, 'time', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))

    ALLOCATE(Lon(ii))
    ALLOCATE(Lat(jj))
    ALLOCATE(Lev(kk))
    ALLOCATE(TMP(ii,jj,kk,tt))    
    ALLOCATE(Q(ii,jj,kk,tt))
    ALLOCATE(P(ii,jj,kk,tt))
    ALLOCATE(Qs(ii,jj,kk,tt))
    ALLOCATE(CRH(ii,jj,tt))

    CALL check(nf90_inq_varid(ncid, 'lon', varid))
    CALL check(nf90_get_var(ncid, varid, Lon))
    CALL check(nf90_inq_varid(ncid, 'lat', varid))
    CALL check(nf90_get_var(ncid, varid, Lat))
    CALL check(nf90_inq_varid(ncid, 'lev', varid))
    CALL check(nf90_get_var(ncid, varid, Lev))
    CALL check(nf90_inq_varid(ncid, 'T', varid))
    CALL check(nf90_get_var(ncid, varid, TMP))
    CALL check(nf90_inq_varid(ncid, 'Q', varid))
    CALL check(nf90_get_var(ncid, varid, Q))
    CALL check(nf90_inq_varid(ncid, 'PRES', varid))
    CALL check(nf90_get_var(ncid, varid, P))
    CALL check(nf90_close(ncid))

    nth_y = MAXLOC(Lat, 1, Lat <= nth)
    sth_y = MAXLOC(Lat, 1, Lat <= sth)
      
    W_ddyy(:,:,:,d) = SUM(W, 4)/tt

    ! calculate CRH
    DO t = 1, tt
      DO k = 1, kk-1
        DO j = 1, jj
          DO i = 1, ii
            IF (P(i,j,k,t) == 0) THEN
              Qs(i,j,k,t) = Q(i,j,k,t)
            ELSE
              !Qs(i,j,k,t) = eps*es0*EXP(-Lv/Rv*(1/TMP(i,j,k,t)-1/T0))/P(i,j,k,t)
              Qs(i,j,k,t) = eps*es0*EXP(53.49-6808/TMP(i,j,k,t)-5.09*LOG(TMP(i,j,k,t)))/P(i,j,k,t) &
                          & *(P(i,j,k+1,t)-P(i,j,k,t))/g
              Q(i,j,k,t) = Q(i,j,k,t)*(P(i,j,k+1,t)-P(i,j,k,t))/g
            END IF
          END DO ! ends DO i = 1, ii
        END DO ! ends DO j = 1, jj
      END DO ! ends DO k = 1, kk-1
      Qs(:,:,kk,t) = Qs(:,:,k,t)
      Q(:,:,kk,t) = Q(:,:,k,t)
    END DO ! ends DO t = 1, tt
    ! no need to Q/(Ps/g) and Qs/(Ps/g) because Ps/g cancealled each out.
    CRH = SUM(Q, 3)/SUM(Qs, 3)*100.
    DEALLOCATE(Lev)
    DEALLOCATE(TMP)
    DEALLOCATE(Q)
    DEALLOCATE(P)
    DEALLOCATE(Qs)
  
    !DO j = sth_y, nth_y
    !  DO i = 1, ii
    !    n = int(CRH_ddyy(i,j,d)/(cwv_max/nn))
    !    IF (n == 0) THEN
    !      n = 1
    !    ELSE IF (n > nn) THEN
    !      n = nn
    !    END IF
    !    MF(n,:,t) = MF(n,:,t)+(-1/g)*W_ddyy(i,j,:,t) ! mass flux (rho*w = -omega/g) [kgm^(-2)s^(-1)]
    !    CNT(n,t)  = CNT(n,t)+1
    !  END DO ! ends DO i = 1, ii
    !END DO ! ends DO j = 1, jj
    !DO n = 1, nn
    !  IF (CNT(n,t) /= 0. ) THEN
    !    MF(n,:,t) = MF(n,:,t)/CNT(n,t)
    !  END IF
    !  CNT(n,t) = CNT(n,t)/(ii*jj)
    !END DO ! ends DO n = 1, nn

    filename = 'CRH/CRH_'//year//'-'//mth//'-'//day//'.nc'
    CALL check(nf90_create(filename, nf90_clobber, ncid))
    CALL check(nf90_def_dim(ncid, 'lon', ii, i_dimid))
    CALL check(nf90_def_dim(ncid, 'lat', jj, j_dimid))
    CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid))

    dimids = (/ i_dimid, j_dimid, t_dimid /)
    CALL check(nf90_def_var(ncid, "lon", nf90_double,  i_dimid, varid1))
    CALL check(nf90_def_var(ncid, "lat", nf90_double,  j_dimid, varid2))
    CALL check(nf90_def_var(ncid, "cwv", nf90_double,  dimids, varid3))
    CALL check(nf90_enddef(ncid))

    CALL check(nf90_put_var(ncid, varid1, Lon))
    CALL check(nf90_put_var(ncid, varid2, Lat))
    CALL check(nf90_put_var(ncid, varid3, CRH))
    CALL check(nf90_close(ncid))

    PRINT *, 'SUCCESS writing ', filename
    DEALLOCATE(Lon)
    DEALLOCATE(Lat)
    DEALLOCATE(CRH)
  END DO ! ends DO d = 1, 365
END DO ! ends DO y = 1, yy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS
SUBROUTINE day2date(day, date)
IMPLICIT NONE
INTEGER, INTENT(in) :: day
INTEGER, DIMENSION(2), INTENT(out) :: date
INTEGER :: i=1
INTEGER, DIMENSION(12), PARAMETER :: days=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

date(1) = i; date(2) = day ! not right
DO WHILE (date(2)-days(i) > 0)
  date(2) = date(2)-days(i)
  date(1) = date(1)+1
  i = i+1
END DO
i=1 ! reset to 1 for the the next day

RETURN
END SUBROUTINE day2date
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE trans_coor(CRH, W, ii, jj, kk, nn, sth_y, nth_y, cwv_max, MF, CNT)
IMPLICIT NONE
INTEGER, INTENT(in) :: ii, jj, kk, nn, sth_y, nth_y
REAL, INTENT(in) :: cwv_max
REAL, DIMENSION(ii,jj), INTENT(in) :: CRH
REAL, DIMENSION(ii,jj,kk), INTENT(in) :: W
REAL, DIMENSION(nn,kk), INTENT(out) :: MF
REAL, DIMENSION(nn), INTENT(out) :: CNT
REAL, PARAMETER :: g=9.8
INTEGER :: i, j, n

MF = 0; CNT = 0
DO j = sth_y, nth_y
   DO i = 1, ii
     n = int(CRH(i,j)/(cwv_max/nn))
     IF (n == 0) THEN
       n = 1
     ELSE IF (n > nn) THEN
       n = nn
     END IF
     MF(n,:) = MF(n,:)+(-1/g)*W(i,j,:) ! mass flux (rho*w=-omega/g) [kgm^(-2)s^(-1)]
     CNT(n)  = CNT(n)+1
   END DO ! ends DO i = 1, ii
END DO ! ends DO j = 1, jj
DO n = 1, nn
  IF (CNT(n) /= 0. ) THEN
    MF(n,:) = MF(n,:)/CNT(n)
  END IF
  CNT(n) = CNT(n)
END DO ! ends DO n = 1, nn

RETURN
END SUBROUTINE trans_coor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sort(Tosort, Sorted, Indices, n)
IMPLICIT NONE
INTEGER , INTENT(in) :: n ! size of data array
REAL, DIMENSION(n), INTENT(in) :: Tosort
REAL, DIMENSION(n), INTENT(out) :: Sorted
INTEGER, DIMENSION(n), INTENT(inout) :: Indices
REAL :: minimum ! temporary variable for swapping
INTEGER :: i,j,k ! counter for do loops

!Sort the data
Sorted = Tosort
DO i = 1, n-1
  k = i
  minimum = Sorted(i)
  ! find the minimum value in sorted(i) to sorted(n)
  ! and store it temporarily in minimum
  DO j = i+1, n
    IF (Sorted(j) < minimum) THEN
      k = j
      minimum = Sorted(k)
    END IF 
  END DO
  ! swap the minimum value with sorted(i)
  Sorted(k) = Sorted(i)
  Sorted(i) = minimum
  Indices(k) = Indices(i)
  Indices(i) = k
ENDDO
RETURN
END SUBROUTINE sort
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE check(STATUS)
                INTEGER, INTENT( in) :: STATUS
                IF (STATUS /= nf90_noerr) THEN
                        PRINT*, trim(nf90_strerror(STATUS))
                        STOP 2
                END IF
        END SUBROUTINE
END PROGRAM crh_olr
