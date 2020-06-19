PROGRAM  mass_flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/04/16                                   !
! Thisi program aims to                                   !
! 1) calculate mass flux on CWV-pressure coordinate        !
! Inputs:  SPCAM outputs                                  !
! Outputs: mass flux (nc)                                 !
! To do: spatial range, climatology, tracl where the second updraft comes from
!        first mean flux then transform coordinate or
! first transform coordinate then mean the flux?
! huge difference of values between +/- mass flux => samller smaple size for +
! MF? how to solve this?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf
IMPLICIT NONE
INTEGER :: ncid, dimid, varid, varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8
INTEGER :: k_dimid, n_dimid, j_dimid, d_dimid, m_dimid
INTEGER, DIMENSION(3) :: dimids3d
INTEGER, DIMENSION(2) :: dimids2d
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=50) :: path ='/data/dadm1/model_output/SPCAM/CPL64'
CHARACTER(len=50) :: filename
CHARACTER(len=15) :: dimnameii, dimnamejj, dimnamekk, dimnamett
REAL, PARAMETER ::eps=0.622,  es0=611., Lv=2.5*10.**6., Rv=461., T0=273.15, cwv_max = 100., g = 9.8, sth = -30., nth = 30.
CHARACTER(len=2) :: year, mth, day
INTEGER, DIMENSION(2) :: date
INTEGER, PARAMETER :: dd=365, mm=12,  nn=100, init_yr=1, finl_yr=10
INTEGER, DIMENSION(12), PARAMETER :: days=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
INTEGER, PARAMETER :: yy=finl_yr-init_yr+1
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
 
DO n = 1, nn
  CWV(n) = cwv_max/nn*n
END DO

CALL system_clock ( t1, clock_rate, clock_max )
DO y = 1, yy
  WRITE(year, '(I2.2)') init_yr+y-1
  PRINT*, 'year = ', year
  DO d = 1, dd
    CALL day2date(d,date)
    WRITE(mth, '(I2.2)') date(1)
    WRITE(day, '(I2.2)') date(2)
    PRINT*, year,'/',mth,'/',day

    filename = '/CPL64.cam.h0.00'//trim(year)//'-'//trim(mth)//'-'//trim(day)//'-00000.nc'
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
    ALLOCATE(W(ii,jj,kk,tt))
    ALLOCATE(TMP(ii,jj,kk,tt))    
    ALLOCATE(Q(ii,jj,kk,tt))
    ALLOCATE(P(ii,jj,kk,tt))
    ALLOCATE(Qs(ii,jj,kk,tt))
    ALLOCATE(CRH(ii,jj,tt))
    IF (d == 1) THEN
      ALLOCATE(CRH_ddyy(ii,jj,dd))
      ALLOCATE(W_ddyy(ii,jj,kk,dd))
      ALLOCATE(MF_ddyy(nn,kk,dd))!; MF_ddyy = 0
    END IF
    IF (y == 1 .AND. d == 1) THEN
      ALLOCATE(CRH_dd(ii,jj,dd))
      ALLOCATE(CRH_mm(ii,jj,mm))
      ALLOCATE(CRH_clim(ii,jj))
      ALLOCATE(W_dd(ii,jj,kk,dd))
      ALLOCATE(W_mm(ii,jj,kk,mm))
      ALLOCATE(W_clim(ii,jj,kk))
      ALLOCATE(MF_dd(nn,kk,dd)); MF_dd = 0
      ALLOCATE(MF_mm(nn,kk,mm))!; MF_mm = 0
      ALLOCATE(MF_clim(nn,kk))!; MF_clim = 0
    END IF

    CALL check(nf90_inq_varid(ncid, 'lon', varid))
    CALL check(nf90_get_var(ncid, varid, Lon))
    CALL check(nf90_inq_varid(ncid, 'lat', varid))
    CALL check(nf90_get_var(ncid, varid, Lat))
    CALL check(nf90_inq_varid(ncid, 'lev', varid))
    CALL check(nf90_get_var(ncid, varid, Lev))
    CALL check(nf90_inq_varid(ncid, 'OMEGA', varid))
    CALL check(nf90_get_var(ncid, varid, W))
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
      DO k = 1, kk
        DO j = 1, jj
          DO i = 1, ii
            IF (P(i,j,k,t) == 0) THEN
              Qs(i,j,k,t) = Q(i,j,k,t)
            ELSE
              !Qs(i,j,k,t) = eps*es0*EXP(-Lv/Rv*(1/TMP(i,j,k,t)-1/T0))/P(i,j,k,t)
              Qs(i,j,k,t) = eps*es0*EXP(53.49-6808/TMP(i,j,k,t)-5.09*LOG(TMP(i,j,k,t)))/P(i,j,k,t)
            END IF
          END DO
        END DO
      END DO
    END DO
    CRH = SUM(Q, 3)/SUM(Qs, 3)*100.
    CRH_ddyy(:,:,d) = SUM(CRH, 3)/tt
    CALL trans_coor(CRH_ddyy(:,:,d), W_ddyy(:,:,:,d), ii, jj, kk, nn, sth_y, nth_y, cwv_max, MF_ddyy(:,:,d), CNT_ddyy(:,d))
    MF_dd(:,:,d) = MF_dd(:,:,d)+MF_ddyy(:,:,d)/yy
    CNT_dd(:,d) = CNT_dd(:,d)+CNT_ddyy(:,d) ! 14:17
    MF_mm(:,:,date(1)) = MF_mm(:,:,date(1))+MF_ddyy(:,:,d)/(yy*days(date(1)))
    CNT_mm(:,date(1)) = CNT_mm(:,date(1))+CNT_ddyy(:,d)
    MF_clim = MF_clim+MF_ddyy(:,:,d)/(yy*dd)
    CNT_clim = CNT_clim+CNT_ddyy(:,d)
    DEALLOCATE(Lon)
    DEALLOCATE(Lat)
    IF(d /= dd) THEN
      DEALLOCATE(Lev)
    END IF
    DEALLOCATE(W)
    DEALLOCATE(TMP)
    DEALLOCATE(Q)
    DEALLOCATE(P)
    DEALLOCATE(Qs)
    DEALLOCATE(CRH)
  END DO ! ends DO d = 1, 365
  
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

  !print*,'here1'
  !filename = 'MF/mass_flux_'//trim(year)//'.nc'
  !CALL check(nf90_create(filename, nf90_clobber, ncid))
  !CALL check(nf90_def_dim(ncid, 'lev', kk, k_dimid))
  !CALL check(nf90_def_dim(ncid, 'cwv', nn, n_dimid))
  !CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid))

  !print*, 'here2'
  !dimids = (/ n_dimid, k_dimid, t_dimid /)
  !CALL check(nf90_def_var(ncid, "cwv", nf90_double,  n_dimid, varid1))
  !CALL check(nf90_def_var(ncid, "lev", nf90_double,  k_dimid, varid2))
  !CALL check(nf90_def_var(ncid, "mass_flux", nf90_double,  dimids, varid3))
  !CALL check(nf90_def_var(ncid, "freq", nf90_double,  n_dimid, varid4))
  !CALL check(nf90_enddef(ncid))

  !print*, 'here3'
  !CALL check(nf90_put_var(ncid, varid1, CWV))
  !print*, 'here4'
  !CALL check(nf90_put_var(ncid, varid2, Lev))
  !print*, 'here5'
  !CALL check(nf90_put_var(ncid, varid3, MF_ddyy))
  !print*, 'here6'
  !CALL check(nf90_put_var(ncid, varid4, CNT_ddyy))
  !CALL check(nf90_close(ncid))

  !PRINT *, 'SUCCESS writing ', filename
!print*,'1'
!    DEALLOCATE(Lon)
!print*,'2'
!    DEALLOCATE(Lat)
!print*,'3'
!    DEALLOCATE(Lev)
!print*,'4'
!    DEALLOCATE(W)
!print*,'5'
!    DEALLOCATE(CRH)
!print*,'6'
!    DEALLOCATE(CNT)
!print*,'7'
!    DEALLOCATE(MF)
!print*,'b'
  !END DO ! ends DO d = 1, 365
  IF(y /= yy) THEN
    DEALLOCATE(Lev)
  END IF
  DEALLOCATE(CRH_ddyy)
  DEALLOCATE(W_ddyy)
  DEALLOCATE(MF_ddyy)
END DO ! ends DO y = 1, yy
filename = 'MF/mass_flux_climatology.nc'
CALL check(nf90_create(filename, nf90_clobber, ncid))
CALL check(nf90_def_dim(ncid, 'lev', kk, k_dimid))
CALL check(nf90_def_dim(ncid, 'cwv', nn, n_dimid))
CALL check(nf90_def_dim(ncid, 'day', dd, d_dimid))
CALL check(nf90_def_dim(ncid, 'month', mm, m_dimid))
print*, 'here2'
CALL check(nf90_def_var(ncid, "cwv", nf90_double,  n_dimid, varid1))
CALL check(nf90_def_var(ncid, "lev", nf90_double,  k_dimid, varid2))
dimids3d = (/ n_dimid, k_dimid, d_dimid /)
CALL check(nf90_def_var(ncid, "mass_flux_dd", nf90_double, dimids3d, varid3))
dimids2d = (/ n_dimid, d_dimid /)
CALL check(nf90_def_var(ncid, "freq_dd", nf90_double, dimids2d, varid4))
dimids3d = (/ n_dimid, k_dimid, m_dimid /)
CALL check(nf90_def_var(ncid, "mass_flux_mm", nf90_double,  dimids3d, varid5))
dimids2d = (/ n_dimid, m_dimid /)
CALL check(nf90_def_var(ncid, "freq_mm", nf90_double,  dimids2d, varid6))
dimids2d = (/ n_dimid, k_dimid /)
CALL check(nf90_def_var(ncid, "mass_flux_clim", nf90_double,  dimids2d, varid7))
CALL check(nf90_def_var(ncid, "freq_clim", nf90_double,  n_dimid, varid8))
CALL check(nf90_enddef(ncid))

print*, 'here3'
CALL check(nf90_put_var(ncid, varid1, CWV))
print*, 'here4'
CALL check(nf90_put_var(ncid, varid2, Lev))
print*, 'here5'
CALL check(nf90_put_var(ncid, varid3, MF_dd))
CALL check(nf90_put_var(ncid, varid4, CNT_dd))
CALL check(nf90_put_var(ncid, varid5, MF_mm))
print*, 'here6'
CALL check(nf90_put_var(ncid, varid6, CNT_mm))
CALL check(nf90_put_var(ncid, varid7, MF_clim))
CALL check(nf90_put_var(ncid, varid8, CNT_clim))
CALL check(nf90_close(ncid))

PRINT *, 'SUCCESS writing ', filename
DEALLOCATE(CRH_dd)
DEALLOCATE(CRH_mm)
DEALLOCATE(CRH_clim)
DEALLOCATE(W_dd)
DEALLOCATE(W_mm)
DEALLOCATE(W_clim)
DEALLOCATE(MF_dd)
DEALLOCATE(MF_mm)
DEALLOCATE(MF_clim)
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
END PROGRAM mass_flux
