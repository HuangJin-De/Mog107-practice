PROGRAM  crh_olr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/04/16                                   !
! Thisi program aims to                                   !
! 1) plot the frequency map of CRH vs. OLR                !
! Inputs:  SPCAM outputs, CRH nc files                    !
! Outputs: crh_olr_frequency.nc                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf
IMPLICIT NONE
INTEGER :: ncid, dimid, varid 
INTEGER :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, varid9, varid10, varid11, varid12, varid13, varid14
INTEGER :: i_dimid, n_dimid, j_dimid, d_dimid, m_dimid
INTEGER, DIMENSION(2) :: dimids1, dimids2
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=50) :: path ='/data/dadm1/model_output/SPCAM/CPL64'
CHARACTER(len=50) :: filename
CHARACTER(len=15) :: dimnameii, dimnamejj, dimnamekk, dimnamett
REAL, PARAMETER :: olr_min=-1, sth=-30., nth=30., crh_bin=1., olr_bin=0.01
CHARACTER(len=2) :: year, mth, day
INTEGER, DIMENSION(2) :: date
INTEGER, PARAMETER :: dd=365, nn=100, mm=175, init_yr=1, finl_yr=1
INTEGER, DIMENSION(12), PARAMETER :: days=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
INTEGER, PARAMETER :: yy=finl_yr-init_yr+1
INTEGER :: sth_y, nth_y
INTEGER :: i, ii, iii, j, jj, jjj, k=13, kk, t, tt, d, y, n, m, cnt=0, cnt_a=0, cnt_b=0, cnt_c=0, cnt_d=0
REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lon, Lev
REAL, DIMENSION(:,:,:), ALLOCATABLE :: CRH
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: OLR
REAL, DIMENSION(:,:,:), ALLOCATABLE :: W
REAL, DIMENSION(nn) :: Bin_crh, PDF_crh=0, pdf_crh_olr=0
REAL, DIMENSION(mm) :: Bin_olr, PDF_olr=0, pdf_w_olr=0
REAL, DIMENSION(nn,mm) :: FRQ=0., W_2pdf=0.
REAL, DIMENSION(nn) :: CNT_clim
REAL, DIMENSION(:,:), ALLOCATABLE :: MAP_A, MAP_B, MAP_C, MAP_D

DO n = 1, nn
  Bin_crh(n) = crh_bin*(n-1)
END DO
DO m = 1, mm
  Bin_olr(m) = olr_min+olr_bin*(m-1)
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
    ALLOCATE(OLR(ii,jj,kk,tt))
    ALLOCATE(W(ii,jj,tt))

    CALL check(nf90_inq_varid(ncid, 'lon', varid))
    CALL check(nf90_get_var(ncid, varid, Lon))
    CALL check(nf90_inq_varid(ncid, 'lat', varid))
    CALL check(nf90_get_var(ncid, varid, Lat))
    CALL check(nf90_inq_varid(ncid, 'OMEGA', varid))
    CALL check(nf90_get_var(ncid, varid, OLR))
    CALL check(nf90_inq_varid(ncid, 'FLNT', varid))
    CALL check(nf90_get_var(ncid, varid, W))
    CALL check(nf90_close(ncid))

    nth_y = MAXLOC(Lat, 1, Lat <= nth)
    sth_y = MAXLOC(Lat, 1, Lat <= sth)
      
    filename = '../hw06/CRH/CRH_'//year//'-'//mth//'-'//day//'.nc'
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    ALLOCATE(CRH(ii,jj,tt))
    CALL check(nf90_inq_varid(ncid, 'cwv', varid))
    CALL check(nf90_get_var(ncid, varid, CRH))
    CALL check(nf90_close(ncid))

    IF (y == 1 .AND. d == 1) THEN
      ALLOCATE(MAP_A(ii,jj)); MAP_A = 0.
      ALLOCATE(MAP_B(ii,jj)); MAP_B = 0.
      ALLOCATE(MAP_C(ii,jj)); MAP_C = 0.
      ALLOCATE(MAP_D(ii,jj)); MAP_D = 0.
    END IF
    ! calculate the CRH_OLR frequency matrix
    DO t = 1, tt
      DO j = 1, jj
        DO i = 1, ii
          IF (j >= sth_y .AND. j <= nth_y) THEN
            m = INT((OLR(i,j,k,t)-olr_min)/olr_bin)+1
            n = INT(CRH(i,j,t)/crh_bin)+1
            IF ((m >= 1 .AND. m <= mm) .AND. (n <= nn .AND. n >= 1)) THEN
              FRQ(n,m) = FRQ(n,m)+1.
              W_2pdf(n,m) = W_2pdf(n,m)+W(i,j,t)
              pdf_crh_olr(n) = pdf_crh_olr(n)+W(i,j,t)
              pdf_w_olr(m) = pdf_w_olr(m)+W(i,j,t)
              PDF_crh(n) = PDF_crh(n)+1.
              PDF_olr(m) = PDF_olr(m)+1.
              cnt = cnt+1
              IF (Bin_olr(m) > 0. .AND. Bin_crh(n) < 60.) THEN
                MAP_A(i,j) = MAP_A(i,j)+1.
                cnt_a = cnt_a+1
              ELSE IF (Bin_olr(m) > 0. .AND. Bin_crh(n) > 60.) THEN
                MAP_B(i,j) = MAP_B(i,j)+1.
                cnt_b = cnt_b+1
              ELSE IF (Bin_olr(m) < 0. .AND. Bin_crh(n) > 60.) THEN
                MAP_C(i,j) = MAP_C(i,j)+1.
                cnt_c = cnt_c+1
              ELSE IF (Bin_olr(m) < 0. .AND. Bin_crh(n) < 60.) THEN
                MAP_D(i,j) = MAP_D(i,j)+1.
                cnt_d = cnt_d+1
              END IF
            ELSE
              PRINT*, 'n,m = ', n, m, '. one exceeds number of bins.'
            END IF

          END IF
        END DO ! ends DO i = 1, ii
      END DO ! ends DO j = 1, jj
    END DO ! ends DO t = 1, tt
    IF (y /= yy .OR. d /= dd) THEN
      DEALLOCATE(Lon)
      DEALLOCATE(Lat)
    END IF
    DEALLOCATE(OLR)
    DEALLOCATE(CRH)
    DEALLOCATE(W)
  END DO ! ends DO d = 1, 365  
END DO ! ends DO y = 1, yy
print*,'A'
DO m = 1,mm
  DO n = 1,nn
    IF (FRQ(n,m) > 0) THEN
      W_2pdf(n,m) = W_2pdf(n,m)/FRQ(n,m)
    END IF
  END DO
END DO
DO n = 1, nn
  IF(PDF_crh(n) > 0) THEN
    pdf_crh_olr(n) = pdf_crh_olr(n)/PDF_crh(n)
  END IF
END DO
DO m = 1, mm
  IF(PDF_crh(m) > 0) THEN
    pdf_w_olr(m) = pdf_w_olr(m)/PDF_olr(m)
  END IF
END DO
print*,'B'
FRQ = FRQ/cnt
MAP_A = MAP_A/cnt_a
MAP_B = MAP_B/cnt_b
MAP_C = MAP_C/cnt_c
MAP_D = MAP_D/cnt_d
PDF_crh = PDF_crh/SUM(PDF_crh)
PDF_olr = PDF_olr/SUM(PDF_olr)

filename = 'CRH_OLR_freq.nc'
CALL check(nf90_create(filename, nf90_clobber, ncid))
CALL check(nf90_def_dim(ncid, 'lon', ii, i_dimid))
CALL check(nf90_def_dim(ncid, 'lat', jj, j_dimid))
CALL check(nf90_def_dim(ncid, 'crh', nn, n_dimid))
CALL check(nf90_def_dim(ncid, 'olr', mm, m_dimid))
print*, 'here2'
dimids1 = (/ n_dimid, m_dimid /)
dimids2 = (/ i_dimid, j_dimid /)
CALL check(nf90_def_var(ncid, "lon", nf90_double, i_dimid, varid1))
CALL check(nf90_def_var(ncid, "lat", nf90_double, j_dimid, varid2))
CALL check(nf90_def_var(ncid, "crh", nf90_double, n_dimid, varid3))
CALL check(nf90_def_var(ncid, "olr", nf90_double, m_dimid, varid4))
CALL check(nf90_def_var(ncid, "frequency", nf90_double,  dimids1, varid5))
CALL check(nf90_def_var(ncid, "map_a", nf90_double,  dimids2, varid6))
CALL check(nf90_def_var(ncid, "map_b", nf90_double,  dimids2, varid7))
CALL check(nf90_def_var(ncid, "map_c", nf90_double,  dimids2, varid8))
CALL check(nf90_def_var(ncid, "pdf_crh", nf90_double,  n_dimid, varid9))
CALL check(nf90_def_var(ncid, "pdf_olr", nf90_double,  m_dimid, varid10))
CALL check(nf90_def_var(ncid, "omega_2pdf", nf90_double,  dimids1, varid11))
CALL check(nf90_def_var(ncid, "pdf_crh_olr", nf90_double,  n_dimid, varid12))
CALL check(nf90_def_var(ncid, "pdf_w_olr", nf90_double,  m_dimid, varid13))
CALL check(nf90_def_var(ncid, "map_d", nf90_double,  dimids2, varid14))
CALL check(nf90_enddef(ncid))

print*, 'here3'
CALL check(nf90_put_var(ncid, varid1, Lon))
CALL check(nf90_put_var(ncid, varid2, Lat))
CALL check(nf90_put_var(ncid, varid3, Bin_crh))
CALL check(nf90_put_var(ncid, varid4, Bin_olr))
CALL check(nf90_put_var(ncid, varid5, FRQ))
CALL check(nf90_put_var(ncid, varid6, MAP_A))
CALL check(nf90_put_var(ncid, varid7, MAP_B))
CALL check(nf90_put_var(ncid, varid8, MAP_C))
CALL check(nf90_put_var(ncid, varid9, PDF_crh))
CALL check(nf90_put_var(ncid, varid10, PDF_olr))
CALL check(nf90_put_var(ncid, varid11, W_2pdf))
CALL check(nf90_put_var(ncid, varid12, pdf_crh_olr))
CALL check(nf90_put_var(ncid, varid13, pdf_w_olr))
CALL check(nf90_put_var(ncid, varid14, MAP_D))
CALL check(nf90_close(ncid))

PRINT *, 'SUCCESS writing ', filename
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
