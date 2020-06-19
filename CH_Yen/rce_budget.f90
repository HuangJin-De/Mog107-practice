PROGRAM rce_budget
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/06/11                                   !
! Thisi program aims to                                   !
! 1) calculate area-mean of radiation, LH, SH in tropical Pacific        !
! Inputs:  TERRA-AQUA, TRMM, ERA                                  !
! Outputs: RCE budget (W/m^2) daily data                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf
IMPLICIT NONE
INTEGER :: i, ii, j, jj, k, kk=3, d, dd, t, tt, y, cnt
INTEGER :: ncid, dimid, varid, varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8
INTEGER :: k_dimid, i_dimid, j_dimid, d_dimid, m_dimid, l_dimid, nofill
INTEGER, DIMENSION(3) :: dimids3d
INTEGER, DIMENSION(2) :: dimids2d
INTEGER :: t1, t2, clock_rate, clock_max
CHARACTER(len=80) :: path(3)
CHARACTER(len=80) :: input, output
CHARACTER(len=15) :: dimnameii, dimnamejj, dimnamekk, dimnamedd, dimnamett
CHARACTER :: var_name*20
REAL, PARAMETER :: Lv=2.5*10.**6., Rv=461., rho_w=1000., g=9.8
REAL, DIMENSION(3) ::  sth, nth, wst, est
REAL ::  fillvalue
REAL, PARAMETER :: crh_max = 100., pre_max = 101000., pre_min = 5000.
REAL, DIMENSION(:), ALLOCATABLE :: Lon, Lat, Lon_LH, Lat_LH, Lon_SH, Lat_SH
REAL, DIMENSION(:,:,:), ALLOCATABLE :: toa_net, sfc_sw_up, sfc_sw_dn, sfc_lw_up, sfc_lw_dn, net_rad, LH, SH
REAL, DIMENSION(:,:), ALLOCATABLE :: pcp, rad_ave, LH_ave, SH_ave, net_flux_ave 
CHARACTER(len=4) :: year, mth, day
INTEGER, DIMENSION(2) :: date
INTEGER, PARAMETER :: mm=12,  nn=100, ll=50, init_yr=2001, finl_yr=2015
INTEGER, DIMENSION(12), PARAMETER :: days=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
INTEGER, PARAMETER :: yy=finl_yr-init_yr+1
INTEGER :: sth_y, nth_y, wst_x, est_x

path(1) = '/data/dadm1/obs/CERES_SYN1deg-Day_Terra-Aqua-MODIS/'
path(2) = '/data/dadm1/obs/TRMM/TRMM3B42/'
path(3) = '/data/dadm1/reanalysis/ECMWF/ITM/daily/SH/'

sth = (/-30., -30., -10./)
nth = (/30., 30., 10./)
wst = (/90., 160., 90./)
est = (/240., 240., 240./)

DO y = 1, yy
  !DO k = 1, kk
    WRITE(year, '(I4.2)') init_yr+y-1
    ! Terra-Aqua radiation data
    input = 'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_'//TRIM(year)//'.nc'
    PRINT*, TRIM(path(1))//TRIM(input)
    CALL check(nf90_open(TRIM(path(1))//TRIM(input), nf90_nowrite, ncid))
    CALL check(nf90_inq_dimid(ncid, 'lon', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
    CALL check(nf90_inq_dimid(ncid, 'lat', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
    CALL check(nf90_inq_dimid(ncid, 'time', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamedd, dd))

    ALLOCATE(Lon(ii))
    ALLOCATE(Lat(jj))
    ALLOCATE(toa_net(ii,jj,dd))
    ALLOCATE(sfc_sw_up(ii,jj,dd))
    ALLOCATE(sfc_sw_dn(ii,jj,dd))
    ALLOCATE(sfc_lw_up(ii,jj,dd))
    ALLOCATE(sfc_lw_dn(ii,jj,dd))
    ALLOCATE(net_rad(ii,jj,dd))
    ALLOCATE(rad_ave(dd,kk)); rad_ave = 0.

    CALL check(nf90_inq_varid(ncid, 'lon', varid))
    CALL check(nf90_get_var(ncid, varid, Lon))
    CALL check(nf90_inq_varid(ncid, 'lat', varid))
    CALL check(nf90_get_var(ncid, varid, Lat))
    CALL check(nf90_inq_varid(ncid, 'toa_net_all_daily', varid))
    CALL check(nf90_get_var(ncid, varid, toa_net))
    CALL check(nf90_inq_varid(ncid, 'ini_sfc_sw_up_all_daily', varid))
    CALL check(nf90_get_var(ncid, varid, sfc_sw_up))
    CALL check(nf90_inq_varid(ncid, 'ini_sfc_sw_down_all_daily', varid))
    CALL check(nf90_get_var(ncid, varid, sfc_sw_dn))
    CALL check(nf90_inq_varid(ncid, 'ini_sfc_lw_up_all_daily', varid))
    CALL check(nf90_get_var(ncid, varid, sfc_lw_up))
    CALL check(nf90_inq_varid(ncid, 'ini_sfc_lw_down_all_daily', varid))
    CALL check(nf90_get_var(ncid, varid, sfc_lw_dn))
    CALL check(nf90_close(ncid))
  DO k = 1, kk
    nth_y = MAXLOC(Lat, 1, Lat <= nth(k))
    sth_y = MAXLOC(Lat, 1, Lat <= sth(k))
    wst_x = MAXLOC(Lon, 1, Lon <= wst(k))
    est_x = MAXLOC(Lon, 1, Lon <= est(k))
    DO d = 1, dd
      DO j = sth_y, nth_y
        DO i = wst_x, est_X
          rad_ave(d,k) = rad_ave(d,k)+toa_net(i,j,d)+sfc_sw_up(i,j,d)-sfc_sw_dn(i,j,d)+sfc_lw_up(i,j,d)-sfc_lw_dn(i,j,d)
        END DO ! ends DO i = wst_x, est_X
      END DO ! ends DO j = sth_y, nth_y
      rad_ave(d,k) = rad_ave(d,k)/(est_x-wst_x)/(nth_y-sth_y)
    END DO ! ends DO d = 1, dd
  END DO ! ends DO k = 1, kk
    net_rad = toa_net+sfc_sw_up-sfc_sw_dn+sfc_lw_up-sfc_lw_dn
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    input = '3B42.'//TRIM(year)//'.3hr.nc'
    PRINT*, TRIM(path(2))//TRIM(input)
    CALL check(nf90_open(TRIM(path(2))//TRIM(input), nf90_nowrite, ncid))
    CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
    CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
    CALL check(nf90_inq_dimid(ncid, 'time', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))

    ALLOCATE(Lon_LH(ii))
    ALLOCATE(Lat_LH(jj))
    ALLOCATE(pcp(ii,jj))
    ALLOCATE(LH(ii,jj,tt/8)); LH = 0
    ALLOCATE(LH_ave(tt/8,kk)); LH_ave = 0

    CALL check(nf90_inq_varid(ncid, 'longitude', varid))
    CALL check(nf90_get_var(ncid, varid, Lon_LH))
    CALL check(nf90_inq_varid(ncid, 'latitude', varid))
    CALL check(nf90_get_var(ncid, varid, Lat_LH))
    !CALL check(nf90_inq_varid(ncid, 'pcp', varid))
    !CALL check(nf90_get_var(ncid, varid, pcp, (1,1,t), (ii,jj,1)))
  !DO k = 1, kk
  !PRINT*,'k=',k
  !  nth_y = MAXLOC(Lat_LH, 1, Lat_LH <= nth(k))
  !  sth_y = MAXLOC(Lat_LH, 1, Lat_LH <= sth(k))
  !  wst_x = MAXLOC(Lon_LH, 1, Lon_LH <= wst(k))
  !  est_x = MAXLOC(Lon_LH, 1, Lon_LH <= est(k)-360.)
    d = 1
    DO t = 1, tt
      CALL check(nf90_inq_varid(ncid, 'pcp', varid))
      CALL check(nf90_inq_var_fill(ncid, varid, nofill, fillvalue))
      CALL check(nf90_get_var(ncid, varid, pcp, (/1,1,t/), (/ii,jj,1/)))
      LH(:,:,d) = LH(:,:,d)+pcp*1E-3*3*rho_w*Lv/86400. ! to get the unit: W/m^2
      IF (MOD(t,8) == 0) THEN
        DO k = 1, kk
          PRINT*,'k=',k
          nth_y = MAXLOC(Lat_LH, 1, Lat_LH <= nth(k))
          sth_y = MAXLOC(Lat_LH, 1, Lat_LH <= sth(k))
          wst_x = MAXLOC(Lon_LH, 1, Lon_LH <= wst(k))
          est_x = MAXLOC(Lon_LH, 1, Lon_LH <= est(k)-360.)
        cnt = 0
        DO j = sth_y, nth_y
          DO i = 1, est_x
            IF (LH(i,j,d) < 0) THEN
              CYCLE
            END IF
            LH_ave(d,k) = LH_ave(d,k)+LH(i,j,d)
            cnt = cnt+1
          END DO
          DO i = wst_x, ii
            IF (LH(i,j,d) < 0) THEN
              CYCLE
            END IF
            LH_ave(d,k) = LH_ave(d,k)+LH(i,j,d)
            cnt = cnt+1
          END DO
        END DO
        LH_ave(d,k) = LH_ave(d,k)/cnt
        END DO ! ends DO k = 1, kk
        d = d+1
print*,d
      END IF
    END DO ! ends DO t = 1, tt, 8
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    input = 'daily_interim_SH_'//TRIM(year)//'.nc'
    PRINT*, TRIM(path(3))//TRIM(input)
    CALL check(nf90_open(TRIM(path(3))//TRIM(input), nf90_nowrite, ncid))
    CALL check(nf90_inq_dimid(ncid, 'longitude', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnameii, ii))
    CALL check(nf90_inq_dimid(ncid, 'latitude', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamejj, jj))
    CALL check(nf90_inq_dimid(ncid, 'time', dimid))
    CALL check(nf90_inquire_dimension(ncid, dimid, dimnamett, tt))

    ALLOCATE(Lon_SH(ii))
    ALLOCATE(Lat_SH(jj))
    ALLOCATE(SH(ii,jj,dd))
    ALLOCATE(SH_ave(dd,kk)); SH_ave = 0

    CALL check(nf90_inq_varid(ncid, 'longitude', varid))
    CALL check(nf90_get_var(ncid, varid, Lon_SH))
    CALL check(nf90_inq_varid(ncid, 'latitude', varid))
    CALL check(nf90_get_var(ncid, varid, Lat_SH))
    CALL check(nf90_inq_varid(ncid, 'sshf', varid))
    CALL check(nf90_get_var(ncid, varid, SH))
  DO k = 1, kk  
    nth_y = MAXLOC(Lat_SH, 1, Lat_SH >= nth(k))
    sth_y = MAXLOC(Lat_SH, 1, Lat_SH <= sth(k))
    wst_x = MAXLOC(Lon_SH, 1, Lon_SH <= wst(k))
    est_x = MAXLOC(Lon_SH, 1, Lon_SH <= est(k))

    DO d = 1, dd
      DO j = nth_y, sth_y
        DO i = wst_x, est_X
          SH_ave(d,k) = SH_ave(d,k)+SH(i,j,d)
        END DO ! ends DO i = wst_x, est_X
      END DO ! ends DO j = sth_y, nth_y
      SH_ave(d,k) = -1.*SH_ave(d,k)/(est_x-wst_x)/(sth_y-nth_y)/86400.
    END DO ! ends DO d = 1, dd
  END DO ! ends DO k = 1, kk
    ALLOCATE(net_flux_ave(dd,kk))
    net_flux_ave = rad_ave+LH_ave+SH_ave

    output = 'rce_budget_'//year//'.nc'
    CALL check(nf90_create(output, nf90_clobber, ncid))
    CALL check(nf90_def_dim(ncid, 'lon', ii, i_dimid))
    CALL check(nf90_def_dim(ncid, 'lat', jj, j_dimid))
    CALL check(nf90_def_dim(ncid, 'day', dd, d_dimid))
    CALL check(nf90_def_dim(ncid, 'domain', kk, k_dimid))

    dimids2d = (/d_dimid, k_dimid/)
    dimids3d = (/i_dimid, j_dimid, d_dimid/)
    CALL check(nf90_def_var(ncid, 'lon', nf90_double, i_dimid, varid1))
    CALL check(nf90_def_var(ncid, 'lat', nf90_double, j_dimid, varid2))
    CALL check(nf90_def_var(ncid, 'net_rad', nf90_double, dimids3d, varid3))
    CALL check(nf90_def_var(ncid, 'rad_ave', nf90_double, dimids2d, varid4))
    CALL check(nf90_def_var(ncid, 'LH_ave', nf90_double, dimids2d, varid5))
    CALL check(nf90_def_var(ncid, 'SH_ave', nf90_double, dimids2d, varid6))
    CALL check(nf90_def_var(ncid, 'net_flux_ave', nf90_double, dimids2d, varid7))
    CALL check(nf90_enddef(ncid))

    CALL check(nf90_put_var(ncid, varid1, lon))
    CALL check(nf90_put_var(ncid, varid2, lat))
    CALL check(nf90_put_var(ncid, varid3, net_rad))
    CALL check(nf90_put_var(ncid, varid4, rad_ave))
    CALL check(nf90_put_var(ncid, varid5, LH_ave))
    CALL check(nf90_put_var(ncid, varid6, SH_ave))
    CALL check(nf90_put_var(ncid, varid7, net_flux_ave))
    CALL check(nf90_close(ncid))
    PRINT *, 'SUCCESS writing ', output
    DEALLOCATE(Lon)
    DEALLOCATE(Lat)
    DEALLOCATE(toa_net)
    DEALLOCATE(sfc_sw_up)
    DEALLOCATE(sfc_sw_dn)
    DEALLOCATE(sfc_lw_up)
    DEALLOCATE(sfc_lw_dn)
    DEALLOCATE(net_rad)
    DEALLOCATE(rad_ave)
 
    DEALLOCATE(Lon_LH)
    DEALLOCATE(Lat_LH)
    DEALLOCATE(pcp)
    DEALLOCATE(LH)
    DEALLOCATE(LH_ave)

    DEALLOCATE(Lon_SH)
    DEALLOCATE(Lat_SH)
    DEALLOCATE(SH)
    DEALLOCATE(SH_ave)

    DEALLOCATE(net_flux_ave)
END DO ! ends DO y = 1, yy
CONTAINS
SUBROUTINE check(STATUS)
        INTEGER, INTENT( in) :: STATUS
        IF (STATUS /= nf90_noerr) THEN
                PRINT*, trim(nf90_strerror(STATUS))
                STOP 2
        END IF
END SUBROUTINE
END PROGRAM rce_budget
