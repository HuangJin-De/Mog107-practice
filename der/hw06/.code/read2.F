PROGRAM readdata
USE netcdf
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da
INTEGER, DIMENSION(12) :: mn_da 
INTEGER :: i,j,k,m,n,o,t,t_out,jj,kk
REAL, PARAMETER :: tripi=4*atan(1.)
REAL :: dum1,dum2,dum3,dum4,dum5
REAL :: dum6,dum7,dum8,dum9
REAL*8 :: day_tot
REAL :: rate, te(81)
REAL :: polysvp1,qv_sat
INTEGER :: idum1,idum2,idum3,idum4,idum5
REAL, DIMENSION(101) :: h_re
REAL, DIMENSION(100) :: hc_re, warm
REAL, DIMENSION(100) :: cth
REAL, DIMENSION(37) :: p, con_t
REAL*8, DIMENSION(37) :: hm_m, con_m, hm_m0, con_m0, hm_m10, con_m10
REAL, DIMENSION(32) :: lat, f_area
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: obj
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sur
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: spa
REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: var
INTEGER :: ierr,ncid1,varid1,ncid2,varid2,ncid3,varid3,ncid4,varid4
INTEGER :: dimid, access
INTEGER, DIMENSION(4) :: startmat, endmat
CHARACTER(300) :: path, filename, region, ncname, varname, run, path2
CHARACTER(10) :: int2str

! count
INTEGER, PARAMETER :: ncri=100
REAL, DIMENSION(ncri+1) :: cri, cri2
REAL, DIMENSION(10+1) :: cri_p
REAL, DIMENSION(ncri,10) :: pdf, w_i, mf, sf, qv_i, rho_i, con, meanhm
REAL, DIMENSION(ncri,10) :: pdf_c, w_c, mf_c, sf_c, qv_c, rho_c, con_c
REAL, DIMENSION(ncri,10) :: pdf_0dbz, pdf_10dbz, meanhm0, meanhm10
REAL, DIMENSION(ncri,10) :: pdf_10dbzc, hm_c
REAL, DIMENSION(ncri,10,3) :: pdf_0dbzc

! compile command
! ifort -free read.F
! -I/data/cloud/.local/include
! -L/data/cloud/.local/lib -lnetcdff -lnetcdf

path="/data/dadm1/model_output/SPCAM/"
run="sed_run_name"

yr=0001
filename=trim(path)//trim(run)//"/"//trim(run)//".cam.h0."//trim(int2str(yr,4))//"-01-01-00000.nc"
ierr = nf90_open(trim(filename),nf90_nowrite,ncid1)
ierr = nf90_inq_varid(ncid1,"lat",varid1)
ierr = nf90_get_var(ncid1,varid1,lat,start=(/ 33 /),count=(/ 32 /))
f_area=cos(lat/180.*tripi)
ierr = nf90_close(ncid1)

path="/data/cloud/der/hw06/for_your_obj/"
path2="/data/cloud/der/hw06/for_you/"
run="sed_run_name"

DO i=1,ncri+1
  dum1=0.
  dum2=1.
  dum4=(dum2-dum1)/real(ncri)
  cri2(i)=dum1+real(i-1)*dum4
ENDDO

! 900 800 700 600 500 400 300 200 100 50
cri_p=(/ 1000.,850.,750.,650.,&
         550.,450.,350.,250.,150.,75.,25. /)

ALLOCATE(spa(144,32,5))

spa=0.
day_tot=0.

DO yr=0001,0010
mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
nt=sum(mn_da)
day_tot=day_tot+real(nt)

ALLOCATE(var(144,32,26,nt,6),sur(144,32,nt),obj(144,32,3,nt))

filename=trim(path2)//trim(run)//"_"//trim(int2str(yr,4))//".dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*32*26*sum(mn_da)*5)
READ(10,REC=1) var(:,:,:,:,1:5)
CLOSE(10)

var(:,:,:,:,6)=var(:,:,:,:,2)+(var(:,:,:,:,4)*9.8+&
               2.5E+6*var(:,:,:,:,3))/1004.5

filename=trim(path)//"object_col_"//trim(run)//"_"//trim(int2str(yr,4))//".dat"
OPEN(30,FILE=trim(filename),ACCESS="direct",RECL=144*32*3*nt)
READ(30,REC=1) obj
CLOSE(30)

DO t=1,nt
DO j=1,32
DO i=1,144

dum6=0.
dum7=0.
DO k=26,1,-1

  dum1=var(i,j,k,t,2)
  dum2=var(i,j,k,t,1)/287.047/dum1

  IF (k<26) THEN
    dum6=dum6-dum2*0.5*(var(i,j,k,t,3)+var(i,j,k+1,t,3))*&
                       (var(i,j,k+1,t,4)-var(i,j,k,t,4))
    IF (var(i,j,k,t,1)>10000.) &
    dum7=dum7-dum2*0.5*(qv_sat(var(i,j,k,t,2),var(i,j,k,t,1),0)&
                       +qv_sat(var(i,j,k+1,t,2),var(i,j,k+1,t,1),0))*&
                       (var(i,j,k+1,t,4)-var(i,j,k,t,4))
  ENDIF  

ENDDO

!write(*,*) dum6, dum7

dum6=dum6/dum7

!write(*,*) dum6

IF (obj(i,j,1,t)>0.) spa(i,j,2)=spa(i,j,2)+dble(1.)
IF (obj(i,j,2,t)>0.) spa(i,j,3)=spa(i,j,3)+dble(1.)
IF (obj(i,j,3,t)>0.) spa(i,j,4)=spa(i,j,4)+dble(1.)
IF (dum6<0.6) THEN
  spa(i,j,1)=spa(i,j,1)+dble(1.)
  IF (obj(i,j,1,t)>0.) spa(i,j,5)=spa(i,j,5)+dble(1.)
ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO n=1,ncri
  IF (dum6>=cri2(n) .AND. dum6<cri2(n+1))THEN
  DO k=24,1,-1
    dum3=var(i,j,k,t,1)/100.
    DO kk=1,10
      IF (dum3<=cri_p(kk) .and. dum3>cri_p(kk+1)) THEN
        dum1=var(i,j,k,t,2)
        dum2=var(i,j,k,t,1)/287.047/dum1
        con_c(n,kk)=con_c(n,kk)+f_area(j)
        mf_c(n,kk)=mf_c(n,kk)+(-9.8*dum2*var(i,j,k,t,5))*dum2/dum4*f_area(j)
        qv_c(n,kk)=qv_c(n,kk)+var(i,j,k,t,3)*dum2/dum4*f_area(j)
        rho_c(n,kk)=rho_c(n,kk)+dum2/dum4*f_area(j)
        hm_c(n,kk)=hm_c(n,kk)+var(i,j,k,t,6)*dum2/dum4*f_area(j)   
        IF (obj(i,j,1,t)>0.) pdf_0dbzc(n,kk,1)=pdf_0dbzc(n,kk,1)+1.
        IF (obj(i,j,2,t)>0.) pdf_0dbzc(n,kk,2)=pdf_0dbzc(n,kk,2)+1.
        IF (obj(i,j,3,t)>0.) pdf_0dbzc(n,kk,3)=pdf_0dbzc(n,kk,3)+1.
        GOTO 654
      ENDIF
    ENDDO
    654 CONTINUE
  ENDDO
  ENDIF
ENDDO

ENDDO
ENDDO
ENDDO

DEALLOCATE(var,sur,obj)
WRITE(*,*) yr
ENDDO

!hm_m=hm_m/con_m
!hm_m0=hm_m0/con_m0
!hm_m10=hm_m10/con_m10

DO k=1,10
  !meanhm(:,k)=real(hm_m(k))
  !meanhm0(:,k)=real(hm_m0(k))
  !meanhm10(:,k)=real(hm_m10(k))
  DO n=1,ncri
    !sf(n,k)=sum(mf(1:n,k),1)*dum3
    sf_c(n,k)=sum(mf_c(1:n,k),1)*dum4
  ENDDO

!con_t(k)=sum(con(:,k))
!pdf_0dbz(:,k)=pdf_0dbz(:,k)/con(:,k)
!pdf_10dbz(:,k)=pdf_10dbz(:,k)/con(:,k)
!con(:,k)=con(:,k)/con_t(k)
!sf(:,k)=sf(:,k)/con_t(k)
!qv_i(:,k)=qv_i(:,k)/con_t(k)
!rho_i(:,k)=rho_i(:,k)/con_t(k)
!mf(:,k)=mf(:,k)/con_t(k)

con_t(k)=sum(con_c(:,k))
pdf_10dbz(:,k)=pdf_10dbzc(:,k)/con_t(k)
con_c(:,k)=con_c(:,k)/con_t(k)
sf_c(:,k)=sf_c(:,k)/con_t(k)
qv_c(:,k)=qv_c(:,k)/con_t(k)
rho_c(:,k)=rho_c(:,k)/con_t(k)
mf_c(:,k)=mf_c(:,k)/con_t(k)
hm_c(:,k)=hm_c(:,k)/con_t(k)
ENDDO

!WHERE (con==0.) con=-1.
!w_i=mf/rho_i
!qv_i=qv_i/rho_i

WHERE (con_c==0.) con_c=-1.
w_c=mf_c/rho_c
qv_c=qv_c/rho_c
hm_c=hm_c/rho_c

spa(:,:,1)=spa(:,:,1)/dble(day_tot)

filename=trim(run)//"_isen_data.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=ncri*10*12)
WRITE(10,REC=1) con_c, mf_c, sf_c, w_c, qv_c, hm_c, rho_c, pdf_0dbzc,&
                pdf_10dbzc, pdf_10dbz

CLOSE(10)

filename=trim(run)//"_map_data.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*32*5)
WRITE(10,REC=1) real(spa)
CLOSE(10)

END PROGRAM readdata

!===========================================================================================

FUNCTION int2str(i,n) result(str)

INTEGER, INTENT(in) :: i, n
CHARACTER(10) :: str

IF (n==2) WRITE(str,122) i
IF (n==4) WRITE(str,123) i
122 FORMAT(I10.2)
123 FORMAT(I10.4)
str=trim(adjustl(str))

END FUNCTION int2str

!===========================================================================================

 real function qv_sat(t_atm,p_atm,i_wrt)

!------------------------------------------------------------------------------------
! Calls polysvp1 to obtain the saturation vapor pressure, and then computes
! and returns the saturation mixing ratio, with respect to either liquid or ice,
! depending on value of 'i_wrt'
!------------------------------------------------------------------------------------

 implicit none

 !Calling parameters:
 real    :: t_atm  !temperature [K]
 real    :: p_atm  !pressure    [Pa]
 integer :: i_wrt  !index, 0 = w.r.t. liquid, 1 = w.r.t. ice

 !Local variables:
 real    :: e_pres         !saturation vapor pressure [Pa]
 real    :: ep_2
 real    :: polysvp1
 !------------------

 ep_2   = 0.622
 e_pres = polysvp1(t_atm,i_wrt)
 qv_sat = ep_2*e_pres/max(1.e-3,(p_atm-e_pres))

 return
 end function qv_sat
!===========================================================================================

 real function polysvp1(T,i_type)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      implicit none

      real    :: DUM,T
      integer :: i_type

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4
! (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

!-------------------------------------------

      if (i_type.EQ.1 .and. T.lt.273.15) then
! ICE

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*&
                    (a4i+dt*(a5i+dt*(a6i+dt*(a7i+&
                    a8i*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                 &
!          log10(273.16/T)+0.876793*(1.-T/273.16)+                        &
!          log10(6.1071))*100.


      elseif (i_type.EQ.0 .or. T.ge.273.15) then
! LIQUID

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*&
                    (a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                         &
!             5.02808*log10(373.16/T)-                                    &
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                  &
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                &
!             log10(1013.246))*100.

         endif

 end function polysvp1
!===========================================================================================


