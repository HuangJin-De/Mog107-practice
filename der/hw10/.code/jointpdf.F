PROGRAM readdata
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da, is,ie,js,je,ii,jj,kk
INTEGER, DIMENSION(12) :: mn_da 
INTEGER :: i,j,k,m,n,o,t,t_out,jjj,nn,nnn,cs,enso
REAL :: tripi
REAL*8 :: ddum1,ddum2
REAL :: dum1,dum2,dum3,dum4,dum5
REAL :: dum6,dum7,dum8,dum9,dum0
REAL, DIMENSION(3) :: day_tot
REAL :: dum13, dum14
REAL :: polysvp1,qv_sat
INTEGER :: idum1,idum2,idum3,idum4,idum5
INTEGER :: idum6,idum7,idum8
REAL, DIMENSION(101) :: h_re
REAL, DIMENSION(100) :: hc_re, warm
REAL, DIMENSION(100) :: cth
REAL, DIMENSION(37) :: p
REAL, DIMENSION(72) :: lat, f_area
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: obj, obj2
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sur, olr
INTEGER, DIMENSION(8,2) :: oni, oni_tem
REAL, DIMENSION(144,72,6) :: spa
REAL, DIMENSION(144,72) :: lsm
REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: var, temp
INTEGER :: ierr,ncid1,varid1,ncid2,varid2,ncid3,varid3,ncid4,varid4
INTEGER :: dimid, access
INTEGER, DIMENSION(4) :: startmat, endmat
CHARACTER(300) :: path, filename, run, model, ncname, varname
CHARACTER(10) :: int2str

! count
INTEGER, PARAMETER :: ncri1=100
REAL, DIMENSION(ncri1+1) :: cri1
REAL*8, DIMENSION(ncri1,37,9) :: crh_ver
REAL*8, DIMENSION(ncri1,3) :: crh

! compile command
! ifort -free read.F
! -I/data/cloud/.local/include
! -L/data/cloud/.local/lib -lnetcdff -lnetcdf

tripi=4*atan(1.)

path="/data/cloud/der/hw09/for_you/"

lat(1)=-88.75
DO i=2,72
  lat(i)=lat(i-1)+2.5 
ENDDO
f_area=cos(lat/180.*tripi)

p=(/ 1000,975,950,925,900,875,850,825,800,&
     775,750,700,650,600,550,500,450,400,350,&
     300,250,225,200,175,150,125,100,70,50,&
     30,20,10,70,50,30,20,10  /)
p=p*100.

oni(:,1)=(/ 1983,1987,1988,1992,1995,1998,2003,2010 /)
oni(:,2)=(/ 1985,1989,1996,1999,2000,2007,2011,2012 /)

DO i=1,ncri1+1
  dum1=0.
  dum2=1.
  dum4=(dum2-dum1)/real(ncri1)
  cri1(i)=dum1+real(i-1)*dum4
ENDDO

crh=0.
crh_ver=0.

day_tot=0.

filename=trim(path)//"daily_interim_lsm_2p5d.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*72)
READ(10,REC=1) lsm
CLOSE(10)
lsm=lsm*100

DO yr=yr_s,yr_e
mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
nt=sum(mn_da)

ALLOCATE(var(144,72,37,5,nt),sur(144,72,nt),olr(144,72,nt))

filename=trim(path)//"daily_interim_"//trim(int2str(yr))//"_2p5d.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*72*37*4+144*72)
DO i=1,nt
  READ(10,REC=i) var(:,:,:,1:4,i), sur(:,:,i)
  var(:,:,:,1,i)=var(:,:,:,1,i)
  var(:,:,:,5,i)=1004.5*var(:,:,:,2,i)&
                +var(:,:,:,1,i)&
                +2.5E+6*var(:,:,:,3,i)
ENDDO
CLOSE(10)

filename=trim(path)//"daily_olr_"//trim(int2str(yr))//"_2p5d.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*72)
DO i=1,nt
  READ(10,REC=i) olr(:,:,i)
ENDDO
CLOSE(10)

t_loop: DO t=1,nt
  IF (sum(mn_da(1:2),1)<t .AND. t<=sum(mn_da(1:11),1)) GOTO 155
  day_tot(1)=day_tot(1)+1
  enso=-1
  oni_tem=oni
  IF (t>sum(mn_da(1:11),1)) oni_tem=oni-1
  DO i=1,8
    IF (yr==oni_tem(i,1)) enso=1
    IF (yr==oni_tem(i,2)) enso=2
  ENDDO
  IF (enso>0) day_tot(enso+1)=day_tot(enso+1)+1
y_loop: DO j=1,72

IF (abs(lat(j))>30.) GOTO 746

x_loop: DO i=1,144

!IF (lsm(i,j)>50.) GOTO 456

dum6=0.
dum7=0.
kk=0
idum3=0
integrate_cwv: DO k=1,37-1
  IF (sur(i,j,t)<p(k)) GOTO 852
  dum1=var(i,j,k,2,t)
  dum2=p(k)/287.047/dum1      

  dum6=dum6+dum2*0.5*(var(i,j,k,3,t)+var(i,j,k+1,3,t))*&
                     (var(i,j,k+1,1,t)-var(i,j,k,1,t))/9.8
  n=0
  m=0
  IF (var(i,j,k,2,t)<273.15) n=1
  IF (var(i,j,k+1,2,t)<273.15) m=1
  IF (p(k)>10000.) dum7=dum7+dum2*0.5*(qv_sat(var(i,j,k,2,t),p(k),n)&
                       +qv_sat(var(i,j,k+1,2,t),p(k+1),m))*&
                       (var(i,j,k+1,1,t)-var(i,j,k,1,t))/9.8

  852 CONTINUE
ENDDO integrate_cwv

!write(*,*) dum6, dum7

dum6=dum6/dum7 
IF (dum6>1.) dum6=1.
!spa(i,j,t)=dum6
!dum7=olr(i,j,t)

IF (dum6<=0.6) spa(i,j,1)=spa(i,j,1)+1
IF (dum6<=0.6 .AND. enso>0) spa(i,j,enso+1)=spa(i,j,enso+1)+1
IF (olr(i,j,t)<150.) spa(i,j,4)=spa(i,j,4)+1
IF (olr(i,j,t)<150. .AND.  enso>0) spa(i,j,enso+4)=spa(i,j,enso+4)+1

DO n=1,ncri1
  IF (dum6>=cri1(n) .AND. dum6<cri1(n+1))THEN
  crh(n,1)=crh(n,1)+f_area(j)
  IF (enso>0) crh(n,enso+1)=crh(n,enso+1)+f_area(j)
  DO k=1,37
    IF (sur(i,j,t)<p(k)) GOTO 853
    dum1=var(i,j,k,2,t)
    dum2=p(k)/287.047/dum1
    crh_ver(n,k,1)=crh_ver(n,k,1)+dble((-var(i,j,k,4,t)/9.8)/dum4*f_area(j))
    crh_ver(n,k,2)=crh_ver(n,k,2)+dble(var(i,j,k,5,t)*dum2/dum4*f_area(j))
    crh_ver(n,k,3)=crh_ver(n,k,3)+dble(dum2/dum4*f_area(j))
    IF (enso>0) THEN
      crh_ver(n,k,3*enso+1)=crh_ver(n,k,3*enso+1)+dble((-var(i,j,k,4,t)/9.8)/dum4*f_area(j))
      crh_ver(n,k,3*enso+2)=crh_ver(n,k,3*enso+2)+dble(var(i,j,k,5,t)*dum2/dum4*f_area(j))
      crh_ver(n,k,3*enso+3)=crh_ver(n,k,3*enso+3)+dble(dum2/dum4*f_area(j))
    ENDIF
    853 CONTINUE
  ENDDO
  GOTO 456
  ENDIF
ENDDO
456 CONTINUE 

ENDDO x_loop
746 CONTINUE
ENDDO y_loop
155 CONTINUE
ENDDO t_loop

DEALLOCATE(var,sur,olr)
WRITE(*,*) yr
ENDDO

DO enso=0,2
  ddum1=sum(crh(:,enso+1))
  crh(:,enso+1)=crh(:,enso+1)/ddum1
  DO i=1,3
  DO k=1,37
    crh_ver(:,k,3*enso+i)=crh_ver(:,k,3*enso+i)/ddum1
  ENDDO
  ENDDO
  crh_ver(:,:,3*enso+2)=crh_ver(:,:,3*enso+2)/crh_ver(:,:,3*enso+3)
ENDDO


filename="isen_data.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=ncri1)
n=1
DO enso=0,2
  DO i=1,1
    WRITE(10,REC=n) real(crh(:,enso+i))
    n=n+1
  ENDDO
  DO i=1,2
  DO k=1,37
    WRITE(10,REC=n) real(crh_ver(:,k,3*enso+i))
    n=n+1
  ENDDO
  ENDDO
ENDDO
CLOSE(10)

DO enso=0,2
  spa(:,:,enso+1)=spa(:,:,enso+1)/real(day_tot(enso+1))
ENDDO
DO enso=0,2
  spa(:,:,enso+4)=spa(:,:,enso+4)/real(day_tot(enso+1))
ENDDO

WHERE(spa==0) spa=-1
filename="map_data.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*72)
DO enso=0,2
  WRITE(10,REC=enso+1) spa(:,:,enso+1)
ENDDO
DO enso=0,2
  WRITE(10,REC=enso+4) spa(:,:,enso+4)
ENDDO
CLOSE(10)

END PROGRAM readdata

!===========================================================================================

FUNCTION int2str(i) result(str)

INTEGER, INTENT(in) :: i
CHARACTER(10) :: str

WRITE(str,123) i
123 FORMAT(I10)
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


