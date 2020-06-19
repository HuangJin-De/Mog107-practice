PROGRAM cal
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da, is,ie,js,je,ii,jj,kk
INTEGER :: yr1, mn1, da1, hr1
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
INTEGER, DIMENSION(8,2) :: oni, oni_tem
REAL, DIMENSION(144,72,3) :: spa
CHARACTER(300) :: path, filename, run, model, ncname, varname
CHARACTER(10) :: int2str

tripi=4*atan(1.)
path="/data/cloud/der/hw15/data/"

oni(:,1)=(/ 1983,1987,1988,1992,1995,1998,2003,2010 /)
oni(:,2)=(/ 1985,1989,1996,1999,2000,2007,2011,2012 /)


spa=0.

DO yr=1979,2015
mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
IF (mod(yr,4)==0) mn_da(2)=mn_da(2)+1
nt=sum(mn_da)
day_tot=day_tot+real(nt)

DO i=1,5

WRITE(filename,111) trim(path),"BT_",i,".txt"
111 FORMAT(A,A,I1,A)

OPEN(11,FILE=trim(filename))

DO
READ(11,184,IOSTAT=o) run, yr1, mn1, da1, hr1, dum2, dum1
184 FORMAT(A9,I5,3I3,2F7.1)

IF (o/=0) GOTO 888
IF (yr1/=yr) GOTO 777
IF (hr1/=0) GOTO 777
IF (abs(dum1)>30.) GOTO 777

! Julian day
t=sum(mn_da(1:mn1-1),1)+da1

IF (sum(mn_da(1:2),1)<t .AND. t<=sum(mn_da(1:11),1)) GOTO 777
day_tot(1)=day_tot(1)+1
enso=-1
oni_tem=oni
IF (t>sum(mn_da(1:11),1)) oni_tem=oni-1
DO ii=1,8
  IF (yr==oni_tem(ii,1)) enso=1
  IF (yr==oni_tem(ii,2)) enso=2
ENDDO
IF (enso>0) day_tot(enso+1)=day_tot(enso+1)+1

! lon lat
dum2=dum2+1.25
IF (dum2<0.) dum2=dum2+360.

dum2=dum2/2.5+1.
dum1=(dum1+90.)/2.5+1.
ii=int(dum2)
jj=int(dum1)


spa(ii,jj,1)=spa(ii,jj,1)+1.
IF (enso>0) spa(ii,jj,enso+1)=spa(ii,jj,enso+1)+1.

777 CONTINUE
ENDDO
888 CONTINUE

CLOSE(11)

ENDDO

WRITE(*,*) yr
ENDDO


!DO enso=0,2
!  spa(:,:,enso+1)=spa(:,:,enso+1)/real(day_tot(enso+1))
!ENDDO

spa(:,:,1)=spa(:,:,1)/37.
spa(:,:,2:3)=spa(:,:,2:3)/8.


WHERE(spa==0) spa=-1
filename="map_data.dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*72)
DO enso=0,2
  WRITE(10,REC=enso+1) spa(:,:,enso+1)
ENDDO
CLOSE(10)

END PROGRAM cal
