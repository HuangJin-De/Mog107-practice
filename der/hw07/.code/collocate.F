PROGRAM readdata
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da
INTEGER :: myid,nproc,ts,te,dnt
INTEGER, DIMENSION(12) :: mn_da 
INTEGER :: i,j,k,m,n,o,t,t_out
INTEGER :: ii,jj,is,ie,js,je
REAL, PARAMETER :: tripi=4*atan(1.)
REAL :: dum1,dum2,dum3,dum4,dum5
REAL :: dum6,dum7,dum8,dum9
REAL :: polysvp1,qv_sat
INTEGER :: idum1,idum2,idum3,idum4,idum5
INTEGER :: idum6,idum7,idum8
REAL, DIMENSION(100) :: cth
REAL, DIMENSION(480) :: lon 
REAL, DIMENSION(81) :: lat
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: sur
INTEGER :: ierr,ncid1,varid1,ncid2,varid2,ncid3,varid3,ncid4,varid4
INTEGER :: dimid, access
INTEGER, DIMENSION(4) :: startmat, endmat
CHARACTER(300) :: path, filename, region, ncname, varname
CHARACTER(10) :: int2str

! compile command
! ifort -free read.F
! -I/data/cloud/.local/include
! -L/data/cloud/.local/lib -lnetcdff -lnetcdf

call mpi_start(myid,nproc)
nt=10
dnt=nt/nproc
ts=myid*dnt+1
te=(myid+1)*dnt
ts=ts+2005
te=te+2005

write(*,*) myid, ts, te

path="/data/dadm1/obs/CloudSat_object/"

DO yr=ts,te
mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
IF (mod(yr,4)==0) mn_da(2)=mn_da(2)+1
nt=sum(mn_da)

ALLOCATE(sur(144,24,2,nt))
sur=0.

filename=trim(path)//"object_"//trim(int2str(yr))//".txt"
OPEN(11,file=trim(filename))

DO
READ(11,555,iostat=idum1) idum2,dum1,dum2,idum8,idum6,idum7,idum3,dum3,idum4,dum4,idum5
555 FORMAT(5X,I5,6X,2F9.2,3I9,5X,2(I3,F7.0),30X,I6)
IF (idum1/=0) GOTO 735

dum8=real(idum6)
dum9=real(idum7)

IF (abs(dum1)>30.) GOTO 841

dum2=dum2+1.25
IF (dum2<0.) dum2=dum2+360.

dum2=dum2/2.5+1.
dum1=(dum1+30.)/2.5+1.
i=int(dum2)
j=int(dum1)

dum2=real(idum5)/2./277.77777778*1.4
idum5=int(dum2)+1

IF (dum8>=6000. .AND. dum9<=2000.)THEN
  sur(i,j,1,idum2)=sur(i,j,1,idum2)+1.
ELSE
  sur(i,j,2,idum2)=sur(i,j,2,idum2)+1.
ENDIF

IF (idum5>=1) THEN

DO jj=j-idum5,j+idum5
DO ii=i-idum5,i+idum5
  is=ii
  js=jj
  IF (is>144) is=ii-144
  IF (is<1) is=ii+144
  IF (js>24) GOTO 152
  IF (js<1) GOTO 152

  dum1=sqrt(real(ii-i)**2+real(jj-j)**2)
  IF (dum1<=dum2) THEN
    IF (dum8>=6000. .AND. dum9<=2000.)THEN
      sur(is,js,1,idum2)=sur(is,js,1,idum2)+1.
    ELSE
      sur(is,js,2,idum2)=sur(is,js,2,idum2)+1.
    ENDIF
  ENDIF 
  152 CONTINUE
ENDDO
ENDDO

ENDIF

841 CONTINUE

ENDDO
735 CONTINUE

filename="for_you/object_col_"//trim(int2str(yr))//".dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*24*2*nt)
WRITE(10,REC=1) sur
CLOSE(10)
CLOSE(11)

WRITE(*,*) yr
DEALLOCATE(sur)
ENDDO

call mpi_end()

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

subroutine mpi_start(myid,nproc)
implicit none
include 'mpif.h'

integer :: ierr, nproc, myid

CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)

end subroutine mpi_start

subroutine mpi_end()
implicit none
include 'mpif.h'

integer :: ierr

CALL MPI_FINALIZE(ierr)
end subroutine mpi_end

