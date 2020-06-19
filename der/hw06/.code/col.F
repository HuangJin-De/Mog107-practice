PROGRAM readdata
USE netcdf
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da, nn, tt
INTEGER :: myid, dnt, nproc
INTEGER, DIMENSION(12) :: mn_da 
INTEGER :: i,j,k,m,n,o,t,t_out
INTEGER :: ii,jj,is,ie,js,je,ts,te
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
CHARACTER(300) :: path, filename, region, ncname, varname, run
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

write(*,*) myid, ts, te

path="/data/cloud/der/hw06/for_your_obj/"
run="sed_run_name"

!ts=1
!te=10

DO yr=ts,te
mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
nt=sum(mn_da)

ALLOCATE(sur(144,32,3,nt))

filename=trim(path)//"object_cloud_"//trim(run)//"_"//trim(int2str(yr,4))//".txt"
OPEN(11,file=trim(filename))

DO
READ(11,555,iostat=idum1) tt, ii, jj, nn, dum8, dum9
555 FORMAT(4I15,2F15.5)
IF (idum1/=0) GOTO 735

!IF (nn<10) GOTO 876

sur(ii,jj,3,tt)=1.
IF (dum8>=6000. .AND. dum9<=2000.) THEN
  sur(ii,jj,1,tt)=1.
ELSE
  sur(ii,jj,2,tt)=1.
ENDIF

876 CONTINUE

ENDDO
735 CONTINUE

write(*,*) "fine"

filename="for_your_obj/object_col_"//trim(run)//"_"//trim(int2str(yr,4))//".dat"
OPEN(10,FILE=trim(filename),ACCESS="direct",RECL=144*32*3*nt)
WRITE(10,REC=1) sur
CLOSE(10)
CLOSE(11)

WRITE(*,*) yr
DEALLOCATE(sur)
ENDDO

call mpi_end() 

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

