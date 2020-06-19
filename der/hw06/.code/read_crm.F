PROGRAM readdata
USE netcdf
IMPLICIT NONE

INTEGER :: nx, ny, nz, nt
INTEGER :: yr, mn, da, myid,nproc, ts,te, dnt
INTEGER, DIMENSION(12) :: mn_da
REAL :: tarray(2), time, dtime
INTEGER :: i,j,k,m,n,o,t,t_out, ii,jj
INTEGER*8 :: nn
INTEGER, DIMENSION(49,24) :: ig,jg
REAL, PARAMETER :: tripi=4*atan(1.)
REAL :: dum1,dum2,dum3,dum4,dum5
REAL :: dum6,dum7,dum8,dum9
REAL :: rate
REAL :: polysvp1,qv_sat
INTEGER :: idum1,idum2,idum3,idum4,idum5
REAL, DIMENSION(101) :: h_re
REAL, DIMENSION(100) :: hc_re, warm
REAL, DIMENSION(100) :: cth
REAL, DIMENSION(32) :: latg
REAL, DIMENSION(144) :: long
REAL, DIMENSION(24) :: lat
REAL, DIMENSION(49) :: lon
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: obj
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sur, z, zt
REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: var, temp
REAL, DIMENSION(64,24) :: vol
INTEGER, DIMENSION(64,24) :: flag_qc
INTEGER, DIMENSION(64*24) :: flagcloud, top, bottom, tap
REAL, DIMENSION(64*24) :: sizecloud
INTEGER :: ierr,ncid1,varid1,ncid2,varid2,ncid3,varid3,ncid4,varid4
INTEGER :: dimid, access
INTEGER, DIMENSION(6) :: startmat, endmat
CHARACTER(300) :: path, run, filename, region, ncname, varname
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

path="/data/dadm1/model_output/SPCAM/"
run="sed_run_name"

ALLOCATE(z(144,32,24),zt(144,32,24))

yr=0001
filename=trim(path)//trim(run)//"/"//trim(run)//".cam.h0."//trim(int2str(yr,4))//"-01-01-00000.nc"
ierr = nf90_open(trim(filename),nf90_nowrite,ncid1)
ierr = nf90_inq_varid(ncid1,"lon",varid1)
ierr = nf90_get_var(ncid1,varid1,long,start=(/ 1 /),count=(/ 144 /))
ierr = nf90_inq_varid(ncid1,"lat",varid1)
ierr = nf90_get_var(ncid1,varid1,latg,start=(/ 33 /),count=(/ 32 /))
ierr = nf90_inq_varid(ncid1,"Z3",varid1)
ierr = nf90_get_var(ncid1,varid1,zt,start=(/ 1, 33, 3, 1 /),count=(/ 144, 32, 24, 1 /))
ierr = nf90_close(ncid1)

DO k=1,24
  z(:,:,25-k)=zt(:,:,k)-zt(:,:,24)
ENDDO

yr=0001
filename=trim(path)//trim(run)//"/"//trim(run)//".cam.h1."//trim(int2str(yr,4))//"-01-01-00000.nc"
ierr = nf90_open(trim(filename),nf90_nowrite,ncid1)
ierr = nf90_inq_varid(ncid1,"LON_60e_to_180e",varid1)
ierr = nf90_get_var(ncid1,varid1,lon,start=(/ 1 /),count=(/ 49 /))
ierr = nf90_inq_varid(ncid1,"LAT_15s_to_30n",varid1)
ierr = nf90_get_var(ncid1,varid1,lat,start=(/ 1 /),count=(/ 24 /))
ierr = nf90_close(ncid1)

DO j=1,24
DO i=1,49
  DO jj=1,32
    IF (lat(j)==latg(jj)) jg(i,j)=jj
  ENDDO
  DO ii=1,144
    IF (lon(i)==long(ii)) ig(i,j)=ii
  ENDDO
ENDDO
ENDDO

vol=1.

!ts=10
!te=10

year_loop: DO yr=ts,te

mn_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
nt=sum(mn_da)

startmat=(/ 1,1,1,1,1,1 /)
endmat=(/ 49,24,64,1,24,48 /)

ALLOCATE(var(49,24,64,24,48),temp(49,24,64,24,48))

filename="for_your_obj/object_cloud_"//trim(run)//"_"//trim(int2str(yr,4))//".txt"
OPEN(10,FILE=trim(filename))

nt=1
DO mn=1,12
DO da=1,mn_da(mn)
!time=dtime(tarray)

filename=trim(path)//trim(run)//"/"//trim(run)//".cam.h1."//trim(int2str(yr,4))&
       //"-"//trim(int2str(mn,2))//"-"//trim(int2str(da,2))//"-00000.nc"

ierr = nf90_open(trim(filename),nf90_nowrite,ncid1)

ierr = nf90_inq_varid(ncid1,"CRM_QC_LON_60e_to_180e_LAT_15s_to_30n",varid1)
ierr = nf90_get_var(ncid1,varid1,var,start=startmat,&
                   count=endmat)
ierr = nf90_close(ncid1)

!time=dtime(tarray)
DO j=1,24
DO i=1,49

DO t=1,48

  call six_connected_2d(64,24,var(i,j,:,:,t),vol,1.E-5 &
                   ,flag_qc,flagcloud,sizecloud,top,bottom)
  
  DO ii=1,64*24
    IF (flagcloud(ii) == -1) GOTO 987 
    WRITE(10,531) nt, ig(i,j), jg(i,j), int(sizecloud(ii))&
                    , z(ig(i,j),jg(i,j),top(ii))&
                    , z(ig(i,j),jg(i,j),bottom(ii))
    531 FORMAT(4I15,2F15.5)
  ENDDO
  987 CONTINUE

ENDDO

ENDDO
ENDDO
!time=dtime(tarray)
!write(*,*) time

nt=nt+1
ENDDO
WRITE(*,*) mn, myid
ENDDO

CLOSE(10)

DEALLOCATE(var)
!WRITE(*,*) yr

ENDDO year_loop


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

SUBROUTINE six_connected_2d(nx,ny,data_array,volume_array,criteria &
                        ,cloudflag_array,cloudflag_seq,cloudsize_seq&
                        ,top,bottom)
IMPLICIT NONE

! input
INTEGER, INTENT(in) :: nx, ny
REAL   , INTENT(in) :: criteria     ! criteria for defining objects
REAL   , DIMENSION(nx,ny), INTENT(in) :: data_array ! data array
REAL   , DIMENSION(nx,ny), INTENT(in) :: volume_array ! volume array

! output
INTEGER, DIMENSION(nx,ny), INTENT(inout) :: cloudflag_array ! cloud flag array
INTEGER, DIMENSION(nx*ny), INTENT(inout) :: cloudflag_seq ! cloud flag
REAL   , DIMENSION(nx*ny), INTENT(inout) :: cloudsize_seq ! cloud size
INTEGER, DIMENSION(nx*ny), INTENT(inout) :: top, bottom

! local variables
REAL :: con_size
INTEGER :: ct, cb
INTEGER :: i, j
INTEGER :: i_end, j_end
INTEGER :: nex, ney, myx, myy
INTEGER :: label, waiting_num, dir
INTEGER, DIMENSION(8), PARAMETER :: &
         move_x_3d(8) =  (/ -1,  0, +1,  0, -1, -1, +1, +1 /), &
         move_y_3d(8) =  (/  0, -1,  0, +1, -1, +1, -1, +1 /)
INTEGER, DIMENSION(2,nx*ny) :: waiting_grid
LOGICAL, DIMENSION(nx,ny) :: search_flag

! initialization
i_end = nx
j_end = ny

label           = 1
cloudflag_array = -1
cloudflag_seq   = -1
cloudsize_seq   = -1.
search_flag     = .FALSE.


! searching loop
j_loop: DO j=1,j_end
i_loop: DO i=1,i_end

  grid_searched_or_not: IF (.not.(search_flag(i,j)))THEN

  search_flag(i,j) = .TRUE.
  waiting_num = 1
  con_size    = 0.
  ct          = 1
  cb          = j_end

  waiting_grid(1,waiting_num) = i
  waiting_grid(2,waiting_num) = j

  searching_loop: DO

    IF (waiting_num < 1) GOTO 888

    myx = waiting_grid(1,waiting_num)
    myy = waiting_grid(2,waiting_num)

    waiting_num = waiting_num - 1

    IF (data_array(myx,myy) >= criteria )THEN
      cloudflag_array(myx,myy) = label
      con_size = con_size + 1. 
      ct = max(ct,myy)
      cb = min(cb,myy)

      check_neighbor: DO dir=1,8
        nex = myx + move_x_3d(dir)
        ney = myy + move_y_3d(dir)

        IF ( nex==0       ) nex = i_end
        IF ( nex==i_end+1 ) nex = 1

        IF ( nex>=1 .AND. nex<=i_end .AND. &
             ney>=1 .AND. ney<=j_end )THEN

          ne_search_or_not: IF (.not.(search_flag(nex,ney)))THEN

            waiting_num = waiting_num +1
            waiting_grid(1,waiting_num) = nex
            waiting_grid(2,waiting_num) = ney

            search_flag(nex,ney) = .TRUE.
          ENDIF ne_search_or_not
        ENDIF

      ENDDO check_neighbor

    ENDIF

  ENDDO searching_loop
  888 CONTINUE

  IF (con_size /= 0)THEN
    cloudflag_seq(label) = label
    cloudsize_seq(label) = con_size
    top(label) = ct
    bottom(label) = cb
    label = label + 1
  ENDIF
  ENDIF grid_searched_or_not
ENDDO i_loop
ENDDO j_loop

END SUBROUTINE six_connected_2d



