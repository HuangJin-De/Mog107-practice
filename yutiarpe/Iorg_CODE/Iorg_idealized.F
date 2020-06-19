PROGRAM TRMM_IorgM
!
!   Program Procedures:
!
!    Main program: 
!       O  Read data
!       O  Use Iorg subroutine to get Iorg, object counts, average prp.
!       O  Write output
!       O  The above procedures are executed for each time step
!       O  This program use mpi, each cpu calculates different time step
!       O  This program is built to read TRMM 3-hourly data for a specific year, 
!          the default domain is 160E~230E, 30S~30N with no subdomains
!
!    Parameter description:
!  
!    NX, NY     : grid counts of the input data in x, y direction
!    NX1, NY1   : grid counts of the subdomain in x, y direction, equal to NX, NY 
!                 if you don't have subdomain
!    int_lon,   
!    int_lat    : grid length in degree
!    LIMFAC     : a factor used to cut off array size
!    tval       : object threshold [mm/hr]
!    tmax       : the required maximum value for an object [mm/hr]
!    r          : earth radius [m]
!    pi     　  : pi
!    d2r        : transfer degree to radiance
!    NT         : total time steps in the input data
!    start_lon,
!    start_lat  : If there are subdomains, the starting longitude/latitude for
!                 each subdomain will be different.
!    org        : array for saving Iorg
!    obj        : array for saving object counts
!    prec       : array for saving average precipitation in the domain 
!    pcp        : array for input data
USE netcdf
IMPLICIT NONE
INTEGER,PARAMETER                 :: NX=280, NY=240, NX1=280, NY1=240
REAL,PARAMETER                    :: int_lon=0.25, int_lat=0.25
INTEGER,PARAMETER                 :: LIMFAC=2, CLDN=2
REAL,PARAMETER                    :: tval=1, tmax=2.5
REAL,PARAMETER                    :: r=6371000,pi=3.14159265359, d2r=pi/180
!==========================================================================
INTEGER                           :: NT, i, j, l
INTEGER                           :: ncid, varid, err, myid, nproc, ff
REAL                              :: start_lon, start_lat
REAL,DIMENSION(NX/NX1,NY/NY1)     :: org, obj, prec
CHARACTER(len=100)                :: path,filen
CHARACTER(len=4)                  :: fyear,leap(4)
REAL,DIMENSION(NX,NY)             :: pcp
DATA leap/'2012','2008','2004','2000'/
call mpi_start(myid,nproc)
 ff=myid+1
call getarg(1,fyear)

 path='/data/dadm1/obs/TRMM/TRMM3B42/3B42.'
 filen=trim(path)//fyear//'.3hr.nc'
 if (any(leap .eq. fyear)) then
  NT=2928
 else
  NT=2920
 endif

do l=(ff-1)*NT/nproc+1,ff*NT/nproc,8 !MPI
 !Read file
 err = nf90_open(trim(filen),NF90_NOWRITE,ncid)
 IF (err/=nf90_noerr) WRITE(*,*) "open fail"
 err = nf90_inq_varid(ncid,"pcp",varid)
 IF (err/=nf90_noerr) WRITE(*,*) "inq_varid fail"
 err = nf90_get_var(ncid,varid,pcp(201:280,:),start=(/ 1361,81,l /),count=(/80,NY,1/))
 err = nf90_get_var(ncid,varid,pcp(1:200,:),start=(/ 1,81,l /),count=(/200,NY,1/))
 ! data in the .nc file is 180W~180E
 err = nf90_close(ncid)
 !Subdomains
 org=0.
 obj=0.
 start_lat=-29.875
 do j=1,NY/NY1
  start_lon=160.125
  do i=1,NX/NX1
   call Iorg(pcp((i-1)*NX1+1:i*NX1,(j-1)*NY1+1:j*NY1), start_lon, start_lat, org(i,j), &
             obj(i,j), prec(i,j))
   start_lon=start_lon+int_lon*NX1
  enddo
  start_lat=start_lat+int_lat*NY1
 enddo
 !write output
 open(unit=200,FILE="Iorg_TRMM_"//fyear//".dat",FORM='UNFORMATTED',STATUS='UNKNOWN',access='direct',&
      recl=3*NX/NX1*NY/NY1)
 write(200,rec=(l-1)/8+1) org(:,:),obj(:,:),prec(:,:)
 close(200)
enddo

call mpi_end()
CONTAINS

SUBROUTINE Iorg(pcp, start_lon, start_lat, org, obj ,prec)
REAL,INTENT(in)                     :: pcp(NX1,NY1)
REAL,INTENT(in)                     :: start_lon, start_lat
REAL,INTENT(out)                    :: org, obj, prec
REAL                                :: dum, rand, area
REAL                                :: clonm, clonn, clatm, clatn, radm, radn, mind
INTEGER                             :: i, j, l, n, m, p, label, k, iter
INTEGER                             :: sortloc(1), opsr, opsl, opsu, opsd
INTEGER                             :: lenx, leny, dum1, putloc(2)
INTEGER                             :: puti, putj, rad100
INTEGER,DIMENSION(NX1,NY1)          :: label_data, filt, filtn, &
                                       edge
INTEGER,DIMENSION(NX1*NY1/LIMFAC)   :: sort
REAL,DIMENSION(NX1,NY1)             :: refill_grid, lon, lat, nobound,weight,pcp2
REAL,DIMENSION(NX1*NY1,CLDN)        :: cld_inc_size
REAL,DIMENSION(NX1,NY1)             :: denoma, dista
REAL,DIMENSION(NX1*NY1/LIMFAC)      :: obsr
REAL,DIMENSION(NX1*NY1/LIMFAC)      :: nncdfx,nncdfy
!==========================================================================
!There are missing values in pcp array, but they doesn't affect objects
!identification.
pcp2=pcp
where (pcp .le. 0) pcp2=0
do j=1,NY1
 denoma(:,j)=cos((start_lat+(j-1)*int_lat)*d2r)
 lat(:,j)=start_lat+(j-1)*int_lat
enddo
do i=1,NX1
 lon(i,:)=start_lon+(i-1)*int_lon
enddo
!Calculate the area (m**2) of the domain
area=sum(denoma*110.574*int_lat*111.32*int_lon)*1000000.!km**2->m**2

!Calculate the meridional weighted average precipitation intensity of the domain
weight=pcp2*denoma
prec=sum(weight)/sum(denoma)


 !pcp          :input field for finding objects
 !label_data   :output field of object index
 !cld_inc_size :information of the object
 !label        :label-1 equals total object counts
 !denoma       :input field for calculating object sizes
 call find_cloud ( pcp(:,:) , label_data , cld_inc_size , label , denoma)
 !clear boundary and MAX----------------------------------
 !call clear_nmax if you set tmax
 call clear_nmax (label_data, cld_inc_size(:,2))
 !call clear_bound if you want to remove objects at the domain boundaries
 !call clear_bound (label_data)
 nobound=label_data
 !call find_cloud again if you call either clear_nmax or clear_bound
 call find_cloud (  nobound , label_data , &
                    cld_inc_size , label , denoma)

 !leave the subroutine if the object counts is less than 3
 obj=label-1
 if (label-1 .lt. 3) then
 org=0
 nncdfx=0
 nncdfy=0
 return
 endif

 filt=0
 where(label_data .gt. 0) filt=1

!m loop = find Nearest neighbor distance for each objects
!obsr : the array for saving NN distances of every objects
 do m=1,label-1
  mind=99999999999
  filtn=filt
  edge=0
  where (label_data .eq. m) 
   filtn=0
   edge=1
  endwhere

  do j=1,NY1
   do i=1,NX1
    if (edge(i,j) .ne. 1) cycle
    dista=asin(sqrt(sin(d2r*(lat-lat(i,j))/2)**2 + cos(d2r*lat(i,j))&
          *cos(d2r*lat) * sin(d2r*(lon-lon(i,j))/2)**2)) * 2 * r !m
    dum=minval(dista,mask=filtn .eq. 1)
    mind=min(mind,dum)
   enddo
  enddo
    obsr(m)=mind
 enddo!m

!sort the obsr array, the shortest NN distance will be obsr(1)
call hpsort(label-1,obsr(1:label-1))

nncdfx=0
nncdfy=0


!nncdfy=observed cdf of NN distances
!nncdfx=idealized cdf of NN distances
do i=1,label-1
 nncdfy(i)=float(i)/float(label-1) !observed
 nncdfx(i)=1.-exp(-sqrt(float(NX1*NY1))*(obsr(i)**2)/area) !idealized 
enddo

!org=Iorg
   org=0
   org=org+nncdfx(1)*nncdfy(1)/2
   do i=1,label-2
    org=org+(nncdfy(i)+nncdfy(i+1))*(nncdfx(i+1)-nncdfx(i))/2
   enddo 
    org=org+1-nncdfx(label-1)

end subroutine
    subroutine find_cloud( varin  , label_data , cld_inc_size ,label ,  denoma)
        real    , intent (in)    :: varin(NX1,NY1),denoma(NX1,NY1)
        real    , intent (out)   :: cld_inc_size(NX1*NY1,CLDN)
        Integer         , intent (out)   :: label_data(NX1,NY1)
        Integer         , intent (out)   :: label
        real                        :: inc, pmax ! add
        integer                     :: x , y ,  stack_ptr , CONNECT=4
        integer                     :: cu_move , nex , ney , myx , myy
        integer                     :: search_flag(NX1,NY1)
        integer                     :: stackX(NX1*NY1) , stackY(NX1*NY1)
        integer                     :: move_x_2d(4)  = (/ -1,  0, +1,  0 /),&
                                       move_y_2d(4)  = (/  0, -1,  0, +1 /)
        label = 1
        label_data(:,:)  = 0
        search_flag(:,:) = 0
        stack_ptr = 0
        stackX(:) = 0
        stackY(:) = 0
        do y = 1 , NY1
            do x =  1 , NX1
                    if ( search_flag(x,y) .EQ. 1 ) cycle
                    stackX(1) = x
                    stackY(1) = y
                    search_flag(x,y) = 1
                    stack_ptr = 1
                    inc = 0.
                    pmax = 0.
                    do while ( stack_ptr .GE. 1 )
                        myx =  stackX(stack_ptr)
                        myy =  stackY(stack_ptr)
                        stack_ptr = stack_ptr - 1
                        if ( varin(myx,myy)    .LT. tval ) cycle 
                        label_data(myx,myy) = label
                        inc = inc + 110.574*int_lat &
                                  *  111.32*int_lon * denoma(myx,myy)
                        pmax=max(pmax,varin(myx,myy))
                        do cu_move = 1 , CONNECT
                                nex = myx + move_x_2d( cu_move )
                                ney = myy + move_y_2d( cu_move )
!x-periodic                   if (nex .gt. NX) nex=nex-NX
!x-periodic                   if (nex .lt. 1)  nex=NX-nex
!x-periodic                   if (ney .GT. NY  ) cycle
!x-periodic                   if (ney .LT. 1   ) cycle
                            if ( (nex .gt. NX1) .or. (ney .gt. NY1)) cycle
                            if ( (nex .lt.  1) .or. (ney .lt.  1)) cycle

                            if ( search_flag(nex,ney) .EQ. 1 ) cycle
                            stackX(stack_ptr + 1) = nex
                            stackY(stack_ptr + 1) = ney
                            search_flag(nex,ney) = 1
                            stack_ptr = stack_ptr + 1
                        end do
                    end do
                     if ( sqrt(inc) .GE. 1E-8  ) then
                         cld_inc_size( label , 1) = sqrt(inc/pi)
                         cld_inc_size( label,  2) = pmax
                        label  = label + 1
                        inc    = 0.0
                    endif
            enddo
        enddo
    end subroutine


    subroutine clear_nmax( label_data , cld_inc_size)
        INTEGER , INTENT (INOUT) :: label_data(NX1,NY1)
        REAL    , INTENT (IN)    :: cld_inc_size(NX1*NY1,1)
        INTEGER                  :: x,y,temp

        do y = 1, NY1
         do x = 1, NX1
          if ( label_data(x,y) .EQ. 0 ) cycle
          if (cld_inc_size( label_data(x,y) , 1 ) .lt. tmax) then
           temp=label_data(x,y)
           where (label_data .eq. temp)  label_data=0
          endif
         enddo
        enddo
    end subroutine
    subroutine clear_bound( label_data )
      INTEGER,INTENT(INOUT) :: label_data(NX1,NY1)
      INTEGER               :: i,j,num
        do i=1,NX1
         if (label_data(i,1) .gt. 0) then
          num=label_data(i,1)
          where (label_data .eq. num) label_data=0
         endif
        enddo
        do i=1,NX1
         if (label_data(i,NY1) .gt. 0) then
          num=label_data(i,NY1)
          where (label_data .eq. num) label_data=0
         endif
        enddo
        do j=1,NY1
         if (label_data(1,j) .gt. 0) then
          num=label_data(1,j)
          where (label_data .eq. num) label_data=0
         endif
        enddo
        do j=1,NY1
         if (label_data(NX1,j) .gt. 0) then
          num=label_data(NX1,j)
          where (label_data .eq. num) label_data=0
         endif
        enddo
     end subroutine
subroutine HPSORT(N,RA)
IMPLICIT NONE
integer::N,IR,J,L,I
real ::RA(N),RRA
L=N/2+1
IR=N
10 continue
 if (L>1)then
  L=L-1
  RRA=RA(L)
 else
  RRA=RA(IR)
  RA(IR)=RA(1)
  IR=IR-1
  if(IR.eq.1)then
   RA(1)=RRA
   return
  end if
 endif
 I=L
 J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1)) J=J+1
  endif
  if(RRA < RA(J))then
   RA(I)=RA(J)
   I=J; J=J+J
  else
   J=IR+1
  endif
  goto 20
  endif
  RA(I)=RRA
  goto 10
end subroutine
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
end program
