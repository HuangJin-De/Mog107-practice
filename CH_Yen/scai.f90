PROGRAM scai
USE netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C.-H. Yen, 2020/04/12                                   !
! This program aims to                                    !
! 1) calculate SCAI within given boxes                    !
! Inputs:  cloud objects data (.txt)                      !
! Outputs: SCAI spatial and temporal distribution (.nc)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: i, ii, j, jj, t, tt, k, n, y, yy, h, w, t1, t2, idnum, eof, row, imin,imax,jmin,jmax
REAL :: radius, ci, cj, N2
INTEGER :: i_dimid, j_dimid, t_dimid
INTEGER, DIMENSION(3) :: dimids
INTEGER, DIMENSION(5) :: varids
INTEGER :: init_yr, finl_yr !1998!2015
CHARACTER(len=30) :: inpfile, outfile, ind_name='scai'
CHARACTER(len=4) :: year
REAL, DIMENSION(:), ALLOCATABLE :: Lat, Lon
REAL, PARAMETER :: wst=-180., est=180., sth=-50., nth=50., reso=0.25, dx=10., dy=10., LL=dx*100 ! km ! width of box
INTEGER, PARAMETER :: l = int(dx/reso)+int(mod(dx/reso,2.))-1 ! number of pixels within the width of box
                                     ! mod(dx/reso,2)-1 is to ensure l is odd
INTEGER, PARAMETER :: Nmax=(l-1)**2!, LL=dx*100 ! km
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Num
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Num2
REAL, DIMENSION(:,:,:), ALLOCATABLE :: D0, IND ! SCAI
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: CLD
!print*,Nmax
! calculate the dimensions of N, D0 given the box width
ii = int((est-wst)/reso)
jj = int((nth-sth)/reso)

tt = 366 ! daily data
ALLOCATE(Lon(ii))
ALLOCATE(Lat(jj))
ALLOCATE(Num(ii,jj,tt)); Num = 0
ALLOCATE(Num2(ii,jj,tt))
ALLOCATE(D0(ii,jj,tt)); D0 = 0
ALLOCATE(IND(ii,jj,tt)); IND = 0
ALLOCATE(CLD(ii,jj,Nmax,2)); CLD = 0

DO i = 1, ii
  Lon(i) = wst+i*reso
END DO
DO j = i, jj
  Lat(j) = sth+j*reso
END DO

PRINT*, 'first year:'
READ(*,*) init_yr
PRINT*, 'last year:'
READ(*,*) finl_yr
yy = finl_yr-init_yr+1

DO y = 1, yy
  WRITE(year, '(I4)') init_yr+y-1
  PRINT*, 'year=', year
  inpfile = 'new_cld_obj/cloud_'//year//'.txt'
  OPEN(10, FILE=inpfile, FORM='FORMATTED', STATUS='UNKNOWN')!, ACCESS='SEQUENTIAL')
  eof = 0
  t1  = 1; t = 1; row=1; imin=1;imax=ii;jmin=1;jmax=jj
  DO WHILE (eof == 0)
      READ(10, *, IOSTAT=eof) t2, idnum, radius, ci, cj
!print*, t2, idnum, radius,ci,cj
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (t2 == t1+8) THEN
        !print*,'finished collections of Num, CLD for time', t1 
        ! calculate SCAI for eiach box
        DO j = 1, jj
          DO i = 1, ii
            IF (Num(i,j,t) > 1) THEN
            ! only calculate the index of box with cloud objects
              N2 = float(Num(i,j,t)*(Num(i,j,t)-1)/2)
              D0(i,j,t) = 1
              DO n = 1, Num(i,j,t)-1
                DO k = n+1, Num(i,j,t)
                  D0(i,j,t) = D0(i,j,t)*(111.*sqrt((CLD(i,j,n,1)-CLD(i,j,k,1))**2.+(CLD(i,j,n,2)-CLD(i,j,k,2))**2.))**(1/N2)
                END DO ! ends DO k = n, N(i,j) 
              END DO ! ends DO n = 1, N(i,j)

!if(D0(i,j,t) > 1110) then
!print*,'--------------------------------------'
!D0(i,j,t)=1
!DO n = 1, Num(i,j,t)-1
!                DO k = n+1, Num(i,j,t)
!                  D0(i,j,t) =D0(i,j,t)*(111.*sqrt((CLD(i,j,n,1)-CLD(i,j,k,1))**2.+(CLD(i,j,n,2)-CLD(i,j,k,2))**2.))**(1./float(Num(i,j,t)))
!                print*,D0(i,j,t), float(Num(i,j,t)) !CLD(i,j,1:Num(i,j,t),1),CLD(i,j,1:Num(i,j,t),2)
!                  END DO ! ends DO k = n, N(i,j)
!              END DO
!endif

              D0(i,j,t) = D0(i,j,t)*111.
              IND(i,j,t) = float(Num(i,j,t))/float(Nmax)*D0(i,j,t)/LL*1000.
            END IF
          END DO ! ends DO i = 1, ii
        END DO ! ends DO j = 1, jj
        CLD = 0. ! clean for the collection of the next time
        t1 = t2
        t = t+1
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      i = int((ci-wst)/reso)+1
      j = int((cj-sth)/reso)+1 ! identify the pixel (i,j) which the object belongs to 
      ! add the data of the new object to the collection'
row=0
      DO h = j-(l-1)/2, j+(l-1)/2
        DO w = i-(l-1)/2, i+(l-1)/2
row=row+1
          IF (h > 0 .AND. w > 0 .AND. w <= ii .AND. h <= jj) THEN ! ensure not exceeds the edge of the domain
            Num(w,h,t) = Num(w,h,t)+1
            CLD(w,h,Num(w,h,t),1) = ci
            CLD(w,h,Num(w,h,t),2) = cj
          END IF
        END DO ! ends DO w = i-(l-1)/2, i+(l-1)/2
      END DO ! ends DO h = j-(l-1)/2, j+(l-1)/2
!print*,l,row ! 4/16 17:55
    IF (eof > 0) THEN
      PRINT*, 'Erorr occured while reading '//inpfile
      EXIT
    ELSE IF (eof < 0) THEN
      PRINT*, 'Successfully finished reading '//inpfile
      EXIT
    END IF
  END DO
  CLOSE(10) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Num2 = FLOAT(Num)
  print*,'here1'
  outfile = trim(ind_name)//'/'//trim(ind_name)//'_'//year//'.nc'
  CALL check(nf90_create(outfile, nf90_clobber, ncid))

  ! Define the dimensions. NetCDF will hand back an ID for each.
  CALL check(nf90_def_dim(ncid, 'latitude', jj, j_dimid)) ! not represent
  CALL check(nf90_def_dim(ncid, 'longitude', ii, i_dimid))
  CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid))

  dimids = (/ i_dimid, j_dimid, t_dimid /)
  print*, 'here2'
  ! Define the variable.
  CALL check(nf90_def_var(ncid, "longitude", nf90_double,  i_dimid, varids(1)))
  CALL check(nf90_def_var(ncid, "latitude", nf90_double,  j_dimid, varids(2)))
  CALL check(nf90_def_var(ncid, "number", nf90_double,  dimids, varids(3)))
  CALL check(nf90_def_var(ncid, "mean_distance", nf90_double,  dimids, varids(4)))
  CALL check(nf90_def_var(ncid, ind_name, nf90_double,  dimids, varids(5)))
  CALL check(nf90_enddef(ncid))
  print*, 'here3'
  CALL check(nf90_put_var(ncid, varids(1), Lon))
  print*,'here4'
  CALL check(nf90_put_var(ncid, varids(2), Lat))
  print*,'here5'
  CALL check(nf90_put_var(ncid, varids(3), Num2))
  CALL check(nf90_put_var(ncid, varids(4), D0))
  print*,'here6'
  CALL check(nf90_put_var(ncid, varids(5), IND))
  CALL check(nf90_close(ncid))
  PRINT *, 'SUCCESS writing ', outfile
END DO ! ends DO y = 1, yy

DEALLOCATE(Num)
DEALLOCATE(Num2)
DEALLOCATE(D0)
DEALLOCATE(IND)
DEALLOCATE(CLD)

CONTAINS
        SUBROUTINE check(STATUS)
                INTEGER, INTENT( in) :: STATUS
                IF (STATUS /= nf90_noerr) THEN
                        PRINT*, trim(nf90_strerror(STATUS))
                        STOP 2
                END IF
        END SUBROUTINE

END PROGRAM scai
