Program aggindex
integer,parameter ::varn=5
real ::dxy=0.25
integer,parameter :: max_num_of_cloud=10000
real, dimension(max_num_of_cloud,varn) :: property, prop
real :: nowt,dis
real :: COP, SCAI, MCAI, D1, D2
real :: Nmax,L
integer :: endid, year,totalcloud
integer :: i, j, k, indt, indbx, indby, indrec,indyear
real,dimension(2) :: bound_x,bound_y
character(len=100) :: casen,datapath,dumc,outpath
real,dimension(varn) :: firstprop
real, dimension(100) :: boundsx, boundsy
real :: resolution,nboundx, nboundy
real,allocatable,dimension(:,:,:) :: output_arr
integer, dimension(4),parameter :: leaps=(/2000,2004,2008,2012/)
integer,dimension(18) :: years
years=(/(i,i=1998,2015)/)
datapath='/data/cloud/yuyu/hw2/TRMMmask'
outpath='/data/cloud/yuyu/hw2/TRMMagg'
casen='3hr'
resolution=30
nboundx=360/resolution
nboundy=(30-(-30))/resolution
allocate(output_arr(int(nboundx),int(nboundy),3))
boundsx(:nboundx+1)=(/(j*10,j=0,int(nboundx))/)
boundsy(:nboundy+1)=(/(-30+j*10, j=0,int(nboundy))/)
print*, nboundx,nboundy

do indyear=1,18
year=years(indyear)
nt=2920
if (any(year.eq.leaps)) nt=2928

write(dumc,'(A,I2,A,I2,A,I4,A)') trim(outpath)//'/'//trim(casen)//'/',int(resolution),'/TRMM_AGGindex_'//trim(casen)//'_',int(resolution),'_',year,'.dat'
open(66,file=dumc,form='unformatted',access='direct',recl=3*nboundx*nboundy)
indrec=1
write(dumc,'(A,I4,A)') trim(datapath)//'/TRMM_'//trim(casen)//'_mask'//'/TRMM_'//trim(casen)//'_',year,'_prop.txt'
open(10,file=dumc,status='old',form='formatted')
do indt=1,nt
totalcloud=0
! read prop.txt-------------------------------
! prop(:,varns)--> id, radius, loc_x, loc_y, precipitation
if(indt==1)then
 j=1
else
 property(1,:)=firstprop(:)
 j=2
endif
nowt=indt
do while(nowt .eq. indt)
 read(10,'(6F)',end=100) nowt,property(j,:)
 !write(*,'(6F)') nowt,id(j),radius(j),loc_x(j),loc_y(j),preci(j)
 !if(prop(j,3)<bound_x(2) .and. prop(j,3)>bound_x(1) .and. prop(j,4)<bound_y(2) .and. prop(j,4)>bound_y(1) .and. (nowt==indt))then
 j=j+1
enddo
100 continue

totalcloud=j-2
firstprop(:)=property(totalcloud+1,:)
property(totalcloud+1,:)=0.
!print*,totalcloud
!print*,(property(j,:),j=1,totalcloud)
!---------------------------------------------


do indby=1,int(nboundy)
do indbx=1,int(nboundx)

!bound_x=(/120,130/)
!bound_y=(/0,10/)
bound_x=(/boundsx(indbx), boundsx(indbx+1)/)
bound_y=(/boundsy(indby), boundsy(indby+1)/)
!print*,bound_x,bound_y



! maskout cloud in specific area
endid=0
do j=1,totalcloud
 if(property(j,3)<bound_x(2) .and. property(j,3)>bound_x(1) .and. &
    property(j,4)<bound_y(2) .and. property(j,4)>bound_y(1))then
    endid=endid+1
    prop(endid,:)=property(j,:)
 endif
enddo
!print*, endid
!print*, (prop(j,:),j=1,endid)


!---------------------------------------------
! calculate aggregation indexs
L=(bound_x(2)-bound_x(1))*100.
Nmax=(bound_x(2)-bound_x(1))/dxy*(bound_y(2)-bound_y(1))/dxy
k=0
D1=1
D2=0
COP=0
if(endid<=1)then
COP=-999000000
SCAI=-999000000
MCAI=-999000000
goto 999
endif

do i=1,endid
 do j=i,endid
 if(i==j) cycle
 call cloud_distance(prop(i,3),prop(i,4),prop(j,3),prop(j,4),dis)
 D1=D1*dis/L
 if((prop(i,2)+prop(j,2))<=dis)then
  D2=D2+dis-(prop(i,2)+prop(j,2))
  COP=COP+(prop(i,2)+prop(j,2))/dis
 else
  D2=D2+0.
  COP=COP+1.
 endif
 k=k+1
 enddo !j
enddo !i
!print*, Nmax, endid
D1=L*D1**(2/(endid*(endid-1)))
D2=D2/((endid*(endid-1))/2)
COP=COP/((endid*(endid-1))/2)
SCAI=endid/Nmax*D1/L*1000
MCAI=endid/Nmax*D2/L*1000

999 continue ! from line 87
!print*,year, nt, COP, SCAI, MCAI, endid,Nmax,k
output_arr(indbx,indby,1)=SCAI
output_arr(indbx,indby,2)=MCAI
output_arr(indbx,indby,3)=COP

print*,year,indt,indbx,indby

enddo !indbx
enddo !indby
write(66,rec=indt) output_arr


enddo !indt
enddo !indyear

!contains
endprogram aggindex

subroutine cloud_distance(lo1,la1,lo2,la2,dis)
real :: lon1,lat1,lon2,lat2
real,intent(in) :: lo1,la1,lo2,la2
real,intent(out) :: dis
real,parameter :: deg2rad=3.1415/180.
! prop(:,varns)--> id, radius, loc_x, loc_y, precipitation

lon1=lo1*deg2rad
lon2=lo2*deg2rad
lat1=la1*deg2rad
lat2=la2*deg2rad
dlon=lon2-lon1
dlat=lat2-lat1
dis=sin(dlat/2)**2+COS(lat1)*COS(lat2)*sin(dlon/2)**2
dis=2*ASIN(sqrt(dis))*6371

end subroutine


