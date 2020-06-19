PROGRAM find_peak
IMPLICIT NONE

INTEGER :: i,j,k,m,n,o
REAL :: dum1, dum2, dum3, dum4
REAL :: dum5, dum6, dum7, dum8
REAL, DIMENSION(100,3) :: crh, era
REAL, DIMENSION(100,6,3) :: crh_var, era_var
CHARACTER(30) model(24)
CHARACTER(200) path

path="/data/der0318/work/cmip6/DATA/"

model( 1)="ACCESS-CM2"
model( 2)="AWI-ESM-1-1-LR"
model( 3)="BCC-CSM2-MR"
model( 4)="BCC-ESM1"
model( 5)="CanESM5"
model( 6)="CESM2"
model( 7)="CESM2-FV2"
model( 8)="CESM2-WACCM"
model( 9)="CESM2-WACCM-FV2"
model(10)="EC-Earth3"
model(11)="GFDL-CM4"
model(12)="GFDL-ESM4"
model(13)="INM-CM4-8"
model(14)="INM-CM5-0"
model(15)="IPSL-CM6A-LR"
model(16)="MIROC6"
model(17)="MPI-ESM-1-2-H"
model(18)="MPI-ESM1-2-HR"
model(19)="MPI-ESM1-2-LR"
model(20)="MRI-ESM2-0"
model(21)="NorESM2-LM"
model(22)="NorESM2-MM"
model(23)="TaiESM1"
model(24)="ERA-I"

OPEN(10,FILE=trim(path)//trim(model(24))//"_pdf_data.dat",ACCESS="direct",RECL=100*19)
READ(10,REC=1) era(:,1), era_var(:,:,1)
READ(10,REC=2) era(:,2), era_var(:,:,2)
READ(10,REC=3) era(:,3), era_var(:,:,3)
CLOSE(10)

dum3=real(maxloc(era(:,2),1))
dum4=real(maxloc(era(:,3),1))

dum5=sum(era(:,2),1)
dum6=sum(era(1:60,2),1)
dum7=sum(era(61:100,3),1)

dum8=real(minloc(abs(era(41:80,2)-era(41:80,3)),1))+40.

OPEN(11,FILE=trim(path)//"multimodel_crhpeak.dat",ACCESS="direct",RECL=6)
m=1
WRITE(11,REC=m) dum4-dum3, maxval(era(:,2),1)-maxval(era(:,3),1)&
               ,dum5, dum6, dum7, dum8
m=m+1

DO n=1,23

OPEN(10,FILE=trim(path)//trim(model(n))//"_pdf_data.dat",ACCESS="direct",RECL=100*19)
READ(10,REC=1) crh(:,1), crh_var(:,:,1)
READ(10,REC=2) crh(:,2), crh_var(:,:,2)
READ(10,REC=3) crh(:,3), crh_var(:,:,3)
CLOSE(10)

dum1=real(maxloc(crh(:,2),1))
dum2=real(maxloc(crh(:,3),1))

dum5=sum(crh(:,2),1)
dum6=sum(crh(1:60,2),1)
dum7=sum(crh(61:100,3),1)

dum8=real(minloc(abs(crh(41:80,2)-crh(41:80,3)),1))+40.

WRITE(11,REC=m) dum2-dum1, maxval(crh(:,2),1)-maxval(crh(:,3),1)&
               ,dum5, dum6, dum7, dum8
m=m+1
!write(*,*) (model(n)), dum1, dum2

ENDDO

END PROGRAM find_peak
