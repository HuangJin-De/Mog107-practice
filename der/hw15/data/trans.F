PROGRAM trans
IMPLICIT NONE

INTEGER :: i,j,k,m,n,o
INTEGER :: yr, mo, da, hr, mn
REAL :: dum1,dum2,dum3
CHARACTER(100) :: temp, nametc


OPEN(10,FILE="BT.txt")

OPEN(11,FILE="BT_1.txt")
OPEN(12,FILE="BT_2.txt")
OPEN(13,FILE="BT_3.txt")

DO

777 CONTINUE
READ(10,184,IOSTAT=o) temp
184 FORMAT(100A)

IF (o/=0) GOTO 888

IF (temp(1:3)=="bwp") WRITE(11,184) trim(temp)
IF (temp(1:3)=="bsh") WRITE(12,184) trim(temp)
IF (temp(1:3)=="bio") WRITE(13,184) trim(temp)

ENDDO
888 CONTINUE

CLOSE(10)

CLOSE(11)
CLOSE(12)
CLOSE(13)

END PROGRAM trans
