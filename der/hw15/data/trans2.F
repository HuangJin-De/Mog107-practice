PROGRAM trans
IMPLICIT NONE

INTEGER :: i,j,k,m,n,o
INTEGER :: yr, mo, da, hr, mn
REAL :: dum1,dum2,dum3
CHARACTER(100) :: temp, nametc


OPEN(10,FILE="hurdat2-nepac-1979-2015.txt")

OPEN(11,FILE="BT_5.txt")

DO

777 CONTINUE
READ(10,184,IOSTAT=o) temp
184 FORMAT(100A)

IF (o/=0) GOTO 888

IF (temp(1:1)=="E" .OR. temp(1:1)=="C") THEN
  nametc="cep"//temp(3:8)
  write(*,*)  nametc
  GOTO 777
ENDIF

READ(temp(1:4),*) yr
READ(temp(5:6),*) mo
READ(temp(7:8),*) da
READ(temp(11:12),*) hr
READ(temp(13:14),*) mn
READ(temp(23:27),*) dum1
READ(temp(31:35),*) dum2

IF(temp(28:28)=="S") dum1=-dum1
IF(temp(36:36)=="W") dum2=-dum2

write(11,185) trim(nametc),yr, mo, da, hr, dum2, dum1
185 FORMAT(A,I5.4,3I3.2,2F7.1)

ENDDO
888 CONTINUE

CLOSE(10)

END PROGRAM trans
