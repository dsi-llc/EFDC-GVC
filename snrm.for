C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION SNRM(N,SX,ITOL)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C----------------------------------------------------------------------C
C
      INCLUDE 'EFDC.PAR'
      DIMENSION SX(LCM)
      IF(ITOL.LE.3)THEN
        SNRM=0.
        DO 11 I=1,N
          SNRM=SNRM+SX(I)**2
11      CONTINUE
        SNRM=SQRT(SNRM)
      ELSE
        ISAMAX=1
        DO 12 I=1,N
          IF(ABS(SX(I)).GT.ABS(SX(ISAMAX))) ISAMAX=I
12      CONTINUE
        SNRM=ABS(SX(ISAMAX))
      ENDIF
      RETURN
      END
