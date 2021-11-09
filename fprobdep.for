C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FPROBDEP(TAUD,TAUB)
C
C**********************************************************************C
C
C **  FPROBDEP CALCULATES PROBABILITY OF DEPOSITION USING PROBABILITY 
C **  INTEGRAL
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
C**********************************************************************C
C
C      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'EFDC.PAR'
C
C **  EVALUATION ASSUMES TAUB.GT.TAUD
C
      YVAL=2.04*LOG(0.25*((TAUB/TAUD)-1.)*EXP(1.27*TAUD))
	INEG=0
	IF(YVAL.LT.0.0)THEN
        INEG=1
        YVAL=ABS(YVAL)
      ENDIF
	XVAL=1.0/(1.0+0.3327*YVAL)
      POLYX=XVAL*(0.4632-0.1202*XVAL+0.9373*XVAL*XVAL)
      EXPY=-0.5*YVAL*YVAL
      FUNY=0.3989*EXP(EXPY)
      TMPVAL=1.0-FUNY*POLYX
      IF(INEG.EQ.1)TMPVAL=1.0-TMPVAL
	FPROBDEP=1.0-TMPVAL
C
C **  NOTES
C     0.3989=1/SQRT(2*PI)
C
C**********************************************************************C
C
      RETURN
      END
