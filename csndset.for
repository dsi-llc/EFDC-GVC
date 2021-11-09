C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION CSNDSET(SND,SDEN,IOPT)
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
C **  CALCULATES HINDERED SETTLING CORRECTION FOR CLASS NS NONCOHESIVE 
C **  SEDIMENT 
C
      INCLUDE 'EFDC.PAR'
C
      ROPT=FLOAT(IOPT)
      CSNDSET=(1.-SDEN*SND)**ROPT
C
      RETURN
      END 
