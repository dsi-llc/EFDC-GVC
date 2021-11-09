C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQZERO4
C
C**********************************************************************C
C
C M. MORTON  02 JUN 1999
C INITIALIZES THE BENTHIC FLUX ARRAYS TO 0.0
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DO LL=2,LA
C
C ZERO THE BENTHIC FLUX ARRAYS:
C
        BFO2SUM(LL)  = 0.0
        BFNH4SUM(LL) = 0.0
        BFNO3SUM(LL) = 0.0
        BFPO4SUM(LL) = 0.0
        BFSADSUM(LL) = 0.0
        BFCODSUM(LL) = 0.0
        BFSMTSUM(LL) = 0.0
        BFBSTSUM(LL) = 0.0
      ENDDO
      TIMEBF = 0.0
      NBFCNT = 0
C
      RETURN
      END
