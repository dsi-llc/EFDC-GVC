C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE HYDZERO
C
C**********************************************************************C
C
C M. MORTON  07 JUN 1999
C INITIALIZES THE BINARY HYDRODYNAMIC ARRAYS TO 0.0
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
      DO K=1,KC
        DO LL=2,LA
C
C ZERO THE HYDRODYNAMIC VARIABLE ARRAYS:
C
          SELSUM(LL)  = 0.0
          DEPSUM(LL) = 0.0
          VXXSUM(LL) = 0.0
          VYYSUM(LL) = 0.0
          QXXSUM(LL) = 0.0
          QYYSUM(LL) = 0.0
        ENDDO
      ENDDO
      TIMEHYD = 0.0
      NHYCNT  = 0
C
      RETURN
      END
