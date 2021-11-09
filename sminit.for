C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SMINIT
C
C**********************************************************************C
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
CXH      INSMICI=40
CXH      INSMRST=40
CXH      ISMORST=45
CXH      ISMOZB=46
C
      SMTSNAME(1) = 'SOM'
      SMTSNAME(2) = 'SIM'
      SMTSNAME(3) = 'SBF'
C
      DO L=2,LA
        SMHYST(L)=.FALSE.
      ENDDO
C
      RETURN
      END
