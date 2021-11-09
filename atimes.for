C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE ATIMES(PCGM,PTMPM)
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION PCGM(LCM),PTMPM(LCM)
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PTMPM(L)=CC(L)*PCGM(L)+CS(L)*PCGM(LS)+CW(L)*PCGM(L-1)
     &       +CE(L)*PCGM(L+1)+CN(L)*PCGM(LN)
      ENDDO
      PTMPM(1)=0.
      PTMPM(LC)=0.
C
      RETURN
      END
