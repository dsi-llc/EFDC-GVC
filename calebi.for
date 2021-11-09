C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEBI
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
C **  CALEBI CALCULATES THE EXTERNAL BUOYANCY INTEGRALS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C     DIMENSION CH(LCM,KCM)
C
C**********************************************************************C
C
      DO L=2,LA
      BI1(L)=0.
      BI2(L)=0.
      CH(L,KC)=DZC(KC)*B(L,KC)
      BE(L)=GP*DZC(KC)*B(L,KC)
      ENDDO
C
      IF(KC.GT.1)THEN
      DO K=KS,1,-1
      DO L=2,LA
      CH(L,K)=CH(L,K+1)+DZC(K)*B(L,K)
      BE(L)=BE(L)+GP*DZC(K)*B(L,K)
      ENDDO
      ENDDO
	ENDIF
C
      DO K=1,KC
      DO L=2,LA
      BI1(L)=BI1(L)+GP*DZC(K)*(CH(L,K)-0.5*DZC(K)*B(L,K))
      BI2(L)=BI2(L)+GP*DZC(K)*(CH(L,K)+Z(K-1)*B(L,K))
      ENDDO
      ENDDO
C
C**********************************************************************C
C
      RETURN
      END
