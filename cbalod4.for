C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CBALOD4
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
C **  SUBROUTINES CBALOD CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  CALCULATE MOMENTUM AND ENERGY DISSIPATION
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
C     UUEOUTO=UUEOUTO+SPB(L)*SPB(L+1)*DXYU(L)
C    &      *(U(L,1)*TBX(L)-U(L,KC)*TSX(L))
C     VVEOUTO=VVEOUTO+SPB(L)*SPB(LN)*DXYV(L)
C    &      *(V(L,1)*TBY(L)-V(L,KC)*TSY(L))
      UUEOUTO=UUEOUTO+0.5*SPB(L)*DXYP(L)*(U(L,1)*TBX(L)
     &       +U(L+1,1)*TBX(L+1)
     &       -U(L,KC)*TSX(L)-U(L+1,KC)*TSX(L+1))
      VVEOUTO=VVEOUTO+0.5*SPB(L)*DXYP(L)*(V(L,1)*TBY(L)+V(LN,1)*TBX(LN)
     &       -V(L,KC)*TSY(L)-V(LN,KC)*TSX(LN))
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      DUTMP=0.5*( U(L,K+1)+U(L+1,K+1)-U(L,K)-U(L+1,K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUTO=UUEOUTO+SPB(L)*2.0*DXYP(L)*AV(L,K)
     &      *( DUTMP*DUTMP )/(DZC(K+1)+DZC(K))
      VVEOUTO=VVEOUTO+SPB(L)*2.0*DXYP(L)*AV(L,K)
     &      *( DVTMP*DVTMP )/(DZC(K+1)+DZC(K))
      BBEOUTO=BBEOUTO+SCB(L)*DXYP(L)*HP(L)
     &      *GP*AB(L,K)*(B(L,K+1)-B(L,K))
      ENDDO
      ENDDO
C 
C**********************************************************************C
C
      RETURN
      END
