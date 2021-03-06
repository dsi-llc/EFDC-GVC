C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T4
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
C 05/01/2002        john hamrick       05/01/2002       john hamrick
C  subroutine added for 2 time-level balances including sed,snd,tox
C----------------------------------------------------------------------C
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      END IF
C
C**********************************************************************C
C
C **  CALCULATE MOMENTUM AND ENERGY DISSIPATION
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
C     UUEOUT=UUEOUT+SPB(L)*SPB(L+1)*DXYU(L)
C    &      *(U(L,1)*TBX(L)-U(L,KC)*TSX(L))
C     VVEOUT=VVEOUT+SPB(L)*SPB(LN)*DXYV(L)
C    &      *(V(L,1)*TBY(L)-V(L,KC)*TSY(L))
      UUEOUT=UUEOUT+0.5*DELT*SPB(L)*DXYP(L)*(U(L,1)*TBX(L)
     &       +U(L+1,1)*TBX(L+1)-U(L,KC)*TSX(L)-U(L+1,KC)*TSX(L+1))
      VVEOUT=VVEOUT+0.5*DELT*SPB(L)*DXYP(L)*(V(L,1)*TBY(L)
     &       +V(LN,1)*TBX(LN)-V(L,KC)*TSY(L)-V(LN,KC)*TSX(LN))
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      DUTMP=0.5*( U(L,K+1)+U(L+1,K+1)-U(L,K)-U(L+1,K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUT=UUEOUT+DELT*SPB(L)*2.0*DXYP(L)*AV(L,K)
     &      *( DUTMP*DUTMP )/(DZC(K+1)+DZC(K))
      VVEOUT=VVEOUT+DELT*SPB(L)*2.0*DXYP(L)*AV(L,K)
     &      *( DVTMP*DVTMP )/(DZC(K+1)+DZC(K))
      BBEOUT=BBEOUT+DELT*SCB(L)*DXYP(L)*HP(L)
     &      *GP*AB(L,K)*(B(L,K+1)-B(L,K))
      ENDDO
      ENDDO
C 
C**********************************************************************C
C
      RETURN
      END
