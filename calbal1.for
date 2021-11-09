C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALBAL1
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
      IF(NBAL.GT.1) RETURN
C
C**********************************************************************C
C
C **  INITIALIZE VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC ENERGY
C **  AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES
C
C----------------------------------------------------------------------C
C
      VOLBEG=0.
      SALBEG=0.
      DYEBEG=0.
      UMOBEG=0.
      VMOBEG=0.
      UUEBEG=0.
      VVEBEG=0.
      PPEBEG=0.
      BBEBEG=0.
C
      VOLOUT=0.
      SALOUT=0.
      DYEOUT=0.
      UMOOUT=0.
      VMOOUT=0.
      UUEOUT=0.
      VVEOUT=0.
      PPEOUT=0.
      BBEOUT=0.
C
      DO L=2,LA
      LN=LNC(L)
      VOLBEG=VOLBEG+SPB(L)*DXYP(L)*HP(L)
      UMOBEG=UMOBEG+SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)
     &                                 +DYIU(L+1)*HUI(L+1)*UHDYE(L+1))
      VMOBEG=VMOBEG+SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)
     &                                 +DXIV(LN)*HVI(LN)*VHDXE(LN))
      PPEBEG=PPEBEG+SPB(L)*0.5*DXYP(L)
     &             *(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
      ENDDO
C
      AMOBEG=SQRT(UMOBEG*UMOBEG+VMOBEG*VMOBEG)
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      SALBEG=SALBEG+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(K)
      DYEBEG=DYEBEG+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
C     UUEBEG=UUEBEG+SPB(L)*0.25*(DXYU(L)*HU(L)*U(L,K)*U(L,K)
C    &      +DXYU(L+1)*HU(L+1)*U(L+1,K)*U(L+1,K))*DZC(K)
C     VVEBEG=VVEBEG+SPB(L)*0.25*(DXYV(L)*HV(L)*V(L,K)*V(L,K)
C    &      +DXYV(LN)*HV(LN)*V(LN,K)*V(LN,K))*DZC(K)
      UUEBEG=UUEBEG+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     &      *( (U(L,K)+U(L+1,K))*(U(L,K)+U(L+1,K)) )
      VVEBEG=VVEBEG+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     &      *( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEBEG=BBEBEG+SPB(L)*GP*DXYP(L)*HP(L)*DZC(K)*( BELV(L) 
     &      +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
      RETURN
      END
