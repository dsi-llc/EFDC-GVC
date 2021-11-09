C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T1
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
      IF(NBAL.GT.1) RETURN
C
C**********************************************************************C
C
C **  INITIALIZE VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC ENERGY
C **  AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES
C
C----------------------------------------------------------------------C
C
cSO   OPEN(1,FILE='BAL2T.INT')
C
      VOLBEG=0.
      VOLBEG2T=0.
      WVOLBEG2T=0.
      BVOLBEG2T=0.
      SALBEG=0.
      DYEBEG=0.
      UMOBEG=0.
      VMOBEG=0.
      UUEBEG=0.
      VVEBEG=0.
      PPEBEG=0.
      BBEBEG=0.
C
      DYEBEG2T=0.0
      DO NS=1,NSED
        SEDBEG2T(NS)=0.0
        SEDBEG2TW(NS)=0.0
        SEDBEG2TB(NS)=0.0
      ENDDO
      DO NS=1,NSND
        SNDBEG2T(NS)=0.0
        SNDBEG2TW(NS)=0.0
        SNDBEG2TB(NS)=0.0
      ENDDO
      DO NT=1,NTOX
        TOXBEG2T(NT)=0.0
        TOXBEG2TW(NT)=0.0
        TOXBEG2TB(NT)=0.0
      ENDDO
C
      VOLOUT=0.
      WVOLOUT=0.
      BVOLOUT=0.
      SALOUT=0.
      DYEOUT=0.
      UMOOUT=0.
      VMOOUT=0.
      UUEOUT=0.
      VVEOUT=0.
      PPEOUT=0.
      BBEOUT=0.
      VOLMORPH2T=0.
C
      DYEOUT2T=0.0
      DO NS=1,NSED
        SEDOUT2T(NS)=0.0
        SEDFLUX2T(NS)=0.0
      ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS)=0.0
        SNDFLUX2T(NS)=0.0
        SBLOUT2T(NS)=0.0
        SNDFBL2T(NS)=0.0
      ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT)=0.0
        TOXFLUXW2T(NT)=0.0
        TOXFLUXB2T(NT)=0.0
        TADFLUX2T(NT)=0.0
        TOXBLB2T(NT)=0.0
        TOXFBL2T(NT)=0.0
      ENDDO
C
      DO L=2,LA
      LN=LNC(L)
      VOLBEG=VOLBEG+SPB(L)*DXYP(L)*HP(L)
      VOLBEG2T=VOLBEG2T+SPB(L)*DXYP(L)*HP(L)
      WVOLBEG2T=WVOLBEG2T+SPB(L)*DXYP(L)*HP(L)
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
      DYEBEG=DYEBEG+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
      DYEBEG2T=DYEBEG2T+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
      SALBEG=SALBEG+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(K)
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
C
      DO NS=1,NSED
        DO K=1,KC
        DO L=2,LA
          SEDBEG2T(NS)=SEDBEG2T(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SED(L,K,NS)
          SEDBEG2TW(NS)=SEDBEG2TW(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SED(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NS=1,NSND
        DO K=1,KC
        DO L=2,LA
          SNDBEG2T(NS)=SNDBEG2T(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SND(L,K,NS)
          SNDBEG2TW(NS)=SNDBEG2TW(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SND(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NT=1,NTOX
        DO K=1,KC
        DO L=2,LA
          TOXBEG2T(NT)=TOXBEG2T(NT)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*TOX(L,K,NT)
          TOXBEG2TW(NT)=TOXBEG2TW(NT)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*TOX(L,K,NT)
        ENDDO
        ENDDO
      ENDDO
C
      DO NS=1,NSED
        DO L=2,LA
        DO K=1,KBT(L)
          SEDBEG2T(NS)=SEDBEG2T(NS)+SCB(L)*DXYP(L)*SEDB(L,K,NS)
          SEDBEG2TB(NS)=SEDBEG2TB(NS)+SCB(L)*DXYP(L)*SEDB(L,K,NS)
        ENDDO
        ENDDO
       ENDDO
      DO NS=1,NSND
        DO L=2,LA
        DO K=1,KBT(L)
          SNDBEG2T(NS)=SNDBEG2T(NS)+SCB(L)*DXYP(L)*SNDB(L,K,NS)
          SNDBEG2TB(NS)=SNDBEG2TB(NS)+SCB(L)*DXYP(L)*SNDB(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NT=1,NTOX
        DO L=2,LA
        DO K=1,KBT(L)
          TOXBEG2T(NT)=TOXBEG2T(NT)+SCB(L)*DXYP(L)*TOXB(L,K,NT)
          TOXBEG2TB(NT)=TOXBEG2TB(NT)+SCB(L)*DXYP(L)*TOXB(L,K,NT)
        ENDDO
        ENDDO
      ENDDO
C
      KK=0
c      WRITE(1,*)KK,VOLBEG2T,WVOLBEG2T
      DO L=2,LA
      DO K=1,KBT(L)
      VOLBEG2T=VOLBEG2T+SPB(L)*DXYP(L)*HBED(L,K)
      BVOLBEG2T=BVOLBEG2T+SPB(L)*DXYP(L)*HBED(L,K)
cjmh     WVOLBEG2T=WVOLBEG2T+SPB(L)*DXYP(L)*PORBED(L,K)*HBED(L,K)
c      WRITE(1,*)KK,VOLBEG2T,WVOLBEG2T,PORBED(L,K)
      ENDDO
      ENDDO
C
      CLOSE(1)
C
C**********************************************************************C
C
      RETURN
      END
