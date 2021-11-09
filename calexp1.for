C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP1 (ISTL)
C
C **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
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
C
C**********************************************************************C
C
       DELT=DT2
       DELTD2=DT
      IF(ISTL.EQ.2)THEN
       DELT=DT
       DELTD2=0.5*DT
      ENDIF
C
      DELTI=1./DELT
C
C**********************************************************************C
C
C **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
C
C----------------------------------------------------------------------C
C
      FCAXE(1)=0.
      FCAYE(1)=0.
      FXE(1)=0.
      FYE(1)=0.
C
      FCAXE(LC)=0.
      FCAYE(LC)=0.
      FXE(LC)=0.
      FYE(LC)=0.
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        FCAXE(L)=0.
        FCAYE(L)=0.
        FXE(L)=0.
        FYE(L)=0.
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2) GOTO 200
      IF(ISTL.EQ.3)THEN
       IF(ISCDMA.EQ.0) GOTO 300
       IF(ISCDMA.EQ.1) GOTO 400
       IF(ISCDMA.EQ.2) GOTO 350
      ENDIF
C
C**********************************************************************C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
C
  200 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         UHDY2(L,K)=UHDY1(L,K)+UHDY(L,K)
         VHDX2(L,K)=VHDX1(L,K)+VHDX(L,K)
         W2(L,K)=W1(L,K)+W(L,K)
        ENDDO
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1S(L,K)=UHDY2(LSC(L),K)
         TVAR1W(L,K)=VHDX2(L-1   ,K)
         TVAR1E(L,K)=UHDY2(L+1   ,K)
         TVAR1N(L,K)=VHDX2(LNC(L),K)
         TVAR2S(L,K)=U1(LSC(L),K)
         TVAR2W(L,K)=V1(L-1   ,K)
         TVAR2E(L,K)=U1(L+1   ,K)
         TVAR2N(L,K)=V1(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL          
         UHC=0.25*(UHDY2(L,K)+TVAR1S(L,K))
         UHB=0.25*(UHDY2(L,K)+TVAR1E(L,K))
         VHC=0.25*(VHDX2(L,K)+TVAR1W(L,K))
         VHB=0.25*(VHDX2(L,K)+TVAR1N(L,K))
         FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
     &            +MIN(UHB,0.)*TVAR2E(L,K)
         FVHU(L,K)=MAX(VHC,0.)*TVAR2S(L,K)
     &            +MIN(VHC,0.)*U1(L,K)
         FUHV(L,K)=MAX(UHC,0.)*TVAR2W(L,K)
     &            +MIN(UHC,0.)*V1(L,K)
         FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
     &            +MIN(VHB,0.)*TVAR2N(L,K)
        ENDDO
       ENDDO
      ENDDO
C
C      DO K=1,KC
C      DO L=2,LA
C      LN=LNC(L)
C      LS=LSC(L)               
C      UHC=0.25*(UHDY2(L,K)+UHDY2(LS,K))
C      UHB=0.25*(UHDY2(L,K)+UHDY2(L+1,K))
C      VHC=0.25*(VHDX2(L,K)+VHDX2(L-1,K))
C      VHB=0.25*(VHDX2(L,K)+VHDX2(LN,K))
C
C      FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
C     &         +MIN(UHB,0.)*U1(L+1,K)
C      FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
C     &         +MIN(VHC,0.)*U1(L,K)
C      FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
C     &         +MIN(UHC,0.)*V1(L,K)
C      FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
C     &         +MIN(VHB,0.)*V1(LN,K)
C      ENDDO
C      ENDDO
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1S(L,K)=W2(LSC(L),K)
         TVAR1W(L,K)=W2(L-1   ,K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         WU=0.25*DXYU(L)*(W2(L,K)+TVAR1W(L,K))
         WV=0.25*DXYV(L)*(W2(L,K)+TVAR1S(L,K))
         FWU(L,K)=MAX(WU,0.)*U1(L,K)
     &           +MIN(WU,0.)*U1(L,K+1)
         FWV(L,K)=MAX(WV,0.)*V1(L,K)
     &           +MIN(WV,0.)*V1(L,K+1)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
C     WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
C     FWU(L,K)=MAX(WU,0.)*U1(L,K)
C    &        +MIN(WU,0.)*U1(L,K+1)
C     FWV(L,K)=MAX(WV,0.)*V1(L,K)
C    &        +MIN(WV,0.)*V1(L,K+1)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      GOTO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1S(L,K)=UHDY(LSC(L),K)
         TVAR1W(L,K)=VHDX(L-1   ,K)
         TVAR1E(L,K)=UHDY(L+1   ,K)
         TVAR1N(L,K)=VHDX(LNC(L),K)
         TVAR2S(L,K)=U1(LSC(L),K)
         TVAR2W(L,K)=V1(L-1   ,K)
         TVAR2E(L,K)=U1(L+1   ,K)
         TVAR2N(L,K)=V1(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL          
         UHC=0.5*(UHDY(L,K)+TVAR1S(L,K))
         UHB=0.5*(UHDY(L,K)+TVAR1E(L,K))
         VHC=0.5*(VHDX(L,K)+TVAR1W(L,K))
         VHB=0.5*(VHDX(L,K)+TVAR1N(L,K))
         FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
     &            +MIN(UHB,0.)*TVAR2E(L,K)
         FVHU(L,K)=MAX(VHC,0.)*TVAR2S(L,K)
     &            +MIN(VHC,0.)*U1(L,K)
         FUHV(L,K)=MAX(UHC,0.)*TVAR2W(L,K)
     &            +MIN(UHC,0.)*V1(L,K)
         FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
     &            +MIN(VHB,0.)*TVAR2N(L,K)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)               
C     UHC=0.5*(UHDY(L,K)+UHDY(LS,K))
C     UHB=0.5*(UHDY(L,K)+UHDY(L+1,K))
C     VHC=0.5*(VHDX(L,K)+VHDX(L-1,K))
C     VHB=0.5*(VHDX(L,K)+VHDX(LN,K))
C     FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
C    &         +MIN(UHB,0.)*U1(L+1,K)
C     FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
C    &         +MIN(VHC,0.)*U1(L,K)
C     FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
C    &         +MIN(UHC,0.)*V1(L,K)
C     FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
C    &         +MIN(VHB,0.)*V1(LN,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1S(L,K)=W(LSC(L),K)
         TVAR1W(L,K)=W(L-1   ,K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         WU=0.5*DXYU(L)*(W(L,K)+TVAR1W(L,K))
         WV=0.5*DXYV(L)*(W(L,K)+TVAR1S(L,K))
         FWU(L,K)=MAX(WU,0.)*U1(L,K)
     &           +MIN(WU,0.)*U1(L,K+1)
         FWV(L,K)=MAX(WV,0.)*V1(L,K)
     &           +MIN(WV,0.)*V1(L,K+1)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
C     WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
C     FWU(L,K)=MAX(WU,0.)*U1(L,K)
C    &        +MIN(WU,0.)*U1(L,K+1)
C     FWV(L,K)=MAX(WV,0.)*V1(L,K)
C    &        +MIN(WV,0.)*V1(L,K+1)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      GOTO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  FIRST HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
C **  WITH TRANSPORT AT (N-1/2) AND TRANSPORTED FIELD AT (N-1)
C **  SECOND HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE  
C **  WITH TRANSPORT AT (N+1/2) AND TRANSPORTED FIELD AT (N)
C
  350 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      U2(L,K)=U1(L,K)+U(L,K)
      V2(L,K)=V1(L,K)+V(L,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)               
      UHC=0.25*(UHDY(L,K)+UHDY(LS,K))
      UHB=0.25*(UHDY(L,K)+UHDY(L+1,K))
      VHC=0.25*(VHDX(L,K)+VHDX(L-1,K))
      VHB=0.25*(VHDX(L,K)+VHDX(LN,K))
      FUHU(L,K)=MAX(UHB,0.)*U2(L,K)
     &         +MIN(UHB,0.)*U2(L+1,K)
      FVHU(L,K)=MAX(VHC,0.)*U2(LS,K)
     &         +MIN(VHC,0.)*U2(L,K)
      FUHV(L,K)=MAX(UHC,0.)*V2(L-1,K)
     &         +MIN(UHC,0.)*V2(L,K)
      FVHV(L,K)=MAX(VHB,0.)*V2(L,K)
     &         +MIN(VHB,0.)*V2(LN,K)
      ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)               
C     UHC=0.125*(UHDY(L,K)+UHDY(LS,K)+UHDY1(L,K)+UHDY1(LS,K))
C     UHB=0.125*(UHDY(L,K)+UHDY(L+1,K)+UHDY1(L,K)+UHDY1(L+1,K))
C     VHC=0.125*(VHDX(L,K)+VHDX(L-1,K)+VHDX1(L,K)+VHDX1(L-1,K))
C     VHB=0.125*(VHDX(L,K)+VHDX(LN,K)+VHDX1(L,K)+VHDX1(LN,K))
C     FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
C    &         +MIN(UHB,0.)*U1(L+1,K)
C     FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
C    &         +MIN(VHC,0.)*U1(L,K)
C     FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
C    &         +MIN(UHC,0.)*V1(L,K)
C     FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
C    &         +MIN(VHB,0.)*V1(LN,K)
C     ENDDO
C     ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)               
C     UHC=0.125*(3.*UHDY(L,K)+3.*UHDY(LS,K)-UHDY1(L,K)-UHDY1(LS,K))
C     UHB=0.125*(3.*UHDY(L,K)+3.*UHDY(L+1,K)-UHDY1(L,K)-UHDY1(L+1,K))
C     VHC=0.125*(3.*VHDX(L,K)+3.*VHDX(L-1,K)-VHDX1(L,K)-VHDX1(L-1,K))
C     VHB=0.125*(3.*VHDX(L,K)+3.*VHDX(LN,K)-VHDX1(L,K)-VHDX1(LN,K))
C     FUHU(L,K)=FUHU(L,K)+MAX(UHB,0.)*U(L,K)
C    &                   +MIN(UHB,0.)*U(L+1,K)
C     FVHU(L,K)=FVHU(L,K)+MAX(VHC,0.)*U(LS,K)
C    &                   +MIN(VHC,0.)*U(L,K)
C     FUHV(L,K)=FUHV(L,K)+MAX(UHC,0.)*V(L-1,K)
C    &                   +MIN(UHC,0.)*V(L,K)
C     FVHV(L,K)=FVHV(L,K)+MAX(VHB,0.)*V(L,K)
C    &                   +MIN(VHB,0.)*V(LN,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      WU=0.25*DXYU(L)*(W(L,K)+W(L-1,K))
      WV=0.25*DXYV(L)*(W(L,K)+W(LS,K))
      FWU(L,K)=MAX(WU,0.)*U2(L,K)
     &        +MIN(WU,0.)*U2(L,K+1)
      FWV(L,K)=MAX(WV,0.)*V2(L,K)
     &        +MIN(WV,0.)*V2(L,K+1)
      ENDDO
      ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     WU=0.125*DXYU(L)*(W(L,K)+W(L-1,K)+W1(L,K)+W1(L-1,K))
C     WV=0.125*DXYV(L)*(W(L,K)+W(LS,K)+W1(L,K)+W1(LS,K))
C     FWU(L,K)=MAX(WU,0.)*U1(L,K)
C    &        +MIN(WU,0.)*U1(L,K+1)
C     FWV(L,K)=MAX(WV,0.)*V1(L,K)
C    &        +MIN(WV,0.)*V1(L,K+1)
C     ENDDO
C     ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     WU=0.125*DXYU(L)*(3.*W(L,K)+3.*W(L-1,K)-W1(L,K)-W1(L-1,K))
C     WV=0.125*DXYV(L)*(3.*W(L,K)+3.*W(LS,K)-W1(L,K)-W1(LS,K))
C     FWU(L,K)=FWU(L,K)+MAX(WU,0.)*U(L,K)
C    &                 +MIN(WU,0.)*U(L,K+1)
C     FWV(L,K)=FWV(L,K)+MAX(WV,0.)*V(L,K)
C    &                 +MIN(WV,0.)*V(L,K+1)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      GOTO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT 
C **  AT (N) AND TRANSPORTED FIELD AT (N)
C
  400 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)  
      FUHU(L,K)=0.25*(UHDY(L+1,K)+UHDY(L,K))*(U(L+1,K)+U(L,K))
      FVHU(L,K)=0.25*(VHDX(L,K)+VHDX(L-1,K))*(U(L,K)+U(LS,K))
      FUHV(L,K)=0.25*(UHDY(L,K)+UHDY(LS,K))*(V(L,K)+V(L-1,K))
      FVHV(L,K)=0.25*(VHDX(L,K)+VHDX(LN,K))*(V(LN,K)+V(L,K))
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)    
      FWU(L,K)=0.25*DXYU(L)*(W(L,K)+W(L-1,K))*(U(L,K+1)+U(L,K))
      FWV(L,K)=0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
      ENDDO
      ENDDO
C
C**********************************************************************C
C
  500 CONTINUE
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FUHU(L,K)=STCUV(L)*FUHU(L,K)
         FVHV(L,K)=STCUV(L)*FVHV(L,K)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=1,LC
C     FUHU(L,K)=STCUV(L)*FUHU(L,K)
C     FVHV(L,K)=STCUV(L)*FVHV(L,K)
C     ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      IF(ISDCCA.EQ.0)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3E(L)=DYU(L+1   )
        TVAR3N(L)=DXV(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1E(L,K)=U(L+1   ,K)
         TVAR1N(L,K)=V(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(TVAR1N(L,K)+V(L,K))*(TVAR3E(L)-DYU(L))
     &        -0.5*SNLT*(TVAR1E(L,K)+U(L,K))*(TVAR3N(L)-DXV(L)) )*HP(L)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)          
C     CAC(L,K)=( FCORC(L)*DXYP(L)
C    &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C    &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
C     ENDDO
C     ENDDO
C
      ELSE
C
      CFMAX=CF
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3E(L)=DYU(L+1   )
        TVAR3N(L)=DXV(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1E(L,K)=U(L+1   ,K)
         TVAR1N(L,K)=V(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(TVAR1N(L,K)+V(L,K))*(TVAR3E(L)-DYU(L))
     &        -0.5*SNLT*(TVAR1E(L,K)+U(L,K))*(TVAR3N(L)-DXV(L)) )*HP(L)
         CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
         CFMAX=MAX(CFMAX,CFEFF)
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)          
C     CAC(L,K)=( FCORC(L)*DXYP(L)
C    &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C    &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
C     CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
C     CFMAX=MAX(CFMAX,CFEFF)
C     ENDDO
C     ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR2S(L,K)=FVHV(LSC(L),K)
         TVAR2W(L,K)=FUHU(L-1   ,K)
         TVAR2E(L,K)=FUHV(L+1   ,K)
         TVAR2N(L,K)=FVHU(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=SAAX(L)*(FUHU(L,K)-TVAR2W(L,K)+TVAR2N(L,K)-FVHU(L,K))
         FY(L,K)=SAAY(L)*(TVAR2E(L,K)-FUHV(L,K)+FVHV(L,K)-TVAR2S(L,K))
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1S(L,K)=U(LSC(L),K)
         TVAR1W(L,K)=V(L-1   ,K)
C         TVAR1E(L,K)=U(L+1   ,K)
C         TVAR1N(L,K)=V(LNC(L),K)
         TVAR2S(L,K)=CAC(LSC(L),K)
         TVAR2W(L,K)=CAC(L-1   ,K)
         TVAR2E(L,K)=U(LSEC(L),K)
         TVAR2N(L,K)=V(LNWC(L),K)
        ENDDO
       ENDDO
      ENDDO
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FCAX(L,K)=ROLD*FCAX(L,K)
     &      +0.25*RNEW*SCAX(L)*(CAC(L,K)*(TVAR1N(L,K)+V(L,K))
     &                     +TVAR2W(L,K)*(TVAR2N(L,K)+TVAR1W(L,K)))
         FCAY(L,K)=ROLD*FCAY(L,K)
     &      +0.25*RNEW*SCAY(L)*(CAC(L,K)*(TVAR1E(L,K)+U(L,K))
     &                     +TVAR2S(L,K)*(TVAR2E(L,K)+TVAR1S(L,K)))
        ENDDO
       ENDDO
      ENDDO
C  
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)      
C     LNW=LNWC(L)
C     LSE=LSEC(L)
C     FCAX(L,K)=ROLD*FCAX(L,K)
C    &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
C    &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
C     FCAY(L,K)=ROLD*FCAY(L,K)
C    &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
C    &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
C     FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K))
C     FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K))
C     ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISHDMF.GE.1)THEN
C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         TVAR1S(L,K)=FMDVY(LSC(L),K)
         TVAR1W(L,K)=FMDUX(L-1   ,K)
         TVAR1E(L,K)=FMDVX(L+1   ,K)
         TVAR1N(L,K)=FMDUY(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)
     &          -SDX(L)*(FMDUX(L,K)-TVAR1W(L,K)+TVAR1N(L,K)-FMDUY(L,K))
         FY(L,K)=FY(L,K)
     &          -SDY(L)*(TVAR1E(L,K)-FMDVX(L,K)+FMDVY(L,K)-TVAR1S(L,K))
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     LN=LNC(L)
C     FX(L,K)=FX(L,K)
C    &        -SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
C     FY(L,K)=FY(L,K)
C    &        -SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
C     ENDDO
C     ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.2)THEN
C
      IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
       ELSE
        WVFACT=1.0
      ENDIF        
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FXE(L)=FXE(L)+WVFACT*SAAX(L)*FXWAVE(L,K)*DZC(K)
         FYE(L)=FYE(L)+WVFACT*SAAY(L)*FYWAVE(L,K)*DZC(K)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        ENDDO
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
         FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
C **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3S(L)=HMP(LSC(L))
        TVAR3W(L)=HMP(L-1   )
        TVAR3E(L)=HP(L-1   )
        TVAR3N(L)=HP(LSC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1S(L,K)=B(LSC(L),K+1)
         TVAR1W(L,K)=B(L-1   ,K+1)
         TVAR2S(L,K)=B(LSC(L),K)
         TVAR2W(L,K)=B(L-1   ,K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*
     &            ( HU(L)*( (B(L,K+1)-TVAR1W(L,K))*DZC(K+1)
     &                     +(B(L,K)-TVAR2W(L,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+TVAR1W(L,K)-TVAR2W(L,K))*
     &            (HMP(L)-TVAR3W(L)-Z(K)*(HP(L)-TVAR3E(L))) )
         FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*
     &            ( HV(L)*( (B(L,K+1)-TVAR1S(L,K))*DZC(K+1)
     &                     +(B(L,K)-TVAR2S(L,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+TVAR1S(L,K)-TVAR2S(L,K))*
     &            (HMP(L)-TVAR3S(L)-Z(K)*(HP(L)-TVAR3N(L))) )
        ENDDO
       ENDDO
      ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*
C    &            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
C    &                     +(B(L,K)-B(L-1,K))*DZC(K) )
C    &           +(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
C    &            (HMP(L)-HMP(L-1)-Z(K)*(HP(L)-HP(L-1))) )
C     FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*
C    &            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
C    &                     +(B(L,K)-B(LS,K))*DZC(K) )
C    &           +(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
C    &            (HMP(L)-HMP(LS)-Z(K)*(HP(L)-HP(LS))) )
C     ENDDO
C     ENDDO
C
C     IF(N.EQ.1)THEN
C       OPEN(1,FILE='BUOY.DIA',STATUS='UNKNOWN')
C       DO L=2,LA
C        DO K=1,KS
C        TMP3D(K)=SUBO(L)*FBBX(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C        DO K=1,KS
C        TMP3D(K)=SVBO(L)*FBBY(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C       ENDDO
C       CLOSE(1)
C     ENDIF
C
 1111 FORMAT(2I5,2X,8E12.4)       
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        RCDZF=CDZF(K)
        DO L=LF,LL
         DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
         DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
        ENDDO
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.2)THEN
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      RETURN
      END
