C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP3 (ISTL)
C
C **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
C **  EXPERIMENTAL VERSION OF SMOLARKIEWCZ-MARGOLING SCHEME
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
      DO K=1,KC
       TVAR1W(1,K)=0.
       TVAR1S(1,K)=0.
       TVAR1E(1,K)=0.
       TVAR1N(1,K)=0.
       TVAR2W(1,K)=0.
       TVAR2S(1,K)=0.
       TVAR2E(1,K)=0.
       TVAR2N(1,K)=0.
       TVAR1W(LC,K)=0.
       TVAR1S(LC,K)=0.
       TVAR1E(LC,K)=0.
       TVAR1N(LC,K)=0.
       TVAR2W(LC,K)=0.
       TVAR2S(LC,K)=0.
       TVAR2E(LC,K)=0.
       TVAR2N(LC,K)=0.
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       TVAR1W(L,K)=0.
       TVAR1S(L,K)=0.
       TVAR1E(L,K)=0.
       TVAR1N(L,K)=0.
       TVAR2W(L,K)=0.
       TVAR2S(L,K)=0.
       TVAR2E(L,K)=0.
       TVAR2N(L,K)=0.
      ENDDO
      ENDDO
C
      TVAR3W(1)=0.
      TVAR3S(1)=0.
      TVAR3E(1)=0.
      TVAR3N(1)=0.
      TVAR3W(LC)=0.
      TVAR3S(LC)=0.
      TVAR3E(LC)=0.
      TVAR3N(LC)=0.
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3W(L)=0.
        TVAR3S(L)=0.
        TVAR3E(L)=0.
        TVAR3N(L)=0.
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE TIME LEVEL N-1 OR TIME LEVEL N INTERNAL STRESSES
C **  CALCULATE FRACTIONAL STEP WITH STRESSES AND CORIOLIS AND 
C **  CURVATURE ACCELERATIONS
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR2W(L,1)=TBX1(L)
        TVAR2S(L,1)=TBY1(L)
        TVAR2E(L,KC)=TSX1(L)
        TVAR2N(L,KC)=TSY1(L)
       ENDDO
      ENDDO
C
      DO K=1,KS
       TMPVAL=DZIG(K)
       DO L=2,LA
        TVAR2E(L,K)=TMPVAL*(U1(L,K+1)-U1(L,K))/AVUI(L,K)
        TVAR2N(L,K)=TMPVAL*(V1(L,K+1)-V1(L,K))/AVVI(L,K)
       ENDDO
      ENDDO
C
      DO K=2,KC
       DO L=2,LA
        TVAR2W(L,K)=TVAR2E(L,K-1)
        TVAR2S(L,K)=TVAR2N(L,K-1)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=DYIU(L)*UHDY1(L,K)
        TVAR1S(L,K)=DXIV(L)*VHDX1(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       TMPVAL=DELTD2*DZIC(K)
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)+TMPVAL*(TVAR2E(L,K)-TVAR2W(L,K))
        TVAR1S(L,K)=TVAR1S(L,K)+TMPVAL*(TVAR2N(L,K)-TVAR2S(L,K))
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)+DELTD2*DXYIU(L)*FCAX1(L,K)
        TVAR1S(L,K)=TVAR1S(L,K)-DELTD2*DXYIV(L)*FCAY1(L,K)
       ENDDO
      ENDDO
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
      WVFACT=DELTD2*WVFACT      
C
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)+WVFACT*DXYIU(L)*SAAX(L)*FXWAVE(L,K)
        TVAR1S(L,K)=TVAR1S(L,K)+WVFACT*DXYIV(L)*SAAY(L)*FYWAVE(L,K)
       ENDDO
      ENDDO
C
      ENDIF
C
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=SUB(L)*TVAR1W(L,K)
        TVAR1S(L,K)=SVB(L)*TVAR1S(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADD PRESSURE GRADIENT TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR1E(L,KC)=DZC(KC)*B1(L,KC)
       ENDDO
      ENDDO
C
      DO K=KS,1,-1
       DO L=2,LA
        TVAR1E(L,K)=TVAR1E(L,K+1)+DZC(K)*B1(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1E(L,K)=TVAR1E(L,K)-0.5*DZC(K)*B1(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1N(L,K)=TVAR1E(L,K)+0.5*(Z(K)+Z(K-1))*B1(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1N(L,K)=GP*TVAR1N(L,K)
        TVAR1E(L,K)=GP*TVAR1E(L,K)
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3W(L)=DELTD2*SUB(L)*H1U(L)*DXIU(L)
        TVAR3S(L)=DELTD2*SVB(L)*H1V(L)*DYIV(L)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        TVAR1W(L,K)=TVAR1W(L,K)-TVAR3W(L)*( P1(L)-P1(L-1) )              
        TVAR1S(L,K)=TVAR1S(L,K)-TVAR3S(L)*( P1(L)-P1(LS ) )              
       ENDDO
      ENDDO
C
C      SBX(L)=0.5*SBX(L)*DYU(L)
C      SBY(L)=0.5*SBY(L)*DXV(L)
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3W(L)=DELT*SUB(L)*H1U(L)*DXYIU(L)*SBX(L)
        TVAR3S(L)=DELT*SVB(L)*H1V(L)*DXYIV(L)*SBY(L)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        TVAR1W(L,K)=TVAR1W(L,K)-TVAR3W(L)*( 
     &          +H1U(L)*(TVAR1E(L,K)-TVAR1E(L-1,K))              
     &          +0.5*GP*(B1(L,K)+B1(L-1,K))*(BELV(L)-BELV(L-1))              
     &          +0.5*(TVAR1N(L,K)+TVAR1N(L-1,K))*(H1P(L)-H1P(L-1)) )              
        TVAR1S(L,K)=TVAR1S(L,K)-TVAR3S(L)*( 
     &          +H1V(L)*(TVAR1E(L,K)-TVAR1E(LS ,K))              
     &          +0.5*GP*(B1(L,K)+B1(LS ,K))*(BELV(L)-BELV(LS ))              
     &          +0.5*(TVAR1N(L,K)+TVAR1N(LS ,K))*(H1P(L)-H1P(LS )) )              
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADJUST INTERNAL TRANSPORT AT OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
C **  FIRST HALF STEP IS COMPLETE TVAR1W = (UH)*, TVAR1S = (VH)*
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
C **  ADVECTION STEP FOR X MOMENTUM EQUATION
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       UHDY2(1,K)=0.
       VHDX2(1,K)=0.
       UHDY2(LC,K)=0.
       VHDX2(LC,K)=0.
       FUHU(1,K)=0.
       FUHU(LC,K)=0.
       FVHU(1,K)=0.
       FVHU(LC,K)=0.
      ENDDO
C
      DO K=0,KC
       W2(1,K)=0. 
       FWU(1,K)=0. 
       W2(LC,K)=0. 
       FWU(LC,K)=0. 
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        W2(L,0)=0. 
        FWU(L,0)=0. 
       ENDDO
      ENDDO
C     
      DO K=1,KC
      DO L=2,LC
       UHDY2(L,K)=0.
       VHDX2(L,K)=0.
       W2(L,K)=0. 
       FUHU(L,K)=0.
       FVHU(L,K)=0.
       FWU(L,K)=0. 
      ENDDO
      ENDDO
C     
      IF(ISTL.EQ.2) GOTO 200
      GOTO 300
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
C **  UHDY2(L,K) IS DYU*U ON POS X U(L) MOM CV FACE (AT B OR P POINT)
C **  VHDX2(L,K) IS DXV*V ON NEG Y U(L) MOM CV FACE (AT CORNER POINT)
C **  W2(L,K)    IS DXYP*W/H ON TOP U(L) MOM CV FACE
C
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=0.25*STCUV(L)*( DYU(L  )*(U1(L  ,K)+U(L  ,K))
     &                              +DYU(L+1)*(U1(L+1,K)+U(L+1,K)) )
        VHDX2(L,K)=0.25*( DXV(L  )*(V1(L  ,K)+V(L  ,K))
     &                   +DXV(L-1)*(V1(L-1,K)+V(L-1,K)) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        W2(L,K)=0.25*( DXYP(L  )*(W1(L  ,K)/H1P(L  )+W(L  ,K)/HP(L  ))
     &                +DXYP(L-1)*(W1(L-1,K)/H1P(L-1)+W(L-1,K)/HP(L-1)) )
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        FUHU(L,K)=MAX(UHDY2(L,K),0.)*TVAR1W(L     ,K)
     &           +MIN(UHDY2(L,K),0.)*TVAR1W(L+1   ,K)
        FVHU(L,K)=MAX(VHDX2(L,K),0.)*TVAR1W(LSC(L),K)
     &           +MIN(VHDX2(L,K),0.)*TVAR1W(L     ,K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWU(L,K)=MAX(W2(L,K),0.)*TVAR1W(L,K)
     &          +MIN(W2(L,K),0.)*TVAR1W(L,K+1)
       ENDDO
      ENDDO
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
C **  UHDY2(L,K) IS DYU*U ON POS X U(L) MOM CV FACE (AT B OR P POINT)
C **  VHDX2(L,K) IS DXV*V ON NEG Y U(L) MOM CV FACE (AT CORNER POINT)
C **  W2(L,K)    IS DXYP*W/H ON TOP U(L) MOM CV FACE
C
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=0.5*STCUV(L)*( DYU(L)*U1(L,K)+DYU(L+1)*U1(L+1,K) )
        VHDX2(L,K)=0.5*( DXV(L)*V1(L,K)+DXV(L-1)*V1(L-1,K) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        W2(L,K)=0.5*( DXYP(L  )*W1(L  ,K)/H1P(L  )
     &               +DXYP(L-1)*W1(L-1,K)/H1P(L-1) )
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        FUHU(L,K)=MAX(UHDY2(L,K),0.)*TVAR1W(L     ,K)
     &           +MIN(UHDY2(L,K),0.)*TVAR1W(L+1   ,K)
        FVHU(L,K)=MAX(VHDX2(L,K),0.)*TVAR1W(LSC(L),K)
     &           +MIN(VHDX2(L,K),0.)*TVAR1W(L     ,K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWU(L,K)=MAX(W2(L,K),0.)*TVAR1W(L,K)
     &          +MIN(W2(L,K),0.)*TVAR1W(L,K+1)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
  500 CONTINUE
C
C     DO K=1,KC
C      DO L=1,LA
C       FUHU(L,K)=STCUV(L-1)*FUHU(L,K)
C       FVHV(L,K)=STCUV(L)*FVHV(L,K)
C      ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  ADVANCE U ADVECTION A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        FX(L,K)=FUHU(L,K)-FUHU(L-1,K)+FVHU(LNC(L),K)-FVHU(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        FX(L,K)=FX(L,K)+(FWU(L,K)-FWU(L,K-1))*DZIC(K)
       ENDDO
      ENDDO
C
      TMPVAL=DELT*SNLT
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)-TMPVAL*SUB(L)*DXYIU(L)*SAAX(L)*FX(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ANTI-DIFFUSE U ADVECTION FOR A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.LE.3) GOTO 800
C
C----------------------------------------------------------------------C
C
C **  CALCULATE U CONTROL VOLUME DIVERGENCE, TVAR1E 
C
      DO K=1,KC
       DO L=2,LA
        TVAR1E(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=0,KC
       DO L=2,LA
        TVAR2C(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=1,KC
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR1E(L,K)=UHDY2(L     ,K)-UHDY2(L-1,K)
     &             +VHDX2(LNC(L),K)-VHDX2(L  ,K)            
     &             +(W2(L,K)-W2(L,K-1))*RDZIG             
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1E(L,K)=0.5*DELT*DXYIU(L)*TVAR1E(L,K)          
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATED X MOMENTUM EQUTION TERMS U*DUDX, V*DUDY, AND W*DUDZ
C
      DO K=1,KC
       DO L=2,LA
        UUU(L,K)=UHDY2(L,K)*(TVAR1W(L+1,K)-TVAR1W(L     ,K))
        VVV(L,K)=VHDX2(L,K)*(TVAR1W(L  ,K)-TVAR1W(LSC(L),K))
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        WWW(L,K)=RDZIG*W2(L,K)*(TVAR1W(L,K+1)-TVAR1W(L,K))         
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE X MOMENTUM EQUTION ANTIDIFFUSIVE VELOCITIES
C
      VELSMAL=1.E-12
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=ABS(UHDY2(L,K))
     &                  *(TVAR1W(L+1,K)-TVAR1W(L,K))
     &                  /(TVAR1W(L+1,K)+TVAR1W(L,K)+VELSMAL)
        TVAR2S(L,K)=ABS(VHDX2(L,K))
     &                  *(TVAR1W(L,K)-TVAR1W(LSC(L),K))
     &                  /(TVAR1W(L,K)+TVAR1W(LSC(L),K)+VELSMAL)
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=RDZIG*ABS(W2(L,K))
     &                       *(TVAR1W(L,K+1)-TVAR1W(L,K))         
     &                       /(TVAR1W(L,K+1)+TVAR1W(L,K)+VELSMAL)         
       ENDDO
      ENDDO
C
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=TVAR2W(L,K)-0.5*UHDY2(L,K)
     &                             *(TVAR1E(L,K)+TVAR1E(L+1   ,K))
        TVAR2S(L,K)=TVAR2S(L,K)-0.5*VHDX2(L,K)
     &                             *(TVAR1E(L,K)+TVAR1E(LSC(L),K))
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=TVAR2C(L,K)-0.5*W2(L,K)
     &                             *(TVAR1E(L,K+1)+TVAR1E(L,K))         
       ENDDO
      ENDDO
C
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=TVAR2W(L,K)-DELT*UHDY2(L,K)
     &   *( UUU(L,K)
     &     +0.25*(VVV(L,K)+VVV(L+1,K)+VVV(LNC(L),K)+VVV(LNEC(L),K))
     &     +0.25*(WWW(L,K)+WWW(L+1,K)+WWW(L,K+1)+WWW(L+1,K+1)) )
     &     /(DXYU(L)*TVAR1W(L,K)+DXYU(L+1)*TVAR1W(L+1,K)+VELSMAL)
        TVAR2S(L,K)=TVAR2S(L,K)-DELT*VHDX2(L,K)
     &   *( VVV(L,K)
     &   +0.25*(UUU(L,K)+UUU(L-1,K)+UUU(LSC(L),K)+UUU(LSWC(L),K))
     &   +0.25*(WWW(L,K)+WWW(LSC(L),K)+WWW(L,K+1)+WWW(LSC(L),K+1)) )
     &   /(DXYU(L)*TVAR1W(L,K)+DXYU(LSC(L))*TVAR1W(LSC(L),K)+VELSMAL)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=STCUV(L)*SPB(L)*TVAR2W(L,K)
C
        TVAR2S(L,K)=SVB(L)*SVB(L-1)*TVAR2S(L,K)
C
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=TVAR2C(L,K)-DELT*W2(L,K)
     &   *( WWW(L,K)
     &   +0.25*(UUU(L,K)+UUU(L-1,K)+UUU(L,K+1)+UUU(L-1,K+1))
     &   +0.25*(VVV(LNC(L),K)+VVV(L,K)+VVV(LNC(L),K+1)+VVV(L,K+1)) )
     &   /(DXYU(L)*TVAR1W(L,K)+DXYU(L)*TVAR1W(L,K+1)+VELSMAL)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCLUATE U EQUATION ANTIDIFFUSIVE FLUXES
C
      DO K=1,KC
       DO L=2,LA
        FUHU(L,K)=MAX(TVAR2W(L,K),0.)*TVAR1W(L     ,K)
     &           +MIN(TVAR2W(L,K),0.)*TVAR1W(L+1   ,K)
        FVHU(L,K)=MAX(TVAR2S(L,K),0.)*TVAR1W(LSC(L),K)
     &           +MIN(TVAR2S(L,K),0.)*TVAR1W(L     ,K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWU(L,K)=MAX(TVAR2C(L,K),0.)*TVAR1W(L,K)
     &          +MIN(TVAR2C(L,K),0.)*TVAR1W(L,K+1)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADVANCE U ANTIDIFFUSIVE ADVECTION A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       FX(L,K)=FUHU(L,K)-FUHU(L-1,K)+FVHU(LNC(L),K)-FVHU(L,K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
       FX(L,K)=FX(L,K)+(FWU(L,K)-FWU(L,K-1))*DZIC(K)
       ENDDO
      ENDDO
C
      TMPVAL=DELT*SNLT
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)-TMPVAL*SUB(L)*SAAX(L)*DXYIU(L)*FX(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADVECTION STEP FOR Y MOMENTUM EQUATION
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
  800 CONTINUE
C
      DO K=1,KC
       UHDY2(1,K)=0.
       VHDX2(1,K)=0.
       UHDY2(LC,K)=0.
       VHDX2(LC,K)=0.
       FUHV(1,K)=0.
       FUHV(LC,K)=0.
       FVHV(1,K)=0.
       FVHV(LC,K)=0.
      ENDDO
C
      DO K=0,KC
       W2(1,K)=0. 
       FWV(1,K)=0. 
       W2(LC,K)=0. 
       FWV(LC,K)=0. 
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        W2(L,0)=0. 
        FWV(L,0)=0. 
       ENDDO
      ENDDO
C     
      DO K=1,KC
      DO L=2,LC
       UHDY2(L,K)=0.
       VHDX2(L,K)=0.
       W2(L,K)=0. 
       FUHV(L,K)=0.
       FVHV(L,K)=0.
       FWV(L,K)=0. 
      ENDDO
      ENDDO
C     
      IF(ISTL.EQ.2) GOTO 201
      GOTO 301
C
C**********************************************************************C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
C
  201 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  UHDY2(L,K) IS DYU*U ON NEG X V(L) MOM CV FACE (CORNER POINT)
C **  VHDX2(L,K) IS DXV*V ON POS Y V(L) MOM CV FACE (B 0R P POINT)
C **  W2(L,K)    IS DXYP*W/H ON TOP V MOM CV FACE
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        LN=LNC(L)
        UHDY2(L,K)=0.25*( DYU(L )*(U1(L ,K)+U(L ,K))
     &                   +DYU(LS)*(U1(LS,K)+U(LS,K)) )
        VHDX2(L,K)=0.25*STCUV(L)*( DXV(L )*(V1(L ,K)+V(L ,K))
     &                            +DXV(LN)*(V1(LN,K)+V(LN,K)) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        LS=LSC(L)
        W2(L,K)=0.25*( DXYP(L )*(W1(L ,K)/H1P(L )+W(L ,K)/HP(L ))
     &                +DXYP(LS)*(W1(LS,K)/H1P(LS)+W(LS,K)/HP(LS)) )
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        FUHV(L,K)=MAX(UHDY2(L,K),0.)*TVAR1S(L-1   ,K)
     &           +MIN(UHDY2(L,K),0.)*TVAR1S(L     ,K)
        FVHV(L,K)=MAX(VHDX2(L,K),0.)*TVAR1S(L     ,K)
     &           +MIN(VHDX2(L,K),0.)*TVAR1S(LNC(L),K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWV(L,K)=MAX(W2(L,K),0.)*TVAR1S(L,K)
     &          +MIN(W2(L,K),0.)*TVAR1S(L,K+1)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      GOTO 501
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
C
  301 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  UHDY2(L,K) IS DYU*U ON NEG X V(L) MOM CV FACE (CORNER POINT)
C **  VHDX2(L,K) IS DXV*V ON POS Y V(L) MOM CV FACE (B 0R P POINT)
C **  W2(L,K)    IS DXYP*W/H ON TOP V MOM CV FACE
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        LN=LNC(L)
        UHDY2(L,K)=0.5*( DYU(L)*U1(L,K)+DYU(LS)*U1(LS,K) )
        VHDX2(L,K)=0.5*STCUV(L)*( DXV(L)*V1(L,K)+DXV(LN)*V1(LN,K) )
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        LS=LSC(L)
        W2(L,K)=0.5*( DXYP(L )*W1(L ,K)/H1P(L )
     &               +DXYP(LS)*W1(LS,K)/H1P(LS) )
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        FUHV(L,K)=MAX(UHDY2(L,K),0.)*TVAR1S(L-1   ,K)
     &           +MIN(UHDY2(L,K),0.)*TVAR1S(L     ,K)
        FVHV(L,K)=MAX(VHDX2(L,K),0.)*TVAR1S(L     ,K)
     &           +MIN(VHDX2(L,K),0.)*TVAR1S(LNC(L),K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWV(L,K)=MAX(W2(L,K),0.)*TVAR1S(L,K)
     &          +MIN(W2(L,K),0.)*TVAR1S(L,K+1)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
  501 CONTINUE
C
C     DO K=1,KC
C      DO L=1,LA
C       FUHU(L,K)=STCUV(L)*FUHU(L,K)
C       FVHV(L,K)=STCUV(L)*FVHV(L,K)
C      ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  ADVANCE V ADVECTION A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       FY(L,K)=FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LSC(L),K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
       FY(L,K)=FY(L,K)+(FWV(L,K)-FWV(L,K-1))*DZIC(K)
       ENDDO
      ENDDO
C
      TMPVAL=DELT*SNLT
      DO K=1,KC
       DO L=2,LA
        TVAR1S(L,K)=TVAR1S(L,K)-TMPVAL*SVB(L)*DXYIV(L)*SAAY(L)*FY(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ANTI-DIFFUSE V ADVECTION FOR A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.LE.3) GOTO 801
C
C----------------------------------------------------------------------C
C
C **  CALCULATE V CONTROL VOLUME DIVERGENCE, TVAR1N
C
      DO K=1,KC
       DO L=2,LA
        TVAR1N(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=0,KC
       DO L=2,LA
        TVAR2C(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=1,KC
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR1N(L,K)=UHDY2(L+1,K)-UHDY2(L     ,K)
     &             +VHDX2(L  ,K)-VHDX2(LSC(L),K)            
     &             +(W2(L,K)-W2(L,K-1))*RDZIG             
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR1N(L,K)=0.25*DELT*DXYIV(L)*TVAR1N(L,K)          
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATED Y MOMENTUM EQUTION TERMS U*DVDX, V*DVDY, AND W*DVDZ
C
      DO K=1,KC
       DO L=2,LA
        UUU(L,K)=UHDY2(L,K)*(TVAR1S(L     ,K)-TVAR1S(L-1,K))
        VVV(L,K)=VHDX2(L,K)*(TVAR1S(LNC(L),K)-TVAR1S(L  ,K))
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        WWW(L,K)=RDZIG*W2(L,K)*(TVAR1S(L,K+1)-TVAR1S(L,K))
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE Y MOMENTUM EQUTION ANTIDIFFUSIVE VELOCITIES
C
      VELSMAL=1.E-12
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=ABS(UHDY2(L,K))
     &                  *(TVAR1S(L,K)-TVAR1S(L-1,K))
     &                  /(TVAR1S(L,K)+TVAR1S(L-1,K)+VELSMAL)
        TVAR2S(L,K)=ABS(VHDX2(L,K))
     &                  *(TVAR1S(LNC(L),K)-TVAR1S(L,K))
     &                  /(TVAR1S(LNC(L),K)+TVAR1S(L,K)+VELSMAL)
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=RDZIG*ABS(W2(L,K))
     &                       *(TVAR1S(L,K+1)-TVAR1S(L,K))
     &                       /(TVAR1S(L,K+1)+TVAR1S(L,K)+VELSMAL)
       ENDDO
      ENDDO
C
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=TVAR2W(L,K)-0.5*UHDY2(L,K)
     &                             *(TVAR1N(L,K)+TVAR1N(L-1,K))
        TVAR2S(L,K)=TVAR2S(L,K)-0.5*VHDX2(L,K)
     &                             *(TVAR1N(LNC(L),K)+TVAR1N(L,K))
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=TVAR2C(L,K)-0.5*W2(L,K)
     &                             *(TVAR1N(L,K+1)+TVAR1N(L,K))         
       ENDDO
      ENDDO
C
C
      DO K=1,KC
       DO L=2,LA
        TVAR2W(L,K)=TVAR2W(L,K)-DELT*UHDY2(L,K)
     &   *( UUU(L,K)
     &   +0.25*(VVV(L,K)+VVV(L-1,K)+VVV(LSC(L),K)+VVV(LSWC(L),K))
     &   +0.25*(WWW(L,K)+WWW(L-1,K)+WWW(L,K+1)+WWW(L-1,K+1)) )
     &   /(DXYV(L)*TVAR1S(L,K)+DXYV(L-1)*TVAR1S(L-1,K)+VELSMAL)
        TVAR2S(L,K)=TVAR2S(L,K)-DELT*VHDX2(L,K)
     &   *( VVV(L,K)
     &   +0.25*(UUU(L,K)+UUU(L+1,K)+UUU(LNC(L),K)+UUU(LNEC(L),K))
     &   +0.25*(WWW(L,K)+WWW(LNC(L),K)+WWW(L,K+1)+WWW(LNC(L),K+1)) )
     &   /(DXYV(L)*TVAR1S(L,K)+DXYV(LNC(L))*TVAR1S(LNC(L),K)+VELSMAL)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR2S(L,K)=STCUV(L)*SPB(L)*TVAR2S(L,K)
C
        TVAR2W(L,K)=SUB(L)*SUB(LSC(L))*TVAR2W(L,K)
C
       ENDDO
      ENDDO
C
      DO K=1,KS
       RDZIG=DZIC(K)
       DO L=2,LA
        TVAR2C(L,K)=TVAR2C(L,K)-DELT*W2(L,K)
     &   *( WWW(L,K)
     &   +0.25*(UUU(L,K)+UUU(L+1,K)+UUU(L,K+1)+UUU(L+1,K+1))
     &   +0.25*(VVV(LSC(L),K)+VVV(L,K)+VVV(LSC(L),K+1)+VVV(L,K+1)) )
     &   /(DXYV(L)*TVAR1S(L,K)+DXYV(L)*TVAR1S(L,K+1)+VELSMAL)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCLUATE V EQUATION ANTIDIFFUSIVE FLUXES
C
      DO K=1,KC
       DO L=2,LA
        FUHV(L,K)=MAX(TVAR2W(L,K),0.)*TVAR1S(L-1   ,K)
     &           +MIN(TVAR2W(L,K),0.)*TVAR1S(L     ,K)
        FVHV(L,K)=MAX(TVAR2S(L,K),0.)*TVAR1S(L     ,K)
     &           +MIN(TVAR2S(L,K),0.)*TVAR1S(LNC(L),K)
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
        FWV(L,K)=MAX(TVAR2C(L,K),0.)*TVAR1S(L,K)
     &          +MIN(TVAR2C(L,K),0.)*TVAR1S(L,K+1)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADVANCE V ANTIDIFFUSIVE ADVECTION A FULL TIME STEP
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       FY(L,K)=FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LSC(L),K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
       FY(L,K)=FY(L,K)+(FWV(L,K)-FWV(L,K-1))*DZIC(K)
       ENDDO
      ENDDO
C
      TMPVAL=DELT*SNLT
      DO K=1,KC
       DO L=2,LA
        TVAR1S(L,K)=TVAR1S(L,K)-TMPVAL*SVB(L)*DXYIV(L)*SAAY(L)*FY(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  SUBSTRACT CORIOLIS IF 3 TIME LEVEL
C
C----------------------------------------------------------------------C
C
  801 CONTINUE
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
       DO L=2,LA
         TVAR1W(L,K)=TVAR1W(L,K)-DELTD2*SUB(L)*DXYIU(L)*FCAX1(L,K)
         TVAR1S(L,K)=TVAR1S(L,K)+DELTD2*SVB(L)*DXYIV(L)*FCAY1(L,K)
        ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      IF(ISDCCA.EQ.0)THEN
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       ENDDO
      ENDDO
C
      ELSE
C
      CFMAX=CF
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
       CFMAX=MAX(CFMAX,CFEFF)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=0.25*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
       FCAY(L,K)=0.25*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
       DO L=2,LA
         TVAR1W(L,K)=TVAR1W(L,K)+DELT*SUB(L)*DXYIU(L)*FCAX(L,K)
         TVAR1S(L,K)=TVAR1S(L,K)-DELT*SVB(L)*DXYIV(L)*FCAY(L,K)
        ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KC
       DO L=2,LA
         TVAR1W(L,K)=TVAR1W(L,K)+DELTD2*SUB(L)*DXYIU(L)*FCAX(L,K)
         TVAR1S(L,K)=TVAR1S(L,K)-DELTD2*SVB(L)*DXYIV(L)*FCAY(L,K)
        ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
       DO L=2,LA
       FCAX1(L,K)=FCAX(L,K)
       FCAY1(L,K)=FCAY(L,K)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  ADD REMAINING WAVE FORCE
C
C----------------------------------------------------------------------C
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
      WVFACT=DELTD2*WVFACT      
C
      DO K=1,KC
       DO L=2,LA
        TVAR1W(L,K)=TVAR1W(L,K)+WVFACT*SAAX(L)*DXYIU(L)*FXWAVE(L,K)
        TVAR1S(L,K)=TVAR1S(L,K)+WVFACT*SAAY(L)*DXYIV(L)*FYWAVE(L,K)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
C     IF(ISHDMF.GE.1)THEN
C
C     DO K=1,KC
C      DO L=2,LA
C      LS=LSC(L)
C      LN=LNC(L)
C      FX(L,K)=FX(L,K)
C    &        -SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
C      FY(L,K)=FY(L,K)
C    &        -SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
C      ENDDO
C     ENDDO
C
C     ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
        FXE(L)=FXE(L)+DYU(L)*TVAR1W(L,K)*DZC(K)
        FYE(L)=FYE(L)+DXV(L)*TVAR1S(L,K)*DZC(K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADD REMAINING SURFACE AND BOTTOM STRESS
C
C----------------------------------------------------------------------C
C 
      DO L=2,LA
        FXE(L)=FXE(L)+DELTD2*SUB(L)*DYU(L)*(TSX(L)-TBX1(L))
        FYE(L)=FYE(L)+DELTD2*SVB(L)*DXV(L)*(TSY(L)-TBY1(L))
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( 2.*(TVAR1W(L,K+1)-TVAR1W(L,K))*DELTI 
     &           -DXYIU(L)*FBBX(L,K) )
        DV(L,K)=RCDZF*( 2.*(TVAR1S(L,K+1)-TVAR1S(L,K))*DELTI
     &           -DXYIV(L)*FBBY(L,K) )
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
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=SBX(L)*GP*HU(L)*
     &            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(L-1,K))*DZC(K) )
     &           -(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     &            (BELV(L)-BELV(L-1)+Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=SBY(L)*GP*HV(L)*
     &            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(LS,K))*DZC(K) )
     &           -(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     &            (BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS))) )
      ENDDO
      ENDDO
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
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=DU(L,K)+2.*RCDZF*DXYIU(L)*FBBX(L,K)
        DV(L,K)=DV(L,K)+2.*RCDZF*DXYIV(L)*FBBY(L,K)
       ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( 2.*(TVAR1W(L,K+1)-TVAR1W(L,K))*DELTI 
     &           +DXYIU(L)*FBBX(L,K) )
        DV(L,K)=RCDZF*( 2.*(TVAR1S(L,K+1)-TVAR1S(L,K))*DELTI
     &           +DXYIV(L)*FBBY(L,K) )
       ENDDO
      ENDDO
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
