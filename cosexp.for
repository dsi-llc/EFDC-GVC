C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE COSEXP (ISTL)
C
C **  SUBROUTINE COSEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
C **  USING THE COSMIC ADVECTION SCHEME
C **  THIS SUBROUTINE IS CURRENT PRODUCTION VERSION
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 18 OCTOBER 1998
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION FUHJ(LCM,KCM),FVHJ(LCM,KCM),CON1(LCM,KCM)
      DIMENSION DELCX(LCM,KCM),DELCY(LCM,KCM),DELCZ(LCM,KCM),
     &          CONCX(LCM,KCM),CONCY(LCM,KCM),CONCZ(LCM,KCM),
     &          CONCXY(LCM,KCM),CONCXZ(LCM,KCM),CONCYZ(LCM,KCM),
     &          CONCYX(LCM,KCM),CONCZX(LCM,KCM),CONCZY(LCM,KCM),
     &          CONCXYZ(LCM,KCM),CONCYZX(LCM,KCM),CONCZXY(LCM,KCM)
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
C      IF(N.EQ.1)THEN
C        OPEN(1,FILE='MFLUX.DIA')
C        CLOSE(1,STATUS='DELETE')
C      ENDIF
C
C      IF(N.LE.4)THEN
C        OPEN(1,FILE='MFLUX.DIA',POSITION='APPEND')
C      ENDIF
C
C**********************************************************************C
C
C **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
C
C----------------------------------------------------------------------C
C
      FCAXE(1)=0.
      FCAYE(1)=0.
      FCAX1E(1)=0.
      FCAY1E(1)=0.
      FXE(1)=0.
      FYE(1)=0.
C
      FCAXE(LC)=0.
      FCAYE(LC)=0.
      FCAX1E(LC)=0.
      FCAY1E(LC)=0.
      FXE(LC)=0.
      FYE(LC)=0.
C
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=LF+LDM-1
        DO L=LF,LL
          FCAXE(L)=0.
          FCAYE(L)=0.
          FCAX1E(L)=0.
          FCAY1E(L)=0.
          FXE(L)=0.
          FYE(L)=0.
        ENDDO
      ENDDO
C
      DO K=1,KC
        DO L=1,LC
          CON1(L,K)=0.0
        ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  U MOMENTUM FLUX
C
C----------------------------------------------------------------------C
C
C **  COSMIC INITIALIZATION
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=0.5*(U(L,K)+U(L+1,K))
       COSMICY(L,K)=0.5*(V(L,K)+V(L-1,K))
      ENDDO 
      ENDDO
C
      DO L=2,LA
       COSMICZ(L,0 )=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZ(L,K)=0.5*(W(L,K)+W(L-1,K))
        ENDDO 
        ENDDO
      ENDIF
C
      ELSE
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=0.25*(U(L,K)+U(L+1,K)+U1(L,K)+U1(L+1,K))
       COSMICY(L,K)=0.25*(V(L,K)+V(L-1,K)+V1(L,K)+V1(L-1,K))
      ENDDO 
      ENDDO
C
      DO L=2,LA
       COSMICZ(L,0 )=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZ(L,K)=0.25*(W(L,K)+W(L-1,K)+W1(L,K)+W1(L-1,K))
        ENDDO 
        ENDDO
      ENDIF
C
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
       RCOSMICX(L,K)=1.
       TMP=COSMICX(L,K)*COSMICX(L-1,K)
       IF(TMP.LT.0.) RCOSMICX(L,K)=0.
       RCOSMICY(L,K)=1.
       TMP=COSMICY(L,K)*COSMICY(LNC(L),K)
       IF(TMP.LT.0.) RCOSMICY(L,K)=0.
       RCOSMICZ(L,K)=1.
       TMP=COSMICZ(L,K)*COSMICZ(L,K-1)
       IF(TMP.LT.0.) RCOSMICZ(L,K)=0.
      ENDDO 
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=DELT*DXIU(L)*COSMICX(L,K)
       COSMICY(L,K)=2.*DELT*COSMICY(L,K)/(DYU(L)+DYU(LSC(L)))
      ENDDO 
      ENDDO
C
      IF(KC.GE.2.AND.ISTL.EQ.3)THEN
      DO K=1,KS
      DO L=2,LA
       COSMICZ(L,K)=DELT*DZIG(K)*COSMICZ(L,K)/HU(L)
      ENDDO 
      ENDDO
      ENDIF
C
      IF(KC.GE.2.AND.ISTL.EQ.2)THEN
      DO K=1,KS
      DO L=2,LA
       COSMICZ(L,K)=2.*DELT*DZIG(K)*COSMICZ(L,K)/(HU(L)+H1U(L))
      ENDDO 
      ENDDO
      ENDIF
C
      DO L=2,LA
       COSMICZ(L,0)=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CX,CY,CZ
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CON1(L,K)-CON1(L-1   ,K)
       DELCY(L,K)=CON1(L,K)-CON1(LSC(L),K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CON1(L,K+1)-CON1(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCX(L,K)=CON1(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCY(L,K)=CON1(L,K)-RCOSMICY(L,K)*
     &           ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &            +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCZ(L,KC)=CON1(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCZ(L,K)=CON1(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CXY,CXZ
C
      DO K=1,KC    
      DO L=2,LA
       DELCY(L,K)=CONCX(L,K)-CONCX(LSC(L),K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CONCX(L,K+1)-CONCX(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCXY(L,K)=CONCX(L,K)-RCOSMICY(L,K)*
     &            ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &             +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCXZ(L,KC)=CONCX(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCXZ(L,K)=CONCX(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CYZ,CYX
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CONCY(L,K)-CONCY(L-1,K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CONCY(L,K+1)-CONCY(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCYX(L,K)=CONCY(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCZ(L,1)=CONCY(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCYZ(L,1)=CONCY(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCYZ(L,KC)=CONCY(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCYZ(L,K)=CONCY(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CZX,CZY
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CONCZ(L,K)-CONCZ(L-1   ,K)
       DELCY(L,K)=CONCZ(L,K)-CONCZ(LSC(L),K)
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCZX(L,K)=CONCZ(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCZY(L,K)=CONCZ(L,K)-RCOSMICY(L,K)*
     &           ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &            +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
        CONCXYZ(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYZ(L,K)
     &                             +CONCZ(L,K)+CONCZY(L,K) )/6. 
        CONCYZX(L,K)=( 2.*CON1(L,K)+CONCZ(L,K)+CONCZX(L,K)
     &                             +CONCX(L,K)+CONCXZ(L,K) )/6. 
        CONCZXY(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYX(L,K)
     &                             +CONCX(L,K)+CONCXY(L,K) )/6. 
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONCXYZ(L-1,K)
     &         +MIN(UHDY2(L,K),0.)*CONCXYZ(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONCYZX(LSC(L),K)
     &         +MIN(VHDX2(L,K),0.)*CONCYZX(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*CONCZXY(L,K)
     &        +MIN(W2(L,K),0.)*CONCZXY(L,K+1)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C

C**********************************************************************C
C
C **  V MOMENTUM FLUX
C
C----------------------------------------------------------------------C
C
C **  COSMIC INITIALIZATION
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=0.5*(U(L,K)+U(LSC(L),K))
       COSMICY(L,K)=0.5*(V(L,K)+V(LNC(L),K))
      ENDDO 
      ENDDO
C
      DO L=2,LA
       COSMICZ(L,0 )=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZ(L,K)=0.5*(W(L,K)+W(LSC(L),K))
        ENDDO 
        ENDDO
      ENDIF
C
      ELSE
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=0.25*(U(L,K)+U(LSC(L),K)+U1(L,K)+U1(LSC(L),K))
       COSMICY(L,K)=0.25*(V(L,K)+V(LNC(L),K)+V1(L,K)+V1(LNC(L),K))
      ENDDO 
      ENDDO
C
      DO L=2,LA
       COSMICZ(L,0 )=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
      IF(KC.GE.2)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZ(L,K)=0.25*(W(L,K)+W(LSC(L),K)+W1(L,K)+W1(LSC(L),K))
        ENDDO 
        ENDDO
      ENDIF
C
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
       RCOSMICX(L,K)=1.
       TMP=COSMICX(L,K)*COSMICX(L+1,K)
       IF(TMP.LT.0.) RCOSMICX(L,K)=0.
       RCOSMICY(L,K)=1.
       TMP=COSMICY(L,K)*COSMICY(LSC(L),K)
       IF(TMP.LT.0.) RCOSMICY(L,K)=0.
       RCOSMICZ(L,K)=1.
       TMP=COSMICZ(L,K)*COSMICZ(L,K-1)
       IF(TMP.LT.0.) RCOSMICZ(L,K)=0.
      ENDDO 
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       COSMICX(L,K)=2.*DELT*COSMICX(L,K)/(DXV(L)+DXV(L-1))
       COSMICY(L,K)=DELT*DYIV(L)*COSMICY(L,K)
      ENDDO 
      ENDDO
C
      IF(KC.GE.2.AND.ISTL.EQ.3)THEN
      DO K=1,KS
      DO L=2,LA
       COSMICZ(L,K)=DELT*DZIG(K)*COSMICZ(L,K)/HV(L)
      ENDDO 
      ENDDO
      ENDIF
C
      IF(KC.GE.2.AND.ISTL.EQ.2)THEN
      DO K=1,KS
      DO L=2,LA
       COSMICZ(L,K)=2.*DELT*DZIG(K)*COSMICZ(L,K)/(HV(L)+H1V(L))
      ENDDO 
      ENDDO
      ENDIF
C
      DO L=2,LA
       COSMICZ(L,0)=0.
       COSMICZ(L,KC)=0.
      ENDDO 
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CX,CY,CZ
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CON1(L,K)-CON1(L-1   ,K)
       DELCY(L,K)=CON1(L,K)-CON1(LSC(L),K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CON1(L,K+1)-CON1(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCX(L,K)=CON1(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCY(L,K)=CON1(L,K)-RCOSMICY(L,K)*
     &           ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &            +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCZ(L,1)=CON1(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCZ(L,KC)=CON1(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCZ(L,K)=CON1(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CXY,CXZ
C
      DO K=1,KC    
      DO L=2,LA
       DELCY(L,K)=CONCX(L,K)-CONCX(LSC(L),K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CONCX(L,K+1)-CONCX(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCXY(L,K)=CONCX(L,K)-RCOSMICY(L,K)*
     &            ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &             +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCXZ(L,1)=CONCX(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCXZ(L,KC)=CONCX(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCXZ(L,K)=CONCX(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CYZ,CYX
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CONCY(L,K)-CONCY(L-1,K)
      ENDDO
      ENDDO
C
      IF(KC.GE.2)THEN
        DO K=1,KS    
        DO L=2,LA
         DELCZ(L,K)=CONCY(L,K+1)-CONCY(L,K)
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC    
      DO L=2,LA
      CONCYX(L,K)=CONCY(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
         CONCZ(L,1)=CONCY(L,1)
        ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
        DO L=2,LA
         CONCYZ(L,1)=CONCY(L,1)
     &     -RCOSMICZ(L,1)*(MIN(COSMICZ(L,1),0.)*DELCZ(L,1) )
        ENDDO
        DO L=2,LA
         CONCYZ(L,KC)=CONCY(L,KC)
     &     -RCOSMICZ(L,KC)*(MAX(COSMICZ(L,KS),0.)*DELCZ(L,KS) )
        ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
        DO K=2,KS    
        DO L=2,LA
         CONCYZ(L,K)=CONCY(L,K)-RCOSMICZ(L,K)*
     &              ( MAX(COSMICZ(L,K-1),0.)*DELCZ(L,K-1) 
     &               +MIN(COSMICZ(L  ,K),0.)*DELCZ(L,K  ) )
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INTERMEDIATE ADVECTION CALCULATIONS
C **  CZX,CZY
C
      DO K=1,KC    
      DO L=2,LA
       DELCX(L,K)=CONCZ(L,K)-CONCZ(L-1   ,K)
       DELCY(L,K)=CONCZ(L,K)-CONCZ(LSC(L),K)
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCZX(L,K)=CONCZ(L,K)-RCOSMICX(L,K)*
     &           ( MAX(COSMICX(L  ,K),0.)*DELCX(L  ,K)
     &            +MIN(COSMICX(L+1,K),0.)*DELCX(L+1,K) )
      ENDDO
      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      CONCZY(L,K)=CONCZ(L,K)-RCOSMICY(L,K)*
     &           ( MAX(COSMICY(L     ,K),0.)*DELCY(L     ,K)
     &            +MIN(COSMICY(LNC(L),K),0.)*DELCY(LNC(L),K) )
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
        CONCXYZ(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYZ(L,K)
     &                             +CONCZ(L,K)+CONCZY(L,K) )/6. 
        CONCYZX(L,K)=( 2.*CON1(L,K)+CONCZ(L,K)+CONCZX(L,K)
     &                             +CONCX(L,K)+CONCXZ(L,K) )/6. 
        CONCZXY(L,K)=( 2.*CON1(L,K)+CONCY(L,K)+CONCYX(L,K)
     &                             +CONCX(L,K)+CONCXY(L,K) )/6. 
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONCXYZ(L-1,K)
     &         +MIN(UHDY2(L,K),0.)*CONCXYZ(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONCYZX(LSC(L),K)
     &         +MIN(VHDX2(L,K),0.)*CONCYZX(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*CONCZXY(L,K)
     &        +MIN(W2(L,K),0.)*CONCZXY(L,K+1)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C


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
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=0.5*(UHDY1(L,K)+UHDY(L,K))
        VHDX2(L,K)=0.5*(VHDX1(L,K)+VHDX(L,K))
        W2(L,K)=W1(L,K)+W(L,K)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY2(L,K)+UHDY2(LS,K))
       UHB=0.5*(UHDY2(L,K)+UHDY2(L+1,K))
       VHC=0.5*(VHDX2(L,K)+VHDX2(L-1,K))
       VHB=0.5*(VHDX2(L,K)+VHDX2(LN,K))
C
       FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
     &         +MIN(UHB,0.)*U1(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
     &         +MIN(VHC,0.)*U1(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
     &         +MIN(UHC,0.)*V1(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
     &         +MIN(VHB,0.)*V1(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
C
C ADD RETURN FLOW MOMENTUM FLUX
C
      DO NWR=1,NQWR
       IF(NQWRMFU(NWR).GT.0)THEN
         IU=IQWRU(NWR)
         JU=JQWRU(NWR)
         KU=KQWRU(NWR)
         LU=LIJ(IU,JU)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LU)*DZC(KU)*DZC(KU)*BQWRMFU(NWR))
         IF(NQWRMFU(NWR).EQ.1)  FUHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.2)  FVHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.3)  FUHJ(LU+1   ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.4)  FVHJ(LNC(LU),KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.-1) FUHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-2) FVHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-3) FUHJ(LU+1   ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-4) FVHJ(LNC(LU),KU)=QUMF 
       ENDIF
       IF(NQWRMFD(NWR).GT.0)THEN
         ID=IQWRD(NWR)
         JD=JQWRD(NWR)
         KD=KQWRD(NWR)
         LD=LIJ(ID,JD)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))
         IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=QUMF
         IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=-QUMF
C         IF(N.LE.4)THEN
C           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,
C     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)
C         ENDIF
       ENDIF
      ENDDO
C
C ** HARDWIRE FOR PEACH BOTTOM
C
C      DO K=1,KC
C       FVHV(535,K)=700./H1P(535)
C      ENDDO
C
C ** END HARDWIRE FOR PEACH BOTTOM
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)
     &        +MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)
     &        +MIN(WV,0.)*V1(L,K+1)
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
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.5*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.5*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.5*(VHDX(L,K)+VHDX(LN,K))
       FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
     &         +MIN(UHB,0.)*U1(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
     &         +MIN(VHC,0.)*U1(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
     &         +MIN(UHC,0.)*V1(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
     &         +MIN(VHB,0.)*V1(LN,K)
       FUHJ(L,K)=0.
       FVHJ(L,K)=0.
       ENDDO
      ENDDO
C
C ADD RETURN FLOW MOMENTUM FLUX
C
      DO NWR=1,NQWR
       IF(NQWRMFU(NWR).GT.0)THEN
         IU=IQWRU(NWR)
         JU=JQWRU(NWR)
         KU=KQWRU(NWR)
         LU=LIJ(IU,JU)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LU)*DZC(KU)*DZC(KU)*BQWRMFU(NWR))
         IF(NQWRMFU(NWR).EQ.1)  FUHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.2)  FVHJ(LU     ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.3)  FUHJ(LU+1   ,KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.4)  FVHJ(LNC(LU),KU)=-QUMF 
         IF(NQWRMFU(NWR).EQ.-1) FUHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-2) FVHJ(LU     ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-3) FUHJ(LU+1   ,KU)=QUMF 
         IF(NQWRMFU(NWR).EQ.-4) FVHJ(LNC(LU),KU)=QUMF 
       ENDIF
       IF(NQWRMFD(NWR).GT.0)THEN
         ID=IQWRD(NWR)
         JD=JQWRD(NWR)
         KD=KQWRD(NWR)
         LD=LIJ(ID,JD)
         NS=NQWRSERQ(NWR)
         QMF=QWR(NWR)+QWRSERT(NS)
         QUMF=QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))
         IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=QUMF
         IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=QUMF
         IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=-QUMF
         IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=-QUMF
C         IF(N.LE.4)THEN
C           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,
C     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)
C         ENDIF
       ENDIF
      ENDDO
C
C ** HARDWIRE FOR PEACH BOTTOM
C
C      DO K=1,KC
C       FVHV(535,K)=700./H1P(535)
C      ENDDO
C
C ** END HARDWIRE FOR PEACH BOTTOM
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
       WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)
     &        +MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)
     &        +MIN(WV,0.)*V1(L,K+1)
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
      DO K=1,KC
       DO L=1,LA
       FUHU(L,K)=STCUV(L)*FUHU(L,K)
       FVHV(L,K)=STCUV(L)*FVHV(L,K)
       ENDDO
      ENDDO
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
C      IF(N.EQ.2)THEN
C       OPEN(1,FILE='CORC.DIA')
C       CLOSE(2,STATUS='DELETE')
C       OPEN(1,FILE='CORC.DIA')
C       K=1
C       DO L=2,LA
C       LN=LNC(L)          
C       WRITE(1,1111)IL(L),JL(L),LN,V(LN,K),V(L,K),DYU(L+1),DYU(L),
C     &        U(L+1,K),U(L,K),DXV(LN),DXV(L),HP(L),CAC(L,K)
C       ENDDO
C       CLOSE(1)
C      ENDIF
 1111 FORMAT(3I5,10E13.4)
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        FCAX1(L,K)=0.
        FCAY1(L,K)=0.     
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.3)THEN
C      IF(IRVEC.EQ.10.OR.IRVEC.EQ.11)THEN
      IF(IRVEC.EQ.14)THEN
C
      DO L=1,LC
       TVAR3E(L)=0.
       TVAR3N(L)=0.
       TVAR3W(L)=0.
       TVAR3S(L)=0.
       TVAR1S(L,1)=0.
       TVAR2S(L,1)=0.
      ENDDO
C
      DO K=1,KC
      DZTMP=DZC(K) 
      DO L=2,LA
       TVAR3E(L)=TVAR3E(L)+DZTMP*U(L,K)
       TVAR3N(L)=TVAR3N(L)+DZTMP*V(L,K)
       TVAR3W(L)=TVAR3W(L)+DZTMP*U1(L,K)
       TVAR3S(L)=TVAR3S(L)+DZTMP*V1(L,K)
      ENDDO
      ENDDO
C
C       CAC(L,K)=( FCORC(L)*DXYP(L)
C     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
C
C       TVAR1S=CAC @ N
C       TVAR2S=CAC @ N-1
C
C       TVAR3E=U @ N
C       TVAR3N=V @ N
C       TVAR3W=U @ N-1
C       TVAR3S=V @ N-1
C
      DO L=2,LA
       LN=LNC(L)          
       TVAR1S(L,1)=( FCORC(L)*DXYP(L)
     &      +0.5*SNLT*(TVAR3N(LN)+TVAR3N(L))*(DYU(L+1)-DYU(L))
     &      -0.5*SNLT*(TVAR3E(L+1)+TVAR3E(L))*(DXV(LN)-DXV(L)) )*HP(L)
       TVAR2S(L,1)=( FCORC(L)*DXYP(L)
     &      +0.5*SNLT*(TVAR3S(LN)+TVAR3S(L))*(DYU(L+1)-DYU(L))
     &      -0.5*SNLT*(TVAR3W(L+1)+TVAR3W(L))*(DXV(LN)-DXV(L)) )*H1P(L)
      ENDDO
C
C       FCAX(L,K)=ROLD*FCAX(L,K)
C     &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
C     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
C       FCAY(L,K)=ROLD*FCAY(L,K)
C     &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
C     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
C
C
C       TVAR1S=CAC @ N
C       TVAR2S=CAC @ N-1
C
C       TVAR3E=U @ N
C       TVAR3N=V @ N
C       TVAR3W=U @ N-1
C       TVAR3S=V @ N-1
C
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=0.0625*SCAX(L)*
     &            ( TVAR2S(L,1)*(TVAR3S(LN)+TVAR3S(L))
     &             +TVAR2S(L-1,1)*(TVAR3S(LNW)+TVAR3S(L-1)) )
       FCAY1E(L)=0.0625*SCAY(L)*
     &            ( TVAR2S(L,1)*(TVAR3W(L+1)+TVAR3W(L))
     &             +TVAR2S(LS,1)*(TVAR3W(LSE)+TVAR3W(LS)) )
      ENDDO
C
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=FCAX1E(L)-0.125*SCAX(L)*
     &            ( TVAR1S(L,1)*(TVAR3N(LN)+TVAR3N(L))
     &             +TVAR1S(L-1,1)*(TVAR3N(LNW)+TVAR3N(L-1)) )
       FCAY1E(L)=FCAY1E(L)-0.125*SCAY(L)*
     &            ( TVAR1S(L,1)*(TVAR3E(L+1)+TVAR3E(L))
     &             +TVAR1S(LS,1)*(TVAR3E(LSE)+TVAR3E(LS)) )
      ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=ROLD*FCAX(L,K)
     &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
       FCAY(L,K)=ROLD*FCAY(L,K)
     &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K)
     &                 +FUHJ(L,K) )
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K)
     &                 +FVHJ(L,K) )
       ENDDO
      ENDDO
C
C      IF(N.EQ.2)THEN
C       OPEN(1,FILE='CORC.DIA',POSITION='APPEND')
C        L=2
C        K=1
C        LN=LNC(L)
C        LNW=LNWC(L)
C        WRITE(1,1111)IL(L),JL(L),LNW,SCAX(L),CAC(L,K),V(LN,K),V(L,K),
C     &                              CAC(L-1,K),V(LNW,K),V(L-1,K)
C       DO L=2,LA
C       WRITE(1,1111)L,IL(L),JL(L),(FCAX(L,K),K=1,KC)
C       ENDDO
C       DO L=2,LA
C       WRITE(1,1111)L,IL(L),JL(L),(FCAY(L,K),K=1,KC)
C       ENDDO
C      ENDIF
C
C**********************************************************************C
C
C **  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISVEG.GE.1)THEN
C
      DO L=2,LA
        FXVEGE(L)=0.
        FYVEGE(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LW=L-1
        LE=L+1
        LS=LSC(L)
        LN=LNC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        VTMPATU=0.25*(V1(L,K)+V1(LW,K)+V1(LN,K)+V1(LNW,K))
        UTMPATV=0.25*(U1(L,K)+U1(LE,K)+U1(LS,K)+U1(LSE,K))
        UMAGTMP=SQRT( U1(L,K)*U1(L,K)+VTMPATU*VTMPATU )
        VMAGTMP=SQRT( UTMPATV*UTMPATV+V1(L,K)*V1(L,K) )
        FXVEG(L,K)=UMAGTMP*SUB(L)*DXYU(L)*FXVEG(L,K)
        FYVEG(L,K)=VMAGTMP*SVB(L)*DXYV(L)*FYVEG(L,K)
        FXVEGE(L)=FXVEGE(L)+FXVEG(L,K)*DZC(K)
        FYVEGE(L)=FYVEGE(L)+FYVEG(L,K)*DZC(K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        FXVEG(L,K)=FXVEG(L,K)*U1(L,K)
        FYVEG(L,K)=FYVEG(L,K)*V1(L,K)
        FX(L,K)=FX(L,K)+FXVEG(L,K)-FXVEGE(L)*U1(L,K)
        FY(L,K)=FY(L,K)+FYVEG(L,K)-FYVEGE(L)*V1(L,K)
       ENDDO
      ENDDO
C
      DO L=2,LA
        FXVEGE(L)=DXYIU(L)*FXVEGE(L)/H1U(L)
        FYVEGE(L)=DXYIV(L)*FYVEGE(L)/H1V(L)
      ENDDO
C
      ENDIF
C
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISHDMF.GE.1)THEN
C
      DO K=1,KC
       DO L=2,LA
       LS=LSC(L)
       LN=LNC(L)
       FX(L,K)=FX(L,K)
     &        -SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
       FY(L,K)=FY(L,K)
     &        -SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
       ENDDO
      ENDDO
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
      IF(ISWAVE.GE.1)THEN
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
      DO K=1,KC
       DO L=2,LA
         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.1)THEN
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
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*
     &            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(L-1,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     &            (HMP(L)-HMP(L-1)-Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*
     &            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(LS,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     &            (HMP(L)-HMP(LS)-Z(K)*(HP(L)-HP(LS))) )
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
C 1111 FORMAT(2I5,2X,8E12.4)       
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
        DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
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
      IF(N.LE.4)THEN
        CLOSE(1)
      ENDIF
C
 1112 FORMAT('N,NW,NS,I,J,K,NF,H,Q,QU,FUU,FVV=',/,2X,7I5,5E12.4)
C
C**********************************************************************C
C
      RETURN
      END
