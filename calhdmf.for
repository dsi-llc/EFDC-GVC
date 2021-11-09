C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALHDMF
C
C **  SUBROUTINE CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND 
C **  DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING
C **  SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO
C
C     SMAGORINSKY, J., 1993: SOME HISTORICAL REMARKS ON THE USE OF 
C     NONLINEAR VISCOSITIES. IN LARGE EDDY SIMULATION OF COMPLEX 
C     ENGINEERING AND GEOPHYSICAL FLOWS. B. GALPERIN AND S. A. ORSZAG, 
C     EDS. CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE, UK.
C
C     COMMENT CSO INDICATE OLD VERSION
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 01/12/2002        John Hamrick       01/12/2002       John Hamrick
C  replaced zbr() with zbrwall for side wall slip friction
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION SXY2CC(LCM,KCM),SXY2EE(LCM,KCM),SXY2NN(LCM,KCM)
      DIMENSION AHCC(LCM,KCM),AHEE(LCM,KCM),AHNN(LCM,KCM)
      DIMENSION ICORDYU(LCM),ICORDXV(LCM)
C
C**********************************************************************C
C
      AHMAX=AHO
C
      DO K=1,KC
      DO L=1,LC
       DXU1(L,K)=0.
       DYV1(L,K)=0.
       DYU1(L,K)=0.
       DXV1(L,K)=0.
       FMDUX(L,K)=0.
       FMDVY(L,K)=0.
       FMDUY(L,K)=0.
       FMDVX(L,K)=0.
       SXY2CC(L,K)=0.
       SXY2EE(L,K)=0.
       SXY2NN(L,K)=0.
       AHCC(L,K)=0.
       AHEE(L,K)=0.
       AHNN(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCUATE TYPE FLAGS
C
C----------------------------------------------------------------------C
C
C ** ICORDYU
C
      DO L=2,LA
       LS=LSC(L)
       IF(SUBO(L).LT.0.5.AND.SUBO(LS).LT.0.5) ICORDYU(L)=0
       IF(SUBO(L).GT.0.5.AND.SUBO(LS).GT.0.5) ICORDYU(L)=1
       IF(SUBO(L).LT.0.5.AND.SUBO(LS).GT.0.5) ICORDYU(L)=2
       IF(SUBO(L).GT.0.5.AND.SUBO(LS).LT.0.5) ICORDYU(L)=3
      ENDDO
C
C ** ICORDXV
C
      DO L=2,LA
       LW=L-1
       IF(SVBO(L).LT.0.5.AND.SVBO(LW).LT.0.5) ICORDXV(L)=0
       IF(SVBO(L).GT.0.5.AND.SVBO(LW).GT.0.5)THEN
         ICORDXV(L)=1
         IF(SUBO(L).LT.0.5) ICORDXV(L)=3
       ENDIF
       IF(SVBO(L).LT.0.5.AND.SVBO(LW).GT.0.5) ICORDXV(L)=2
       IF(SVBO(L).GT.0.5.AND.SVBO(LW).LT.0.5) ICORDXV(L)=3
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VELOCITY SHEARS - OLD
C
CSO      DO K=1,KC
CSO      DO L=2,LA
CSO      LS=LSC(L)      
CSO      LN=LNC(L)      
CSO      DXU1(L,K)=(U1(L+1,K)-U1(L,K))/DXP(L)
CSO      DYV1(L,K)=(V1(LN,K)-V1(L,K))/DYP(L)
CSOC      DYU1(L,K)=2.*SUB(LS)*(U1(L,K)-U1(LS,K))/(DYU(L)+DYU(LS))
CSOC      DXV1(L,K)=2.*SVB(L-1)*(V1(L,K)-V1(L-1,K))/(DXV(L)+DXV(L-1))
CSO      DYU1(L,K)=2.*SVB(L)*SVB(L-1)*(U1(L,K)-U1(LS,K))/(DYU(L)+DYU(LS))
CSO      DXV1(L,K)=2.*SUB(L)*SUB(LS)*(V1(L,K)-V1(L-1,K))/(DXV(L)+DXV(L-1))
CSOC      IF(SUB(LS).LT.0.5) DYU1(L,K)=2.*U1(L,K)/DYU(L)
CSOC      IF(SVB(L-1).LT.0.5) DXV1(L,K)=2.*V1(L,K)/DXV(L)
CSO      IF(SVB(L).LT.0.5.OR.SVB(L-1).LT.0.5) DYU1(L,K)=2.*U1(L,K)/DYU(L)
CSO      IF(SUB(L).LT.0.5.OR.SUB(LS).LT.0.5)  DXV1(L,K)=2.*V1(L,K)/DXV(L)
CSO      ENDDO
CSO      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VELOCITY SHEARS - NEW
C
C     SXX-SYY DEFINED AT CELL CENTERS AND STORED IN DXU1(L<K)
C
      SLIPCO=0.5/SQRT(AHD)
C     
      DO K=1,KC
      DO L=2,LA
        LN=LNC(L)      
        DXU1(L,K)=(U(L+1,K)-U(L,K))/DXP(L)
        DYV1(L,K)=(V(LN,K)-V(L,K))/DYP(L)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        LN=LNC(L)      
        DXU1(L,K)=DXU1(L,K)-DYV1(L,K)
     &  +0.5*( (V(LN ,K)+V(L,K))*DXDJ(L)
     &        -(U(L+1,K)+U(L,K))*DYDI(L) )/DXYP(L)
      ENDDO
      ENDDO
C
C     DYU1(L,K) DEFINED AT CELL CORNER
C
      DO K=1,KC
      DO L=2,LA
        LS=LSC(L)
        DYU1(L,K)=0.
        IF(ICORDYU(L).EQ.1) DYU1(L,K)=2.*(U(L,K)-U(LS,K))
     &                                  /(DYU(L)+DYU(LS))
        IF(ICORDYU(L).EQ.2.AND.ISHDMF.EQ.2)THEN
C          DY2DZBR=1.+0.5*DYU(LS)/ZBR(LS)
          DY2DZBR=1.+0.5*DYU(LS)/ZBRWALL
          CSDRAG=0.16/((LOG(DY2DZBR))**2)
          SLIPFAC=SLIPCO*SQRT(CSDRAG)
          DYU1(L,K)=-2.*SLIPFAC*U(LS,K)/DYU(LS)
        ENDIF
        IF(ICORDYU(L).EQ.3.AND.ISHDMF.EQ.2)THEN
C          DY2DZBR=1.+0.5*DYU(L)/ZBR(L)
          DY2DZBR=1.+0.5*DYU(L)/ZBRWALL
          CSDRAG=0.16/((LOG(DY2DZBR))**2)
          SLIPFAC=SLIPCO*SQRT(CSDRAG)
          DYU1(L,K)=2.*SLIPFAC*U(L,K)/DYU(L)
        ENDIF
      ENDDO
      ENDDO
C
C     DXV1(L,K) DEFINED AT CELL CORNER
C
      DO K=1,KC
      DO L=2,LA
        LW=L-1
        DXV1(L,K)=0.
        IF(ICORDXV(L).EQ.1) DXV1(L,K)=2.*(V(L,K)-V(LW,K))
     &                                  /(DXV(L)+DXV(LW))
        IF(ICORDXV(L).EQ.2.AND.ISHDMF.EQ.2)THEN
C          DX2DZBR=1.+0.5*DXV(LW)/ZBR(LW)
          DX2DZBR=1.+0.5*DXV(LW)/ZBRWALL
          CSDRAG=0.16/((LOG(DX2DZBR))**2)
          SLIPFAC=SLIPCO*SQRT(CSDRAG)
          DXV1(L,K)=-2.*SLIPFAC*V(LW,K)/DXV(LW)
        ENDIF
        IF(ICORDXV(L).EQ.3.AND.ISHDMF.EQ.2)THEN
C          DX2DZBR=1.+0.5*DXV(L)/ZBR(L)
          DX2DZBR=1.+0.5*DXV(L)/ZBRWALL
          CSDRAG=0.16/((LOG(DX2DZBR))**2)
          SLIPFAC=SLIPCO*SQRT(CSDRAG)
          DXV1(L,K)=2.*SLIPFAC*V(L,K)/DXV(L)
        ENDIF
      ENDDO
      ENDDO
C
C     SXY STORED IN DYU1
C
      DO K=1,KC
      DO L=2,LA
        LS=LSC(L)      
        LSW=LSWC(L)
        TMPVAL=1.+SUBO(L)+SVBO(L)+SUBO(L)*SVBO(L)
        DXYCOR=(DXYP(L)+SUBO(L)*DXYP(L-1)+SVBO(L)*DXYP(LS)
     &                 +SUBO(L)*SVBO(L)*DXYP(LSW))/TMPVAL     
        DYU1(L,K)=DYU1(L,K)+DXV1(L,K)
     &  -0.5*( SUB(L)*SUB(LS )*(V(L,K)+V(L-1,K))*(DYV(L)-DYV(L-1))  
     &        +SVB(L)*SVB(L-1)*(U(L,K)+U(LS ,K))*(DXU(L)-DXU(LS ))
     &       )/DXYCOR  
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        SXY2CC(L,K)=DYU1(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        IF(SUBO(L+1).GT.0.5) SXY2EE(L,K)=DYU1(L+1,K)
        IF(SUBO(L+1).LT.0.5)THEN
          IF(ISHDMF.EQ.1) SXY2EE(L,K)=0.0
          IF(ISHDMF.EQ.2)THEN
C            DX2DZBR=1.+0.5*DXV(L)/ZBR(L)
            DX2DZBR=1.+0.5*DXV(L)/ZBRWALL
            CSDRAG=0.16/((LOG(DX2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            SXY2EE(L,K)=-2.*SLIPFAC*V(L,K)/DXV(L)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        LN=LNC(L)
        IF(SVBO(LN).GT.0.5) SXY2NN(L,K)=DYU1(LN,K)
        IF(SVBO(LN).LT.0.5)THEN
          IF(ISHDMF.EQ.1) SXY2NN(L,K)=0.0
          IF(ISHDMF.EQ.2)THEN
C            DY2DZBR=1.+0.5*DYU(L)/ZBR(L)
            DY2DZBR=1.+0.5*DYU(L)/ZBRWALL
            CSDRAG=0.16/((LOG(DY2DZBR))**2)
            SLIPFAC=SLIPCO*SQRT(CSDRAG)
            SXY2NN(L,K)=-2.*SLIPFAC*U(L,K)/DYU(L)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
C
C
      IF(N.EQ.2)THEN
        OPEN(1,FILE='AHSXY.DIA')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='AHSXY.DIA')
        DO L=2,LA
C        WRITE(1,1112)IL(L),JL(L),(AH(L,K),K=1,KC),(AHC(L,K),K=1,KC)
        WRITE(1,1112)IL(L),JL(L),SXY2CC(L,KC),SXY2EE(L,KC),
     &                                    SXY2NN(L,KC)
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY - OLD
C
CSO      IF(AHD.GT.0.0)THEN
C
CSO      DO K=1,KC
CSO      DO L=2,LA
CSO      LN=LNC(L)
CSO      LS=LSC(L)     
CSO      LNW=LNWC(L)
CSO      LNE=LNEC(L)
CSO      LSW=LSWC(L)
CSO      LSE=LSEC(L)
CSO      AH(L,K)=AHO+AHD*DXP(L)*DYP(L)*SQRT( 2.*DXU1(L,K)*DXU1(L,K)
CSO     &                                   +2.*DYV1(L,K)*DYV1(L,K)
CSO     &  +0.0625*(DYU1(L,K)+DYU1(LN,K)+DYU1(L+1,K)+DYU1(LNE,K)
CSO     &          +DXV1(L,K)+DXV1(LN,K)+DXV1(L+1,K)+DXV1(LNE,K))**2 )
CSO      AHC(L,K)=AHO+AHD*0.0625*( (DXP(L)+DXP(L-1)+DXP(LS)+DXP(LSW))**2 )
CSO     &   *SQRT( 0.125*(DXU1(L,K)+DXU1(L-1,K)+DXU1(LS,K)+DXU1(LSW,K))**2
CSO     &         +0.125*(DYV1(L,K)+DYV1(L-1,K)+DYV1(LS,K)+DYV1(LSW,K))**2
CSO     &         +DYU1(L,K)*DYU1(L,K)+2.*DYU1(L,K)*DXV1(L,K)
CSO     &         +DXV1(L,K)*DXV1(L,K) )
CSO      ENDDO
CSO      ENDDO
C
CSO      ELSE
C
CSO      DO K=1,KC
CSO      DO L=2,LA
CSO       AH(L,K)=AHO
CSO       AHC(L,K)=AHO
CSO      ENDDO
CSO      ENDDO
C
CSO      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY - NEW
C
      IF(AHD.GT.0.0)THEN
C
      IF(N.EQ.2)THEN
        OPEN(1,FILE='AHNN.DIA')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='AHNN.DIA')
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
       LN=LNC(L)     
       LNE=LNEC(L)
       DSP=MIN(DXP(L),DYP(L))
       TMPVAL=AHD*DSP*DSP
       DSQR=DXU1(L,K)*DXU1(L,K)+0.25*(SXY2CC(L,K)*SXY2CC(L,K)
     &                               +SXY2EE(L,K)*SXY2EE(L,K)
     &                               +SXY2NN(L,K)*SXY2NN(L,K)
     &                               +SXY2CC(LNE,K)*SXY2CC(LNE,K))
       AH(L,K)=AHO+TMPVAL*SQRT(DSQR)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       LS=LSC(L)     
       LSW=LSWC(L)
       DXTMP=0.5*(DXV(L)+DXV(L-1))
       IF(SUBO(L).LT.0.5) DXTMP=DXV(L)
       DYTMP=0.5*(DYU(L)+DYU(LS ))
       IF(SVBO(L).LT.0.5) DYTMP=DYU(L)
       DSP=MIN(DXTMP,DYTMP)
       TMPVAL=AHD*DSP*DSP
       DSQR=SXY2CC(L,K)*SXY2CC(L,K)+0.25*(DXU1(L  ,K)*DXU1(L  ,K)
     &                           +SUBO(L)*DXU1(L-1,K)*DXU1(L-1,K)
     &                           +SVBO(L)*DXU1(LS ,K)*DXU1(LS ,K)
     &                   +SUBO(L)*SVBO(L)*DXU1(LSW,K)*DXU1(LSW,K))
       AHC(L,K)=AHO+TMPVAL*SQRT(DSQR)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        AHCC(L,K)=AHC(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        LS=LSC(L)
        LSE=LSEC(L)
        IF(SUB(L+1).GT.0.5) AHEE(L,K)=AHC(L+1,K)
        IF(SUB(L+1).LT.0.5)THEN
          IF(ISHDMF.EQ.1) AHEE(L,K)=0.0
          IF(ISHDMF.EQ.2)THEN
            DXTMP=0.5*(DXV(L+1)+DXV(L))
            IF(SUBO(L+1).LT.0.5) DXTMP=DXV(L)
C            DYTMP=0.5*(DYU(L+1)+DYU(LSW))
C            IF(SVBO(L+1).LT.0.5) DYTMP=DYU(L)
C            DSP=MIN(DXTMP,DYTMP)
            DSP=DXTMP
            TMPVAL=AHD*DSP*DSP
            DSQR=SXY2EE(L,K)*SXY2EE(L,K)
     &                      +0.25*(DXU1(L  ,K)*DXU1(L  ,K)
     &                  +SUBO(L+1)*DXU1(L+1,K)*DXU1(L+1,K)
     &                  +SVBO(L  )*DXU1(LS ,K)*DXU1(LS ,K)
     &        +SUBO(LSE)*SVBO(L+1)*DXU1(LSE,K)*DXU1(LSE,K))
            AHEE(L,K)=AHO+TMPVAL*SQRT(DSQR)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        LN=LNC(L)
        LNW=LNWC(L)
        IF(SVBO(LN).GT.0.5) AHNN(L,K)=AHC(LN,K)
        IF(SVBO(LN).LT.0.5)THEN
          IF(ISHDMF.EQ.1) AHNN(L,K)=0.0
          IF(ISHDMF.EQ.2)THEN
C            DXTMP=0.5*(DXV(LN)+DXV(LNW))
C            IF(SUBO(LN).LT.0.5) DXTMP=DXV(LN)
            DYTMP=0.5*(DYU(LN)+DYU(L))
            IF(SVBO(LN).LT.0.5) DYTMP=DYU(L)
C            DSP=MIN(DXTMP,DYTMP)
            DSP=DYTMP
            TMPVAL=AHD*DSP*DSP
            DSQR=SXY2NN(L,K)*SXY2NN(L,K)
     &                      +0.25*(DXU1(L  ,K)*DXU1(L  ,K)
     &                   +SVBO(LN)*DXU1(LN ,K)*DXU1(LN ,K)
     &                   +SUBO(L )*DXU1(L-1,K)*DXU1(L-1,K)
     &          +SUBO(LN)*SVBO(LNW)*DXU1(LNW,K)*DXU1(LNW,K))
            IF(K.EQ.KC)THEN
            IF(N.EQ.2)THEN
              TMP1=+0.25*(DXU1(L  ,K)*DXU1(L  ,K)
     &                   +SVBO(LN)*DXU1(LN ,K)*DXU1(LN ,K)
     &                   +SUBO(L )*DXU1(L-1,K)*DXU1(L-1,K)
     &          +SUBO(LN)*SVBO(LNW)*DXU1(LNW,K)*DXU1(LNW,K))
              TMP1=SQRT(TMP1)
              WRITE(1,1112)IL(L),JL(L),SXY2NN(L,KC),TMP1,TMPVAL
            ENDIF
            ENDIF
            AHNN(L,K)=AHO+TMPVAL*SQRT(DSQR)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KC
      DO L=2,LA
       AH(L,K)=AHO
       AHC(L,K)=AHO
       AHCC(L,K)=AHO
       AHEE(L,K)=AHO
       AHNN(L,K)=AHO
      ENDDO
      ENDDO
C
      IF(N.EQ.2)THEN
        CLOSE(1)
      ENDIF
C
      ENDIF
C
      IF(N.EQ.2)THEN
        OPEN(1,FILE='AHDIFF.DIA')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='AHDIFF.DIA')
        DO L=2,LA
C        WRITE(1,1112)IL(L),JL(L),(AH(L,K),K=1,KC),(AHC(L,K),K=1,KC)
        WRITE(1,1112)IL(L),JL(L),AH(L,KC),AHCC(L,KC),AHEE(L,KC),
     &                                    AHNN(L,KC),SXY2NN(L,KC)
        ENDDO
        CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL DIFFUSION DUE TO WAVE BREAKING
C
      IF(ISWAVE.GE.2)THEN
      IF(WVLSH.GT.0.0.OR.WVLSX.GT.0.0)THEN
C
       IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
        ELSE
         WVFACT=1.0
       ENDIF      
C
C     DO K=1,KC
C     DO L=2,LA
C     DTMP=(WVFACT*WVDISP(L,K)/HMP(L))**0.3333
C     RLTMP=WVLSH*( HMP(L)**1.3333)+WVLSX*( DXYP(L)**0.6667 )
C     AH(L,K)=AH(L,K)+DTMP*RLTMP
C     ENDDO
C     ENDDO
C
C     DO K=2,KC
C     DO L=2,LA
C     LS=LSC(L)     
C     LSW=LSWC(L)
C     DTMP=0.25*WVFACT*(WVDISP(L,K)+WVDISP(L-1,K)
C    &          +WVDISP(LS,K)+WVDISP(LSW,K))
C     HCTMP=0.25*(HMP(L)+HMP(L-1)
C    &          +HMP(LS)+HMP(LSW))
C     DTMP=(DTMP/HCTMP)**0.3333
C     RLTMP=WVLSH*(HCTMP**1.3333)+WVLSX*( DXYP(L)**0.6667 )
C     AHC(L,K)=AHC(L,K)+DTMP*RLTMP
C     ENDDO
C     ENDDO
C
C
      AHWVX=WVFACT*WVLSX*WVPRD*WVPRD
C
      DO K=1,KC
      DO L=2,LA
      DTMPH=(WVFACT*WVDISP(L,K))**0.3333
      DTMPX=WVDISP(L,K)/HMP(L)
      AH(L,K)=AH(L,K)+WVLSH*DTMPH*HMP(L)+AHWVX*DTMPX
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LSW=LSWC(L)
      DTMP=0.25*(WVDISP(L,K)+WVDISP(L-1,K)
     &          +WVDISP(LS,K)+WVDISP(LSW,K))
      DTMPH=(WVFACT*DTMP)**0.3333
      DTMPX=DTMP/HMC(L)
      AHC(L,K)=AHC(L,K)+WVLSH*DTMPH*HMC(L)+AHWVX*DTMPX
      ENDDO
      ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE DIFFUSIVE MOMENTUM FLUXES - OLD
C
CSO      DO K=1,KC
CSO      DO L=2,LA
CSO      LN=LNC(L)
CSO      LS=LSC(L)      
CSO      FMDUX(L,K)=2.*DYP(L)*H1P(L)*AH(L,K)*DXU1(L,K)
CSO      FMDUY(L,K)=0.5*(DXU(L)+DXU(LS))*H1C(L)*AHC(L,K)
CSO     &              *(DYU1(L,K)+DXV1(L,K))
CSO      FMDVY(L,K)=2.*DXP(L)*H1P(L)*AH(L,K)*DYV1(L,K)
CSO      FMDVX(L,K)=0.5*(DYV(L)+DYV(L-1))*H1C(L)*AHC(L,K)
CSO     &              *(DYU1(L,K)+DXV1(L,K))
CSO      ENDDO
CSO      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE DIFFUSIVE MOMENTUM FLUXES - NEW
C
C      DO K=1,KC
C      DO L=2,LA
C       LS=LSC(L)      
C       FMDUX(L,K)= DYP(L)*HP(L)*AH(L,K)*DXU1(L,K)
C       FMDVY(L,K)=-DXP(L)*HP(L)*AH(L,K)*DXU1(L,K)
C       FMDUY(L,K)=0.5*(DXU(L)+DXU(LS ))*H1C(L)*AHC(L,K)*DYU1(L,K)
C       FMDVX(L,K)=0.5*(DYV(L)+DYV(L-1))*H1C(L)*AHC(L,K)*DYU1(L,K)
C      ENDDO
C      ENDDO
C
C **  ADJUST CORNER FLUXES FOR SIDE WALL BOUNDARY CELLS
C
      CWALL=0.0
      IF(ISHDMF.EQ.2) CWALL=0.01
C
C      DO K=1,KC
C      DO L=2,LA
C        IF(SVBO(L).LT.0.5.OR.SVBO(L-1).LT.0.5)THEN
C        FMDUY(L,K)=DXU(L)*H1C(L)*CWALL*ABS(U(L,K))*U(L,K)
C        ENDIF
C      ENDDO
C      ENDDO
C
C      DO K=1,KC
C      DO L=2,LA
C        LS=LSC(L)      
C        IF(SUBO(L).LT.0.5.OR.SUBO(LS).LT.0.5)THEN
C        FMDVX(L,K)=DYV(L)*H1C(L)*CWALL*ABS(V(L,K))*V(L,K)
C        ENDIF
C      ENDDO
C      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        LN=LNC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        LS=LSC(L) 
        TMPVAL=1.+SUBO(L)+SVBO(LN)+SUBO(LN)*SVBO(LNW)
        H1CN=(H1P(L)+SUBO(L)*H1P(L-1)+SVBO(LN)*H1P(LN)
     &               +SUBO(LN)*SVBO(LNW)*H1P(LNW))/TMPVAL
        TMPVAL=1.+SUBO(L+1)+SVBO(L)+SUBO(LSE)*SVBO(L+1)
        H1CE=(H1P(L)+SUBO(L+1)*H1P(L+1)+SVBO(L)*H1P(LS)
     &               +SUBO(LSE)*SVBO(L+1)*H1P(LSE))/TMPVAL
        FMDUX(L,K)=DYP(L  )*HP(L  )*AH(L  ,K)*DXU1(L  ,K)
     &    -SUBO(L)*DYP(L-1)*HP(L-1)*AH(L-1,K)*DXU1(L-1,K)
        IF(SUBO(LN).GT.0.5) DXULN=DXU(LN)
        IF(SUBO(LN).LT.0.5)THEN
          DDXDDDY=2.*(DXU(L)-DXU(LS))/(DYU(L)+DYU(LS))
          DXULN=DXU(L)+0.5*DDXDDDY*DYU(L)
        ENDIF              
        IF(SUBO(LS).GT.0.5) DXULS=DXU(LS)
        IF(SUBO(LS).LT.0.5)THEN
          DDXDDDY=2.*(DXU(LN)-DXU(L))/(DYU(LN)+DYU(L))
          DXULS=DXU(L)-0.5*DDXDDDY*DYU(L)
        ENDIF              
        FMDUY(L,K)=0.5*(DXULN+DXU(L))*H1CN*AHNN(L,K)*SXY2NN(L,K)
     &            -0.5*(DXU(L)+DXULS)*H1C(L)*AHCC(L,K)*SXY2CC(L,K)    
        FMDVY(L,K)=-DXP(L )*HP(L )*AH(L ,K)*DXU1(L ,K)
     &     +SVBO(L)*DXP(LS)*HP(LS)*AH(LS,K)*DXU1(LS,K)
        IF(SVBO(L+1).GT.0.5) DYVLP=DYV(L+1)
        IF(SVBO(L+1).LT.0.5)THEN
          DDYDDDX=2.*(DYV(L)-DYV(L-1))/(DXV(L)+DXV(L-1))
          DYVLP=DYV(L)+0.5*DDYDDDX*DXV(L)
        ENDIF              
        IF(SVBO(L-1).GT.0.5) DYVLM=DYV(L-1)
        IF(SVBO(L-1).LT.0.5)THEN
          DDYDDDX=2.*(DYV(L+1)-DYV(L))/(DXV(L+1)+DXV(L))
          DYVLM=DXV(L)-0.5*DDYDDDX*DXV(L)
        ENDIF              
        FMDVX(L,K)=0.5*(DYVLP+DYV(L))*H1CE*AHEE(L,K)*SXY2EE(L,K)
     &            -0.5*(DYV(L)+DYVLM)*H1C(L)*AHCC(L,K)*SXY2CC(L,K)
      ENDDO
      ENDDO

      IF(N.EQ.2)THEN
      OPEN(1,FILE='AHD2.DIA')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='AHD2.DIA')
      DO L=2,LA
      I=IL(L)
      J=JL(L)
      DO K=1,KC
       WRITE(1,1111)N,I,J,K,FMDUX(L,K),FMDVY(L,K),FMDUY(L,K),
     &              FMDVX(L,K),AHC(L,K),DYU1(L,K)
      ENDDO
      ENDDO
      CLOSE(1)
      ENDIF
C
 1111 FORMAT(4I5,6E13.4)
 1112 FORMAT(2I5,8E13.4)
C
C**********************************************************************C
C
      RETURN
      END
