C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV9C (ISTL)
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
C ** SUBROUTINE CALPUV9 CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOWS WITH PROVISIONS FOR WETTING
C ** AND DRYING OF CELLS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION QCHANUT(NCHANM),QCHANVT(NCHANM)
      DIMENSION QSUMTMP(LCM)
C      DIMENSION XSMLC(LCM-2),XLSMLC(LCM-2),XUSMLC(LCM-2)
C      DIMENSION WSMLC(LCM*LCM)
C      DIMENSION IWSMLC(2*LCM)
C      DIMENSION ASMLC(1,LCM)
C      DIMENSION BSMLC(1)
C
C**********************************************************************C
C
C      WRITE(6,6000)N
C 6000 FORMAT(' CALLED CALPUV9, N = ',I10)
C
      IF(N.EQ.1.AND.ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(N.EQ.2.AND.ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
C**********************************************************************C
C
       DELT=DT2
       DELTD2=DT
      IF(ISTL.EQ.2)THEN
       DELT=DT
       DELTD2=0.5*DT
      ENDIF
      DELTI=1./DELT
C
C**********************************************************************C
C
C **  SET SWITCHES FOR DRYING AND WETTING  
C
C----------------------------------------------------------------------C
C     
      ITERHP=0
      NCORDRY=0
      ICORDRY=0
      DO L=1,LC
      ISCDRY(L)=0
      ENDDO
C
C**********************************************************************C
C
C **  INITIALIZE SUBGRID SCALE CHANNEL INTERACTIONS
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          QCHANUT(NMD)=QCHANUN(NMD)
          QCHANVT(NMD)=QCHANVN(NMD)
        ENDDO
        IF(ISTL.EQ.3)THEN
        DO NMD=1,MDCHH
          QCHANUN(NMD)=QCHANU(NMD)
          QCHANVN(NMD)=QCHANV(NMD)
        ENDDO
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
C
      IF(BSC.GT.1.E-6) CALL CALEBI
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL PRESSURE GRADIENTS 
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C **  SNLPX=SNLPX*GID2*DYU & SNLPY=SNLPY*GID2*DXV
C
C----------------------------------------------------------------------C
C
      IF(BSC.GT.1.E-6)THEN 
      DO L=2,LA
      LS=LSC(L)      
      FPGXE(L)=ROLD*FPGXE(L)+RNEW*(
     &           -SBX(L)*HU(L)*((BI2(L)+BI2(L-1))*(HP(L)-HP(L-1))
     &           +2.*HU(L)*(BI1(L)-BI1(L-1))
     &           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1))) )
      FPGYE(L)=ROLD*FPGYE(L)+RNEW*(
     &           -SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
     &           +2.*HV(L)*(BI1(L)-BI1(LS))
     &           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS))) )
      ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
        DO L=2,LA
        HUTMP(L)=0.5*(HU(L)+H1U(L))
        HVTMP(L)=0.5*(HV(L)+H1V(L))
        ENDDO
       ELSE
        DO L=2,LA
        HUTMP(L)=HU(L)
        HVTMP(L)=HV(L)
        ENDDO
      ENDIF
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDY1E(L)
     &         -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P1(L)-P1(L-1))
     &         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-RITB1*TBX1(L))
     &         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      FVHDXE(L)=VHDX1E(L)
     &         -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P1(L)-P1(LS ))
     &         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-RITB1*TBY1(L))
     &         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ENDDO
C
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='FUV.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='FUV1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='FUV2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
C
C**********************************************************************C
C
C **  SET IMPLICIT BOTTOM AND VEGETATION DRAG AS APPROPRIATE
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      RCX(L)=1./(1.+FBODYFXI(L))
      RCY(L)=1./(1.+FBODYFYI(L))
      ENDDO
C    
      RCX(1)=0.
      RCY(1)=0.
      RCX(LC)=0.
      RCY(LC)=0.
C
C * SINGLE LAYER NO VEGETATION
C
      IF(KC.EQ.1.AND.ISVEG.EQ.0)THEN
        DO L=2,LA
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * SINGLE LAYER WITH VEGETATION
C
      IF(KC.EQ.1.AND.ISVEG.GE.1)THEN
        DO L=2,LA
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1))
     &    +DELT*FXVEGE(L) )
C     &    +RITB*DELT*FXVEGE(L) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
     &    +DELT*FYVEGE(L) )
C     &    +RITB*DELT*FYVEGE(L) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C
C * MULTIPLE LAYERS NO VEGETATION
C
      IF(KC.GT.1.AND.ISVEG.EQ.0)THEN
        DO L=2,LA
	  TMPX=1.0
	  TMPY=1.0
	  IF(UHE(L).NE.0.0) TMPX=U(L,1)*HU(L)/UHE(L)
	  IF(VHE(L).NE.0.0) TMPY=V(L,1)*HV(L)/VHE(L)
        RCX(L)=1./( 1.+FBODYFXI(L)
     &  +TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &  +TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * MULTIPLE LAYERS WITH VEGETATION
C
      IF(KC.GT.1.AND.ISVEG.GE.1)THEN
        DO L=2,LA
	  TMPX=1.0
	  TMPY=1.0
	  IF(UHE(L).NE.0.0) TMPX=U(L,1)*HU(L)/UHE(L)
	  IF(VHE(L).NE.0.0) TMPY=V(L,1)*HV(L)/VHE(L)
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1))
     &    +DELT*FXVEGE(L) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
     &    +DELT*FYVEGE(L) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * MULTIPLE LAYERS WITH VEGETATION  ORIGINAL
C
CORG      IF(KC.GT.1.AND.ISVEG.GE.1)THEN
CORG        DO L=2,LA
C        RCX(L)=1./( 1.+RITB*DELT*FXVEGE(L) )
C        RCY(L)=1./( 1.+RITB*DELT*FYVEGE(L) )
CORG        RCX(L)=1./( 1.+DELT*FXVEGE(L) )
CORG        RCY(L)=1./( 1.+DELT*FYVEGE(L) )
CORG        FUHDYE(L)=FUHDYE(L)*RCX(L)
CORG        FVHDXE(L)=FVHDXE(L)*RCY(L)
CORG        ENDDO
CORG      ENDIF
C
C**********************************************************************C
C
C **  RESET BOUNDARY CONDITIONS SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
        SUB(L)=SUBO(L)
        SVB(L)=SVBO(L)
        SBX(L)=SBXO(L)
        SBY(L)=SBYO(L)
        SUB(L+1)=SUBO(L+1)
        SVB(LN)=SVBO(LN)
        SBX(L+1)=SBXO(L+1)
        SBY(LN)=SBYO(LN)
      ENDDO
C
      DO L=1,LC
      FP(L)=0.
      FP1(L)=0.
      ENDDO
C
C**********************************************************************C
C
C **  BLOCK DRY CELLS WITH NO POTENTIAL FOR WETTING
C
C----------------------------------------------------------------------C
C
C      DO L=2,LA
C       LN=LNC(L)
C       LS=LSC(L)
C       IF(H1P(L).LT.HWET)THEN
C         IF(SUBO(L).GT.0.5)THEN
C           IF(P(L-1).LT.P(L).OR.HP(L-1).LT.HDRY) SUB(L)=0.0
C           IF(P(L-1).LT.P(L)) SUB(L)=0.0
C         ENDIF
C         IF(SUBO(L+1).GT.0.5)THEN
C           IF(P(L+1).LT.P(L).OR.HP(L+1).LT.HDRY) SUB(L+1)=0.0
C           IF(P(L+1).LT.P(L)) SUB(L+1)=0.0
C         ENDIF
C         IF(SVBO(L).GT.0.5)THEN
C           IF(P(LS).LT.P(L).OR.HP(LS).LT.HDRY) SVB(L)=0.0
C           IF(P(LS).LT.P(L)) SVB(L)=0.0
C         ENDIF
C         IF(SVBO(LN).GT.0.5)THEN
C           IF(P(LN).LT.P(L).OR.HP(LN).LT.HDRY) SVB(LN)=0.0
C           IF(P(LN).LT.P(L)) SVB(LN)=0.0
C         ENDIF
C       ENDIF
C      ENDDO
C
C**********************************************************************C
C
C **  SET OPEN BOUNDARY SURFACE ELEVATIONS 
C
C----------------------------------------------------------------------C
C
      IF(ISDYNSTP.EQ.0)THEN
        TN=DT*FLOAT(N)+TCON*TBEGIN
      ELSE
        TN=TIMESEC
      ENDIF
C
      DO M=1,MTIDE
      TM=MOD(TN,TCP(M))
      TM=PI2*TM/TCP(M)
      CCCOS(M)=COS(TM)
      SSSIN(M)=SIN(TM)
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBW
      L=LPBW(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      FP1(L)=PSERT(NPSERW(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
      ENDDO
      CET=0.5*DELTD2*G*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/TMP
       CE(L)=-CET
       FP1(L)=CET*(2.*FP1(L)
     & -SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/TMP
      ELSE
       FP1(L+1)=CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.      
      FP1(L)=PSERT(NPSERE(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
      ENDDO
      CWT=0.5*DELTD2*G*HRUO(L  )*RCX(L)*HUTMP(L  )
      IF(ISPBE(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/TMP
       CW(L)=-CWT
       FP1(L)=CWT*(2.*FP1(L)
     & +SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/TMP
      ELSE
       FP1(L-1)=CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LN=LNC(L)
      FP1(L)=PSERT(NPSERS(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
      ENDDO
      CNT=0.5*DELTD2*G*HRVO(LN )*RCY(LN)*HVTMP(LN )
      IF(ISPBS(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HMV(LN))*DYIV(LN)                  ! *** HMV Fixed at IC
       CC(L)=CNT*(1.+TMP)/TMP
       CN(L)=-CNT
       FP1(L)=CNT*(2.*FP1(L)                               
     & -SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/TMP    ! *** HMV Fixed at IC 
      ELSE
       FP1(LN)=CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LS=LSC(L)
      FP1(L)=PSERT(NPSERN(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
      ENDDO
      CST=0.5*DELTD2*G*HRVO(L  )*RCY(L)*HVTMP(L  )
      IF(ISPBN(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HMV(L))*DYIV(L)                    ! *** HMV Fixed at IC (ORIGINAL)
       ! TMP=DELTD2*SQRT(G*HV(L))*DYIV(L)                   ! *** Variable HV
       CC(L)=CST*(1.+TMP)/TMP
       CS(L)=-CST
       FP1(L)=CST*(2.*FP1(L)                                ! *** HMV Fixed at IC (ORIGINAL)
     & +SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/TMP      
       !FP1(L)=CST*(2.*FP1(L)                               ! *** Variable HV
       !    & +SQRT(G*HV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/TMP
      ELSE
       FP1(LS)=CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  ADJUST VOLUME SOURCE AND SINKS
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.EQ.0)THEN
C
      DO L=2,LA
      IF(QSUME(L).LE.0.)THEN
        IF(H1P(L).LE.HDRY)THEN
          QSUMTMP(L)=0.
         ELSE
          QSUMTMP(L)=-(H1P(L)-HDRY)*DXYP(L)*DELTI
          QSUMTMP(L)=MAX(QSUMTMP(L),QSUME(L))
        ENDIF
       ELSE
        QSUMTMP(L)=QSUME(L)
      ENDIF
      ENDDO
C 
      DO L=2,LA
       DIFQVOL=QSUME(L)-QSUMTMP(L)
       DO K=1,KC
       QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
       ENDDO
       QSUME(L)=QSUMTMP(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C 
C **  ADJUST SOURCES AND SINKS ESTIMATING SURFACE AND GROUNDWATER
C **  AVAILABLE FOR EVAPOTRANSPIRATON AND INFILTRATION
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.GE.1)THEN
C
      DO L=2,LA
      RIFTR(L)=0.
      EVAPSW(L)=0.
      EVAPGW(L)=0.
      IF(H1P(L).GT.HDRY)THEN
C       APPLY MAXIMUM ET
        IF(EVAPCVT.LT.0.)THEN 
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/     
     &              (1.+0.00412*TEM(L,KC))))            
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)       
        ENDIF                                            
        EVAPSW(L)=EVAPT(L)*DXYP(L)
        RIFTR(L)=0.
C       CALCULATE DEPTH OF ACTIVE GROUNDWATER ELEV BELOW SURFACE
        DTAGW=BELV(L)-AGWELV(L)
        IF(DTAGW.GT.0.0)THEN
C         INFLITRATION CAN OCCUR, CALCULATE LIMITING RATE TO BRING
C         GW ELEV TO SOIL SURFACE
          RIFTRL=RNPOR*DTAGW*DELTI
C         SET RIFTRL TO MIN OF LIMITING RATE OR ACTUAL RATE
          RIFTRL=MIN(RIFTRM,RIFTRL)
C         ESTIMATE RATE BASED ON AVAILABLE SURFACE WATER 
          RAVAIL=(H1P(L)-HDRY)*DELTI-EVAPT(L)
C         SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
          RIFTRL=MIN(RAVAIL,RIFTRL)
C         CONVERT TO VOLUME FLOW UNITS
          RIFTR(L)=RIFTRL*DXYP(L)         
        ENDIF
C       ADJUST VOLUME OUTFLOWS OF WET CELLS
        IF(QSUME(L).LT.0.0)THEN
          QSUMIET=RIFTR(L)+EVAPSW(L)
          QEAVAIL=DXYP(L)*(H1P(L)-HDRY)*DELTI-QSUMIET
          QEAVAIL=MAX(QEAVAIL,0.0)
          QEAVAIL=-QEAVAIL
          QSUMTMP(L)=MAX(QSUME(L),QEAVAIL)
         ELSE
          QSUMTMP(L)=QSUME(L)
        ENDIF         
       ELSE
        RIFTR(L)=0.
        EVAPSW(L)=0.
        QSUMTMP(L)=MAX(QSUME(L),0.0)       
      ENDIF
      ENDDO
C
      DO L=2,LA
      DIFQVOL=QSUME(L)-QSUMTMP(L)
      DO K=1,KC
      QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
      ENDDO
      QSUME(L)=QSUMTMP(L)
      ENDDO
C       
      ENDIF
C
C**********************************************************************C
C
C **  SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR 
C **  NON BOUNDARY POINTS
C **  HRU=HMU*DYU/DXU & HRV=HMV*DXV/DYV 
C **  DXYIP=1/(DXP*DYP)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*P1(L)
     &      -0.5*G*(UHDY1E(L+1)-UHDY1E(L)
     &             +VHDX1E(LN )-VHDX1E(L)) )
      ENDDO
C
C**********************************************************************C
C
C **  ADVANCE EXTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        UHDY2E(L)=UHDY1E(L)
        VHDX2E(L)=VHDX1E(L)
        UHDY1E(L)=UHDYE(L)
        VHDX1E(L)=VHDXE(L)
        U1V(L)=UV(L)
        V1U(L)=VU(L)
C       DELP=P(L)-P1(L)
        P1(L)=P(L)
C       P(L)=P(L)+DELP
        H1U(L)=HU(L)
        H1V(L)=HV(L)
        H1UI(L)=HUI(L)
        H1VI(L)=HVI(L)
        H2P(L)=H1P(L)
        H1P(L)=HP(L)
        AGWELV2(L)=AGWELV1(L)
        AGWELV1(L)=AGWELV(L)
        ENDDO
      ENDIF
C
CC     DO L=2,LA
CC      PAM(L)=P(L)
CC      ENDDO
C
C**********************************************************************C
C
C **  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
C **  HOST-GUEST CHANNAL INTERACTION FOR NON BOUNDARY POINTS 
C
C **  REENTER AT 1000 FOR WETTING-DRYING CORRECTION AND CHANNEL 
C **  INTERACTION
C
C----------------------------------------------------------------------C
C
 1000 CONTINUE
C
      DO L=2,LA
      LN=LNC(L)
      FP(L)=FP1(L)-0.5*G*SPB(L)*
     &      ( SUB(L+1)*FUHDYE(L+1)-SUB(L)*FUHDYE(L)
     &       +SVB(LN )*FVHDXE(LN )-SVB(L)*FVHDXE(L)
     &       -2.0*QSUME(L) )
CC      P(L)=0.
      ENDDO
C
      IF(ISGWIE.GE.1)THEN
        DO L=2,LA
        FP(L)=FP(L)-G*SPB(L)*(RIFTR(L)+EVAPSW(L))
        ENDDO
      ENDIF
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='FP.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='FP1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='FP2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      LN=LNC(L)
      C1=-0.5*DELTD2*G*SPB(L)
      CS(L)=C1*SVB(L  )*HRVO(L  )*RCY(L  )*HVTMP(L  )
C    &     +(1.-SPB(L))*CS(L)
      CW(L)=C1*SUB(L  )*HRUO(L  )*RCX(L  )*HUTMP(L  )
C    &     +(1.-SPB(L))*CW(L)
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
C    &     +(1.-SPB(L))*CE(L)
      CN(L)=C1*SVB(LN )*HRVO(LN )*RCY(LN )*HVTMP(LN )
C    &     +(1.-SPB(L))*CN(L)
      CC(L)=SPB(L)*(DELTI*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
C    &     +(1.-SPB(L))*CC(L)
      ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  INSERT IMPLICT SUB-GRID SCALE CHANNEL INTERACTIONS
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
C
        IF(ISTL.EQ.3)THEN
C
        DO NMD=1,MDCHH
          CCCCHU(NMD)=0.
          CCCCHV(NMD)=0.
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
          GBELVTMP=G*BELV(LHOST)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            IACTIVE=0
            IF(P(LCHNU).GT.P(LHOST).AND.H2P(LCHNU).GT.HWET)THEN
              HCHNCOR=HWET
              IF(ISCDRY(LCHNU).EQ.0) IACTIVE=1
            ENDIF
            IF(P(LHOST).GT.P(LCHNU).AND.H2P(LHOST).GT.HDRY)THEN
              HCHNCOR=HDRY
              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            IF(IACTIVE.EQ.1)THEN
              WCHAN=DXP(LCHNU)
C              RLCHN=0.5*DYP(LCHNU)+0.25*DYP(LHOST)
C              HCHAN=0.5*DYP(LCHNU)*HP(LCHNU)+0.25*DYP(LHOST)*HP(LHOST)
              RLCHN=0.5*DYP(LCHNU)+CHANLEN(NMD)
              HCHAN=0.5*DYP(LCHNU)*H2P(LCHNU)+CHANLEN(NMD)*H2P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBX(LCHNU)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*QCHANUT(NMD))
                HCHANT=HCHAN-HCHNCOR
                HCHANT=MAX(HCHANT,0.0)
                CCCCHV(NMD)=DELT*HCHANT/RLCHN
              ENDIF
            ENDIF
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            IACTIVE=0
            IF(P(LCHNV).GT.P(LHOST).AND.H2P(LCHNV).GT.HWET)THEN
              HCHNCOR=HWET
              IF(ISCDRY(LCHNV).EQ.0) IACTIVE=1
            ENDIF
            IF(P(LHOST).GT.P(LCHNV).AND.H2P(LHOST).GT.HDRY)THEN
              HCHNCOR=HDRY
              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            IF(IACTIVE.EQ.1)THEN
              WCHAN=DYP(LCHNV)
C              RLCHN=0.5*DXP(LCHNV)+0.25*DXP(LHOST)
C              HCHAN=0.5*DXP(LCHNV)*HP(LCHNV)+0.25*DXP(LHOST)*HP(LHOST)
              RLCHN=0.5*DXP(LCHNV)+CHANLEN(NMD)
              HCHAN=0.5*DXP(LCHNV)*H2P(LCHNV)+CHANLEN(NMD)*H2P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBY(LCHNV)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*QCHANVT(NMD))
                HCHANT=HCHAN-HCHNCOR
                HCHANT=MAX(HCHANT,0.0)
                CCCCHV(NMD)=DELT*HCHANT/RLCHN
              ENDIF
            ENDIF
          ENDIF
          CC(LHOST)=CC(LHOST)+G*CCCCHU(NMD)*CCCCHV(NMD)
          CCCCHH(NMD)=CCCCHU(NMD)*CCCCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            CC(LCHNU)=CC(LCHNU)+G*CCCCHU(NMD)*CCCCHV(NMD)
            FP(LHOST)=FP(LHOST)+G*CCCCHU(NMD)*QCHANUT(NMD)
            FP(LCHNU)=FP(LCHNU)-G*CCCCHU(NMD)*QCHANUT(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            CC(LCHNV)=CC(LCHNV)+G*CCCCHU(NMD)*CCCCHV(NMD)
            FP(LHOST)=FP(LHOST)+G*CCCCHU(NMD)*QCHANVT(NMD)
            FP(LCHNV)=FP(LCHNV)-G*CCCCHU(NMD)*QCHANVT(NMD)
          ENDIF
        ENDDO
C
        ELSE
C
        DO NMD=1,MDCHH
          CCCCHU(NMD)=0.
          CCCCHV(NMD)=0.
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
          GBELVTMP=G*BELV(LHOST)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            IACTIVE=0
            IF(P(LCHNU).GT.P(LHOST).AND.H1P(LCHNU).GT.HWET)THEN
              HCHNCOR=HWET
              IF(ISCDRY(LCHNU).EQ.0) IACTIVE=1
            ENDIF
            IF(P(LHOST).GT.P(LCHNU).AND.H1P(LHOST).GT.HDRY)THEN
              HCHNCOR=HDRY
              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            IF(IACTIVE.EQ.1)THEN
              WCHAN=DXP(LCHNU)
C              RLCHN=0.5*DYP(LCHNU)+0.25*DYP(LHOST)
C              HCHAN=0.5*DYP(LCHNU)*HP(LCHNU)+0.25*DYP(LHOST)*HP(LHOST)
              RLCHN=0.5*DYP(LCHNU)+CHANLEN(NMD)
              HCHAN=0.5*DYP(LCHNU)*H1P(LCHNU)+CHANLEN(NMD)*H1P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBX(LCHNU)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*QCHANUT(NMD))
                HCHANT=HCHAN-HCHNCOR
                HCHANT=MAX(HCHANT,0.0)
                CCCCHV(NMD)=DELT*HCHANT/RLCHN
              ENDIF
            ENDIF
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            IACTIVE=0
            IF(P(LCHNV).GT.P(LHOST).AND.H1P(LCHNV).GT.HWET)THEN
              HCHNCOR=HWET
              IF(ISCDRY(LCHNV).EQ.0) IACTIVE=1
            ENDIF
            IF(P(LHOST).GT.P(LCHNV).AND.H1P(LHOST).GT.HDRY)THEN
              HCHNCOR=HDRY
              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            IF(IACTIVE.EQ.1)THEN
              WCHAN=DYP(LCHNV)
C              RLCHN=0.5*DXP(LCHNV)+0.25*DXP(LHOST)
C              HCHAN=0.5*DXP(LCHNV)*HP(LCHNV)+0.25*DXP(LHOST)*HP(LHOST)
              RLCHN=0.5*DXP(LCHNV)+CHANLEN(NMD)
              HCHAN=0.5*DXP(LCHNV)*H1P(LCHNV)+CHANLEN(NMD)*H1P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBY(LCHNV)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*QCHANVT(NMD))
                HCHANT=HCHAN-HCHNCOR
                HCHANT=MAX(HCHANT,0.0)
                CCCCHV(NMD)=DELT*HCHANT/RLCHN
              ENDIF
            ENDIF
          ENDIF
          CC(LHOST)=CC(LHOST)+G*CCCCHU(NMD)*CCCCHV(NMD)
          CCCCHH(NMD)=CCCCHU(NMD)*CCCCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            CC(LCHNU)=CC(LCHNU)+G*CCCCHU(NMD)*CCCCHV(NMD)
            FP(LHOST)=FP(LHOST)+G*CCCCHU(NMD)*QCHANUT(NMD)
            FP(LCHNU)=FP(LCHNU)-G*CCCCHU(NMD)*QCHANUT(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            CC(LCHNV)=CC(LCHNV)+G*CCCCHU(NMD)*CCCCHV(NMD)
            FP(LHOST)=FP(LHOST)+G*CCCCHU(NMD)*QCHANVT(NMD)
            FP(LCHNV)=FP(LCHNV)-G*CCCCHU(NMD)*QCHANVT(NMD)
          ENDIF
        ENDDO
C
        ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
C
C----------------------------------------------------------------------C
C
      CCMNM=1.E+18
      DO L=2,LA
        CCMNM=MIN(CCMNM,CC(L))
        FPTMP(L)=FP(L)
      ENDDO
C
      CCMNMI=1./CCMNM
C     WRITE(6,666) CCMNM,CCMNMI
C 666 FORMAT(' CCMIN, CCMINI = ',E12.4,1X,E12.4)
C
      DO LL=1,NPBW
      IF(ISPBW(LL).EQ.0)THEN
        L=LPBW(LL)
        CW(L+1)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBE
      IF(ISPBE(LL).EQ.0)THEN
        L=LPBE(LL)
        CE(L-1)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBS
      IF(ISPBS(LL).EQ.0)THEN
        L=LPBS(LL)
        LN=LNC(L)
        CS(LN)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBN
      IF(ISPBN(LL).EQ.0)THEN
        L=LPBN(LL)
        LS=LSC(L)
        CN(LS)=0.
      ENDIF
      ENDDO
C
      CC(1)=1.
      CC(LC)=1.
C
C----------------------------------------------------------------------C
C
C **  SCALE BY MINIMUM DIAGONAL
C
      IF(IRVEC.EQ.9)THEN
C
      DO L=2,LA
      CCS(L)=CS(L)*CCMNMI
      CCW(L)=CW(L)*CCMNMI
      CCE(L)=CE(L)*CCMNMI
      CCN(L)=CN(L)*CCMNMI
      CCC(L)=CC(L)*CCMNMI
      FPTMP(L)=FPTMP(L)*CCMNMI
      CCCI(L)=1./CCC(L)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          CCCCHH(NMD)=CCCCHU(NMD)*CCMNMI
        ENDDO
      ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  SCALE TO NORMAL FORM - NOTE NORMAL FORM CURRENTLY CANNOT BE USED
C **  WITH SUB-GRID SCALE CHANNEL MODEL
C
C      IF(IRVEC.EQ.99)THEN
C
C      DO L=2,LA
C        CCS(L)=CS(L)/SQRT( CC(L)*CC(LSC(L)) )
C        CCW(L)=CW(L)/SQRT( CC(L)*CC(L-1   ) )
C        CCE(L)=CE(L)/SQRT( CC(L)*CC(L+1   ) )
C        CCN(L)=CN(L)/SQRT( CC(L)*CC(LNC(L)) )
C        CCC(L)=1.
C        FPTMP(L)=FPTMP(L)/SQRT( CC(L) )
C        P(L)=P(L)*SQRT( CC(L) )
C        CCCI(L)=1.
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
C **  CALL EQUATION SOLVER
C
      IF(MDCHH.EQ.0) CALL CONGRAD (ISTL)
      IF(MDCHH.GE.1) CALL CONGRADC (ISTL)
C
C----------------------------------------------------------------------C
C
C BEGIN CONSTRAINED MINIMIZATION SOLUTION
C
C      DO L=2,LA
C       XSMLC(L-1)=P(L)
CJH       XLSMLC(L-1)=G*(BELV(L)+HDRY)
C       XLSMLC(L-1)=G*BELV(L)
C       XUSMLC(L-1)=1.E4
C      ENDDO
C
C      IF(IRVEC.EQ.99)THEN
C        DO L=2,LA
C          XLSMLC(L-1)=XLSMLC(L-1)*SQRT( CC(L) )
C          XUSMLC(L-1)=XUSMLC(L-1)*SQRT( CC(L) )
C        ENDDO
C      ENDIF
C
C      DO L=2,LA
C       XSMLC(L-1)=MAX(XSMLC(L-1),XLSMLC(L-1))
C      ENDDO
C
C      NTMP=LC-2
C      MTMP=0
C      MEQ=0
CJH     ASMLC=ASMLC(1,1)
C      LDA=1
CJH     BSMLC=BSMLC(1)
C      IPRINT=0
C      NFMAX=ITERM
CJH     IWSMLC=IWSM(LIW)
C      LIWSM=4+MTMP+2*NTMP
C      LWSM=3+MTMP+NTMP*(NTMP+16)
CJH     WSMLC=IWSM(LIW)
C      CALL SMLC01(NTMP,MTMP,MEQ,ASMLC,LDA,BSMLC,XLSMLC,
C     &     XUSMLC,XSMLC,RSQM,IPRINT,NFMAX,IWSMLC,LIWSM,WSMLC,LWSM)
C      DO L=2,LA
C       P(L)=XSMLC(L-1)
C      ENDDO
C
C END CONSTRAINED MINIMIZATION SOLUTION
C
C----------------------------------------------------------------------C
C
C **  REVERSE SCALE IF NORMAL FORM SOLUTION IS USED
C
C      IF(IRVEC.EQ.99)THEN
C
C      DO L=2,LA
C      P(L)=P(L)/SQRT( CC(L) )
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
C ** DIAGNOSTICS
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
C        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
        WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
        ENDDO
        IF(MDCHH.GE.1)THEN
          DO NMD=1,MDCHH
            WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                                    CCCCHV(NMD)
          ENDDO
        ENDIF
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='EQCOEF1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          SURFTMP=GI*P(L)
C          WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
          WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
          ENDDO
          IF(MDCHH.GE.1)THEN
            DO NMD=1,MDCHH
             WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                                     CCCCHV(NMD)
            ENDDO
          ENDIF
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='EQCOEF2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          SURFTMP=GI*P(L)
C          WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
          WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
          ENDDO
          IF(MDCHH.GE.1)THEN
            DO NMD=1,MDCHH
             WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                                     CCCCHV(NMD)
            ENDDO
          ENDIF
          CLOSE(1)
        ENDIF
      ENDIF
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='EQTERM1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='EQTERM2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
 1001 FORMAT(2I5,10(1X,E12.4))
 1002 FORMAT(3I4,10(1X,E9.2))
C
C**********************************************************************C
C
C **  CALCULATE UHEX AND VHEX AND TOTAL DEPTHS AT TIME LEVEL (N+1)
C **  HRU=SUB*DYU/DXU & HRV=SVB*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)      
      UTMP=SUB(L)*( FUHDYE(L)
     &            -DELTD2*HRUO(L)*RCX(L)*HUTMP(L)*(P(L)-P(L-1)) )
      VTMP=SVB(L)*( FVHDXE(L)
     &            -DELTD2*HRVO(L)*RCY(L)*HVTMP(L)*(P(L)-P(LS )) )
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE NEW SUB-GRID SCALE CHANNEL EXCHANGE FLOWS
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            QCHANU(NMD)=CCCCHU(NMD)*QCHANUT(NMD)
     &               -CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNU))
            QCHANV(NMD)=0.
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            QCHANV(NMD)=CCCCHU(NMD)*QCHANVT(NMD)
     &               -CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNV))
            QCHANU(NMD)=0.
          ENDIF
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL 
C **  TRANSPORTS AT (N+1)
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     &       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP
     &                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        ENDDO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP
     &                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  ADD CHANNEL INTERACTION EXCHANGES
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        IF(MDCHTYP(NMD).EQ.1)THEN
          HP(LMDCHH(NMD))=HP(LMDCHH(NMD))
     &                   +DELT*DXYIP(LMDCHH(NMD))*QCHANU(NMD)
          HP(LMDCHU(NMD))=HP(LMDCHU(NMD))
     &                   -DELT*DXYIP(LMDCHU(NMD))*QCHANU(NMD)
        ENDIF            
        IF(MDCHTYP(NMD).EQ.2)THEN
          HP(LMDCHH(NMD))=HP(LMDCHH(NMD))
     &                   +DELT*DXYIP(LMDCHH(NMD))*QCHANV(NMD)
          HP(LMDCHV(NMD))=HP(LMDCHV(NMD))
     &                   -DELT*DXYIP(LMDCHV(NMD))*QCHANV(NMD)
        ENDIF            
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      ENDDO
C
      DO L=2,LA
      LS=LSC(L)      
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LS )*DYP(LS)*HP(LS ))
     &           /(DXV(L)*DYV(L))
      ENDDO
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      ENDDO
C
C
C**********************************************************************C
C
C **  CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.GT.0)THEN
      OPEN(1,FILE='DRYWET.LOG',POSITION='APPEND',STATUS='UNKNOWN')
C
      ICORDRY=0
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
C
      IF(HP(L).LE.HDRY)THEN
        IF(ISCDRY(L).EQ.0)THEN
          ISCDRY(L)=1
          ICORDRY=1
C         WRITE(1,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(6,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(8,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
        ENDIF
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
      ENDIF
C
      ENDDO
C
      IF(ICORDRY.EQ.1)THEN
        NCORDRY=NCORDRY+1
        GOTO 1000
      ENDIF
C
      CLOSE(1)
      ENDIF
C
cSO   WRITE(6,6960)NCORDRY
CTMP      WRITE(8,6960)NCORDRY
C
 6960 FORMAT(' NCORDRY =', I5)
 6961 FORMAT(' UNSTABLE, NCORDRY =', I5)
C
 9999 CONTINUE
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      ENDDO
C
      DO L=2,LA
      LS=LSC(L)      
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LS )*DYP(LS )*HP(LS ))   ! 2014-01 BUG FIX
     &           /(DXV(L)*DYV(L))
      ENDDO
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      ENDDO
C
C**********************************************************************C
C
C **  PERFORM UPDATE ON GROUNDWATER ELEVATION
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.GE.1)THEN 
C
        DO L=2,LA
        QSUM(L,KC)=QSUM(L,KC)-EVAPSW(L)
        QSUM(L,1 )=QSUM(L,1 )-RIFTR(L)
        ENDDO
C
C       INFILTRATION STEP
C
        RNPORI=1./RNPOR
        IF(ISTL.EQ.3)THEN
          DO L=2,LA
          AGWELV(L)=AGWELV2(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          ENDDO
         ELSE
          DO L=2,LA
          AGWELV(L)=AGWELV1(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          ENDDO
        ENDIF
        DO L=2,LA
        AGWELV(L)=MIN(AGWELV(L),BELV(L))
        ENDDO
C
C       ET STEP
C
        DO L=2,LA
        IF(EVAPCVT.LT.0.)THEN                    
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/ 
     &              (1.+0.00412*TEM(L,KC))))      
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)       
        ENDIF                                    
        ETGWTMP=EVAPT(L)-EVAPSW(L)*DXYIP(L)
        ETGWTMP=MAX(ETGWTMP,0.0)
        ETGWAVL=RNPOR*DELTI*(AGWELV(L)-BELAGW(L))
        ETGWAVL=MAX(ETGWAVL,0.0)
        ETGWTMP=MIN(ETGWTMP,ETGWAVL)
        EVAPGW(L)=ETGWTMP*DXYP(L)
        ENDDO
        DO L=2,LA
        AGWELV(L)=AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)
        ENDDO
        DO L=2,LA
        AGWELV(L)=MAX(AGWELV(L),BELAGW(L))
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CHECK FOR NEGATIVE DEPTHS
C
C----------------------------------------------------------------------C
C
      IF(ISNEGH.GE.1)THEN
      INEGFLG=0
C
      DO L=2,LA
      IF(HP(L).LT.0.)THEN
      INEGFLG=1
      LN=LNC(L)
      LS=LSC(L)
      WRITE (6,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
      WRITE (8,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
      ENDIF
      ENDDO
C
C      DO L=2,LA
C      IF(HU(L).LT.0.)THEN
C      INEGFLG=1
C      LN=LNC(L)
C      LS=LSC(L)
C      WRITE (6,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (6,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (6,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (6,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (6,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      WRITE (8,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (8,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (8,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (8,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (8,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      ENDIF
C      ENDDO
C
C
C      DO L=2,LA
C      IF(HV(L).LT.0.)THEN
C      INEGFLG=1
C      LN=LNC(L)
C      LS=LSC(L)
C      WRITE (6,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (6,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (6,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (6,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (6,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      WRITE (8,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (8,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (8,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (8,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (8,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      ENDIF
C      ENDDO
C
      ENDIF
C
      IF(ISNEGH.EQ.2)THEN
      IF(INEGFLG.EQ.1)THEN
C
        IF(MDCHH.GT.0)THEN
        DO NMD=1,MDCHH
          WRITE(8,8000)
          LHOST=LMDCHH(NMD)
          IHOST=IL(LHOST)        
          JHOST=JL(LHOST)        
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
          ICHNU=IL(LCHNU)
          JCHNU=JL(LCHNU)
          WRITE(8,8001)NMD,MDCHTYP(NMD),ICHNU,JCHNU,ISCDRY(LCHNU),
     &                 P(LCHNU),HP(LCHNU),QCHANU(NMD)
          WRITE(8,8002)IHOST,JHOST,P(LHOST),ISCDRY(LHOST),
     &                 HP(LHOST),QCHANUN(NMD)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
          ICHNV=IL(LCHNV)
          JCHNV=JL(LCHNV)
          WRITE(8,8001)NMD,MDCHTYP(NMD),ICHNV,JCHNV,ISCDRY(LCHNV),
     &                 P(LCHNV),HP(LCHNV),QCHANV(NMD)
          WRITE(8,8002)IHOST,JHOST,ISCDRY(LHOST),P(LHOST),
     &                 HP(LHOST),QCHANVN(NMD)
          ENDIF
        ENDDO
        ENDIF
C
        CALL RESTOUT(1)
C
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &               FPTMP(L),SURFTMP
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='CFLMAX.OUT')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='CFLMAX.OUT')
        DO L=2,LA
         WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
         WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
         WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
         WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
        ENDDO
        CLOSE(1)
C
        STOP
C
      ENDIF
      ENDIF
C
 1991 FORMAT(2I5,12F8.3)
 1992 FORMAT(10X,12F8.3)
 6060 FORMAT('  NEG DEPTH AT I,J =',3I4,'  HP,H1P,H2P =',3(2X,E12.4))
 6061 FORMAT('  NEG DEPTH AT I,J =',3I4,'  HUW,H1UW =',2(2X,E12.4))
 6062 FORMAT('  NEG DEPTH AT I,J =',3I4,'  HUE,H1UE =',2(2X,E12.4))
 6063 FORMAT('  NEG DEPTH AT I,J =',3I4,'  HVS,H1VS =',2(2X,E12.4))
 6064 FORMAT('  NEG DEPTH AT I,J =',3I4,'  HVN,H1VN =',2(2X,E12.4))
 8001 FORMAT(5I5,3E13.4)
 8002 FORMAT(10X,3I5,3E13.4)
 8000 FORMAT(' NMD MTYP   I    J   IDRY      P         H            Q')
C
C**********************************************************************C
C
 6910 FORMAT('  DRYING AT N,I,J =',I10,2I6,'  HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6911 FORMAT('  DRY W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6912 FORMAT('  DRY E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6913 FORMAT('  DRY S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6914 FORMAT('  DRY N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6920 FORMAT('  WETTING AT N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6921 FORMAT('  WET S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6922 FORMAT('  WET W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6923 FORMAT('  WET E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6924 FORMAT('  WET N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6930 FORMAT('  WET BY VOL  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6940 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6941 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HUE,HP,H1P ='
     &         ,3(2X,E12.4))
 6942 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HUW,HP,H1P ='
     &         ,3(2X,E12.4))
 6943 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HVS,HP,H1P ='
     &         ,3(2X,E12.4))
 6944 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HVN,HP,H1P ='
     &         ,3(2X,E12.4))
 6945 FORMAT('  RESOLVE NEG,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6950 FORMAT('  RESOLVE, NEG DEP N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
C
C**********************************************************************C
C
C **  CALCULATE THE EXTERNAL DIVERGENCE 
C
C----------------------------------------------------------------------C
C
      IF(ISDIVEX.EQ.1)THEN
C
      DIVEXMX=0.
      DIVEXMN=1000000.
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO L=2,LA
      IF(SPB(L).NE.0)THEN
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI
     &     +0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     &     +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))-QSUME(L)
     &     +RIFTR(L)+EVAPSW(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      ENDIF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      ENDIF
      ENDIF
      ENDDO
C
      ELSE
C
      DO L=2,LA
      IF(SPB(L).NE.0)THEN
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI
     &     +0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L)
     &     +RIFTR(L)+EVAPSW(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      ENDIF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      ENDIF
      ENDIF
      ENDDO
C
      ENDIF
C
      IMAX=IL(LMAX)
      JMAX=JL(LMAX)
      IMIN=IL(LMIN)
      JMIN=JL(LMIN)
C
      WRITE(6,6628)DIVEXMX,IMAX,JMAX
      WRITE(6,6629)DIVEXMN,IMIN,JMIN
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C----------------------------------------------------------------------C
C
  566 FORMAT('  I=',I5,3X,'J=',I5,3X,'HP=',F12.4)
 6628 FORMAT('  DIVEXMX=',E13.5,5X,2I10)
 6629 FORMAT('  DIVEXMN=',E13.5,5X,2I10)
C
C**********************************************************************C
C
C **  UPDATE ZERO DIMENSION VOLUME BALANCE
C 
C----------------------------------------------------------------------C
C
      IF(ISDRY.GE.1.AND.ISTL.EQ.3)THEN
        VOLADD=0.
        DO L=2,LA
        IF(SPB(L).NE.0)THEN
          VOLADD=VOLADD+QSUME(L)-RIFTR(L)-EVAPSW(L)
        ENDIF
        ENDDO
        VOLADD=VOLADD*DT
        VOLZERD=VOLZERD+VOLADD
        VETZERD=VETZERD+VOLADD+DT*EVAPSW(L)
      ENDIF
 
 5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E13.5))
C                         
C**********************************************************************C
C
      RETURN
      END
