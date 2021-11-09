C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUVB (ISTL)
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
C ** SUBROUTINE CALPUVB CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOWS USING AN IMPLICIT SCHEME WITH
C ** TWO-TIME LEVEL CONTINUITY FOR BOTH THE 3TL AND 2TL STEPS
C ** IRVEC=11
C
C    FIRST PASS 3TL HAS FORM
C
C    A**2-1      A        A**2+2A+1     U
C      A       A**2-1     A**2+2A+1     V  
C     A+1       A+1          A-1        P
C
C    SECOND PASS 3TL HAS FORM
C
C    A**2-1     A**2+2A+1   A**2+2A+1     U
C   A**2+2A+1    A**2-1     A**2+2A+1     V  
C     A+1         A+1         A-1         P
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION QSUMTMP(LCM)
C
C**********************************************************************C
C
C      WRITE(6,6000)N
C 6000 FORMAT(' CALLED CALPUV9, N = ',I10)
C
      IF(N.EQ.1.AND.ISDSOLV.GT.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF1.OUT',STATUS='UNKNOWN')
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
C **  CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
C
      CALL CALEBI
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL PRESSURE GRADIENTS 
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C **  SNLPX=SNLPX*GID2*DYU & SNLPY=SNLPY*GID2*DXV
C
C----------------------------------------------------------------------C
C
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
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDY1E(L)
     &         -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P1(L)-P1(L-1))
     &         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-TBX1(L))
     &         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      FVHDXE(L)=VHDX1E(L)
     &         -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P1(L)-P1(LS ))
     &         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-TBY1(L))
     &         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDY1E(L)
     &         -0.5*DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P1(L)-P1(L-1))
     &         -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P(L)-P(L-1))
     &         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-TBX1(L))
     &         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      FVHDXE(L)=VHDX1E(L)
     &         -0.5*DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P1(L)-P1(LS ))
     &         -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P(L)-P(LS ))
     &         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-TBY1(L))
     &         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ENDDO
C
      ENDIF
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
      IF(ISTL.EQ.2)THEN
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
      CET=-0.5*DELTD2*G*HRUO(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/(1.-TMP)
       CE(L)=CET
       FP1(L)=CET*(4.*FP1(L)
     & -2.*SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/(1.-TMP)
      ELSE
       FP1(L+1)=-CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ELSE
C
      DO LL=1,NPBW
      L=LPBW(LL)
       CC(L)=2.*DELTI*DXYP(L)
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
      CET=-0.25*DELTD2*G*HRUO(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/(1.-TMP)
       CE(L)=CET
       FP1(L)=CET*(4.*FP1(L)
     & -2.*SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/(1.-TMP)
      ELSE
       FP1(L+1)=-CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
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
      CWT=-0.5*DELTD2*G*HRUO(L  )*HUTMP(L  )
      IF(ISPBE(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/(1.-TMP)
       CW(L)=CWT
       FP1(L)=CWT*(4.*FP1(L)
     & +2.*SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/(1.-TMP)
      ELSE
       FP1(L-1)=-CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ELSE
C
      DO LL=1,NPBE
      L=LPBE(LL)
       CC(L)=2.*DELTI*DXYP(L)
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
      CWT=-0.25*DELTD2*G*HRUO(L  )*HUTMP(L  )
      IF(ISPBE(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/(1.-TMP)
       CW(L)=CWT
       FP1(L)=CWT*(4.*FP1(L)
     & +2.*SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/(1.-TMP)
      ELSE
       FP1(L-1)=-CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
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
      CNT=-0.5*DELTD2*G*HRVO(LN )*HVTMP(LN )
      IF(ISPBS(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
       CC(L)=CNT*(1.+TMP)/(1.-TMP)
       CN(L)=CNT
       FP1(L)=CNT*(4.*FP1(L)
     & -2.*SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/(1.-TMP)
      ELSE
       FP1(LN)=-CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ELSE
C
      DO LL=1,NPBS
      L=LPBS(LL)
       CC(L)=2.*DELTI*DXYP(L)
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
      CNT=-0.25*DELTD2*G*HRVO(LN )*HVTMP(LN )
      IF(ISPBS(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
       CC(L)=CNT*(1.+TMP)/(1.-TMP)
       CN(L)=CNT
       FP1(L)=CNT*(4.*FP1(L)
     & -2.*SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/(1.-TMP)
      ELSE
       FP1(LN)=-CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
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
      CST=-0.5*DELTD2*G*HRVO(L  )*HVTMP(L  )
      IF(ISPBN(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
       CC(L)=CST*(1.+TMP)/(1.-TMP)
       CS(L)=CST
       FP1(L)=CST*(4.*FP1(L)
     & +2.*SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/(1.-TMP)
      ELSE
       FP1(LS)=-CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ELSE
C
      DO LL=1,NPBN
      L=LPBN(LL)
       CC(L)=2.*DELTI*DXYP(L)
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
      CST=-0.25*DELTD2*G*HRVO(L  )*HVTMP(L  )
      IF(ISPBN(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
       CC(L)=CST*(1.+TMP)/(1.-TMP)
       CS(L)=CST
       FP1(L)=CST*(4.*FP1(L)
     & +2.*SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/(1.-TMP)
      ELSE
       FP1(LS)=-CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
      ENDIF
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
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      LN=LNC(L)
      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*P1(L)
     &      -0.5*G*(UHDY1E(L+1)-UHDY1E(L)
     &             +VHDX1E(LN )-VHDX1E(L)) )
      ENDDO
C 
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      FP1(L)=FP1(L)+SPB(L)*( 2.*DELTI*DXYP(L)*P(L)
     &      -0.5*G*(UHDYE(L+1)-UHDYE(L)
     &             +VHDXE(LN )-VHDXE(L)) )
      ENDDO
C
      ENDIF
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
C----------------------------------------------------------------------C
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
      CCMNM=1.E+18
C
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      LN=LNC(L)
      C1=-0.5*DELTD2*G*SPB(L)
      CS(L)=C1*SVB(L  )*HRVO(L  )*HVTMP(L  )
C    &     +(1.-SPB(L))*CS(L)
      CW(L)=C1*SUB(L  )*HRUO(L  )*HUTMP(L  )
C    &     +(1.-SPB(L))*CW(L)
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*HUTMP(L+1)
C    &     +(1.-SPB(L))*CE(L)
      CN(L)=C1*SVB(LN )*HRVO(LN )*HVTMP(LN )
C    &     +(1.-SPB(L))*CN(L)
      CC(L)=SPB(L)*(DELTI*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
C    &     +(1.-SPB(L))*CC(L)
      ENDIF
      CCMNM=MIN(CCMNM,CC(L))
      FPTMP(L)=FP(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      LN=LNC(L)
      C1=-0.25*DELTD2*G*SPB(L)
      CS(L)=C1*SVB(L  )*HRVO(L  )*HVTMP(L  )
C    &     +(1.-SPB(L))*CS(L)
      CW(L)=C1*SUB(L  )*HRUO(L  )*HUTMP(L  )
C    &     +(1.-SPB(L))*CW(L)
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*HUTMP(L+1)
C    &     +(1.-SPB(L))*CE(L)
      CN(L)=C1*SVB(LN )*HRVO(LN )*HVTMP(LN )
C    &     +(1.-SPB(L))*CN(L)
      CC(L)=SPB(L)*(2.*DELTI*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
C    &     +(1.-SPB(L))*CC(L)
      ENDIF
      CCMNM=MIN(CCMNM,CC(L))
      FPTMP(L)=FP(L)
      ENDDO
C
      ENDIF
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
C **  SCALE BY MINIMUM DIAGONAL
C
C      DO L=2,LA
C      CCS(L)=CS(L)*CCMNMI
C      CCW(L)=CW(L)*CCMNMI
C      CCE(L)=CE(L)*CCMNMI
C      CCN(L)=CN(L)*CCMNMI
C      CCC(L)=CC(L)*CCMNMI
C      FPTMP(L)=FPTMP(L)*CCMNMI
C      CCCI(L)=1./CCC(L)
C      ENDDO
C
C
C **  SCALE TO NORMAL FORM
C
      DO L=2,LA
      CCS(L)=CS(L)/SQRT( CC(L)*CC(LSC(L)) )
      CCW(L)=CW(L)/SQRT( CC(L)*CC(L-1   ) )
      CCE(L)=CE(L)/SQRT( CC(L)*CC(L+1   ) )
      CCN(L)=CN(L)/SQRT( CC(L)*CC(LNC(L)) )
      CCC(L)=1.
      FPTMP(L)=FPTMP(L)/SQRT( CC(L) )
      CCCI(L)=1.
      ENDDO
C
      CALL CONGRAD (ISTL)
C
      DO L=2,LA
      P(L)=P(L)/SQRT( CC(L) )
      ENDDO
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &               FPTMP(L),SURFTMP
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='EQCOEF1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
         SURFTMP=GI*P(L)
         WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &               FPTMP(L),SURFTMP
          ENDDO
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
      ENDIF
C
 1001 FORMAT(2I5,7(1X,E12.4))
 1002 FORMAT(3I4,8(1X,E9.2))
C
C**********************************************************************C
C
C **  CALCULATE UHEX AND VHEX AND TOTAL DEPTHS AT TIME LEVEL (N+1)
C **  HRU=SUB*DYU/DXU & HRV=SVB*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      LS=LSC(L)      
      UTMP=SUB(L)*( FUHDYE(L)
     &            -DELTD2*HRUO(L)*HUTMP(L)*(P(L)-P(L-1)) )
      VTMP=SVB(L)*( FVHDXE(L)
     &            -DELTD2*HRVO(L)*HVTMP(L)*(P(L)-P(LS )) )
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LS=LSC(L)      
      UTMP=SUB(L)*( FUHDYE(L)
     &            -0.5*DELTD2*HRUO(L)*HUTMP(L)*(P(L)-P(L-1)) )
      VTMP=SVB(L)*( FVHDXE(L)
     &            -0.5*DELTD2*HRVO(L)*HVTMP(L)*(P(L)-P(LS )) )
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
C** ADD FILTER
      UHDY1E(L)=UHDY1E(L)+FILT3TL*(UHDYE(L)-2.*UHDY1E(L)+UHDY2E(L))
      VHDX1E(L)=VHDX1E(L)+FILT3TL*(VHDXE(L)-2.*VHDX1E(L)+VHDX2E(L))
C** ADD FILTER
      ENDDO
C
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
        HPPTMP=H1P(L)+0.5*DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
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
      HU(L)=0.5*(HP(L)+HP(L-1))
      HV(L)=0.5*(HP(L)+HP(LS))
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
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
      ENDIF
      ENDDO
C
      DO L=2,LA
      IF(HU(L).LT.0.)THEN
      INEGFLG=1
      LN=LNC(L)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
      ENDIF
      ENDDO
C
C
      DO L=2,LA
      IF(HV(L).LT.0.)THEN
      INEGFLG=1
      LN=LNC(L)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
      ENDIF
      ENDDO
C
      ENDIF
C
      IF(ISNEGH.EQ.2)THEN
      IF(INEGFLG.EQ.1)THEN
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
        STOP
C
      ENDIF
      ENDIF
C
 6060 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HP,H1P,H2P =',3(2X,E12.4))
 6061 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUW,H1UW =',2(2X,E12.4))
 6062 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUE,H1UE =',2(2X,E12.4))
 6063 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVS,H1VS =',2(2X,E12.4))
 6064 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVN,H1VN =',2(2X,E12.4))
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
      DIVEX=2.*SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI
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
