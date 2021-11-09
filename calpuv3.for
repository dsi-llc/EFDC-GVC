C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV3 (ISTL)
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
C ** SUBROUTINE CALPUV3 CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOW WITHOUT WETTING AND DRYING
C ** USING THE SMOLARKIEWICZ AND MARGOLIN SCHEME
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(ISDSOLV.EQ.1)THEN
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
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=FXE(L)-SUB(L)*DELTD2*DXIU(L)*FPGXE(L)
      FVHDXE(L)=FYE(L)-SVB(L)*DELTD2*DYIV(L)*FPGYE(L)
      ENDDO
C
      ENDIF
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
      FPGXE(L)= -SBX(L)*HU(L)*((BI2(L)+BI2(L-1))*(HP(L)-HP(L-1))
     &           +2.*HU(L)*(BI1(L)-BI1(L-1))
     &           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1))) 
      FPGYE(L)= -SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
     &           +2.*HV(L)*(BI1(L)-BI1(LS))
     &           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS)))
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=FUHDYE(L)+SUB(L)*DELT*DXIU(L)*FPGXE(L)
      FVHDXE(L)=FVHDXE(L)+SVB(L)*DELT*DYIV(L)*FPGYE(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=FXE(L)+SUB(L)*DELTD2*DXIU(L)*FPGXE(L)
      FVHDXE(L)=FYE(L)+SVB(L)*DELTD2*DYIV(L)*FPGYE(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  RESET RELAXATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO L=2,LA
      LN=LNC(L)
      C1=-G*DT*DT*SPB(L)*DXYIP(L)
      CS(L)=C1*SVB(L)*HRVO(L)*(2.*HV(L)-H1V(L))
      CW(L)=C1*SUB(L)*HRUO(L)*(2.*HU(L)-H1U(L))
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*(2.*HU(L+1)-H1U(L+1))
      CN(L)=C1*SVB(LN)*HRVO(LN)*(2.*HV(LN)-H1V(LN))
      CC(L)=1.-CS(L)-CW(L)-CE(L)-CN(L)
      CCI(L)=1./CC(L)
      ENDDO
C
      NCORDRY=0
C
      IF(ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),FP(L)
        ENDDO
        CLOSE(1)
      ENDIF
C
      IF(ISDSOLV.EQ.1.AND.N.EQ.1)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
      ENDIF
      IF(ISDSOLV.EQ.1.AND.N.EQ.2)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
         ENDDO
        CLOSE(1)
      ENDIF
 1001 FORMAT(2I5,6(1X,E12.4))
 1002 FORMAT(3I4,8(1X,E9.2))
C
      DO LR=1,NRC
      L=LRC(LR)
      CSR(LR)=CS(L)*CCI(L)
      CWR(LR)=CW(L)*CCI(L)
      CER(LR)=CE(L)*CCI(L)
      CNR(LR)=CN(L)*CCI(L)
      ENDDO
C
      DO LB=1,NBC
      L=LBC(LB)
      CSB(LB)=CS(L)*CCI(L)
      CWB(LB)=CW(L)*CCI(L)
      CEB(LB)=CE(L)*CCI(L)
      CNB(LB)=CN(L)*CCI(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      C1=-G*DT*DT*SPB(L)*DXYIP(L)
      C1=0.25*C1
      CCS(L)=C1*SVB(L)*HRVO(L)*HV(L)
      CCW(L)=C1*SUB(L)*HRUO(L)*HU(L)
      CCE(L)=C1*SUB(L+1)*HRUO(L+1)*HU(L+1)
      CCN(L)=C1*SVB(LN)*HRVO(LN)*HV(LN)
      CCC(L)=1.-CCS(L)-CCW(L)-CCE(L)-CCN(L)
      CCCI(L)=1./CCC(L)
      ENDDO
C
C
      IF(ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &               FP(L)
        ENDDO
        CLOSE(1)
      ENDIF
C
      IF(ISDSOLV.EQ.1.AND.N.EQ.1)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
      ENDIF
      IF(ISDSOLV.EQ.1.AND.N.EQ.2)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
         ENDDO
        CLOSE(1)
      ENDIF
C
      DO LR=1,NRC
      L=LRC(LR)
      CCSR(LR)=CCS(L)*CCCI(L)
      CCWR(LR)=CCW(L)*CCCI(L)
      CCER(LR)=CCE(L)*CCCI(L)
      CCNR(LR)=CCN(L)*CCCI(L)
      ENDDO
C
      DO LB=1,NBC
      L=LBC(LB)
      CCSB(LB)=CCS(L)*CCCI(L)
      CCWB(LB)=CCW(L)*CCCI(L)
      CCEB(LB)=CCE(L)*CCCI(L)
      CCNB(LB)=CCN(L)*CCCI(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  SET RHS OF EXTERNAL CONTINUITY EQUATION
C **  HRU=HMU*DYU/DXU & HRV=HMV*DXV/DYV 
C **  DXYIP=1/(DXP*DYP)
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      LN=LNC(L)
      FP(L)=SPB(L)*(P1(L)
     &      -G*DELTD2*DXYIP(L)*(UHDY1E(L+1)-UHDY1E(L)+FUHDYE(L+1)
     &      -FUHDYE(L)+VHDX1E(LN)-VHDX1E(L)+FVHDXE(LN)-FVHDXE(L)))
     &      *CCCI(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      FP(L)=SPB(L)*(P1(L)
     &      -G*DELTD2*DXYIP(L)*(UHDY1E(L+1)-UHDY1E(L)+FUHDYE(L+1)
     &      -FUHDYE(L)+VHDX1E(LN)-VHDX1E(L)+FVHDXE(LN)-FVHDXE(L)))
     &      *CCI(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  SET TIDAL ELEVATION AND VOLUME SOURCE BOUNDARY CONDITIONS
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
      FP(L)=PSERT(NPSERW(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
      ENDDO
      IF(ISPBW(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
       CET=(1.-TMP)/(1.+TMP)
       CE(L)=CET
       CCE(L)=CET
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CER(LR)=CET
        CCER(LR)=CET
       ELSE
        LB=LLBC(L)
        CEB(LB)=CET
        CCEB(LB)=CET
       ENDIF
       FP(L)=(4.*FP(L)-2.*SQRT(G/HMU(L+1))*FUHDYE(L+1)*DYIU(L+1))/
     &       (1.+TMP)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)      
      FP(L)=PSERT(NPSERE(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
      ENDDO
      IF(ISPBE(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
       CWT=(1.-TMP)/(1.+TMP)
       CW(L)=CWT
       CCW(L)=CWT
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CWR(LR)=CWT
        CCWR(LR)=CWT
       ELSE
        LB=LLBC(L)
        CWB(LB)=CWT
        CCWB(LB)=CWT
       ENDIF
       FP(L)=(4.*FP(L)+2.*SQRT(G/HMU(L))*FUHDYE(L)*DYIU(L))/(1.+TMP)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      FP(L)=PSERT(NPSERS(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
      ENDDO
      IF(ISPBS(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
       CNT=(1.-TMP)/(1.+TMP)
       CN(L)=CNT
       CCN(L)=CNT
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CNR(LR)=CNT
        CCNR(LR)=CNT
       ELSE
        LB=LLBC(L)
        CNB(LB)=CNT
        CCNB(LB)=CNT
       ENDIF
       FP(L)=(4.*FP(L)-2.*SQRT(G/HMV(LN))*FVHDXE(LN)*DXIV(LN))/(1.+TMP)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
      LS=LSC(L)
      FP(L)=PSERT(NPSERN(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
      ENDDO
      IF(ISPBN(LL).GE.1)THEN
       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
       CST=(1.-TMP)/(1.+TMP)
       CS(L)=CST
       CCS(L)=CST
       IF(ISRED(L).EQ.1)THEN
        LR=LLRC(L)
        CSR(LR)=CST
        CCSR(LR)=CST
       ELSE
        LB=LLBC(L)
        CSB(LB)=CST
        CCSB(LB)=CST
       ENDIF
       FP(L)=(4.*FP(L)+2.*SQRT(G/HMV(L))*FVHDXE(L)*DXIV(L))/(1.+TMP)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(N.GT.NTSNB)THEN
      IF(ISTL.EQ.2)THEN
        DO L=2,LA
        FP(L)=FP(L)+DELT*G*SPB(L)*QSUME(L)*DXYIP(L)*CCCI(L)
        ENDDO
       ELSE
        DO L=2,LA
        FP(L)=FP(L)+DELT*G*SPB(L)*QSUME(L)*DXYIP(L)*CCI(L)
        ENDDO
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO LR=1,NRC
      L=LRC(LR)
      FPR(LR)=FP(L)
      ENDDO
      DO LB=1,NBC
      L=LBC(LB)
      FPB(LB)=FP(L)
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
        DELP=P(L)-P1(L)
        P1(L)=P(L)
        P(L)=P(L)+DELP
        H1U(L)=HU(L)
        H1V(L)=HV(L)
        H1UI(L)=HUI(L)
        H1VI(L)=HVI(L)
        H2P(L)=H1P(L)
        H1P(L)=HP(L)
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  SOLVE EQUATIONS FOR P 
C
C----------------------------------------------------------------------C
C
      DO LR=1,NRC
      L=LRC(LR)
      PRED(LR)=P(L)
      ENDDO
C
      DO LB=1,NBC
      L=LBC(LB)
      PBLK(LB)=P(L)
      ENDDO
C
      IF(ISTL.EQ.3)THEN
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      FPR(L)=FPR(L)-PRED(L)-CSR(L)*PBLK(LS)-CWR(L)*PBLK(LW)
     &                     -CER(L)*PBLK(LE)-CNR(L)*PBLK(LN)
      ENDDO
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      FPB(L)=FPB(L)-PBLK(L)-CSB(L)*PRED(LS)-CWB(L)*PRED(LW)
     &                     -CEB(L)*PRED(LE)-CNB(L)*PRED(LN)
      ENDDO
C
      ELSE
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      FPR(L)=FPR(L)-PRED(L)-CCSR(L)*PBLK(LS)-CCWR(L)*PBLK(LW)
     &                     -CCER(L)*PBLK(LE)-CCNR(L)*PBLK(LN)
      ENDDO
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      FPB(L)=FPB(L)-PBLK(L)-CCSB(L)*PRED(LS)-CCWB(L)*PRED(LW)
     &                     -CCEB(L)*PRED(LE)-CCNB(L)*PRED(LN)
      ENDDO
C
      ENDIF 
C
      DO L=2,LA
      PAM(L)=P(L)
      P(L)=0.
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(IRVEC.EQ.0) CALL RELAX (ISTL)
      IF(IRVEC.EQ.1) CALL RELAXV (ISTL)
      IF(IRVEC.EQ.2) CALL CONGRAD (ISTL)
      IF(IRVEC.EQ.3) CALL CGRS (ISTL)
      IREPRX=0
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=P(L)+PAM(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE UHEX AND VHEX AND TOTAL DEPTHS AT TIME LEVEL (N+1)
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)      
      UTMP=SUB(L)*( FUHDYE(L)-0.5*DELTD2*HRUO(L)*(P(L)-P(L-1))*
     &     (GI*(P(L)+P(L-1))-BELV(L)-BELV(L-1)) )
      VTMP=SVB(L)*( FVHDXE(L)-0.5*DELTD2*HRVO(L)*(P(L)-P(LS))*
     &     (GI*(P(L)+P(LS))-BELV(L)-BELV(LS)) )
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     &       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        ENDDO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      ENDDO
C
      DO L=2,LA
      LS=LSC(L)      
      HU(L)=0.5*GI*(P(L)+P(L-1))-0.5*(BELV(L)+BELV(L-1))
      HV(L)=0.5*GI*(P(L)+P(LS))-0.5*(BELV(L)+BELV(LS))
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
C **  CHECK FOR NEGATIVE DEPTHS
C
C----------------------------------------------------------------------C
C
      IF(ISNEGH.EQ.1)THEN
C
      DO L=2,LA
      IF(HP(L).LT.0.)THEN
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
 6060 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HP,H1P,H2P =',3(2X,E12.4))
 6061 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUW,H1UW =',2(2X,E12.4))
 6062 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUE,H1UE =',2(2X,E12.4))
 6063 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVS,H1VS =',2(2X,E12.4))
 6064 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVN,H1VN =',2(2X,E12.4))
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
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI
     &     +0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     &     +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))-QSUME(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      ENDIF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      ENDIF
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI
     &     +0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      ENDIF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
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
 6628 FORMAT('  DIVEXMX=',E12.4,5X,2I10)
 6629 FORMAT('  DIVEXMN=',E12.4,5X,2I10)
C
C**********************************************************************C
C
      RETURN
      END
