C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV (ISTL)
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
C ** SUBROUTINE CALPUV CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE FOR FLOWS WITH A RIGID LID CONDITION AT THE SURFACE
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
     &           -(BE(L)+BE(L-1))*(HMP(L)-HMP(L-1)))
     &           -SNLPX(L)*(P(L)+P(L-1))*(P(L)-P(L-1)) )
      FPGYE(L)=ROLD*FPGYE(L)+RNEW*(
     &           -SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
     &           +2.*HV(L)*(BI1(L)-BI1(LS))
     &           -(BE(L)+BE(LS))*(HMP(L)-HMP(LS)))
     &           -SNLPY(L)*(P(L)+P(LS))*(P(L)-P(LS)) )
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDY1E(L)-DELTD2*HRU(L)*(P1(L)-P1(L-1))
     &         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-TBX1(L))
     &         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      FVHDXE(L)=VHDX1E(L)-DELTD2*HRV(L)*(P1(L)-P1(LS))
     &         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-TBY1(L))
     &         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ENDDO
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
C **  SET PRESSURE AND VOLUMN SOURCE BOUNDARY CONDITIONS
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
        CE(L)=CET
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
C
      IF(ISTL.EQ.2)THEN
C
      DO L=2,LA
      FP(L)=FP(L)+DELT*G*SPB(L)*QSUME(L)*DXYIP(L)*CCCI(L)
      ENDDO
C
      ELSE
C
      DO L=2,LA
      FP(L)=FP(L)+DELT*G*SPB(L)*QSUME(L)*DXYIP(L)*CCI(L)
      ENDDO
C
      ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO LR=1,NRC
      L=LRC(LR)
      FPR(LR)=FP(L)
      ENDDO
C
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
C
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
C
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
      UTMP=FUHDYE(L)-DELTD2*HRU(L)*(P(L)-P(L-1))
      VTMP=FVHDXE(L)-DELTD2*HRV(L)*(P(L)-P(LS))
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO L=2,LA
      LN=LNC(L)
      HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
     &     -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     &     +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
      HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
      ENDDO
C
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     &     -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
      HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
      ENDDO
C
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
      WRITE (6,6066)L,IL(L),JL(L),HP(L)
      ENDIF
      ENDDO
C
      ENDIF
C
 6066 FORMAT('  NEG DEPTH AT L,I,J,HP=',(3I6),2X,E12.4)
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
 6628 FORMAT('  DIVEXMX=',E13.5,5X,2I10)
 6629 FORMAT('  DIVEXMN=',E13.5,5X,2I10)
C
C**********************************************************************C
C
      RETURN
      END
