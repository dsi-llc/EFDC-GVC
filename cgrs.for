C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CGRS (ISTL)
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
C **  SUBROUTINE CGRS SOLVES THE EXTERNAL MODE BY A RED-BLACK
C **  ORDERED REDUCED SYSTEM CONJUGATEC GRADIENT SCHEME
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C     DIMENSION CGSR(LCM),CGWR(LCM),CGER(LCM),CGNR(LCM),
C    &          CGSB(LCM),CGWB(LCM),CGEB(LCM),CGNB(LCM),
C    &          PCG(LCM),RCG(LCM),APCG(LCM),PBTMP(LCM),FPRT(LCM)
CDHP  DIMENSION PSOUT(LCM),PWEST(LCM),PEAST(LCM),PNORT(LCM),RSDT(LCM)
C
C**********************************************************************C
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C**********************************************************************C
C
      IF(ISTL.EQ.2)THEN
       DO L=1,LC
       CGSR(L)=CCSR(L)
       CGWR(L)=CCWR(L)
       CGER(L)=CCER(L)
       CGNR(L)=CCNR(L)
       CGSB(L)=CCSB(L)
       CGWB(L)=CCWB(L)
       CGEB(L)=CCEB(L)
       CGNB(L)=CCNB(L)
       ENDDO
      ELSE
       DO L=1,LC
       CGSR(L)=CSR(L)
       CGWR(L)=CWR(L)
       CGER(L)=CER(L)
       CGNR(L)=CNR(L)
       CGSB(L)=CSB(L)
       CGWB(L)=CWB(L)
       CGEB(L)=CEB(L)
       CGNB(L)=CNB(L)
       ENDDO
      ENDIF
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      FPRT(L)=FPR(L)-CGSR(L)*FPB(LS)-CGWR(L)*FPB(LW)
     &              -CGER(L)*FPB(LE)-CGNR(L)*FPB(LN)
      ENDDO
C
C**********************************************************************C
C
C **  LOAD P INTO PRED AND PBLK
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
C**********************************************************************C
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBTMP(L)=CGSB(L)*PRED(LS)+CGWB(L)*PRED(LW)
     &        +CGEB(L)*PRED(LE)+CGNB(L)*PRED(LN)
      ENDDO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)-CGSR(L)*PBTMP(LS)-CGWR(L)*PBTMP(LW)
     &           -CGER(L)*PBTMP(LE)-CGNR(L)*PBTMP(LN)-FPRT(L)
      RCG(L)=-RSD
      PCG(L)=-RSD
      ENDDO
C
      ITER=0
C
C**********************************************************************C
C
  100 CONTINUE
C
      ITER=ITER+1
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBTMP(L)=CGSB(L)*PCG(LS)+CGWB(L)*PCG(LW)
     &        +CGEB(L)*PCG(LE)+CGNB(L)*PCG(LN)
      ENDDO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      APCG(L)=PCG(L)-CGSR(L)*PBTMP(LS)-CGWR(L)*PBTMP(LW)
     &              -CGER(L)*PBTMP(LE)-CGNR(L)*PBTMP(LN)
      ENDDO
C
      PAPCG=0.
      RPCG=0.
C
      DO L=1,NRC
      PAPCG=PAPCG+APCG(L)*PCG(L)
      RPCG=RPCG+RCG(L)*PCG(L)
      ENDDO
C
      ALPHA=RPCG/PAPCG
C
      DO L=1,NRC
      PRED(L)=PRED(L)+ALPHA*PCG(L)
      ENDDO
C
C**********************************************************************C
C
      RSQ=0.
C  
      DO L=1,NRC
      RCG(L)=RCG(L)-ALPHA*APCG(L)
      ENDDO
C
      DO L=1,NRC
      RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
C
      IF(RSQ .LE. RSQM) GOTO 200
C
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
       STOP
      ENDIF
C  
      BETA=0.
C
      DO L=1,NRC
      BETA=BETA+RCG(L)*APCG(L)
      ENDDO
C
      BETA=-BETA/PAPCG
C
      DO L=1,NRC
      PCG(L)=RCG(L)+BETA*PCG(L)
      ENDDO
C
      GOTO 100
C
C**********************************************************************C
C
C **  CALCULATE PBLK AND FINAL RESIDUAL
C
  200 CONTINUE
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBLK(L)=FPB(L)-CGSB(L)*PRED(LS)-CGWB(L)*PRED(LW)
     &              -CGEB(L)*PRED(LE)-CGNB(L)*PRED(LN)
      ENDDO
C
      RSQ=0.
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CGSB(L)*PRED(LS)+CGWB(L)*PRED(LW)
     &           +CGEB(L)*PRED(LE)+CGNB(L)*PRED(LN)-FPB(L)
      RSQ=RSQ+RSD*RSD
      ENDDO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CGSR(L)*PBLK(LS)+CGWR(L)*PBLK(LW)
     &           +CGER(L)*PBLK(LE)+CGNR(L)*PBLK(LN)-FPR(L)
      RSQ=RSQ+RSD*RSD
      ENDDO
C
C
C**********************************************************************C
C
C **  LOAD PRED AND PBLK INTO P
C
      DO LR=1,NRC
      L=LRC(LR)
      P(L)=PRED(LR)
      ENDDO
C
      DO LB=1,NBC
      L=LBC(LB)
      P(L)=PBLK(LB)
      ENDDO
C
      IF(ISCRAY.EQ.0)THEN
        TCGRS=TCGRS+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCGRS=TCGRS+T2TMP-T1TMP
        WTCGRS=WTCGRS+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
  600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
      RETURN
      END
