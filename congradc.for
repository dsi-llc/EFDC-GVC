C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CONGRADC (ISTL)
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
C **  SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE 
C **  GRADIENT SCHEME    
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION PNORTH(LCM),PSOUTH(LCM),TMPCG(LCM)
C
C     DIMENSION CGS(LCM),CGW(LCM),CGE(LCM),CGN(LCM),
C    &          PCG(LCM),RCG(LCM),APCG(LCM)
C
C**********************************************************************C
C
C     IF(ISTL.EQ.2)THEN
C      DO L=2,LA
C      CIT=1./CCC(L)
C      CGS(L)=CIT*CCS(L)
C      CGW(L)=CIT*CCW(L)
C      CGE(L)=CIT*CCE(L)
C      CGN(L)=CIT*CCN(L)
C      ENDDO
C     ELSE
C      DO L=2,LA
C      CIT=1./CC(L)
C      CGS(L)=CIT*CS(L)
C      CGW(L)=CIT*CW(L)
C      CGE(L)=CIT*CE(L)
C      CGN(L)=CIT*CN(L)
C      ENDDO
C     ENDIF
C
C
C     IF(ISTL.EQ.2)THEN
C      DO L=2,LA
C      CG(L)=CCC(L)
C      CGS(L)=CCS(L)
C      CGW(L)=CCW(L)
C      CGE(L)=CCE(L)
C      CGN(L)=CCN(L)
C      ENDDO
C     ELSE
C      DO L=2,LA
C      CG(L)=CC(L)
C      CGS(L)=CS(L)
C      CGW(L)=CW(L)
C      CGE(L)=CE(L)
C      CGN(L)=CN(L)
C      ENDDO
C     ENDIF
C
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)
C     FPTMP(L)=FP(L)-CG(L)*PAM(L)-CGS(L)*PAM(LS)-CGW(L)*PAM(L-1)
C    &        -CGE(L)*PAM(L+1)-CGN(L)*PAM(LN)
C     ENDDO
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
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
C
      DO L=2,LA
C       LN=LNC(L)
C       LS=LSC(L)
C       RSD=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
C    &        +CCN(L)*PNORTH(L)-FPTMP(L)
C       RCG(L)=-RSD
        RCG(L)=FPTMP(L)-CCC(L)*P(L)-CCN(L)*PNORTH(L)-CCS(L)*PSOUTH(L)
     &        -CCW(L)*P(L-1)-CCE(L)*P(L+1)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            RCG(LCHNU)=RCG(LCHNU)+CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)+CCCCHH(NMD)*P(LCHNU)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            RCG(LCHNV)=RCG(LCHNV)+CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)+CCCCHH(NMD)*P(LCHNV)
          ENDIF
        ENDDO
      ENDIF
C
      DO L=2,LA
        PCG(L)=RCG(L)*CCCI(L)
      ENDDO
C
      RPCG=0.
      DO L=2,LA
      RPCG=RPCG+RCG(L)*PCG(L)
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
      DO L=2,LA
	  PNORTH(L)=PCG(LNC(L))
	  PSOUTH(L)=PCG(LSC(L))
	ENDDO
C
      DO L=2,LA
        APCG(L)=CCC(L)*PCG(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)
     &       +CCW(L)*PCG(L-1)+CCE(L)*PCG(L+1)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            APCG(LCHNU)=APCG(LCHNU)+CCCCHH(NMD)*PCG(LHOST)
            APCG(LHOST)=APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNU)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            APCG(LCHNV)=APCG(LCHNV)+CCCCHH(NMD)*PCG(LHOST)
            APCG(LHOST)=APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNV)
          ENDIF
        ENDDO
      ENDIF
C
      PAPCG=0.
COLD     RPCG=0.
C
      DO L=2,LA
      PAPCG=PAPCG+APCG(L)*PCG(L)
COLD     RPCG=RPCG+RCG(L)*PCG(L)
      ENDDO
C
      ALPHA=RPCG/PAPCG
C
      DO L=2,LA
      P(L)=P(L)+ALPHA*PCG(L)
      ENDDO
C
C**********************************************************************C
C
      DO L=2,LA
        RCG(L)=RCG(L)-ALPHA*APCG(L)
      ENDDO
c
      DO L=2,LA
        TMPCG(L)=CCCI(L)*RCG(L)
      ENDDO
C   
      RPCGN=0.
      RSQ=0.
      DO L=2,LA
        RPCGN=RPCGN+RCG(L)*TMPCG(L)
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
C
      IF(RSQ .LE. RSQM) GOTO 200
C
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
	 DO L=2,LA
	   WRITE(8,800)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &                 FPTMP(L)
       ENDDO
       STOP
      ENDIF
C  
COLD     BETA=0.
C
COLD     DO L=2,LA
COLD     BETA=BETA+RCG(L)*APCG(L)
COLD     ENDDO
C
COLD     BETA=-BETA/PAPCG
      BETA=RPCGN/RPCG
      RPCG=RPCGN
C
      DO L=2,LA
COLD     PCG(L)=RCG(L)+BETA*PCG(L)
      PCG(L)=TMPCG(L)+BETA*PCG(L)
      ENDDO
C
      GOTO 100
C
  600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
C**********************************************************************C
C
C ** CALCULATE FINAL RESIDUAL
C
  200 CONTINUE
C
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
C
      RSQ=0.
C
      DO L=2,LA
      RCG(L)=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)
     &        +CCW(L)*P(L-1)+CCE(L)*P(L+1)-FPTMP(L)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            RCG(LCHNU)=RCG(LCHNU)-CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)-CCCCHH(NMD)*P(LCHNU)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            RCG(LCHNV)=RCG(LCHNV)-CCCCHH(NMD)*P(LHOST)
            RCG(LHOST)=RCG(LHOST)-CCCCHH(NMD)*P(LCHNV)
          ENDIF
        ENDDO
      ENDIF
C
      DO L=2,LA
        RCG(L)=RCG(L)*CCCI(L)
      ENDDO
      DO L=2,LA
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
C
      IF(ISCRAY.EQ.0)THEN
        TCONG=TCONG+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCONG=TCONG+T2TMP-T1TMP
        WTCONG=WTCONG+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
  800 FORMAT(2I6,6E13.4)
C
      RETURN
      END
