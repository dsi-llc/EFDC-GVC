C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CONGRAD (ISTL)
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
C**********************************************************************C
C
      IF(LA.EQ.2)THEN
	  P(L)=FPTMP(L)/CCC(L)
	  RETURN
	ENDIF
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
      DO L=2,LA
	  PNORTH(L)=P(LNC(L))
	  PSOUTH(L)=P(LSC(L))
	ENDDO
C
      DO L=2,LA
        RCG(L)=FPTMP(L)-CCC(L)*P(L)-CCN(L)*PNORTH(L)-CCS(L)*PSOUTH(L)
     &        -CCW(L)*P(L-1)-CCE(L)*P(L+1)
      ENDDO
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
      PAPCG=0.
C
      DO L=2,LA
        PAPCG=PAPCG+APCG(L)*PCG(L)
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
        DO L=1,LC
	    CDIADOM=CCC(L)+CCE(L)+CCN(L)+CCS(L)+CCW(L)
          WRITE(8,808)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &                CDIADOM,FPTMP(L),HU(L),HV(L)
        END DO
  	  CLOSE(8)
	  CALL RESTOUT(1)
        STOP
      ENDIF
C  
      BETA=RPCGN/RPCG
      RPCG=RPCGN
C
      DO L=2,LA
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
C
      DO L=2,LA
      RCG(L)=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)
     &        +CCW(L)*P(L-1)+CCE(L)*P(L+1)-FPTMP(L)
      ENDDO
C
      DO L=2,LA
        RCG(L)=RCG(L)*CCCI(L)
      ENDDO
      DO L=2,LA
        RSQ=RSQ+RCG(L)*RCG(L)
      ENDDO
C
C**********************************************************************C
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
  800 FORMAT(I5,8E13.4)
  808 FORMAT(2I5,9E13.4)
C
      RETURN
      END
