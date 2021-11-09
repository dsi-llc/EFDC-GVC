C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE VCONGRAD(ISTL)
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
C **  SUBROUTINE VCONGRAD SOLVES THE EXTERNAL MODE LINEAR SYSTEM 
C **  BY A PRECONDITIONED CONJUGATE GRADIENT SCHEME USING 
C **  ALTIVEC VECTOR PROCESSING
C
C **  THE LINEAR SYSTEM P IS SPECIFIED ON A 5 POINT STENCIL BY
C
C     CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)+CCN(L)*P(LN)
C     =FPTMP(L)
C
C     WHERE LN=LNC(L) AND LS=LSC(L) ARE LOOK UP TABLES FOR THE NORTH 
C     AND SOUTH STENCIL POINTS.
C
C     THE ACTUAL SIZE OF THE SYSTEM IS LC-2, WITH VARIABLES OCCUPING
C     POSITIONS (2:LC-1) IN THE ARRAYS WHICH MUST BE DIMENSIONED 
C     AT LCM = OR > THAN LC.  THE LOCATIONS 1 AND LC ARE DUMMIES SUCH
C     THAT L-1 AND L+1 ARE DEFINED AND ARE DEFINED AND CAN BE SET TO
C     ZERO.  VECTOR ALTIVEC OPERATIONS MUST THE CARRIED OUT ON THE 
C     FULL STORAGE DIMENSION, LCM WHICH IS AN INTEGER MULTIPE OF 4, OF 
C     THE ARRAYS FOR PROPER DIMENSIONING OF THE C VECTOR SUBROUTINE 
C     ARGUMENTS AS VECTOR FLOATS.  
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C ** HAMRICK'S C VECTOR SUBROUTINES MUST BE LINKED WITH EXECUTABLE
C ** FOR VECTOR EXECUTION ON MAC G4.  FOR PC'S, A DUMMY SUBROUTINE
C ** CVMULAD2 IS INCLUDED IN THE SOURCE CODE
C
C     CVMULAD2(VAL(N),A,B,C,D) => [A]=[B]+[C]*[D]
C
C     [A],[B],[C] AND [D] ARE 1D REAL*4 ARRAYS OF DIMENSION (1:LCM)
C
C**********************************************************************C
C
C **  START TIMER
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
C **  INITIALIZATION
C
      IVECDIG=0
      LCM4=LCM/4
C
      P(1)=0.0
      PW(1)=0.0
      PE(1)=0.0
      PS(1)=0.0
      PN(1)=0.0
      RCG(1)=0.0
      PCG(1)=0.0
      APCG(1)=0.0
      RA4(1)=0.0
      RTMP(1)=0.0
      RTMP1(1)=0.0
      RNULL(1)=0.0
      FPTMP(1)=0.0
      PCGW(1)=0.0
      PCGE(1)=0.0
      PCGS(1)=0.0
      PCGN(1)=0.0
      CCC(1)=0.0
      CCW(1)=0.0
      CCE(1)=0.0
      CCS(1)=0.0
      CCN(1)=0.0
      CCCI(1)=0.0
C
      DO L=LC,LCM
        P(L)=0.0
        PW(L)=0.0
        PE(L)=0.0
        PS(L)=0.0
        PN(L)=0.0
        RCG(L)=0.0
        PCG(L)=0.0
        APCG(L)=0.0
        RA4(L)=0.0
        RTMP(L)=0.0
        RTMP1(L)=0.0
        RNULL(L)=0.0
        FPTMP(L)=0.0
        PCGW(L)=0.0
        PCGE(L)=0.0
        PCGS(L)=0.0
        PCGN(L)=0.0
        CCC(L)=0.0
        CCW(L)=0.0
        CCE(L)=0.0
        CCS(L)=0.0
        CCN(L)=0.0
        CCCI(L)=0.0
      ENDDO
C
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        PW(L)=P(L-1)
        PE(L)=P(L+1)
        PS(L)=P(LS)
        PN(L)=P(LN)
      ENDDO
C
C**********************************************************************C
C
C     NOTE: LA=LC-1
C
C     DO L=2,LA
C       LN=LNC(L)
C       LS=LSC(L)
C       RCG(L)=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
C    &     +CCN(L)*P(LN)-FPTMP(L)
C     ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 0'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,CCW(L),CCC(L),CCE(L),CCS(L),CCCI(L),CCN(L),P(L)
        END DO
      ENDIF
C
      DO L=2,LA
        RCG(L)=-FPTMP(L)
        RNULL(L)=0.0
        RTMP(L)=0.0
        RTMP1(L)=0.0
      END DO
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 1'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RCG(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCC,P)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 2'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RTMP(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),RCG,RTMP,CCS,PS)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 3'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RCG(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCW,PW)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 4'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RTMP(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),RCG,RTMP,CCE,PE)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 5'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RCG(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RCG,CCN,PN)
C
C     DO L=2,LA
C       RCG(L)=-RCG(L)
C       PCG(L)=RCG(L)*CCCI(L)
C     ENDDO
C
      DO L=2,LA
        RCG(L)=-RTMP(L)
      END DO
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 6'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RCG(L)
        END DO
      ENDIF
C
cval      CALL CVMULAD2(VAL(LCM),PCG,RNULL,RCG,CCCI)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 7'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,PCG(L)
        END DO
      ENDIF
C
C     RPCG=0.
C     DO L=2,LA
C       RPCG=RPCG+RCG(L)*RCG(L)*CCCI(L)
C     ENDDO
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,RCG,CCCI)
cval      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RCG,RTMP)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 8'
C
      RPCG=0.
      DO L=2,LA
        RPCG=RPCG+RTMP1(L)
      ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RPCG
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
        LN=LNC(L)
        LS=LSC(L)
        PCGS(L)=PCG(LS)
        PCGW(L)=PCG(L-1)
        PCGE(L)=PCG(L+1)
        PCGN(L)=PCG(LN)
      ENDDO
C
C     DO L=2,LA
C       LN=LNC(L)
C       LS=LSC(L)
C       APCG(L)=CCC(L)*PCG(L)+CCS(L)*PCG(LS)+CCW(L)*PCG(L-1)
C    &         +CCE(L)*PCG(L+1)+CCN(L)*PCG(LN)
C     ENDDO
C
cval      CALL CVMULAD2(VAL(LCM),APCG,RNULL,CCC,PCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 9'
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,APCG,CCS,PCGS)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 10'
C
cval      CALL CVMULAD2(VAL(LCM),APCG,RTMP,CCW,PCGW)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 11'
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,APCG,CCE,PCGE)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 12'
C
cval      CALL CVMULAD2(VAL(LCM),APCG,RTMP,CCN,PCGN)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 13'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,APCG(L)
        END DO
      END IF
C
C     PAPCG=0.
C     DO L=2,LA
C       PAPCG=PAPCG+APCG(L)*PCG(L)
C     END DO
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,APCG,PCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 14'
C
      PAPCG=0.
      DO L=2,LA
        PAPCG=PAPCG+RTMP(L)
        RTMP1(L)=P(L)
      ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,800)LA,PAPCG
C
      ALPHA=RPCG/PAPCG
C
C     DO L=2,LA
C       P(L)=P(L)+ALPHA*PCG(L)
C     ENDDO
C
      DO L=2,LA
        RA4(L)=ALPHA
      END DO
C
cval      CALL CVMULAD2(VAL(LCM),P,RTMP1,RA4,PCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 15'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,P(L)
        END DO
      END IF
C
C**********************************************************************C
C
C     DO L=2,LA
C       RCG(L)=RCG(L)-ALPHA*APCG(L)
C     ENDDO
C
      DO L=2,LA
        RA4(L)=-ALPHA
        RTMP(L)=RCG(L)
      END DO
C
cval      CALL CVMULAD2(VAL(LCM),RCG,RTMP,RA4,APCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 16'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,RCG(L)
        END DO
      END IF
C
C     RPCGN=0.
C     DO L=2,LA
C       RPCGN=RPCGN+RCG(L)*RCG(L)*CCCI(L)
C     ENDDO
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,RCG,CCCI)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 17'
C
cval      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RTMP,RCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 18'
C
      RPCGN=0.
      DO L=2,LA
        RPCGN=RPCGN+RTMP1(L)
        RTMP(L)=RCG(L)
      ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RPCGN
C
C     RSQ=0.
C     DO L=2,LA
C       RSQ=RSQ+RCG(L)*RCG(L)
C     ENDDO
C
cval      CALL CVMULAD2(VAL(LCM),RTMP1,RNULL,RTMP,RCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 19'
C
      RSQ=0.
      DO L=2,LA
        RSQ=RSQ+RTMP1(L)
      ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RSQ
C
      IF(RSQ.LE.RSQM)GOTO 200
C
      IF(ITER.GE.ITERM)THEN
        WRITE(6,600)
        STOP
      ENDIF
C  
      BETA=RPCGN/RPCG
      RPCG=RPCGN
C
C     DO L=2,LA
C       PCG(L)=CCCI(L)*RCG(L)+BETA*PCG(L)
C     ENDDO
C
      DO L=2,LA
        RA4(L)=BETA
      ENDDO
C
cval      CALL CVMULAD2(VAL(LCM),RTMP,RNULL,PCG,RA4)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 20'
C
cval      CALL CVMULAD2(VAL(LCM),PCG,RTMP,CCCI,RCG)
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 21'
      IF(IVECDIG.EQ.1)THEN
        DO L=1,LC
          WRITE(8,800)L,PCG(L)
        ENDDO
      ENDIF
C
      GOTO 100
C
  600 FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
C**********************************************************************C
C
C ** CALCULATE FINAL RESIDUAL
C
  200 CONTINUE
C
      RSQ=0.
C
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     $        +CCN(L)*P(LN)-FPTMP(L)
        RSD=RSD*CCCI(L)
        RSQ=RSQ+RSD*RSD
      ENDDO
C
      IF(IVECDIG.EQ.1)WRITE(8,*)'VCONGRAD 22'
      IF(IVECDIG.EQ.1)WRITE(8,800)LA,RSQ
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
C
C**********************************************************************C
C
      RETURN
      END
