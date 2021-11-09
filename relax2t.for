C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RELAX2T
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 02/15/2002        John Hamrick       02/15/2002       John Hamrick
C  added this subroutine relax2t
C----------------------------------------------------------------------C
C
C **  SUBROUTINE RELAX SOLVES THE FINITE DIFFERENCE FORM
C **  OF A PSEUDO HEMHOLTZ EQUATION
C **
C **              CS(L)*P(LS)+CW(L)*P(L-1)                
C **              +CC(L)*P(L)+CE(L)*P(L+1)                
C **              +CN(L)*P(LN) = FP(L)                    
C **                                                      
C **  BY SUCCESSIVE OVER RELAXATION USING A RED-BLACK ORDERING 
C **  WITH CONVERGENCE MEASURED BY A GLOBAL SQUARE ERROR RSQ.  
C **  NON-CONVERGENCE IS SIGNALED WHEN THE ITERATIONS EXCEED A 
C **  MAXIMUM.                                                 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      RJ2=RP
C
C      PAVG=0.0
C      DO L=2,LA
C      PAVG=PAVG+P(L)
C      ENDDO
C	PAVG=PAVG/FLOAT(LA-1)
C
C      DO L=2,LA
C        FPTMP(L)=FPTMP(L)-PAVG*(CCC(L)+CCS(L)+CCW(L)+CCE(L)+CCN(L))
C        P(L)=P(L)-PAVG
C      ENDDO
C
      FPSQ=0.
      DO L=2,LA
        FPSQ=FPSQ+FPTMP(L)*FPTMP(L)
      ENDDO
C
      ITER=1
C
  200 CONTINUE
      RSQ=0.
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      IF(ITER.EQ.1) RPT=1.0
	IF(ITER.GT.1) RPT=1.0/(1.0-0.25*RJ2*RPT)
C
      DO L=2,LA
      K=IL(L)+JL(L)
      IVAL=MOD(K,2)
      IF(IVAL.EQ.0)THEN
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     &        +CCN(L)*P(LN)-FPTMP(L)
        P(L)=P(L)-RPT*RSD/CCC(L)
        RSQ=RSQ+RSD*RSD
      ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
C
      IF(ITER.EQ.1) RPT=1.0/(1.0-0.5*RJ2)
	IF(ITER.GT.1) RPT=1.0/(1.0-0.25*RJ2*RPT)
C
      DO L=2,LA
      K=IL(L)+JL(L)
      IVAL=MOD(K,2)
      IF(IVAL.NE.0)THEN
        LN=LNC(L)
        LS=LSC(L)
        RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     &        +CCN(L)*P(LN)-FPTMP(L)
        P(L)=P(L)-RPT*RSD/CCC(L)
        RSQ=RSQ+RSD*RSD
      ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      RSQ=SQRT(RSQ)/SQRT(FPSQ)
      IF(RSQ .LE. RSQM) GOTO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
       WRITE(6,601)RSQ
       WRITE(8,600)
       WRITE(8,601)RSQ
       DO L=2,LA
         WRITE(8,800)L,CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),P(L),FPTMP(L)
       ENDDO

       STOP
      ENDIF
C  
      ITER=ITER+1
      GOTO 200
C
  400 CONTINUE
C
C      DO L=2,LA
C	 P(L)=P(L)+PAVG
C      END DO
C
C**********************************************************************C
C
  600 FORMAT(' MAX ITERATIONS EXCEEDED IN EXTERNAL SOLUTION, relax2t')
  601 FORMAT(' RSQ = ',E14.5)
  800 FORMAT(I6,7E13.5)
C
C**********************************************************************C
C
      RETURN
      END
