C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RELAX (ISTL)
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
C **  SELECT TIME LEVEL
C
      IF(ISTL.EQ.3) GOTO 250
C
C**********************************************************************C
C
C **  BEGIN TWO TIME LEVEL (CRANK-NICHOLSON) SOR PROCEDURE
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
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CCSR(L)*PBLK(LS)+CCWR(L)*PBLK(LW)+CCER(L)*PBLK(LE)
     &   +CCNR(L)*PBLK(LN)-FPR(L)
      PRED(L)=PRED(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      ENDDO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CCSB(L)*PRED(LS)+CCWB(L)*PRED(LW)+CCEB(L)*PRED(LE)
     &   +CCNB(L)*PRED(LN)-FPB(L)
      PBLK(L)=PBLK(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      ENDDO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      IF(RSQ .LE. RSQM) GOTO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
       STOP
      ENDIF
C  
      ITER=ITER+1
      GOTO 200
C     
C**********************************************************************C
C**********************************************************************C
C
  250 CONTINUE
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN THREE TIME LEVEL (LEAP-FROG) SOR PROCEDURE
C
      ITER=1
C
  300 CONTINUE
      RSQ=0.
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CSR(L)*PBLK(LS)+CWR(L)*PBLK(LW)+CER(L)*PBLK(LE)
     &   +CNR(L)*PBLK(LN)-FPR(L)
      PRED(L)=PRED(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      ENDDO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CSB(L)*PRED(LS)+CWB(L)*PRED(LW)+CEB(L)*PRED(LE)
     &   +CNB(L)*PRED(LN)-FPB(L)
      PBLK(L)=PBLK(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      ENDDO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      IF(RSQ .LE. RSQM) GOTO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITERM)THEN
       WRITE(6,600)
       STOP
      ENDIF
C  
      ITER=ITER+1
      GOTO 300
C     
C**********************************************************************C
C
C **  LOAD PRED AND PBLK INTO P
C
  400 CONTINUE
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
C**********************************************************************C
C
  600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
C**********************************************************************C
C
      RETURN
      END
