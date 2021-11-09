C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RELAXV (ISTL)
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
C **  GATHER-SCATTER VECTORIZABLE VERSION
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
C     DIMENSION RSDR(LCM),RSDB(LCM)
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
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      DO L=1,NRC
      RSDR(L)=PRED(L)-FPR(L)
      ENDDO
C
      DO L=1,NRC
      LS=LBSRC(L)
      RSDR(L)=RSDR(L)+CCSR(L)*PBLK(LS)
      ENDDO
C
      DO L=1,NRC
      LW=LBWRC(L)
      RSDR(L)=RSDR(L)+CCWR(L)*PBLK(LW)
      ENDDO
C
      DO L=1,NRC
      LE=LBERC(L)
      RSDR(L)=RSDR(L)+CCER(L)*PBLK(LE)
      ENDDO
C
      DO L=1,NRC
      LN=LBNRC(L)
      RSDR(L)=RSDR(L)+CCNR(L)*PBLK(LN)
      ENDDO
C
      DO L=1,NRC
      PRED(L)=PRED(L)-RP*RSDR(L)
      ENDDO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      RSDB(L)=PBLK(L)-FPB(L)
      ENDDO
C
      DO L=1,NBC
      LS=LRSBC(L)
      RSDB(L)=RSDB(L)+CCSB(L)*PRED(LS)
      ENDDO
C
      DO L=1,NBC
      LW=LRWBC(L)
      RSDB(L)=RSDB(L)+CCWB(L)*PRED(LW)
      ENDDO
C
      DO L=1,NBC
      LE=LREBC(L)
      RSDB(L)=RSDB(L)+CCEB(L)*PRED(LE)
      ENDDO
C
      DO L=1,NBC
      LN=LRNBC(L)
      RSDB(L)=RSDB(L)+CCNB(L)*PRED(LN)
      ENDDO
C
      DO L=1,NBC
      PBLK(L)=PBLK(L)-RP*RSDB(L)
      ENDDO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      RSQ=0.
      DO L=1,NRC
      RSQ=RSQ+RSDR(L)*RSDR(L)
      ENDDO
      DO L=1,NBC
      RSQ=RSQ+RSDB(L)*RSDB(L)
      ENDDO
C
      IF(RSQ .LE. RSQM) GOTO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITERM) STOP
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
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      DO L=1,NRC
      RSDR(L)=PRED(L)-FPR(L)
      ENDDO
C
      DO L=1,NRC
      LS=LBSRC(L)
      RSDR(L)=RSDR(L)+CSR(L)*PBLK(LS)
      ENDDO
C
      DO L=1,NRC
      LW=LBWRC(L)
      RSDR(L)=RSDR(L)+CWR(L)*PBLK(LW)
      ENDDO
C
      DO L=1,NRC
      LE=LBERC(L)
      RSDR(L)=RSDR(L)+CER(L)*PBLK(LE)
      ENDDO
C
      DO L=1,NRC
      LN=LBNRC(L)
      RSDR(L)=RSDR(L)+CNR(L)*PBLK(LN)
      ENDDO
C
      DO L=1,NRC
      PRED(L)=PRED(L)-RP*RSDR(L)
      ENDDO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      RSDB(L)=PBLK(L)-FPB(L)
      ENDDO
C
      DO L=1,NBC
      LS=LRSBC(L)
      RSDB(L)=RSDB(L)+CSB(L)*PRED(LS)
      ENDDO
C
      DO L=1,NBC
      LW=LRWBC(L)
      RSDB(L)=RSDB(L)+CWB(L)*PRED(LW)
      ENDDO
C
      DO L=1,NBC
      LE=LREBC(L)
      RSDB(L)=RSDB(L)+CEB(L)*PRED(LE)
      ENDDO
C
      DO L=1,NBC
      LN=LRNBC(L)
      RSDB(L)=RSDB(L)+CNB(L)*PRED(LN)
      ENDDO
C
      DO L=1,NBC
      PBLK(L)=PBLK(L)-RP*RSDB(L)
      ENDDO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      RSQ=0.
      DO L=1,NRC
      RSQ=RSQ+RSDR(L)*RSDR(L)
      ENDDO
      DO L=1,NBC
      RSQ=RSQ+RSDB(L)*RSDB(L)
      ENDDO
C
      IF(RSQ .LE. RSQM) GOTO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITERM) STOP
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
      RETURN
      END
