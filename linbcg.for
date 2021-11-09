C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LINBCG (N,B,X,ITOL,TOL,ITMAX,ITER,ERR,IRVEC)
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
C                LINBCG (LC,FPTMP,P,ITMP,RSQM,ITERM,ITER,RSQ,IRVEC)
C
      INCLUDE 'EFDC.PAR'
C
      PARAMETER (EPS=1.D-14)
C
CU    USES ATIMES,ASOLVE,SNRM
C
      DIMENSION P(LCM),PP(LCM),R(LCM),RR(LCM),Z(LCM),ZZ(LCM),
     &          B(LCM),X(LCM)
C
      ITER=0
      CALL ATIMES(X,R)
      DO 11 J=1,N
        R(J)=B(J)-R(J)
        RR(J)=R(J)
11    CONTINUE
      IF(IRVEC.GE.50) CALL ATIMES(R,RR)
      ZNRM=1.D0
      IF(ITOL.EQ.1)THEN
        BNRM=SNRM(N,B,ITOL)
      ELSE IF(ITOL.EQ.2)THEN
        CALL ASOLVE(B,Z)
        BNRM=SNRM(N,Z,ITOL)
      ELSE IF(ITOL.EQ.3.OR.ITOL.EQ.4)THEN
        CALL ASOLVE(B,Z)
        BNRM=SNRM(N,Z,ITOL)
        CALL ASOLVE(R,Z)
        ZNRM=SNRM(N,Z,ITOL)
      ELSE
        PAUSE 'ILLEGAL ITOL IN LINBCG'
      ENDIF
      CALL ASOLVE(R,Z)
100   IF(ITER.LE.ITMAX)THEN
        ITER=ITER+1
        ZM1NRM=ZNRM
        CALL ASOLVE(RR,ZZ)
        BKNUM=0.D0
        DO 12 J=1,N
          BKNUM=BKNUM+Z(J)*RR(J)
12      CONTINUE
        IF(ITER.EQ.1)THEN
          DO 13 J=1,N
            P(J)=Z(J)
            PP(J)=ZZ(J)
13        CONTINUE
        ELSE
          BK=BKNUM/BKDEN
          DO 14 J=1,N
            P(J)=BK*P(J)+Z(J)
            PP(J)=BK*PP(J)+ZZ(J)
14        CONTINUE
        ENDIF
        BKDEN=BKNUM
        CALL ATIMES(P,Z)
        AKDEN=0.D0
        DO 15 J=1,N
          AKDEN=AKDEN+Z(J)*PP(J)
15      CONTINUE
        AK=BKNUM/AKDEN
        CALL ATIMES(PP,ZZ)
        DO 16 J=1,N
          X(J)=X(J)+AK*P(J)
          R(J)=R(J)-AK*Z(J)
          RR(J)=RR(J)-AK*ZZ(J)
16      CONTINUE
        CALL ASOLVE(R,Z)
        IF(ITOL.EQ.1.OR.ITOL.EQ.2)THEN
          ZNRM=1.D0
          ERR=SNRM(N,R,ITOL)/BNRM
        ELSE IF(ITOL.EQ.3.OR.ITOL.EQ.4)THEN
          ZNRM=SNRM(N,Z,ITOL)
          IF(ABS(ZM1NRM-ZNRM).GT.EPS*ZNRM)THEN
            DXNRM=ABS(AK)*SNRM(N,P,ITOL)
            ERR=ZNRM/ABS(ZM1NRM-ZNRM)*DXNRM
          ELSE
            ERR=ZNRM/BNRM
            GOTO 100
          ENDIF
          XNRM=SNRM(N,X,ITOL)
          IF(ERR.LE.0.5D0*XNRM)THEN
            ERR=ERR/XNRM
          ELSE
            ERR=ZNRM/BNRM
            GOTO 100
          ENDIF
        ENDIF
C       WRITE (*,*) ' ITER=',ITER,' ERR=',ERR
      IF(ERR.GT.TOL) GOTO 100
      ENDIF
      RETURN
      END
