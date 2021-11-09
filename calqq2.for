C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALQQ2 (ISTL)
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
C **  SUBROUTINE CALQQ2 CALCULATES THE TURBULENT INTENSITY SQUARED AT 
C **  TIME LEVEL (N+1).  THE VALUE OF ISTL INDICATES THE NUMBER OF 
C **  TIME LEVELS INVOLVED.  THIS VERSION USES A SEPARATE ADVECTIVE
C **  TRANSPORT SUBROUTINE CALTRANQ
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF(ISTL.EQ.2)THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      ENDIF
      BSMALL=1.E-12
C
C**********************************************************************C
C
C **  MOVE UHDY2, VHDX2 AND W2 TO QQ TRANSPORT LOCATIONS
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      U2(L,K)=0.5*(U2(L,K)+U2(L,K+1))
      V2(L,K)=0.5*(V2(L,K)+V2(L,K+1))
      UHDY2(L,K)=0.5*(UHDY2(L,K)+UHDY2(L,K+1))
      VHDX2(L,K)=0.5*(VHDX2(L,K)+VHDX2(L,K+1))
      ENDDO
      ENDDO
C
      DO K=0,KS
      DO L=2,LA
      W2(L,K)=0.5*(W2(L,K)+W2(L,K+1))
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ELIMINATE INFLOWS ACROSS OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF(VHDX2(LN,K).GT.0.) VHDX2(LN,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO LL=1,NCBW
      L=LCBW(LL)      
      IF(UHDY2(L+1,K).GT.0.) UHDY2(L+1,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO LL=1,NCBE
      L=LCBE(LL)      
      IF(UHDY2(L,K).LT.0.) UHDY2(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      IF(VHDX2(L,K).LT.0.) VHDX2(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALL ADVECTIVE TRANSPORT SUBROUTINE FOR QQ AND QQL
C
C----------------------------------------------------------------------C
C
      CALL CALTRANQ (ISTL,0,QQ,QQ1)
C
      CALL CALTRANQ (ISTL,0,QQL,QQL1)
C
      DO L=2,LA
      W2(L,0)=0.
      W2(L,KC)=0.
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE PRODUCTION, LOAD BOUNDARY CONDITIONS AND SOLVE 
C **  TRANSPORT EQUATIONS
C
C **  CU1=CUQ, CU2=CUQL, UUU=QQH, VVV=QQLH
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQ=DELT*(AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
     &    +AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &    +AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2)
      UUU(L,K)=QQ(L,K)*HP(L)+2.*PQQ
      VVV(L,K)=QQL1(L,K)*HP(L)+CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
      CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
      CMQTMP=1.-CLQTMP-CUQTMP
     &       +2.*DELT*SQRT(QQ(L,1))/(CTURB*DML(L,1)*HP(L))
      CMQLTMP=1.-CLQTMP-CUQTMP
     &       +DELT*(SQRT(QQ(L,1))/(CTURB*DML(L,1)*HP(L)))*(1.
     &       +CTE4*DML(L,1)*DML(L,1)*FPROX(1))
      EQ=1./CMQTMP
      EQL=1./CMQLTMP
      CU1(L,1)=CUQTMP*EQ
      CU2(L,1)=CUQTMP*EQL
      UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0))*EQ
      VVV(L,1)=VVV(L,1)*EQL
      CUQTMP=-DELT*CDZKKP(KS)*AQ(L,KC)*HPI(L)
      UUU(L,KS)=UUU(L,KS)-CUQTMP*HP(L)*QQ(L,KC)
      ENDDO
C
      DO K=2,KS
      DO L=2,LA
      CLQTMP=-DELT*CDZKK(K)*AQ(L,K)*HPI(L)
      CUQTMP=-DELT*CDZKKP(K)*AQ(L,K+1)*HPI(L)
      CMQTMP=1.-CLQTMP-CUQTMP
     &       +2.*DELT*SQRT(QQ(L,K))/(CTURB*DML(L,K)*HP(L))
      CMQLTMP=1.-CLQTMP-CUQTMP
     &       +DELT*(SQRT(QQ(L,K))/(CTURB*DML(L,K)*HP(L)))*(1.
     &       +CTE4*DML(L,K)*DML(L,K)*FPROX(K))
      EQ=1./(CMQTMP-CLQTMP*CU1(L,K-1))
      EQL=1./(CMQLTMP-CLQTMP*CU2(L,K-1))
      CU1(L,K)=CUQTMP*EQ
      CU2(L,K)=CUQTMP*EQL     
      UUU(L,K)=(UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
      VVV(L,K)=(VVV(L,K)-CLQTMP*VVV(L,K-1))*EQL
      ENDDO
      ENDDO
C
      DO K=KS-1,1,-1
      DO L=2,LA
      UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
      VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      QQ1(L,K)=S2TL*QQ1(L,K)+S3TL*QQ(L,K)
      QQHDH=UUU(L,K)*HPI(L)
      QQ(L,K)=MAX(QQHDH,QQMIN)
      QQ(L,K)=SPB(L)*QQ(L,K)+(1.-SPB(L))*QQMIN
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      QQL1(L,K)=S2TL*QQL1(L,K)+S3TL*QQL(L,K)
      QQHDH=VVV(L,K)*HPI(L)
      QQL(L,K)=MAX(QQHDH,QQLMIN)
      QQL(L,K)=SPB(L)*QQL(L,K)+(1.-SPB(L))*QQLMIN
      DMLTMP=QQL(L,K)/QQ(L,K)
      DMLTMP=MAX(DMLTMP,DMLMIN)
      DELB=B(L,K)-B(L,K+1)
      IF(DELB.GT.0.0.AND.ISLLIM.EQ.2)THEN
       DMLMAX=SQRT(RIQMAX)*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
       DML(L,K)=MIN(DMLMAX,DMLTMP)
	 QQL(L,K)=QQ(L,K)*DML(L,K)
      ELSE
       DML(L,K)=DMLTMP
      ENDIF
      DML(L,K)=SPB(L)*DML(L,K)+(1.-SPB(L))*DMLMIN
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISTOPT(0).GE.2)THEN
      DO K=1,KS
       DO L=2,LA
        DELBSQ=( DZIG(K)*(B(L,K+1)-B(L,K)) )**2
        BBT(L,K)=CTURBB2(L,K)*DML(L,K)*AB(L,K)*DELBSQ/SQRT(QQ(L,K))
       ENDDO
      ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C     DO K=1,KS
C     DO LL=1,NCBS
C     LN=LNC(L)
C     QQ(L,K)=QQ(LN,K)
C     QQL(L,K)=QQL(LN,K)
C     DML(L,K)=DML(LN,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KS
C     DO LL=1,NCBW
C     L=LCBW(LL)
C     QQ(L,K)=QQ(L+1,K)
C     QQL(L,K)=QQL(L+1,K)
C     DML(L,K)=DML(L+1,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KS
C     DO LL=1,NCBE
C     L=LCBE(LL)
C     QQ(L,K)=QQ(L-1,K)
C     QQL(L,K)=QQL(L-1,K)
C     DML(L,K)=DML(L-1,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KS
C     DO LL=1,NCBN
C     L=LCBN(LL)
C     LS=LSC(L)
C     QQ(L,K)=QQ(LS,K)
C     QQL(L,K)=QQL(LS,K)
C     DML(L,K)=DML(LS,K)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
      RETURN
      END
