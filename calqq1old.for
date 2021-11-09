C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALQQ1OLD (ISTL)
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
C **  SUBROUTINE CALQQ CALCULATES THE TURBULENT INTENSITY SQUARED AT 
C **  TIME LEVEL (N+1).  THE VALUE OF ISTL INDICATES THE NUMBER OF 
C **  TIME LEVELS INVOLVED
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
CDIAG      DIMENSION UUUTMP(LCM),EQTMP(LCM)
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
      IF(ISTL.EQ.3)THEN
       IF(ISCDCA(0).EQ.2)THEN
        DO K=1,KS
        DO L=2,LA
        QQ2(L,K)=QQ1(L,K)+QQ(L,K)
        QQL2(L,K)=QQL1(L,K)+QQL(L,K)
        ENDDO
        ENDDO
       ELSE
        DO K=1,KS
        DO L=2,LA
        QQ2(L,K)=QQ1(L,K)+QQ1(L,K)
        QQL2(L,K)=QQL1(L,K)+QQL1(L,K)
        ENDDO
        ENDDO
       ENDIF
      ENDIF
C
      DO K=1,KC
       DO L=1,LC
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FUHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH TRANSPORT
C **  AVERAGED BETWEEN (N) AND (N+1) AND TRANSPORTED FIELD AT (N) OR
C **  TRANSPORT BETWEEN (N-1) AND (N+1) AND TRANSPORTED FIELD AT (N-1)
C **  FOR ISTL EQUAL TO 2 AND 3 RESPECTIVELY 
C
C **  FUHQQ=FUHU, FVHQQ=FVHU, FUHQQL=FUHV, FVHQQL=FVHV
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KC
      DO L=2,LA
        WB=0.5*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)=MAX(WB,0.)*QQ1(L,K-1)
     &         +MIN(WB,0.)*QQ1(L,K)
        FWQQL(L,K)=MAX(WB,0.)*QQL1(L,K-1)
     &          +MIN(WB,0.)*QQL1(L,K)
      ENDDO
      ENDDO
C
      ELSE
C
      IF(ISCDCA(0).EQ.1)THEN
C
      DO K=1,KC
      DO L=2,LA
        WB=0.25*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)=WB*(QQ(L,K-1)+QQ(L,K))
        FWQQL(L,K)=WB*(QQL(L,K-1)+QQL(L,K))
      ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KC
      DO L=2,LA
        WB=0.25*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)=MAX(WB,0.)*QQ2(L,K-1)
     &         +MIN(WB,0.)*QQ2(L,K)
        FWQQL(L,K)=MAX(WB,0.)*QQL2(L,K-1)
     &          +MIN(WB,0.)*QQL2(L,K)
      ENDDO
      ENDDO
C
      ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KS
      DO L=2,LA
        LS=LSC(L)
        UHUW=0.5*(UHDY2(L,K)+UHDY2(L,K+1))
        VHVW=0.5*(VHDX2(L,K)+VHDX2(L,K+1))
        FUHU(L,K)=MAX(UHUW,0.)*QQ1(L-1,K)
     &          +MIN(UHUW,0.)*QQ1(L,K)
        FVHU(L,K)=MAX(VHVW,0.)*QQ1(LS,K)
     &          +MIN(VHVW,0.)*QQ1(L,K)
        FUHV(L,K)=MAX(UHUW,0.)*QQL1(L-1,K)
     &          +MIN(UHUW,0.)*QQL1(L,K)
        FVHV(L,K)=MAX(VHVW,0.)*QQL1(LS,K)
     &          +MIN(VHVW,0.)*QQL1(L,K)
      ENDDO
      ENDDO
C
      ELSE
C
      IF(ISCDCA(0).EQ.1)THEN
C
      DO K=1,KS
      DO L=2,LA
        LS=LSC(L)
        UHUW=0.25*(UHDY2(L,K)+UHDY2(L,K+1))
        VHVW=0.25*(VHDX2(L,K)+VHDX2(L,K+1))
        FUHU(L,K)=UHUW*(QQ(L-1,K)+QQ(L,K))
        FVHU(L,K)=VHVW*(QQ(LS,K)+QQ(L,K))
        FUHV(L,K)=UHUW*(QQL(L-1,K)+QQL(L,K))
        FVHV(L,K)=VHVW*(QQL(LS,K)+QQL(L,K))
      ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KS
      DO L=2,LA
        LS=LSC(L)
        UHUW=0.25*(UHDY2(L,K)+UHDY2(L,K+1))
        VHVW=0.25*(VHDX2(L,K)+VHDX2(L,K+1))
        FUHU(L,K)=MAX(UHUW,0.)*QQ2(L-1,K)
     &         +MIN(UHUW,0.)*QQ2(L,K)
        FVHU(L,K)=MAX(VHVW,0.)*QQ2(LS,K)
     &         +MIN(VHVW,0.)*QQ2(L,K)
        FUHV(L,K)=MAX(UHUW,0.)*QQL2(L-1,K)
     &         +MIN(UHUW,0.)*QQL2(L,K)
        FVHV(L,K)=MAX(VHVW,0.)*QQL2(LS,K)
     &         +MIN(VHVW,0.)*QQL2(L,K)
      ENDDO
      ENDDO
C
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE PRODUCTION, LOAD BOUNDARY CONDITIONS AND SOLVE 
C **  TRANSPORT EQUATIONS
C
C **  FUHQQ=FUHU, FVHQQ=FVHU, FUHQQL=FUHV, FVHQQL=FVHV
C
C **  CU1=CUQ, CU2=CUQL, UUU=QQH, VVV=QQLH
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF(FVHU(LN,K).GT.0)THEN
      FVHU(LN,K)=0.0
      FVHV(LN,K)=0.0
      ENDIF
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(FUHU(L+1,K).GT.0)THEN
      FUHU(L+1,K)=0.0
      FUHV(L+1,K)=0.0
      ENDIF
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(FUHU(L,K).LT.0.)THEN
      FUHU(L,K)=0.0
      FUHV(L,K)=0.0
      ENDIF
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      IF(FVHU(L,K).LT.0.)THEN
      FVHU(L,K)=0.0
      FVHV(L,K)=0.0
      ENDIF
      ENDDO
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISWAVE.LE.1)THEN
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      UUU(L,K)=QQ1(L,K)*H1P(L)
     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)
      VVV(L,K)=QQL1(L,K)*H1P(L)
     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      UUU(L,K)=MAX(UUU(L,K),0.)
      VVV(L,K)=MAX(VVV(L,K),0.)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQB=AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
      PQQU=AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &    +AV(L,K)*DZIGSD4(K)*(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 
CJMH      RFNUM=-PQQB/PQQU
      PQQ=DELT*(PQQB+PQQU)
      UUU(L,K)=UUU(L,K)+2.*PQQ
      VVV(L,K)=VVV(L,K)+CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      UUU(L,K)=QQ1(L,K)*H2P(L)
     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)
      VVV(L,K)=QQL1(L,K)*H2P(L)
     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      UUU(L,K)=MAX(UUU(L,K),0.)
      VVV(L,K)=MAX(VVV(L,K),0.)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      PQQB=AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
      PQQU=AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &    +AV(L,K)*DZIGSD4(K)*(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2
CJMH      RFNUM=-PQQB/PQQU
      PQQ=DELT*(PQQB+PQQU)
      UUU(L,K)=UUU(L,K)+2.*PQQ
      VVV(L,K)=VVV(L,K)+CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISWAVE.GE.2)THEN
C
       IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
        ELSE
         WVFACT=1.0
       ENDIF      
C
      IF(ISTL.EQ.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1W(L,K)=(WVDTKEM(K)*WVDISP(L,K)
     &               +WVDTKEP(K)*WVDISP(L,K+1))
        ENDDO
       ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQ=DELT*( AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
     &    +AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &    +AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2
     &    +WVFACT*TVAR1W(L,K) )
      UUU(L,K)=QQ1(L,K)*H1P(L)
     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)+2.*PQQ
      VVV(L,K)=QQL1(L,K)*H1P(L)
     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
     &      +CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
C
      ELSE
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
         TVAR1W(L,K)=(WVDTKEM(K)*WVDISP(L,K)
     &               +WVDTKEP(K)*WVDISP(L,K+1))
        ENDDO
       ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQ=DELT*(AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
     &    +AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &    +AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2
     &    +WVFACT*TVAR1W(L,K) )
      UUU(L,K)=QQ1(L,K)*H2P(L)
     &      +DELT*(FUHU(L,K)-FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K)
     &      +(FWQQ(L,K)-FWQQ(L,K+1))*DZIG(K))*DXYIP(L)+2.*PQQ
      VVV(L,K)=QQL1(L,K)*H2P(L)
     &      +DELT*(FUHV(L,K)-FUHV(L+1,K)+FVHV(L,K)-FVHV(LN,K)
     &      +(FWQQL(L,K)-FWQQL(L,K+1))*DZIG(K))*DXYIP(L)
     &      +CTE1*DML(L,K)*PQQ
      ENDDO
      ENDDO
C
      ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(KC.LE.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP
     &       +2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP
     &       +DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.
     &       +CTE4*DML(L,1)*DML(L,1)*FPROX(1))
        EQ=1./CMQTMP
CDIAG        EQTMP(L)=EQ
        EQL=1./CMQLTMP
        CU1(L,1)=CUQTMP*EQ
        CU2(L,1)=CUQTMP*EQL
CDIAG        UUUTMP(L)=UUU(L,1)
        UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0)
     &                    -CUQTMP*HP(L)*QQ(L,KC))*EQ
        VVV(L,1)=VVV(L,1)*EQL
       ENDDO
      ENDDO
C
CDIAG      IF(N.EQ.300)THEN
CDIAG        OPEN(1,FILE='TURB.DIA',STATUS='UNKNOWN')
CDIAG        CLOSE(1,STATUS='DELETE')
CDIAG        OPEN(1,FILE='TURB.DIA',STATUS='UNKNOWN')
CDIAG        WRITE(1,110)
CDIAG        DO L=2,LA
CDIAG         QQTMP=UUU(L,1)*HPI(L)
CDIAG         WRITE(1,111) IL(L),JL(L),QQ(L,0),QQTMP,
CDIAG     &                QQ(L,KC),UUUTMP(L),EQTMP(L)
CDIAG        ENDDO
CDIAG        CLOSE(1)
CDIAG      ENDIF
C
      ENDIF
C
      IF(KC.GT.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
        CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
        CMQTMP=1.-CLQTMP-CUQTMP
     &       +2.*DELT*SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L))
        CMQLTMP=1.-CLQTMP-CUQTMP
     &       +DELT*(SQRT(QQ(L,1))/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.
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
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        DO L=LF,LL
         CLQTMP=-DELT*CDZKK(K)*AQ(L,K)*HPI(L)
         CUQTMP=-DELT*CDZKKP(K)*AQ(L,K+1)*HPI(L)
         CMQTMP=1.-CLQTMP-CUQTMP
     &       +2.*DELT*SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L))
         CMQLTMP=1.-CLQTMP-CUQTMP
     &       +DELT*(SQRT(QQ(L,K))/(CTURBB1(L,K)*DML(L,K)*HP(L)))*(1.
     &       +CTE4*DML(L,K)*DML(L,K)*FPROX(K))
         EQ=1./(CMQTMP-CLQTMP*CU1(L,K-1))
         EQL=1./(CMQLTMP-CLQTMP*CU2(L,K-1))
         CU1(L,K)=CUQTMP*EQ
         CU2(L,K)=CUQTMP*EQL     
C        IF(EQ.0) PAUSE
         UUU(L,K)=(UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
         VVV(L,K)=(VVV(L,K)-CLQTMP*VVV(L,K-1))*EQL
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KS-1,1,-1
        DO L=LF,LL
         UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
         VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      QQ1(L,K)=S2TL*QQ1(L,K)+S3TL*QQ(L,K)
      QQHDH=UUU(L,K)*HPI(L)
      QQ(L,K)=MAX(QQHDH,QQMIN)
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
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
C      QQMXSV=-1.E+12
C      QQMNSV=1.E+12
C      QQLMXSV=-1.E+12
C      QQLMNSV=1.E+12
C      
C      DO K=1,KS
C       DO L=2,LA
C        QQ1(L,K)=S2TL*QQ1(L,K)+S3TL*QQ(L,K)
C        QQL1(L,K)=S2TL*QQL1(L,K)+S3TL*QQL(L,K)
C        QQTMP=UUU(L,K)*HPI(L)
C        QQLTMP=VVV(L,K)*HPI(L)
C      QQMXSV=MAX(QQTMP,QQMXSV)
C      QQMNSV=MIN(QQTMP,QQMNSV)
C      IF(QQMNSV.LT.0.)THEN
C        WRITE(8,601)IL(L),JL(L),K,QQMNSV
C      ENDIF
C      QQLMXSV=MAX(QQLTMP,QQLMXSV)
C      QQLMNSV=MIN(QQLTMP,QQLMNSV)
C        QQ(L,K)=MAX(QQTMP,QQMIN)
C        QQL(L,K)=MAX(QQLTMP,QQLMIN)
C        DMLTMP=QQL(L,K)/QQ(L,K)
C        DML(L,K)=MAX(DMLTMP,DMLMIN)
C       ENDDO
C      ENDDO
C
C      IF(ISTOPT(0).EQ.1)THEN
C        DO K=1,KS
C        DO L=2,LA
C          DELB=B(L,K)-B(L,K+1)
C          IF(DELB.GT.0.0)THEN
C            DMLMAX=CTE3*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
C            DML(L,K)=MIN(DMLMAX,DML(L,K))
C          ENDIF
C        ENDDO
C        ENDDO
C      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      QQ(L,K)=QQ(LN,K)
      QQL(L,K)=QQL(LN,K)
      DML(L,K)=DML(LN,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBW
      L=LCBW(LL)
      QQ(L,K)=QQ(L+1,K)
      QQL(L,K)=QQL(L+1,K)
      DML(L,K)=DML(L+1,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBE
      L=LCBE(LL)
      QQ(L,K)=QQ(L-1,K)
      QQL(L,K)=QQL(L-1,K)
      DML(L,K)=DML(L-1,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      QQ(L,K)=QQ(LS,K)
      QQL(L,K)=QQL(LS,K)
      DML(L,K)=DML(LS,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
  110 FORMAT('    I    J   QQ BOT        QQ MID        QQ SURF',
     &       '       PROD+ADV      1./DIAGON')
  111 FORMAT(2I5,5E14.5)
C
  600 FORMAT('N,QX,QN,QLX,QLN,CX,CN=',I5,6E12.4)
  601 FORMAT('NEG QQ I,J,K,QQ=',3I5,E13.5)
C
C----------------------------------------------------------------------C
C
      RETURN
      END
