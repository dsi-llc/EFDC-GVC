C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTRWQGVC (M,NW,CON,CON1)
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
C **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE
C **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
C **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES 
C **  THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION CON(LCM,KCM), CON1(LCM,KCM)
      DIMENSION CONTMX(LCM,KCM),CONTMN(LCM,KCM)
      DIMENSION FQCPAD(LCM,KCM),QSUMPAD(LCM,KCM),QSUMNAD(LCM,KCM)
C     DIMENSION  CH(LCM,KCM), CMAX(LCM,KCM), CMIN(LCM,KCM),
C    &         CONT(LCM,KCM), CON2(LCM,KCM)
C
C**********************************************************************C
C
      MVAR=M
C
      BSMALL=1.0E-6
C
      DELT=DT2
      DELTA=DT2
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      ISTL=3
C
      IF(IS2TIM.GE.1) THEN
	  ISTL=2
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
          DELTA=DT
          DELTD2=0.5*DT
          S3TL=0.0
          S2TL=1.0
          ISUD=0
         ELSE
          DELT=DTDYN
          DELTA=DTDYN
          DELTD2=0.5*DTDYN
          S3TL=0.0
          S2TL=1.0
          ISUD=0
        END IF
	ENDIF
C
      IF(ISLTMT.GE.1) ISUD=1
C
      DO K=1,KC
       DO L=1,LC
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FVHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       ENDDO
      ENDDO
C
C
      DO L=1,LC
        FWU(L,0)=0.
        FWU(L,KC)=0.
      ENDDO

C
C**********************************************************************C
C
      DO K=1,KC
      DO L=1,LC
      CONT(L,K)=0.
      CMAX(L,K)=0.
      CMIN(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATED EXTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
C      CALL CALFQC (ISTL,MVAR,M,CON,CON1)
      CALL CALFQC (ISTL,0,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,
     &                   QSUMNAD) 
C
CWQJH     IF(ISTRAN(M).GE.1)THEN
CWQJH       IF(M.EQ.7) CALL CALFQC (ISTL,M,CON,CON1)
CWQJH       IF(M.EQ.8)THEN
CWQJH         DO K=1,KC 
CWQJH          DO L=2,LA
CWQJH           FQC(L,K)=0.
CWQJH          ENDDO
CWQJH         ENDDO
CWQJH       ENDIF
CWQJH     ENDIF
C
C**********************************************************************C
C
C **  SELECT TRANSPORT OPTION, ISPLIT=1 FOR HORIZONTAL-VERTICAL
C **  OPERATOR SPLITTING
C
C     IF(ISPLIT(M).EQ.1) GOTO 1000
C
C**********************************************************************C     
C**********************************************************************C
C
C **  BEGIN COMBINED ADVECTION SCHEME
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
CWQJH     CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H2WQ(L)
      CONT(L,K)=CON1(L,K)
      ENDDO
      ENDDO
C
CWQJH     WFQC=0.5
C
C----------------------------------------------------------------------C
C
      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
      FUHU(L,K)=MAX(UHDYWQ(L,K),0.)*CONT(L-1,K)
     &         +MIN(UHDYWQ(L,K),0.)*CONT(L,K)
      FVHU(L,K)=MAX(VHDXWQ(L,K),0.)*CONT(LS,K)
     &         +MIN(VHDXWQ(L,K),0.)*CONT(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(WWQ(L,K),0.)*CONT(L,K)
     &        +MIN(WWQ(L,K),0.)*CONT(L,K+1)
      ENDDO
      ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     FWU(L,K)=0.5*WWQ(L,K)*(CON(L,K+1)+CON(L,K))
C     ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX
C
C----------------------------------------------------------------------C
C
C     IF(ISTRAN(M).GE.1) CALL CALDIFF (ISTL,M,CON1)
C
C**********************************************************************C
C
C **  ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      DO K=1,KC
       DO L=2,LA
	   FUHU(L,K)=GVCSCLU(L)*FUHU(L,K)
	   FVHU(L,K)=GVCSCLV(L)*FVHU(L,K)
	 ENDDO
	ENDDO

c
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)      
CWQJH      CH(L,K)=CONT(L,K)*H2WQ(L)
CWQJH     &       +DELT*((WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CONT(L,K)*GVCSCLP(L)*H2WQ(L)
     &       +DELT*((RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &                       +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     &                       +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      ENDDO
      ENDDO
C
      IF(ISFCT(M).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*GVCSCLPI(L)/HWQ(L)+(1.-SCB(L))*CON1(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES FORM M=4
C
      IF(M.EQ.4)THEN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
C
      IF(VHDXWQ(LN,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(GVCSCLV(LN)*VHDXWQ(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
       CLOS(LL,K,M)=CON(L,K)        
       NLOS(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOS(LL,K,M)
       IF(NMNLO.GE.NTSCRS(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOS(LL,K,M)
     &         +(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
C
      IF(UHDYWQ(L+1,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(GVCSCLU(L+1)*UHDYWQ(L+1,K)*CON1(L,K)
     &      -FUHU(L+1,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOW(LL,K,M)
       IF(NMNLO.GE.NTSCRW(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOW(LL,K,M)
     &         +(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)      
C
      IF(UHDYWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &   -GVCSCLU(L)*UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOE(LL,K,M)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOE(LL,K,M)
     &         +(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
C
      IF(VHDXWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &   -GVCSCLV(L)*VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLON(LL,K,M)
       IF(NMNLO.GE.NTSCRN(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLON(LL,K,M)
     &         +(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES FOR WQVARS, M=8
C
      IF(M.EQ.8)THEN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBS
      NSID=IWQOBS(LL,NW)
      L=LIJ( IWQCBS(LL),JWQCBS(LL) )
      LN=LNC(L)
C
      IF(VHDXWQ(LN,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(GVCSCLV(LN)*VHDXWQ(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
       CWQLOS(LL,K,NW)=CON(L,K)        
       NWQLOS(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCS(LL,1,NW)+WTCI(K,2)*WQOBCS(LL,2,NW)
     &    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOS(LL,K,NW)
       IF(NMNLO.GE.NTSCRS(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOS(LL,K,NW)
     &         +(CBT-CWQLOS(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBW
      NSID=IWQOBW(LL,NW)
      L=LIJ( IWQCBW(LL),JWQCBW(LL) )
C
      IF(UHDYWQ(L+1,K).LT.0.)THEN
       CTMP=CON1(L,K)+DELT*(GVCSCLU(L+1)*UHDYWQ(L+1,K)*CON1(L,K)
     &      -FUHU(L+1,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
       CWQLOW(LL,K,NW)=CON(L,K)
       NWQLOW(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCW(LL,1,NW)+WTCI(K,2)*WQOBCW(LL,2,NW)
     &    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOW(LL,K,NW)
       IF(NMNLO.GE.NTSCRW(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOW(LL,K,NW)
     &         +(CBT-CWQLOW(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBE
      NSID=IWQOBE(LL,NW)
      L=LIJ( IWQCBE(LL),JWQCBE(LL) )
C
      IF(UHDYWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &   -GVCSCLU(L)*UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
       CWQLOE(LL,K,NW)=CON(L,K)
       NWQLOE(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCE(LL,1,NW)+WTCI(K,2)*WQOBCE(LL,2,NW)
     &    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOE(LL,K,NW)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOE(LL,K,NW)
     &         +(CBT-CWQLOE(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBN
      NSID=IWQOBN(LL,NW)
      L=LIJ( IWQCBN(LL),JWQCBN(LL) )
      LS=LSC(L)
C
      IF(VHDXWQ(L,K).GT.0.)THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &  -GVCSCLV(L)*VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)*GVCSCLPI(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF(CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
       CWQLON(LL,K,NW)=CON(L,K)
       NWQLON(LL,K,NW)=N
      ELSE
       CBT=WTCI(K,1)*WQOBCN(LL,1,NW)+WTCI(K,2)*WQOBCN(LL,2,NW)
     &    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLON(LL,K,NW)
       IF(NMNLO.GE.NTSCRN(LL))THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLON(LL,K,NW)
     &         +(CBT-CWQLON(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
      IF(ISADAC(MVAR).EQ.0) RETURN
      IF(ISCDCA(MVAR).EQ.1) RETURN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=1,LC
      UUU(L,K)=0.0
      VVV(L,K)=0.0
	WWW(L,K)=0.0
      ENDDO
      ENDDO
C
      DO L=1,LC
	WWW(L,0)=0.0
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      UUU(L,K)=UWQ(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
      VVV(L,K)=VWQ(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      ENDDO
      ENDDO
C
C MODIFIED FOR GVC - JMH 08/04/04
C

      DO K=1,KS
      RDZIG=DZIG(K)
      DO L=2,LA
      WWW(L,K)=WWQ(L,K)*(CON(L,K+1)-CON(L,K))*RDZIG*GVCSCLPI(L)/HWQ(L)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      AUHU=GVCSCLU(L)*ABS(UHDYWQ(L,K))
      AVHV=GVCSCLV(L)*ABS(VHDXWQ(L,K))
c      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
c      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
      UTERM=AUHU*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(CON(L,K)-CON(LS,K))
      IF(ISADAC(MVAR).GE.2)THEN
C
C MODIFIED FOR GVC - JMH 08/04/04
C
        SSCORUE=DELTA*RDZIC*DXYIP(L  )*GVCSCLPI(L  )*HWQI(L  )
     &               *(FQCPAD(L  ,K)-QSUMPAD(L  ,K)*CON(L  ,K))
        SSCORUW=DELTA*RDZIC*DXYIP(L-1)*GVCSCLPI(L-1)*HWQI(L-1)
     &               *(FQCPAD(L-1,K)-QSUMPAD(L-1,K)*CON(L-1,K))
        SSCORVN=DELTA*RDZIC*DXYIP(L  )*GVCSCLPI(L  )*HWQI(L  )
     &               *(FQCPAD(L  ,K)-QSUMPAD(L  ,K)*CON(L  ,K))
        SSCORVS=DELTA*RDZIC*DXYIP(LS )*GVCSCLPI(LS )*HWQI(LS )
     &               *(FQCPAD(LS ,K)-QSUMPAD(LS ,K)*CON(LS ,K))
        SSCORU=MAX(UHDYWQ(L,K),0.0)*SSCORUW+MIN(UHDYWQ(L,K),0.0)*SSCORUE
        SSCORV=MAX(VHDXWQ(L,K),0.0)*SSCORVS+MIN(VHDXWQ(L,K),0.0)*SSCORVN
	  UTERM=UTERM+GVCSCLU(L)*SSCORU
	  VTERM=VTERM+GVCSCLV(L)*SSCORV
	ENDIF
CADD1023
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
C     VTERM=VTERM*(CON(L,K)-CON(LS,K))
C
CADD1031      UTERM=UTERM-0.25*DELTA*UHDYWQ(L,K)*
CADD1031     &      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
CADD1031     &      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
CADD1031     VTERM=VTERM-0.25*DELTA*VHDXWQ(L,K)*
CADD1031     &      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
CADD1031     &      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
C
C MODIFIED FOR GVC - JMH 08/04/04
C
      IF(UHDYWQ(L,K).GE.0.0)THEN
        UTERM=UTERM-0.5*DELTA*GVCSCLU(L)*UHDYWQ(L,K)*
     &      (VVV(LNW,K)+VVV(L-1,K)+WWW(L-1,K)+WWW(L-1,K-1)
     &      +UUU(L,K)+UUU(L-1,K))
      ELSE
        UTERM=UTERM-0.5*DELTA*GVCSCLU(L)*UHDYWQ(L,K)*
     &      (VVV(LN,K)+VVV(L,K)+WWW(L,K)+WWW(L,K-1)
     &      +UUU(L,K)+UUU(L+1,K))
	ENDIF
C
      IF(VHDXWQ(L,K).GE.0.0)THEN
        VTERM=VTERM-0.5*DELTA*GVCSCLV(L)*VHDXWQ(L,K)*
     &      (UUU(LS,K)+UUU(LSE,K)+WWW(LS,K)+WWW(LS,K-1)
     &      +VVV(LS,K)+VVV(L,K))
      ELSE
        VTERM=VTERM-0.5*DELTA*GVCSCLV(L)*VHDXWQ(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+WWW(L,K)+WWW(L,K-1)
     &      +VVV(LN,K)+VVV(L,K))
	ENDIF
C
      IF(ISFCT(MVAR).GE.2)THEN
        FUHU(L,K)=0.5*UTERM
        FVHU(L,K)=0.5*VTERM
        IF(ISFCT(MVAR).EQ.3)THEN
C JMH WHY ARE WE NOT DIVIDING BY 2 HERE
          FUHU(L,K)=UTERM
          FVHU(L,K)=VTERM
        ENDIF
      ELSE
        UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
        VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
C     AUHU=ABS(UHU)
C     AVHV=ABS(VHV)
C     UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
C     VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
        FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)
     &         +MIN(UHU,0.)*CON(L,K)
        FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)
     &         +MIN(VHV,0.)*CON(L,K)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(WWQ(L,K))
C      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)/HWQ(L))*(CON(L,K+1)-CON(L,K))
      WTERM=AWW*(CON(L,K+1)-CON(L,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORWA=DELTA*DZIG(K+1)*GVCSCLPI(L)*HWQI(L)*DXYIP(L)
     &         *(FQCPAD(L,K+1)-QSUMPAD(L,K+1)*CON(L,K+1))
        SSCORWB=DELTA*DZIG(K)*GVCSCLPI(L)*HWQI(L)*DXYIP(L)
     &         *(FQCPAD(L,K  )-QSUMPAD(L,K  )*CON(L,K  ))
        SSCORW=MAX(WWQ(L,K),0.0)*SSCORWB+MIN(WWQ(L,K),0.0)*SSCORWA
	  WTERM=WTERM+SSCORW
      ENDIF
C
CCC1031      WTERM=WTERM-0.25*DELTA*WWQ(L,K)*
CCC1031     &      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
CCC1031     &      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
C
      IF(WWQ(L,K).GE.0.0)THEN
        WTERM=WTERM-0.5*DELTA*WWQ(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+VVV(L,K)+VVV(LN,K)
     &      +WWW(L,K)+WWW(L,K-1))
	ELSE
        WTERM=WTERM-0.5*DELTA*WWQ(L,K)*
     &      (UUU(L+1,K+1)+UUU(L,K+1)+VVV(LN,K+1)+VVV(L,K+1)
     &      +WWW(L,K)+WWW(L,K+1))
	ENDIF
C
      IF(ISFCT(MVAR).GE.2)THEN
        FWU(L,K)=0.5*WTERM
        IF(ISFCT(MVAR).EQ.3)THEN
          FWU(L,K)=WTERM
        ENDIF
      ELSE
        WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
        FWU(L,K)=MAX(WW,0.)*CON(L,K)
     &        +MIN(WW,0.)*CON(L,K+1)
	ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
c  ** ANTIDIFFUSION TURNED OFF FOR SOURCE CELLS
C
      IF(ISADAC(MVAR).EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
C        IF(ABS(QSUM(L,K)).GT.1.E-12)THEN
        IF(QSUMPAD(L,K).GT.0.0)THEN
          LN=LNC(L)
          FUHU(L  ,K)=0.
          FUHU(L+1,K)=0.
          FVHU(L  ,K)=0.
          FVHU(LN ,K)=0.
          FWU(L,K  )=0.
          FWU(L,K-1)=0.
          CONT(L,K)=0.
        ELSE
          CONT(L,K)=1.
        ENDIF
      ENDDO
      ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR OPEN BOUNDARY CELLS
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      FVHU(LN,K)=0.0
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      FUHU(L+1,K)=0.0
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      FUHU(L,K)=0.0
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      FVHU(L,K)=0.0
      ENDDO
C
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
C
      IF(ISFCT(MVAR).EQ.0) GOTO 600
C
C----------------------------------------------------------------------C
C
C **  DETERMINE MAX AND MIN CONCENTRATIONS
C
      DO K=1,KC
      DO L=1,LC
      CONTMX(L,K)=0.0
      CONTMN(L,K)=0.0
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
	IF(LGVCP(L,K))THEN
        CONTMX(L,K)=MAX(CON(L,K),CON2(L,K))
        CONTMN(L,K)=MIN(CON(L,K),CON2(L,K))
	ENDIF
      ENDDO
      ENDDO
C
      DO L=2,LA
	  KBTMP=KGVCP(L)
	  KLTMP=KC-KGVCP(L)+1
	  IF(KLTMP.EQ.1)THEN
          CMAX(L,KC)=CONTMX(L,KC)
          CMIN(L,KC)=CONTMN(L,KC)
	  ENDIF
	  IF(KLTMP.GE.2)THEN
          CMAX(L,KBTMP)=MAX(CONTMX(L,KBTMP),CONTMX(L,KBTMP+1))
          CMAX(L,KC)=MAX(CONTMX(L,KS),CONTMX(L,KC))
          CMIN(L,KBTMP)=MIN(CONTMN(L,KBTMP),CONTMN(L,KBTMP+1))
          CMIN(L,KC)=MIN(CONTMN(L,KS),CONTMN(L,KC))
        ENDIF
      ENDDO
C
	IF(KS.GE.2)THEN
        DO K=2,KS
        DO L=2,LA
	    KBTMPP=KGVCP(L)+1
	    IF(K.GE.KBTMPP)THEN
            CMAXT=MAX(CONTMX(L,K-1),CONTMX(L,K+1))
            CMAX(L,K)=MAX(CONTMX(L,K),CMAXT)
            CMINT=MIN(CONTMN(L,K-1),CONTMN(L,K+1))
            CMIN(L,K)=MIN(CONTMN(L,K),CMINT)
          ENDIF
        ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
	IF(LGVCP(L,K))THEN
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB3D(L,K)*CONTMX(L-1,K) 
      CEMAX=SUB3D(L+1,K)*CONTMX(L+1,K)
      CSMAX=SVB3D(L,K)*CONTMX(LS,K)
      CNMAX=SVB3D(LN,K)*CONTMX(LN,K)
      CMAXT=MAX(CNMAX,CEMAX)
      CMAXT=MAX(CMAXT,CSMAX)
      CMAXT=MAX(CMAXT,CWMAX)
      CMAX(L,K)=MAX(CMAX(L,K),CMAXT)
      CWMIN=SUB3D(L,K)*CONTMN(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB3D(L+1,K)*CONTMN(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB3D(L,K)*CONTMN(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB3D(LN,K)*CONTMN(LN,K)+1.E+6*(1.-SVB(LN))
      CMINT=MIN(CNMIN,CEMIN)
      CMINT=MIN(CMINT,CSMIN)
      CMINT=MIN(CMINT,CWMIN)
      CMIN(L,K)=MIN(CMIN(L,K),CMINT)
	ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  SEPARATE POSITIVE AND NEGATIVE FLUXES PUTTING NEGATIVE FLUXES 
C **  INTO FUHV, FVHV, AND FWV
C
      DO K=1,KC
      DO L=2,LA
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS AND LOAD
C **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      IF(LGVCP(L,K))THEN
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)
     &                               +FVHU(L,K)-FVHV(LN,K))
     &              +RDZIC*(FWU(L,K-1)-FWV(L,K)) )*GVCSCLPI(L)*HWQI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)
     &                               +FVHU(LN,K)-FVHV(L,K))
     &              +RDZIC*(FWU(L,K)-FWV(L,K-1)) )*GVCSCLPI(L)*HWQI(L)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DU(LN,K)=0.
      DV(LN,K)=0.
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DU(L+1,K)=0.
      DV(L+1,K)=0.
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DU(L-1,K)=0.
      DV(L-1,K)=0.
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      DU(LS,K)=0.
      DV(LS,K)=0.
      ENDDO
C
      ENDDO
C
CX      DO NS=1,NQSIJ
CX       L=LQS(NS)
CX       NQSTMP=NQSERQ(NS)
CX       DO K=1,KC
CX        QQQTMP=ABS(QSS(K,NS)+QSERT(K,NQSTMP))
CX        IF(QQQTMP.GE.1.E-12)THEN
CX          DU(L,K)=0.
CX          DV(L,K)=0.
C         DU(L+1,K)=0.
C         DV(LNC(L),K)=0.
CX        ENDIF       
CX       ENDDO
CX      ENDDO
C
CX      DO NCTL=1,NQCTL
CX       IU=IQCTLU(NCTL)
CX       JU=JQCTLU(NCTL)
CX       LU=LIJ(IU,JU)
CX       ID=IQCTLD(NCTL)
CX       JD=JQCTLD(NCTL)
CX       IF(ID.EQ.0.AND.JD.EQ.0)THEN
CX         LD=LC
CX        ELSE
CX         LD=LIJ(ID,JD)
CX       ENDIF
CX       DO K=1,KC
CX        QQQTMP=ABS(QCTLT(K,NCTL))
CX        IF(QQQTMP.GE.1.E-12)THEN
CX          DU(LU,K)=0.
CX          DV(LU,K)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
CX          DU(LD,K)=0.
CX          DV(LD,K)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
CX        ENDIF
CX       ENDDO
CX      ENDDO
C
CX      DO NWR=1,NQWR
CX       IU=IQWRU(NWR)
CX       JU=JQWRU(NWR)
CX       KU=KQWRU(NWR)
CX       ID=IQWRD(NWR)
CX       JD=JQWRD(NWR)
CX       KD=KQWRD(NWR)
CX       LU=LIJ(IU,JU)
CX       LD=LIJ(ID,JD)
CX        NQSTMP=NQWRSERQ(NWR)
CX        QQQTMP=ABS(QWR(NWR)+QWRSERT(NQSTMP))
CX        IF(QQQTMP.GE.1.E-12)THEN
CX          DU(LU,KU)=0.
CX          DV(LU,KU)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
CX          DU(LD,KD)=0.
CX          DV(LD,KD)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
CX        ENDIF
CX      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DU AND DV
C
      DO K=1,KC
      DO L=2,LA
      IF(LGVCP(L,K))THEN
      IF(DU(L,K).GT.0.) DU(L,K)=(CMAX(L,K)-CON(L,K))/(DU(L,K)+BSMALL)
      DU(L,K)=MIN(DU(L,K),1.)
      IF(DV(L,K).GT.0.) DV(L,K)=(CON(L,K)-CMIN(L,K))/(DV(L,K)+BSMALL)
      DV(L,K)=MIN(DV(L,K),1.)
      ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  LIMIT FLUXES
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)
     &         +MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)
     &         +MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MIN(DV(L,K),DU(L,K+1))*FWU(L,K)
     &        +MIN(DU(L,K),DV(L,K+1))*FWV(L,K)
      ENDDO
      ENDDO  
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=SWB3D(L,K)*FWU(L,K)
      ENDDO
      ENDDO  
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTION CALCULATION
C
  600 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      IF(LGVCP(L,K))THEN
      LN=LNC(L)
      CH(L,K)=CON(L,K)*GVCSCLP(L)*HWQ(L)
     &       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     &              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     &             +(FWU(L,K-1)-FWU(L,K))*RDZIC)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      IF(LGVCP(L,K))THEN
      CON(L,K)=SCB(L)*CH(L,K)*GVCSCLPI(L)/HWQ(L)+(1.-SCB(L))*CON(L,K)
      ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  DIAGNOSE FCT SCHEME
C
      IF(ISFCT(M).EQ.99)THEN
      WRITE(6,6010)N
C
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF(CCMAX.GT.0.)THEN
       WRITE(6,6011)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      ENDIF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF(CCMIN.GT.0.)THEN
       WRITE(6,6012)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C
 6010 FORMAT('  FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT('  CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT('  CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
C
C----------------------------------------------------------------------C
C
      RETURN
      END
