C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTRAN1D (ISTL,MVAR,M,CON,CON1)
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
C    &         CONT(LCM,KCM), CON2(LCM,KCM),
C
C**********************************************************************C
C
CDIA      IF(N.EQ.1)THEN
CDIA        OPEN(1,FILE='CALTRAN1D.DIA')
CDIA        CLOSE(1,STATUS='DELETE')
CDIA        OPEN(1,FILE='CALTRAN1D.DIA')
CDIA        WRITE(1,1110)N,ISTL
CDIA       ELSE
CDIA        OPEN(1,FILE='CALTRAN1D.DIA',POSITION='APPEND')
CDIA        WRITE(1,1110)N,ISTL
CDIA      ENDIF
C
 1110 FORMAT(' N,ISTL = ',2I5)
C
      BSMALL=1.0E-6
CBUG    WFQC=0.0
C
      DELT=DT2
      DELTA=DT2
      IF(ISCDCA(MVAR).EQ.2) DELTA=DT
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF(ISTL.NE.3)THEN
       DELT=DT
       DELTA=DT
       DELTD2=0.5*DT
       S3TL=0.0
       S2TL=1.0
       ISUD=1
      ENDIF
C
      DELTA4=0.25*DELTA
C
      IF(ISLTMT.GE.1) ISUD=1
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
       DO L=1,LC
        IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
        IF(LCT(L).EQ.6) AREAOLD(L)=FADYP1(L)
        IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
        IF(LCT(L).EQ.7) AREAOLD(L)=FADXP1(L)
       ENDDO
C
C**********************************************************************C
C 
        DO K=1,KC
         DO L=1,LC
         CONT(L,K)=0.
         CMAX(L,K)=0.
         CMIN(L,K)=0.
         CON1(L,K)=CON(L,K)
         ENDDO
        ENDDO
C
C**********************************************************************C
C
C **  CALCULATED EXTERNAL SOURCES AND SINKS 
C
C----------------------------------------------------------------------C
C
      CALL CALFQC (ISTL,IS2TL,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,
     &                   QSUMNAD) 
C     IF(ISTRAN(M).EQ.1) CALL CALFQC (ISTL,M,CON,CON1)
C     IF(ISTRAN(M).EQ.3) CALL CALFQC (ISTL,M,CON,CON1)
C     IF(M.EQ.4)THEN
C      IF(ISTOPT(4).EQ.2) CALL CALSED2(ISTL,0.5,CON1)
C      NTSWVD=ISWVSD
C      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C     ENDIF
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
C **  ADVECTIVE FLUX CALCULATION
C
      IF(ISTL.EQ.2) GOTO 300
      IF(ISCDCA(MVAR).EQ.0) GOTO 300
      IF(ISCDCA(MVAR).EQ.1) GOTO 400
      IF(ISCDCA(MVAR).EQ.2) GOTO 350
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
CX      IF(ISTL.EQ.2)THEN
C
CX      DO K=1,KC
CX      DO L=2,LA
CX      CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H1P(L)
CX      ENDDO
CX      ENDDO
C
CX      ELSE
C
CX      DO K=1,KC
CX      DO L=2,LA
CX      CONT(L,K)=CON1(L,K)+DELT*0.5*FQC(L,K)*DXYIP(L)/H2P(L)
CX      ENDDO
CX      ENDDO
C
CX      ENDIF
C
CX      WFQC=0.5
C
C----------------------------------------------------------------------C
C
CX      DO K=1,KC    
CX      DO L=2,LA
CX      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONT(L-1,K)
CX     &         +MIN(UHDY2(L,K),0.)*CONT(L,K)
CX      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONT(LSC(L),K)
CX     &         +MIN(VHDX2(L,K),0.)*CONT(L,K)
CX      ENDDO
CX      ENDDO
CX      DO K=1,KS
CX      DO L=2,LA
CX      FWU(L,K)=MAX(W2(L,K),0.)*CONT(L,K)
CX     &        +MIN(W2(L,K),0.)*CONT(L,K+1)
CX      ENDDO
CX      ENDDO
C
C1D      DO K=1,KC    
C1D      DO L=2,LA
C1D      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)
C1D     &         +MIN(UHDY2(L,K),0.)*CON1(L,K)
C1D      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)
C1D     &         +MIN(VHDX2(L,K),0.)*CON1(L,K)
C1D      ENDDO
C1D      ENDDO
C1D      DO K=1,KS
C1D      DO L=2,LA
C1D      FWU(L,K)=MAX(W2(L,K),0.)*CON1(L,K)
C1D     &        +MIN(W2(L,K),0.)*CON1(L,K+1)
C1D      ENDDO
C1D      ENDDO
C
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)
     &         +MIN(UHDY2(L,K),0.)*CON1(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)
     &         +MIN(VHDX2(L,K),0.)*CON1(L,K)
      ENDDO
      ENDDO
C
CDIA      WRITE(1,1111) IL(2),JL(2),K,AREAOLD(2),AREANEW(2),
CDIA     &             CON1(2,1),UHDY2(2,1),UHDY2(3,1),FQC(2,1)
CDIA     &             CON1(2,1),FUHU(2,1),FUHU(3,1),FQC(2,1)
C
      GOTO 500
C
 1111 FORMAT(3I5,6E12.4)
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN  (N-1) AND (N+1) AND ADVECTED FIELD AVERAGED 
C **  BETWEEN AT (N-1) AND (N) IF ISTL 3 ONLY
C
  350 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=0.5*(CON(L,K)+CON1(L,K))
     &         +DELT*0.5*FQC(L,K)*DXYIP(L)/H2P(L)
      ENDDO
      ENDDO
C
CX     WFQC=0.5
C
C----------------------------------------------------------------------C
C
      DO K=1,KC    
      DO L=2,LA
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*CONT(L-1,K)
     &         +MIN(UHDY2(L,K),0.)*CONT(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*CONT(LSC(L),K)
     &         +MIN(VHDX2(L,K),0.)*CONT(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*CONT(L,K)
     &        +MIN(W2(L,K),0.)*CONT(L,K+1)
      ENDDO
      ENDDO
C
      GOTO 500
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT 
C **  AVERAGED BETWEEN (N+1) AND (N-1) AND TRANSPORTED FIELD AT (N)
C
  400 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=CON1(L,K)
      ENDDO
      ENDDO
C
CX     WFQC=1.0
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)    
      FUHU(L,K)=0.5*UHDY2(L,K)*(CON(L,K)+CON(L-1,K))
      FVHU(L,K)=0.5*VHDX2(L,K)*(CON(L,K)+CON(LS,K))
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF(VHDX2(LN,K).LT.0.) FVHU(LN,K)=VHDX2(LN,K)*CON1(LN,K)
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(UHDY2(L+1,K).LT.0.) FUHU(L+1,K)=UHDY2(L+1,K)*CON1(L+1,K)
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(UHDY2(L,K).GT.0.) FUHU(L,K)=UHDY2(L,K)*CON1(L-1,K)
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS =LSC(L)
      IF(VHDX2(L,K).GT.0.) FVHU(L,K)=VHDX2(L,K)*CON1(LS,K)
      ENDDO
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=0.5*W2(L,K)*(CON(L,K+1)+CON(L,K))
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX FOR CENTRAL DIFF
C
      IF(ISHDMF.EQ.2) CALL CALDIFF (ISTL,M,CON1)
C     IF(ISTRAN(M).GE.1) CALL CALDIFF (ISTL,M,CON1)
C
C----------------------------------------------------------------------C
C
C **  ADD SECOND CALL TO EXPERIMENTAL SED CALCULATION FOR NO 
C **  ANTIDIFFUSIVE CORRECTION
C
C     IF(ISCDCA(MVAR).EQ.0)THEN
C     IF(M.EQ.4)THEN
C       IF(ISTRAN(4).LE.2) CALL CALSED2(ISTL,0.5,CON1)
C       IF(ISTRAN(4).EQ.3) CALL CALSED3(ISTL,0.5,CON1)
C     ENDIF
C     ENDIF
C
C**********************************************************************C
C
C **  STANDARD ADVECTION CALCULATION
C
  500 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE
C
C BEGIN IF ON TRANSPORT OPTION CHOICE
C
      IF(ISCDCA(MVAR).EQ.0)THEN
C
C BEGIN IF ON TIME LEVEL CHOICE FOR ISCDAC=0
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H1P(L)
CX     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*AREAOLD(L)
     &       +DELT*( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &                        +FVHU(L,K)-FVHU(LNC(L),K) )*DXYIP(L)
C     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
C     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
C     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
C
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
C  
C ELSE ON TIME LEVEL CHOICE FOR ISCDAC=0
C
      ELSE
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H2P(L)
CX     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*AREAOLD(L)
     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
C
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON(L,K)
       ENDDO
       ENDDO
      ENDIF
C
      ENDIF
C
C ENDIF ON TIME LEVEL CHOICE FOR ISCDAC=0
C
C1D      IF(ISUD.EQ.1)THEN
C1D      DO K=1,KC
C1D      DO L=2,LA
C1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
C1D      ENDDO
C1D      ENDDO
C1D      ENDIF
C
C1D      DO K=1,KC
C1D      DO L=2,LA
C1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
C1D      ENDDO
C1D      ENDDO
C
      DO K=1,KC
      DO L=2,LA
C      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CON(L,K)=SCB(L)*(CH(L,K)/AREANEW(L))+(1.-SCB(L))*CON(L,K)
CDIA      IF(L.EQ.2) WRITE(1,1111) IL(L),JL(L),K,CON(L,K),
CDIA     &             CON1(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
C
C **  ADD REMAINING SEDIMENT SETTLING AND FLUX
C
C     IF(M.EQ.4)THEN
C      IF(ISTOPT(4).EQ.2) CALL CALSED2(ISTL,0.5,CON )
C      NTSWVD=ISWVSD
C      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C     ENDIF
C
C **  IF ISACAC NE 0 DO NOT INCLUDE FQC MASS SOURCES IN UPDATE
C
C ELSE ON TRANSPORT OPTION CHOICE
C
      ELSE
C
C BEGIN IF ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H1P(L)
CX     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*AREAOLD(L)
     &       +DELT*( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K) )*DXYIP(L)
C1D     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
C1D     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
C1D     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
C
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       ENDDO
       ENDDO
      ENDIF
C
C ELSE ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
C
      ELSE
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H2P(L)
CX     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*AREAOLD(L)
     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
C
      IF(ISFCT(MVAR).GE.1)THEN
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON(L,K)
       ENDDO
       ENDDO
      ENDIF
C
      ENDIF
C
C ENDIF ON TIME LEVEL CHOICE FOR ISCDAC.NE.0
C
C1D      IF(ISUD.EQ.1)THEN
C1D      DO K=1,KC
C1D      DO L=2,LA
C1D      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
C1D      ENDDO
C1D      ENDDO
C1D      ENDIF
C
      DO K=1,KC
      DO L=2,LA
C      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      CON(L,K)=SCB(L)*(CH(L,K)/AREANEW(L))+(1.-SCB(L))*CON(L,K)
      CONT(L,K)=0.0
      ENDDO
      ENDDO
C  
      ENDIF
C
C ENDIF ON TRANSPORT OPTION CHOICE
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES
C
  700 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
C
      IF(VHDX2(LN,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)/AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)/AREANEW(L)
     &     +DELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))
     &     -FVHU(LN,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBSTMP=CBS(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBSTMP) CON(L,K)=CBSTMP
       CLOS(LL,K,M)=CON(L,K)
       NLOS(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
C      WRITE(6,6001)N,K,CBT
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
 6001 FORMAT('N,K,CBTS = ',2I10,F12.3)
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
C
      IF(UHDY2(L+1,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)
     &      -FUHU(L+1,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)
     &   +DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)/AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)/AREANEW(L)
     &     +DELT*(0.5*UHDY2(L+1,K)*(CON1(L,K)+CON(L,K))
     &     -FUHU(L+1,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBWTMP=CBW(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBWTMP) CON(L,K)=CBWTMP
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
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
      IF(UHDY2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
CDIA        WRITE(1,3331)N,UHDY2(L,K),CTMP,CON1(L,K),CON(L,K)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &   +0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)/AREANEW(L)+DELT*(FUHU(L,K)
     &   -0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBETMP=CBE(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBETMP) CON(L,K)=CBETMP
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
CDIA       WRITE(1,3332)N,NLOE(LL,K,M),CBETMP,CLOE(LL,K,M)
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOE(LL,K,M)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CON(L,K)=CBT
CDIA        WRITE(1,3333)N,NMNLO,NLOE(LL,K,M),CBT
       ELSE
        CON(L,K)=CLOE(LL,K,M)
     &         +(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
CDIA        WRITE(1,3334)N,NMNLO,NLOE(LL,K,M),CBT,CLOE(LL,K,M)
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
 3331 FORMAT(' W1 ',I5,4F12.4)
 3332 FORMAT(' W2 ',2I5,4F12.4)
 3333 FORMAT(' W3 ',2I5,4F12.4)
 3334 FORMAT(' W4 ',3I5,4F12.4)
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
C
      IF(VHDX2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)/AREANEW(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &   +0.5*(CON1(L,K)-CON(L,K))*AREAOLD(L)/AREANEW(L)+DELT*(FVHU(L,K)
     &   -0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)/AREANEW(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBNTMP=CBN(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBNTMP) CON(L,K)=CBNTMP
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
C      WRITE(6,6002)N,K,CBT
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
 6002 FORMAT('N,K,CBTN = ',2I10,F12.3)
C
CDIA      CLOSE(1)
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
      IF(ISADAC(MVAR).EQ.0) RETURN
      IF(ISCDCA(MVAR).EQ.1) RETURN
C
C**********************************************************************C
C
C **  STANDARD ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      UUU(L,K)=U2(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
      VVV(L,K)=V2(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      ENDDO
      ENDDO
C
      DO K=1,KS
      RDZIG=DZIG(K)
      DO L=2,LA
      WWW(L,K)=W2(L,K)*(CON(L,K+1)-CON(L,K))*HPI(L)*RDZIG
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
      AUHU=ABS(UHDY2(L,K))
      AVHV=ABS(VHDX2(L,K))
C1D      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
C1D      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
C     VTERM=VTERM*(CON(L,K)-CON(LS,K))
      UTERM=UTERM-DELTA4*UHDY2(L,K)*
     &      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
     &      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
      VTERM=VTERM-DELTA4*VHDX2(L,K)*
     &      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
     &      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
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
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(W2(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
      WTERM=WTERM-DELTA4*W2(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
     &      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
      WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
      FWU(L,K)=MAX(WW,0.)*CON(L,K)
     &        +MIN(WW,0.)*CON(L,K+1)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR SOURCE CELLS JMH20MARCH97
C
CLH
CLH   DO K=1,KC
CLH   DO L=2,LA
CLH     CONT(L,K)=1.
CLH     ABSQSUM=ABS(QSUM(L,K))
CLH     IF(ABSQSUM.GT.1.E-12) CONT(L,K)=0.
CLH   ENDDO
CLH   ENDDO
C    
CLH   DO K=1,KC
CLH   DO L=2,LA
CLH       TMPVAL=CONT(L,K)
CLH       LN=LNC(L)
CLH       FUHU(L  ,K)=TMPVAL*FUHU(L  ,K)
CLH       FUHU(L+1,K)=TMPVAL*FUHU(L+1,K)
CLH       FVHU(L  ,K)=TMPVAL*FVHU(L  ,K)
CLH       FVHU(LN ,K)=TMPVAL*FVHU(LN ,K)
CLH       FWU(L,K  )=TMPVAL*FWU(L,K  )
CLH       FWU(L,K-1)=TMPVAL*FWU(L,K-1)
CLH   ENDDO
CLH   ENDDO
C
      DO K=1,KC
CDIR& IVDEP
      DO L=2,LA
        IF(ABS(QSUM(L,K)).GT.1.E-12)THEN
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
C     DO K=1,KC
      FVHU(LN,K)=0.0
C     ENDDO
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
C     DO K=1,KC
      FUHU(L+1,K)=0.0
C     ENDDO
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
C     DO K=1,KC
      FUHU(L,K)=0.0
C     ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
C     DO K=1,KC
      FVHU(L,K)=0.0
C     ENDDO
      ENDDO
C
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
C
      IF(ISFCT(MVAR).EQ.0) GOTO 1100
C
C----------------------------------------------------------------------C
C
C **  DETERMINE MAX AND MIN CONCENTRATIONS
C
      DO K=1,KC
      DO L=2,LA
      CONTMX(L,K)=MAX(CON(L,K),CON2(L,K))
      CONTMN(L,K)=MIN(CON(L,K),CON2(L,K))
      ENDDO
      ENDDO
C
      DO L=2,LA
      CMAX(L,1)=MAX(CONTMX(L,1),CONTMX(L,2))
      CMAX(L,KC)=MAX(CONTMX(L,KS),CONTMX(L,KC))
      CMIN(L,1)=MIN(CONTMN(L,1),CONTMN(L,2))
      CMIN(L,KC)=MIN(CONTMN(L,KS),CONTMN(L,KC))
      ENDDO
C
      DO K=2,KS
      DO L=2,LA
      CMAXT=MAX(CONTMX(L,K-1),CONTMX(L,K+1))
      CMAX(L,K)=MAX(CONTMX(L,K),CMAXT)
      CMINT=MIN(CONTMN(L,K-1),CONTMN(L,K+1))
      CMIN(L,K)=MIN(CONTMN(L,K),CMINT)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB(L)*CONTMX(L-1,K) 
      CEMAX=SUB(L+1)*CONTMX(L+1,K)
      CSMAX=SVB(L)*CONTMX(LS,K)
      CNMAX=SVB(LN)*CONTMX(LN,K)
      CMAXT=MAX(CNMAX,CEMAX)
      CMAXT=MAX(CMAXT,CSMAX)
      CMAXT=MAX(CMAXT,CWMAX)
      CMAX(L,K)=MAX(CMAX(L,K),CMAXT)
      CWMIN=SUB(L)*CONTMN(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB(L+1)*CONTMN(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB(L)*CONTMN(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB(LN)*CONTMN(LN,K)+1.E+6*(1.-SVB(LN))
      CMINT=MIN(CNMIN,CEMIN)
      CMINT=MIN(CMINT,CSMIN)
      CMINT=MIN(CMINT,CWMIN)
      CMIN(L,K)=MIN(CMIN(L,K),CMINT)
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
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)
     &                               +FVHU(L,K)-FVHV(LN,K))
     &                      +RDZIC*(FWU(L,K-1)-FWV(L,K)) )*HPI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)
     &                               +FVHU(LN,K)-FVHV(L,K))
     &                      +RDZIC*(FWU(L,K)-FWV(L,K-1)) )*HPI(L)
      ENDDO
      ENDDO
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
C     DO K=1,KC
      DU(LN,K)=0.
      DV(LN,K)=0.
C     ENDDO
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
C     DO K=1,KC
      DU(L+1,K)=0.
      DV(L+1,K)=0.
C     ENDDO
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
C     DO K=1,KC
      DU(L-1,K)=0.
      DV(L-1,K)=0.
C     ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
C     DO K=1,KC
      DU(LS,K)=0.
      DV(LS,K)=0.
C     ENDDO
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
      IF(DU(L,K).GT.0.) DU(L,K)=(CMAX(L,K)-CON(L,K))/(DU(L,K)+BSMALL)
      DU(L,K)=MIN(DU(L,K),1.)
      IF(DV(L,K).GT.0.) DV(L,K)=(CON(L,K)-CMIN(L,K))/(DV(L,K)+BSMALL)
      DV(L,K)=MIN(DV(L,K),1.)
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
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTION CALCULATION
C
 1100 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  ADD SECOND CALL TO EXPERIMENTAL SED CALCULATION FOR  
C **  ANTIDIFFUSIVE CORRECTION
C !!  THIS IS A TEMPORARY TRIAL IF DOESN'T WORK MAKE CALL 
C !!  FOR ALL CASES AT PREVIOUS LOCATION
C !!  NOT ALSO THAT FQC HAS BEEN MOVED TO HERE FOR ISADAC=1
C !!  IF PROBLEMS ARISE PUT BACK IN AS BEFORE AT OLD LOCATION
C !!  RIGHT AFTER SECOND CALL TO CALSED2
C
C     IF(M.EQ.4) CALL CALSED2(ISTL,0.5,CON)
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CON(L,K)*HP(L)
CX     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON(L,K)*HP(L)
     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      ENDDO
      ENDDO
C
C     CH(L,K)=CONT(L,K)*HP(L)
C    &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
C    &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
C    &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
C
CRAY      DO K=1,KC
CRAY      DO L=2,LA
CRAY      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
CRAY      ENDDO
CRAY      ENDDO
C
C **  ADD REMAINING SEDIMENT SETTLING AND FLUX
C
C     IF(M.EQ.4)THEN
C      IF(ISTOPT(6).EQ.2) CALL CALSED2(ISTL,0.5)
C      NTSWVD=ISWVSD
C      IF(ISTOPT(4).EQ.13) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.14.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.15) CALL CALSED3(ISTL,0.5,SED)
C      IF(ISTOPT(4).EQ.16.AND.N.GE.NTSWVD) CALL CALSED3(ISTL,0.5,SED)
C     ENDIF
C
C**********************************************************************C
C
C **  DIAGNOSE FCT SCHEME
C
C----------------------------------------------------------------------C
C
      IF(ISFCT(MVAR).EQ.2)THEN
      WRITE(6,6110)N
C
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF(CCMAX.GT.0.)THEN
       WRITE(6,6111)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      ENDIF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF(CCMIN.GT.0.)THEN
       WRITE(6,6112)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C
 6110 FORMAT('  FCT DIAGNOSTICS AT N = ',I5)
 6111 FORMAT('  CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6112 FORMAT('  CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
C
C----------------------------------------------------------------------C
C
      RETURN
      END
