C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTRAN (ISTL,IS2TL,MVART,M,CON,CON1)
C
C                                    MVAR=VARIALBE ID 1-8
C                                    M=ARRAY STORAGE ID FOR TIME SERIES
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 03/05/2002        john hamrick       03/05/2002       john hamrick
C  added transport bypass mask, imaskdry for dry cells
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
      ITRANFLOC=0
	IF(MVART.LE.8)MVAR=MVART
	IF(MVART.EQ.9)THEN
	  ITRANFLOC=1
	  MVAR=6
	ENDIF
C
      BSMALL=1.0E-6
CBUG    WFQC=0.0
C
      IF(ISDYNSTP.EQ.0)THEN
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
          ISUD=0
        ENDIF
      ELSE
        DELT=DTDYN
        DELTA=DTDYN
        DELTD2=0.5*DTDYN
        S3TL=0.0
        S2TL=1.0
        ISUD=0
      END IF
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
        FVHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       ENDDO
      ENDDO
C
      DO L=1,LC
        FWU(L,0)=0.
        FWU(L,KC)=0.
      ENDDO

C
      IF(IS2TL.EQ.1)THEN
        ISUD=1
        DO K=1,KC
        DO L=1,LC
          CON1(L,K)=CON(L,K)
        ENDDO
        ENDDO
      ENDIF
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
      CALL CALFQC (ISTL,IS2TL,MVAR,M,CON,CON1,FQCPAD,QSUMPAD,QSUMNAD)
C
      IF(ITRANFLOC.EQ.1)THEN
	  DO K=1,KC
	    DO L=2,LA
            FQC(L,K)=FQC(L,K)*FLOCDIA(L,K)
	    ENDDO
	  ENDDO
	ENDIF
C
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
      IF(IDRYTBP.EQ.0)THEN
C
      DO K=1,KC    
        DO L=2,LA
          FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)
     &         +MIN(UHDY2(L,K),0.)*CON1(L,K)
          FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)
     &         +MIN(VHDX2(L,K),0.)*CON1(L,K)
        ENDDO
      ENDDO
      DO K=1,KS
        DO L=2,LA
          FWU(L,K)=MAX(W2(L,K),0.)*CON1(L,K)
     &        +MIN(W2(L,K),0.)*CON1(L,K+1)
        ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KC    
        DO L=2,LA
C          IF(LMASKDRY(L))THEN
          IF(LMASKDRY(L))THEN
C            FUHU(L,K)=MAX(UHDY2(L,K),0.)*CON1(L-1,K)
C     &         +MIN(UHDY2(L,K),0.)*CON1(L,K)
C            FVHU(L,K)=MAX(VHDX2(L,K),0.)*CON1(LSC(L),K)
C     &         +MIN(VHDX2(L,K),0.)*CON1(L,K)
            FUHU(L,K)=UHDY2(L,K)*CON1(LUPU(L,K),K)
            FVHU(L,K)=VHDX2(L,K)*CON1(LUPV(L,K),K)
          ENDIF
        ENDDO
      ENDDO
	IF(KC.GT.1)THEN
      DO K=1,KS
        DO L=2,LA
          IF(LMASKDRY(L))THEN
C            FWU(L,K)=MAX(W2(L,K),0.)*CON1(L,K)
C     &        +MIN(W2(L,K),0.)*CON1(L,K+1)
            FWU(L,K)=W2(L,K)*CON1(L,KUPW(L,K))
          ENDIF
        ENDDO
      ENDDO
	ENDIF
C
      ENDIF
C
      GOTO 500
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
c
c  GVC DIAGNOSTIC
C
C      IF(N.LE.NTSPTC)THEN
C        DO K=1,KC
C        DO L=2,LA
C	  WRITE(8,811)N,IL(L),JL(L),K,MVAR,CON(L,K),FUHU(L,K),UHDY2(L,K),
C     &              FWU(L,K),W2(L,K)
C        ENDDO
C        ENDDO  
C	ENDIF
c
C  811 FORMAT('CK ADV SIG ',5I6,8E14.6)
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
      IF(IDRYTBP.EQ.0)THEN
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H1P(L)
CX     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
C1023      CH(L,K)=CON1(L,K)*H1P(L)
C1023     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
C1023     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
C1023     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)
     &                        +FUHU(L,K)-FUHU(L+1,K)
     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
C      TMPFAC=1.0/(1.0-0.5*RDZIC*DELT*HPI(L)*QSUMNAD(L,K)*DXYIP(L))
C      CH(L,K)=TMPFAC*CH(L,K)
      ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
	  IF(IMASKDRY(L).EQ.0)
     &    CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
	  IF(IMASKDRY(L).EQ.1)
     &    CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( FQC(L,K) )*DXYIP(L) )
	  IF(IMASKDRY(L).EQ.2)
     &    CH(L,K)=CON1(L,K)*H1P(L)
      ENDDO
      ENDDO
C
      ENDIF
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
c1023      CH(L,K)=CON1(L,K)*H2P(L)
c1023     &       +DELT*( ( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
c1023     &                        +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
c1023     &                        +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      CH(L,K)=CON1(L,K)*H2P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)
     &                +FUHU(L,K)-FUHU(L+1,K)
     &                +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &                +(FWU(L,K-1)-FWU(L,K))*RDZIC )
C      TMPFAC=1.0/(1.0-0.5*DELT*RDZIC*HPI(L)*QSUMNAD(L,K)*DXYIP(L))
C      CH(L,K)=TMPFAC*CH(L,K)
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
      IF(ISUD.EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
      ENDDO
      ENDDO
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
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
      IF(IDRYTBP.EQ.0)THEN
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
CX      CH(L,K)=CONT(L,K)*H1P(L)
CX     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      ENDDO
      ENDDO
C
      ELSE
C
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
      IF(IMASKDRY(L).EQ.0)
     &  CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      IF(IMASKDRY(L).EQ.1)
     &  CH(L,K)=CON1(L,K)*H1P(L)
     &       +DELT*( ( FQC(L,K) )*DXYIP(L) )
      IF(IMASKDRY(L).EQ.2)
     &  CH(L,K)=CON1(L,K)*H1P(L)
      ENDDO
      ENDDO
C
      ENDIF
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
      CH(L,K)=CON1(L,K)*H2P(L)
     &       +DELT*( ( RDZIC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
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
      IF(ISUD.EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
      ENDDO
      ENDDO
      ENDIF
C
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
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
      IF(ITRANFLOC.EQ.0)THEN
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
C
      IF(VHDX2(LN,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)
     &     +DELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))
     &     -FVHU(LN,K))*DXYIP(L)*HPI(L)
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
	ENDIF
C
 6001 FORMAT('N,K,CBTS = ',2I10,F12.3)
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.0)THEN
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
C
      IF(UHDY2(L+1,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)
     &      -FUHU(L+1,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)
     &   +DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)
     &     +DELT*(0.5*UHDY2(L+1,K)*(CON1(L,K)+CON(L,K))
     &     -FUHU(L+1,K))*DXYIP(L)*HPI(L)
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
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.0)THEN
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)      
C
      IF(UHDY2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &      +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)+DELT*(FUHU(L,K)
     &      -0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
       CBETMP=CBE(LL,1,M)+CSERT(1,NSID,M)
       IF(M.EQ.1.AND.CON(L,K).GT.CBETMP) CON(L,K)=CBETMP
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
       IF(ISUD.EQ.1) CON1(L,K)=CON(L,K)
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
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.0)THEN
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
C
      IF(VHDX2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &      +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)+DELT*(FVHU(L,K)
     &      -0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
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
      ENDIF
C
 6002 FORMAT('N,K,CBTN = ',2I10,F12.3)
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.1)THEN
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
C
      IF(VHDX2(LN,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(VHDX2(LN,K)*CON1(L,K)
     &      -FVHU(LN,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)
     &     +DELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))
     &     -FVHU(LN,K))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
      ENDIF
C
      ENDDO
      ENDDO
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.1)THEN
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
C
      IF(UHDY2(L+1,K).LT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(UHDY2(L+1,K)*CON1(L,K)
     &      -FUHU(L+1,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)
     &   +DELT*(UHDY2(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &     +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)
     &     +DELT*(0.5*UHDY2(L+1,K)*(CON1(L,K)+CON(L,K))
     &     -FUHU(L+1,K))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
      ENDIF
C
      ENDDO
      ENDDO
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.1)THEN
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)      
C
      IF(UHDY2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     &      -UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &      +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)+DELT*(FUHU(L,K)
     &      -0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
      ENDIF
C
      ENDDO
      ENDDO
	ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ITRANFLOC.EQ.1)THEN
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
C
      IF(VHDX2(L,K).GT.0.)THEN
       IF(ISTL.EQ.2)THEN
        CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
       ELSE
        IF(ISCDCA(MVAR).NE.2) CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     &      -VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)
        IF(ISCDCA(MVAR).EQ.2) CTMP=0.5*(CON1(L,K)+CON(L,K))
     &      +0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPI(L)+DELT*(FVHU(L,K)
     &      -0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPI(L)
        CON1(L,K)=CON(L,K)
       ENDIF
       CON(L,K)=CTMP
      ENDIF
C
      ENDDO
      ENDDO
      ENDIF
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
      IF(IDRYTBP.EQ.0)THEN
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
      RDZIC=DZIC(K)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
      AUHU=ABS(UHDY2(L,K))
      AVHV=ABS(VHDX2(L,K))
c      UTERM=AUHU*(1.-DELTA*AUHU*DXYIU(L)*HUI(L))*(CON(L,K)-CON(L-1,K))
c      VTERM=AVHV*(1.-DELTA*AVHV*DXYIV(L)*HVI(L))*(CON(L,K)-CON(LS,K))
      UTERM=AUHU*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(CON(L,K)-CON(LS,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORUE=DELTA*RDZIC*DXYIP(L  )*HPI(L  )*(FQCPAD(L  ,K)
     &                                 -QSUMPAD(L  ,K)*CON(L  ,K))
        SSCORUW=DELTA*RDZIC*DXYIP(L-1)*HPI(L-1)*(FQCPAD(L-1,K)
     &                                 -QSUMPAD(L-1,K)*CON(L-1,K))
        SSCORVN=DELTA*RDZIC*DXYIP(L  )*HPI(L  )*(FQCPAD(L  ,K)
     &                                 -QSUMPAD(L  ,K)*CON(L  ,K))
        SSCORVS=DELTA*RDZIC*DXYIP(LS )*HPI(LS )*(FQCPAD(LS ,K)
     &                                 -QSUMPAD(LS ,K)*CON(LS ,K))
        SSCORU=MAX(UHDY2(L,K),0.0)*SSCORUW+MIN(UHDY2(L,K),0.0)*SSCORUE
        SSCORV=MAX(VHDX2(L,K),0.0)*SSCORVS+MIN(VHDX2(L,K),0.0)*SSCORVN
	  UTERM=UTERM+SSCORU
	  VTERM=VTERM+SSCORV
	ENDIF
CADD1023
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
C     VTERM=VTERM*(CON(L,K)-CON(LS,K))
C
C      UTERM=UTERM-DELTA4*UHDY2(L,K)*
C     &      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
C     &      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
C      VTERM=VTERM-DELTA4*VHDX2(L,K)*
C     &      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
C     &      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
C
      IF(UHDY2(L,K).GE.0.0)THEN
        UTERM=UTERM-0.5*DELTA*UHDY2(L,K)*
     &      (VVV(LNW,K)+VVV(L-1,K)+WWW(L-1,K)+WWW(L-1,K-1)
     &      +UUU(L,K)+UUU(L-1,K))
      ELSE
        UTERM=UTERM-0.5*DELTA*UHDY2(L,K)*
     &      (VVV(LN,K)+VVV(L,K)+WWW(L,K)+WWW(L,K-1)
     &      +UUU(L,K)+UUU(L+1,K))
	ENDIF
C
      IF(VHDX2(L,K).GE.0.0)THEN
        VTERM=VTERM-0.5*DELTA*VHDX2(L,K)*
     &      (UUU(LS,K)+UUU(LSE,K)+WWW(LS,K)+WWW(LS,K-1)
     &      +VVV(LS,K)+VVV(L,K))
      ELSE
        VTERM=VTERM-0.5*DELTA*VHDX2(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+WWW(L,K)+WWW(L,K-1)
     &      +VVV(LN,K)+VVV(L,K))
	ENDIF
C
      IF(ISFCT(MVAR).GE.2)THEN
        FUHU(L,K)=0.5*UTERM
        FVHU(L,K)=0.5*VTERM
        IF(ISFCT(MVAR).EQ.3)THEN
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
     &           +MIN(UHU,0.)*CON(L,K)
        FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)
     &           +MIN(VHV,0.)*CON(L,K)
      ENDIF
      ENDDO
      ENDDO
C
      IF(KC.GT.1)THEN
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(W2(L,K))
c      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
      WTERM=AWW*(CON(L,K+1)-CON(L,K))
      IF(ISADAC(MVAR).GE.2)THEN
        SSCORWA=DELTA*DZIG(K+1)*HPI(L)*DXYIP(L)
     &         *(FQCPAD(L,K+1)-QSUMPAD(L,K+1)*CON(L,K+1))
        SSCORWB=DELTA*DZIG(K)*HPI(L)*DXYIP(L)
     &         *(FQCPAD(L,K  )-QSUMPAD(L,K  )*CON(L,K  ))
        SSCORW=MAX(W2(L,K),0.0)*SSCORWB+MIN(W2(L,K),0.0)*SSCORWA
	  WTERM=WTERM+SSCORW
      ENDIF
C
C      WTERM=WTERM-DELTA4*W2(L,K)*
C     &      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
C     &      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
C
      IF(W2(L,K).GE.0.0)THEN
        WTERM=WTERM-0.5*DELTA*W2(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+VVV(L,K)+VVV(LN,K)
     &      +WWW(L,K)+WWW(L,K-1))
	ELSE
        WTERM=WTERM-0.5*DELTA*W2(L,K)*
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
     &          +MIN(WW,0.)*CON(L,K+1)
      ENDIF
      ENDDO
      ENDDO
	ELSE
	DO L=2,LA
	FWU(L,KC)=0.0
	ENDDO
	ENDIF
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
C1023      DO K=1,KC
CDIR& IVDEP
C1023      DO L=2,LA
C1023        IF(ABS(QSUM(L,K)).GT.1.E-12)THEN
C1023          LN=LNC(L)
C1023          FUHU(L  ,K)=0.
C1023          FUHU(L+1,K)=0.
C1023          FVHU(L  ,K)=0.
C1023          FVHU(LN ,K)=0.
C1023          FWU(L,K  )=0.
C1023          FWU(L,K-1)=0.
C1023          CONT(L,K)=0.
C1023        ELSE
C1023          CONT(L,K)=1.
C1023        ENDIF
C1023      ENDDO
C1023      ENDDO
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
      DO L=1,LC
      CONTMX(L,K)=0.0
      CONTMN(L,K)=0.0
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      CONTMX(L,K)=MAX(CON(L,K),CON2(L,K))
      CONTMN(L,K)=MIN(CON(L,K),CON2(L,K))
      ENDDO
      ENDDO
C
      IF(KC.EQ.1)THEN
        DO L=2,LA
          CMAX(L,KC)=CONTMX(L,KC)
          CMIN(L,KC)=CONTMN(L,KC)
        ENDDO
	ELSE
        DO L=2,LA
          CMAX(L,1)=MAX(CONTMX(L,1),CONTMX(L,2))
          CMAX(L,KC)=MAX(CONTMX(L,KS),CONTMX(L,KC))
          CMIN(L,1)=MIN(CONTMN(L,1),CONTMN(L,2))
          CMIN(L,KC)=MIN(CONTMN(L,KS),CONTMN(L,KC))
        ENDDO
	  IF(KS.GE.2)THEN
          DO K=2,KS
          DO L=2,LA
            CMAXT=MAX(CONTMX(L,K-1),CONTMX(L,K+1))
            CMAX(L,K)=MAX(CONTMX(L,K),CMAXT)
            CMINT=MIN(CONTMN(L,K-1),CONTMN(L,K+1))
            CMIN(L,K)=MIN(CONTMN(L,K),CMINT)
          ENDDO
          ENDDO
	  ENDIF
	ENDIF
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
c  GVC DIAGNOSTIC
C
C      IF(N.LE.NTSPTC)THEN
C        DO K=1,KC
C        DO L=2,LA
C	  WRITE(8,801)N,IL(L),JL(L),K,DV(L-1,K),DU(L,K),FUHU(L,K),
C     &         DU(L-1,K),DV(L,K),FUHV(L,K)
C        ENDDO
C        ENDDO  
C	ENDIF
c
C  801 FORMAT('FLUX LIM 1 SIG ',4I6,8E14.6)
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
c
c  GVC DIAGNOSTIC
C
C      IF(N.LE.NTSPTC)THEN
C        DO K=1,KC
C        DO L=2,LA
C	  WRITE(8,801)N,IL(L),JL(L),K,CON(L,K),CH(L,K),
C     &                              FUHU(L,K),DU(L,K),DV(L,K)
C        ENDDO
C        ENDDO  
C	ENDIF
c
C  802 FORMAT('FLUX LIM 2 SIG ',4I6,8E14.6)
C
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
      ENDIF
C
C     ENDIF ON NON-DRYING AND WETTING ANTIDIFFUSION CALCULATION
C
C**********************************************************************C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION WITH DRY BYPASS
C
      IF(IDRYTBP.GT.0)THEN 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
        UUU(L,K)=0.
        VVV(L,K)=0.
        FUHU(L,K)=0.
        FVHU(L,K)=0.
        FUHV(L,K)=0.
        FVHV(L,K)=0.
	END DO
	END DO
C
      DO K=1,KS
      DO L=2,LA
        WWW(L,K)=0.
        FWU(L,K)=0.
        FWV(L,K)=0.
      ENDDO
      ENDDO

      DO K=1,KC
      DO L=2,LA
	IF(LMASKDRY(L))THEN
        LS=LSC(L)
        UUU(L,K)=U2(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
        VVV(L,K)=V2(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KS
      RDZIG=DZIG(K)
      DO L=2,LA
	IF(LMASKDRY(L))THEN
        WWW(L,K)=W2(L,K)*(CON(L,K+1)-CON(L,K))*HPI(L)*RDZIG
      ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
      AUHU=ABS(UHDY2(L,K))
      AVHV=ABS(VHDX2(L,K))
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
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
        LN=LNC(L)
        AWW=ABS(W2(L,K))
        WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
        WTERM=WTERM-DELTA4*W2(L,K)*
     &      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
     &      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
        WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
        FWU(L,K)=MAX(WW,0.)*CON(L,K)
     &        +MIN(WW,0.)*CON(L,K+1)
      ENDIF
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
      DO L=2,LA
	  IF(LMASKDRY(L))THEN
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
      IF(ISFCT(MVAR).EQ.0) GOTO 1101
C
C----------------------------------------------------------------------C
C
C **  DETERMINE MAX AND MIN CONCENTRATIONS
C
      DO K=1,KC
      DO L=2,LA
      CMIN(L,K)=0.
      CMAX(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      CONTMX(L,K)=MAX(CON(L,K),CON2(L,K))
      CONTMN(L,K)=MIN(CON(L,K),CON2(L,K))
      ENDIF
      ENDDO
      ENDDO
C
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      CMAX(L,1)=MAX(CONTMX(L,1),CONTMX(L,2))
      CMAX(L,KC)=MAX(CONTMX(L,KS),CONTMX(L,KC))
      CMIN(L,1)=MIN(CONTMN(L,1),CONTMN(L,2))
      CMIN(L,KC)=MIN(CONTMN(L,KS),CONTMN(L,KC))
      ENDIF
      ENDDO
C
      DO K=2,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      CMAXT=MAX(CONTMX(L,K-1),CONTMX(L,K+1))
      CMAX(L,K)=MAX(CONTMX(L,K),CMAXT)
      CMINT=MIN(CONTMN(L,K-1),CONTMN(L,K+1))
      CMIN(L,K)=MIN(CONTMN(L,K),CMINT)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
	IF(LMASKDRY(L))THEN
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
	IF(LMASKDRY(L))THEN
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      ENDIF
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS AND LOAD
C **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
C
      DO K=1,KC
      DO L=2,LA
      DU(L,K)=0.
      DV(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)
     &                               +FVHU(L,K)-FVHV(LN,K))
     &                      +RDZIC*(FWU(L,K-1)-FWV(L,K)) )*HPI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)
     &                               +FVHU(LN,K)-FVHV(L,K))
     &                      +RDZIC*(FWU(L,K)-FWV(L,K-1)) )*HPI(L)
      ENDIF
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
	IF(LMASKDRY(L))THEN
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
	IF(LMASKDRY(L))THEN
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)
     &         +MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)
     &         +MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
	IF(LMASKDRY(L))THEN
      FWU(L,K)=MIN(DV(L,K),DU(L,K+1))*FWU(L,K)
     &        +MIN(DU(L,K),DV(L,K+1))*FWV(L,K)
      ENDIF
      ENDDO
      ENDDO  
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTION CALCULATION
C
 1101 CONTINUE
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
      IF(LMASKDRY(L))THEN
CX      CH(L,K)=CON(L,K)*HP(L)
CX     &       +DELT*( (WFQC*FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON(L,K)*HP(L)
     &       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
     &               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     &               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      ENDIF
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
      ENDIF
C    
C     ENDIF ON DRY MASKED ANTIDIFFUSION CALCULATION
C
C**********************************************************************C
C
C **  DIAGNOSE FCT SCHEME
C
C----------------------------------------------------------------------C
C
      IF(ISFCT(MVAR).EQ.99)THEN
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
