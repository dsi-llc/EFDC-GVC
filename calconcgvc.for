C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALCONCGVC (ISTL,IS2TL)
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
C 05/01/2002        John Hamrick       05/01/2002       John Hamrick
C  modified calls to calbal and budget subroutines
C  added calls to bal2t2, bal2t3
C----------------------------------------------------------------------C
C
C **  SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND 
C **  SUSPENDED CONSTITUTENTS, INCLUDING SALINITY, TEMPERATURE, DYE AND
C **  AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL 
C **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION EEB(LCM),CCLBTMP(LCM)
	DIMENSION TOXASM(NTXM),SEDASM(NSCM),SNDASM(NSNM)
C
CDHP  DIMENSION CUBTMP(LCM),CMBTMP(LCM),CLBTMP(LCM),EB(LCM),VTMP(LCM),
CDHP &          ABHPI(LCM,KSM)
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF(ISTL.EQ.2)THEN
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
        ELSE
          DELT=DTDYN
        END IF
        S3TL=0.0
        S2TL=1.0
      ENDIF
      DELTD2=DELT
C     DELTD2=0.5*DELT
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T3A
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION EXPLICIT HALF STEP CALCULATION
C
C----------------------------------------------------------------------C
C
C      IF(KC.EQ.1) GOTO 500
C
C      TTMP=SECNDS(0.0)
C
C      K=1
C      RCDZKK=-DELTD2*CDZKK(1)
C      DO L=2,LA
C      CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
C      CCMBTMP=1.-CCUBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCUBTMP*(SAL1(L,K+1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCUBTMP*(TEM1(L,K+1)
C      DU(L,K) =CCMBTMP*DYE1(L,K)+CCUBTMP*(DYE1(L,K+1)
C      DV(L,K) =CCMBTMP*SED1(L,K)+CCUBTMP*(SED1(L,K+1)
C      ENDDO
C
C      DO K=2,KS
C      RCDZKMK=-DELTD2*CDZKMK(K)
C      RCDZKK=-DELTD2*CDZKK(K)
C      DO L=2,LA
C      CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)
C      CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
C      CCMBTMP=1.-CCUBTMP-CCLBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCUBTMP*SAL1(L,K+1)
C     &                          +CCLBTMP*SAL1(L,K-1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCUBTMP*TEM1(L,K+1)
C     &                          +CCLBTMP*TEM1(L,K-1)
C      DU(L,K) =CCMBTMP*DYE1(L,K)+CCUBTMP*DYE1(L,K+1)
C     &                          +CCLBTMP*DYE1(L,K-1)
C      DV(L,K) =CCMBTMP*SED1(L,K)+CCUBTMP*SED1(L,K+1)
C     &                          +CCLBTMP*SED1(L,K-1)
C      ENDDO
C      ENDDO
C
C      K=KC
C      RCDZKMK=-DELTD2*CDZKMK(K)
C      DO L=2,LA
C      CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)
C      CCMBTMP=1.-CCLBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCLBTMP*SAL1(L,K-1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCLBTMP*TEM1(L,K-1)
C      DU(L,K)=CCMBTMP*DYE1(L,K)+CCLBTMP*DYE1(L,K-1)
C      DV(L,K)=CCMBTMP*SED1(L,K)+CCLBTMP*SED1(L,K-1)
C      ENDDO
C
C      TVDIF=TVDIF+SECNDS(TTMP)
C
C**********************************************************************C
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION-COSMIC INITIALIZATION 
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.EQ.0)THEN
      IF(ISCOSMIC.EQ.1)THEN
C
      DO K=1,KC
       RCOSMICX(1 ,K)=0.
       RCOSMICX(LC,K)=0.
       RCOSMICY(1 ,K)=0.
       RCOSMICY(LC,K)=0.
       RCOSMICZ(1 ,K)=0.
       RCOSMICZ(LC,K)=0.
       COSMICXP(1 ,K)=0.
       COSMICXP(LC,K)=0.
       COSMICYP(1 ,K)=0.
       COSMICYP(LC,K)=0.
       COSMICZP(1 ,K)=0.
       COSMICZP(LC,K)=0.
       COSMICXN(1 ,K)=0.
       COSMICXN(LC,K)=0.
       COSMICYN(1 ,K)=0.
       COSMICYN(LC,K)=0.
       COSMICZN(1 ,K)=0.
       COSMICZN(LC,K)=0.
      ENDDO
C
      DO L=1,LC
       COSMICZP(L,0)=0.
       COSMICZP(L,KC)=0.
       COSMICZN(L,0)=0.
       COSMICZN(L,KC)=0.
      ENDDO 
C
      DO K=1,KC
      DO L=2,LA
       RCOSMICX(L,K)=-1.
       TMP=U2(L,K)*U2(L+1,K)
       IF(TMP.LT.0.) RCOSMICX(L,K)=0.
       RCOSMICY(L,K)=-1.
       TMP=V2(L,K)*V2(LNC(L),K)
       IF(TMP.LT.0.) RCOSMICY(L,K)=0.
       RCOSMICZ(L,K)=-1.
       TMP=W2(L,K)*W2(L,K-1)
       IF(TMP.LT.0.) RCOSMICZ(L,K)=0.
      ENDDO 
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       COSMICXP(L,K)=DELT*DXIU(L)*U2(L,K)
       COSMICYP(L,K)=DELT*DYIV(L)*V2(L,K)
      ENDDO 
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
       COSMICXN(L,K)=MIN(COSMICXP(L,K),0.)
       COSMICYN(L,K)=MIN(COSMICYP(L,K),0.)
       COSMICXP(L,K)=MAX(COSMICXP(L,K),0.)
       COSMICYP(L,K)=MAX(COSMICYP(L,K),0.)
      ENDDO 
      ENDDO
C
      IF(KC.GE.2.AND.ISTL.EQ.3)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZP(L,K)=DELT*DZIG(K)*W2(L,K)/H1P(L)
        ENDDO 
        ENDDO
        DO K=1,KS
        DO L=2,LA
         COSMICZN(L,K)=MIN(COSMICZP(L,K),0.)
         COSMICZP(L,K)=MAX(COSMICZP(L,K),0.)
        ENDDO 
        ENDDO
      ENDIF
C
      IF(KC.GE.2.AND.ISTL.EQ.2)THEN
        DO K=1,KS
        DO L=2,LA
         COSMICZP(L,K)=2.*DELT*DZIG(K)*W2(L,K)/(HP(L)+H1P(L))
        ENDDO 
        ENDDO
        DO K=1,KS
        DO L=2,LA
         COSMICZN(L,K)=MIN(COSMICZP(L,K),0.)
         COSMICZP(L,K)=MAX(COSMICZP(L,K),0.)
        ENDDO 
        ENDDO
      ENDIF
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION 
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.EQ.0)THEN
c      IF(ISCOSMIC.EQ.0)THEN
C 
C **  PRESPECIFY THE UPWIND CELLS FOR 3D ADVECTION
C
      DO K=1,KC
      DO L=2,LA
        IF(LMASKDRY(L))THEN
	    IF(UHDY2(L,K).GE.0.0)THEN
            LUPU(L,K)=L-1	    
          ELSE
	      LUPU(L,K)=L
          END IF
	    IF(VHDX2(L,K).GE.0.0)THEN
            LUPV(L,K)=LSC(L)	    
          ELSE
	      LUPV(L,K)=L
          END IF
        END IF
	ENDDO
	ENDDO
C
      IF(KC.GT.1)THEN
      DO K=1,KS
      DO L=2,LA
        IF(LMASKDRY(L))THEN
	    IF( W2(L,K).GE.0.0 .OR. K.EQ.1 )THEN       ! *** DSI
            KUPW(L,K)=K   
          ELSE
	      KUPW(L,K)=K-1
          END IF
        END IF
	ENDDO
	ENDDO
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C                  SUBROUTINE CALTRAN (ISTL,IS2TL,MVAR,M,CON,CON1)
C 
      IF(ISTRAN(1).EQ.1.AND.ISCDCA(1).LT.4)
     &  CALL CALTRANGVC (ISTL,IS2TL,1,1,SAL,SAL1)
      IF(ISTRAN(2).EQ.1.AND.ISCDCA(2).LT.4)
     &  CALL CALTRANGVC (ISTL,IS2TL,2,2,TEM,TEM1)
      IF(ISTRAN(3).EQ.1.AND.ISCDCA(3).LT.4)
     &  CALL CALTRANGVC (ISTL,IS2TL,3,3,DYE,DYE1)
C     IF(ISTRAN(3).EQ.1) CALL CALTRANGVC (ISTL,IS2TL,4,4,SFL,SFL1)
      IF(ISTRAN(5).EQ.1.AND.ISCDCA(5).LT.4)THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         ENDDO
         ENDDO
C
C  SUBROUTINE CALTRAN (ISTL,IS2TL,MVAR,M,CON,CON1)
C 
         CALL CALTRANGVC (ISTL,IS2TL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(6).EQ.1.AND.ISCDCA(6).LT.4)THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         ENDDO
         ENDDO
C
C  SUBROUTINE CALTRAN (ISTL,IS2TL,MVAR,M,CON,CON1)
C 
         CALL CALTRANGVC (ISTL,IS2TL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(7).EQ.1.AND.ISCDCA(7).LT.4)THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         ENDDO
         ENDDO
C
C  SUBROUTINE CALTRAN (ISTL,IS2TL,MVAR,M,CON,CON1)
C 
         CALL CALTRANGVC (ISTL,IS2TL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
C     IF(ISTRAN(1).EQ.1) CALL CALTRAN (ISTL,1,SAL,UUU)
C     IF(ISTRAN(2).EQ.1) CALL CALTRAN (ISTL,2,TEM,VVV)
C     IF(ISTRAN(3).EQ.1) CALL CALTRAN (ISTL,3,DYE,DU)
C     IF(ISTRAN(4).EQ.1) CALL CALTRAN (ISTL,4,SED,DV)
C
C     IF(ISTRAN(1).EQ.2) CALL CALTRANI (ISTL,1,SAL,SAL1)
C     IF(ISTRAN(2).EQ.2) CALL CALTRANI (ISTL,2,TEM,TEM1)
C     IF(ISTRAN(3).EQ.2) CALL CALTRANI (ISTL,3,DYE,DYE1)
C     IF(ISTRAN(4).EQ.2) CALL CALTRANI (ISTL,4,SED,SED1)
C     IF(ISTRAN(1).EQ.2) CALL CALTRANI (ISTL,1,SAL,UUU)
C     IF(ISTRAN(2).EQ.2) CALL CALTRANI (ISTL,2,TEM,VVV)
C     IF(ISTRAN(3).EQ.2) CALL CALTRANI (ISTL,3,DYE,DU)
C     IF(ISTRAN(4).EQ.2) CALL CALTRANI (ISTL,4,SED,DV)
C
      IF(ISCRAY.EQ.0)THEN
        TSADV=TSADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSADV=TSADV+T2TMP-T1TMP
        WTSADV=WTSADV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
c      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  3D COSMIC ADVECTI0N TRANSPORT CALCULATION 
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.EQ.0)THEN
      IF(ISCOSMIC.EQ.1)THEN
C
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C 
      IF(ISTRAN(1).EQ.1.AND.ISCDCA(1).EQ.4)
     &  CALL COSTRANW (ISTL,IS2TL,1,1,SAL,SAL1)
      IF(ISTRAN(2).EQ.1.AND.ISCDCA(2).EQ.4)
     &  CALL COSTRANW (ISTL,IS2TL,2,2,TEM,TEM1)
      IF(ISTRAN(3).EQ.1.AND.ISCDCA(3).EQ.4)
     &  CALL COSTRANW (ISTL,IS2TL,3,3,DYE,DYE1)
c
      IF(ISTRAN(1).EQ.1.AND.ISCDCA(1).EQ.5)
     &  CALL COSTRAN (ISTL,IS2TL,1,1,SAL,SAL1)
      IF(ISTRAN(2).EQ.1.AND.ISCDCA(2).EQ.5)
     &  CALL COSTRAN (ISTL,IS2TL,2,2,TEM,TEM1)
      IF(ISTRAN(3).EQ.1.AND.ISCDCA(3).EQ.5)
     &  CALL COSTRAN (ISTL,IS2TL,3,3,DYE,DYE1)
c
C     IF(ISTRAN(3).EQ.1) CALL COSTRAN (ISTL,IS2TL,4,4,SFL,SFL1)
C
      IF(ISTRAN(5).EQ.1.AND.ISCDCA(5).EQ.4)THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         ENDDO
         ENDDO
         CALL COSTRANW (ISTL,IS2TL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(5).EQ.1.AND.ISCDCA(5).EQ.5)THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         ENDDO
         ENDDO
         CALL COSTRAN (ISTL,IS2TL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(6).EQ.1.AND.ISCDCA(6).EQ.4)THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         ENDDO
         ENDDO
         CALL COSTRANW (ISTL,IS2TL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(6).EQ.1.AND.ISCDCA(6).EQ.5)THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         ENDDO
         ENDDO
         CALL COSTRAN (ISTL,IS2TL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).EQ.1.AND.ISCDCA(7).EQ.4)THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         ENDDO
         ENDDO
         CALL COSTRANW (ISTL,IS2TL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(7).EQ.1.AND.ISCDCA(7).EQ.5)THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         ENDDO
         ENDDO
         CALL COSTRAN (ISTL,IS2TL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TSADV=TSADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSADV=TSADV+T2TMP-T1TMP
        WTSADV=WTSADV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  1D ADVECTI0N TRANSPORT CALCULATION 
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.GE.1)THEN
C
       IF(ISTL.EQ.3)THEN
         DO L=2,LA
          IF(LCT(L).EQ.6) AREAOLD(L)=FADYP2(L)
          IF(LCT(L).EQ.7) AREAOLD(L)=FADXP2(L)
          IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
          IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
         ENDDO
C       ELSE
         DO L=2,LA
          IF(LCT(L).EQ.6) AREAOLD(L)=FADYP1(L)
          IF(LCT(L).EQ.7) AREAOLD(L)=FADXP1(L)
          IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
          IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
         ENDDO
       ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C 
      IF(ISTRAN(1).EQ.1) CALL CALTRAN1D (ISTL,1,1,SAL,SAL1)
      IF(ISTRAN(2).EQ.1) CALL CALTRAN1D (ISTL,2,2,TEM,TEM1)
      IF(ISTRAN(3).EQ.1) CALL CALTRAN1D (ISTL,3,3,DYE,DYE1)
C     IF(ISTRAN(4).EQ.1) CALL CALTRAN1D (ISTL,4,4,SFL,SFL1)
      IF(ISTRAN(5).EQ.1)THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         ENDDO
         ENDDO
         CALL CALTRAN1D (ISTL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(6).EQ.1)THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         ENDDO
         ENDDO
         CALL CALTRAN1D (ISTL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
      IF(ISTRAN(7).EQ.1)THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         ENDDO
         ENDDO
         CALL CALTRAN1D (ISTL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         ENDDO
         ENDDO
        ENDDO
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TSADV=TSADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSADV=TSADV+T2TMP-T1TMP
        WTSADV=WTSADV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
C
      IF(ISTRAN(2).GE.1) CALL CALHEATGVC(ISTL)
C
C**********************************************************************C
C
C **  FULL IMPLICIT DYE AND TOXIC CONTAMINANT DECAY CALCULATION
C
      IF(ISTRAN(3).GE.1)THEN
C
      CDYETMP=1./(1.+DELT*RKDYE)
C     CDYETMP=(1.-DELTD2*RKDYE)/(1.+DELTD2*RKDYE)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         DYE(L,K)=CDYETMP*DYE(L,K)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C
C      IF(ISTRAN(5).GE.1)THEN
C      DO NT=1,NTOX
C      CDYETMP=1./(1.+DELT*RKTOXW(NT))
C     CDYETMP=(1.-DELTD2*RKTOXW(NT))/(1.+DELTD2*RKTOXW(NT))
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO K=1,KC
C        DO L=LF,LL
C         TOX(L,K,NT)=CDYETMP*TOX(L,K,NT)
C        ENDDO
C       ENDDO
C      ENDDO
C      ENDDO
C      ENDIF
C
C**********************************************************************C
C
C **  BOTTOM AND INTERNAL SEDIMENT AND TOXIC CONTAMINAT
C **  SOURCE-SINK CALCULATION
C
CJH5/13/77      IF(ISTRAN(6).GE.1) CALL CALSED(ISTL,1.0)
C
CJH5/13/77           IF(ISTRAN(7).GE.1)THEN
CJH5/13/77           NTSWVD=ISWVSD
CJH5/13/77            IF(ISTOPT(7).LE.1) CALL CALSED3(ISTL,1.0)
CJH5/13/77            IF(ISTOPT(7).EQ.2.AND.N.GE.NTSWVD)  CALL CALSED3(ISTL,1.0)
CJH5/13/77           ENDIF
C
CJH5/13/77           IF(ISTRAN(5).GE.1) CALL CALTOX(ISTL,1.0)
C
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C **  SEDIMENT AND TOXICS SETTLING,DEPOSITION,RESUSPENSION,ETC 
C **  FOR TWO TIME LEVEL SIMULATION
C 
      IF(IS2TIM.GE.1)THEN
C
CX       NTMP=MOD(N,2)
CX       IF(NTMP.EQ.0)THEN
         IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
	     ISEDDTC=ISEDDTC+1
	     IF(ISEDDTC.EQ.1)THEN
	       DTSED=DELT
	     ELSE
	       DTSED=DTSED+DELT
	     ENDIF
           IBALSTDT=0
c	     WRITE(8,888)N,ISEDDTC,ISEDDT,DTSED,DELT
	     IF(ISEDDTC.EQ.ISEDDT)THEN
	       IF(ISHOUSATONIC.EQ.0)THEN
               CALL SSEDTOX(ISTL,IS2TL,1.0)
             ELSE
               CALL SSEDTOXHOUS(ISTL,IS2TL,1.0)
             ENDIF
c	       WRITE(8,889)N,ISEDDTC,ISEDDT,DTSED
             IBALSTDT=1
             ISEDDTC=0
	     ENDIF
         ENDIF
CX       ENDIF
C
      ENDIF
C
C **  SEDIMENT AND TOXICS SETTLING,DEPOSITION,RESUSPENSION,ETC 
C **  FOR THREE TIME LEVEL SIMULATION
C 
      IF(IS2TIM.EQ.0)THEN
C
         IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
           NTSTBCM1=NTSTBC-1
	     IF(NCTBC.EQ.NTSTBCM1)THEN 
             IBALSTDT=0
             DTSED=FLOAT(NTSTBC)*DT
	       IF(ISHOUSATONIC.EQ.0)THEN
               CALL SSEDTOX(ISTL,IS2TL,1.0)
             ELSE
               CALL SSEDTOXHOUS(ISTL,IS2TL,1.0)
             ENDIF
             IBALSTDT=1
           ENDIF
         ENDIF
C
      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TSSED=TSSED+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSSED=TSSED+T2TMP-T1TMP
        WTSSED=WTSSED+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
  888 FORMAT('N,IC,I,DTS,DT = ',3I5,2F12.8)
  889 FORMAT('N,IC,I,DTS = ',3I5,F12.8,12X,'SSEDTOX CALLED')
C
C**********************************************************************C
C
C **  OPTIONAL MASS BALANCE CALCULATION
C
C----------------------------------------------------------------------C
C
C     IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
C
      IF(IS2TIM.EQ.0) THEN
      IF(ISTL.NE.2.AND.ISBAL.GE.1)THEN
         CALL CALBAL2
         CALL CALBAL3
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV2
           CALL CBALEV3
          ELSE
           CALL CBALOD2
           CALL CBALOD3
         ENDIF
      ENDIF
      ENDIF
C
C **  CALLS TO TWO-TIME LEVEL BALANCES
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          CALL BAL2T2
C          CALL BAL2T3(IBALSTDT)
          CALL BAL2T3B(IBALSTDT)
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  SEDIMENT BUDGET CALCULATION    (DLK 10/15)
C
C----------------------------------------------------------------------C
C
C      IF(NCTBC.NE.NTSTBC)THEN
      IF(IS2TIM.EQ.0) THEN
      IF(ISTL.NE.2.AND.ISSBAL.GE.1)THEN
         CALL BUDGET2
         CALL BUDGET3
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL BUDGEV2
C           CALL BUDGEV3
C          ELSE
C           CALL BUDGOD2
C           CALL BUDGOD3
C         ENDIF
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION
C
C----------------------------------------------------------------------C
C
      IF(KC.EQ.1) GOTO 1500
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
      ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      RCDZKK=-DELTD2*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*GVCSCLPI(L)*HPI(L)*SWB3D(L,1)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB(L)=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB(L)
       ENDDO
       IF(ISTRAN(1).GE.1)THEN
         DO L=LF,LL
          SAL(L,1)=SAL(L,1)*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(2).GE.1)THEN
         DO L=LF,LL
          TEM(L,1)=TEM(L,1)*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(3).GE.1)THEN
         DO L=LF,LL
          DYE(L,1)=DYE(L,1)*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(5).GE.1)THEN
         DO NT=1,NTOX      
          DO L=LF,LL
           TOX(L,1,NT)=TOX(L,1,NT)*EEB(L)
          ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(6).GE.1)THEN
         DO NS=1,NSED
          DO L=LF,LL
           SED(L,1,NS)=SED(L,1,NS)*EEB(L)
          ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
         DO L=LF,LL
          SND(L,1,NS)=SND(L,1,NS)*EEB(L)
         ENDDO
        ENDDO
       ENDIF
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELTD2*CDZKMK(K)
        RCDZKK=-DELTD2*CDZKK(K)
        DO L=LF,LL
          CCLBTMP(L)=RCDZKMK*GVCSCLPI(L)*HPI(L)*SWB3D(L,K-1)*AB(L,K-1)
          CCUBTMP=RCDZKK*GVCSCLPI(L)*HPI(L)*SWB3D(L,K)*AB(L,K)
          CCMBTMP=1.-CCLBTMP(L)-CCUBTMP
          IF( ABS(CCMBTMP-CCLBTMP(L)*CU1(L,K-1)) > 1.E-8 )THEN    ! DSI 2014-07  ADDED TO PREVENT SINGLE PRECISION DIVIDE BY ZERO
            EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
            CU1(L,K)=CCUBTMP*EEB(L)
          ELSE
            EEB(L)=1.
            CCLBTMP(L)=0.
            CU1(L,K)=0.
          ENDIF
        ENDDO
        IF(ISTRAN(1).GE.1)THEN
          DO L=LF,LL
           SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)
          ENDDO
        ENDIF
        IF(ISTRAN(2).GE.1)THEN
          DO L=LF,LL
           TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)
          ENDDO
        ENDIF
        IF(ISTRAN(3).GE.1)THEN
          DO L=LF,LL
           DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)
          ENDDO
        ENDIF
        IF(ISTRAN(5).GE.1)THEN
          DO NT=1,NTOX
           DO L=LF,LL
            TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
           ENDDO
          ENDDO
        ENDIF
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
           DO L=LF,LL
            SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)
           ENDDO
          ENDDO
        ENDIF
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
           DO L=LF,LL
            SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)
           ENDDO
          ENDDO
        ENDIF
       ENDDO
      ENDDO
C
      K=KC
      RCDZKMK=-DELTD2*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         CCLBTMP(L)=RCDZKMK*GVCSCLPI(L)*HPI(L)*SWB3D(L,K-1)*AB(L,K-1)
         CCMBTMP=1.-CCLBTMP(L)
         IF( ABS(CCMBTMP-CCLBTMP(L)*CU1(L,K-1)) > 1.E-8 )THEN   ! DSI 2014-07  ADDED TO PREVENT SINGLE PRECISION DIVIDE BY ZERO
           EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
         ELSE
           EEB(L)=1.
           CCLBTMP(L)=0.
         ENDIF
       ENDDO
       IF(ISTRAN(1).GE.1)THEN
         DO L=LF,LL
          SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(2).GE.1)THEN
         DO L=LF,LL
          TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(3).GE.1)THEN
         DO L=LF,LL
          DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)
         ENDDO
       ENDIF
       IF(ISTRAN(5).GE.1)THEN
         DO NT=1,NTOX
          DO L=LF,LL
           TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
          ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(6).GE.1)THEN
         DO NS=1,NSED
          DO L=LF,LL
           SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)
          ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(7).GE.1)THEN
         DO NS=1,NSND
          DO L=LF,LL
           SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)
          ENDDO
         ENDDO
       ENDIF
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        IF(ISTRAN(1).GE.1)THEN
          DO L=LF,LL
           SAL(L,K)=SAL(L,K)-CU1(L,K)*SAL(L,K+1)
          ENDDO
        ENDIF
        IF(ISTRAN(2).GE.1)THEN
          DO L=LF,LL
           TEM(L,K)=TEM(L,K)-CU1(L,K)*TEM(L,K+1)
          ENDDO
        ENDIF
        IF(ISTRAN(3).GE.1)THEN
          DO L=LF,LL
           DYE(L,K)=DYE(L,K)-CU1(L,K)*DYE(L,K+1)
          ENDDO
        ENDIF
        IF(ISTRAN(5).GE.1)THEN
          DO NT=1,NTOX
           DO L=LF,LL
            TOX(L,K,NT)=TOX(L,K,NT)-CU1(L,K)*TOX(L,K+1,NT)
           ENDDO
          ENDDO
        ENDIF
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
           DO L=LF,LL
            SED(L,K,NS)=SED(L,K,NS)-CU1(L,K)*SED(L,K+1,NS)
           ENDDO
          ENDDO
        ENDIF
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
           DO L=LF,LL
            SND(L,K,NS)=SND(L,K,NS)-CU1(L,K)*SND(L,K+1,NS)
           ENDDO
          ENDDO
        ENDIF
       ENDDO
      ENDDO
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=1,KB
      DO NS=1,NSED
       DO L=1,LC
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSND
      DO K=1,KB
       DO L=1,LC
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSED
       DO K=1,KC
        DO L=1,LC
         SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
C
      DO NS=1,NSND
       DO K=1,KC
        DO L=1,LC
         SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)
        ENDDO
       ENDDO
      ENDDO
C
CJH5/13/97      DO NT=1,NTOX
CJH5/13/97       DO K=1,KC
CJH5/13/97        DO L=1,LC
CJH5/13/97         TOXPFT(L,K,NT)=0.
CJH5/13/97        ENDDO
CJH5/13/97       ENDDO
CJH5/13/97      ENDDO
C
CJH5/13/97      DO NT=1,NTOX
CJH5/13/97       DO NS=1,NSED+NSND
CJH5/13/97        DO K=1,KC
CJH5/13/97         DO L=1,LC
CJH5/13/97          TOXPFT(L,K,NT)=TOXPFT(L,K,NT)+TOXPF(L,K,NS,NT)
CJH5/13/97         ENDDO
CJH5/13/97        ENDDO
CJH5/13/97       ENDDO
CJH5/13/97      ENDDO
C
      IF(ISCRAY.EQ.0)THEN
        TVDIF=TVDIF+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TVDIF=TVDIF+T2TMP-T1TMP
        WTVDIF=WTVDIF+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
 1500 CONTINUE
C
C**********************************************************************C
C
C **  DATA ASSIMILATION
C
      IF(NLCDA.GT.0)THEN
C
        SALASM=0.0
	  TEMASM=0.0
	  DYEASM=0.0
	  SFLASM=0.0
	  DO NT=1,NTOX
	   TOXASM(NT)=0.0
	  ENDDO
	  DO NS=1,NSED
	   SEDASM(NS)=0.0
	  ENDDO
	  DO NS=1,NSND
	   SNDASM(NS)=0.0
	  ENDDO
C
	  IWASM=0
C
        IF(N.EQ.1)THEN
	    OPEN(1,FILE='CDATASM.DIA')
          CLOSE(1,STATUS='DELETE')
	    OPEN(1,FILE='CDATASM.DIA')
	    IWASM=1
	    DO NLC=1,NLCDA
            DO NDAYA=1,NTC
	        FSALASM(NDAYA,NLC)=0.
              FVOLASM(NDAYA,NLC)=0.
              FTEMASM(NDAYA,NLC)=0.
	      ENDDO
	    ENDDO
	  ENDIF
C
	  NDAYA=MOD(N,NTSPTC)
	  NDAYA=1+(N-NDAYA)/NTSPTC
C	  WRITE(6,1212)N,NDAYA
C
        IF(N.EQ.NTSPTC)THEN
	    OPEN(1,FILE='CDATASM.DIA',POSITION='APPEND')
	    IWASM=1
	    WRITE(1,1212)N,NDAYA
	  ENDIF
C
        IF(ISCDA(1).GT.0)THEN
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=SAL(L,K)
        NSID=NCSERA(NLC,1)
        IF(IWASM.EQ.1) WRITE(1,1111)N,NLC,ICDA(NLC),JCDA(NLC),NSID,
     &                 CSERT(K,NSID,1),SAL(L,K)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,1).GT.0)THEN
	      FSALASM(NDAYA,NLC)=FSALASM(NDAYA,NLC)+TSCDA*DZC(K)*
     &              DXYP(L)*GVCSCLP(L)*HP(L)*(CSERT(K,NSID,1)-SAL(L,K))
c	      FVOLASM(NDAYA,NLC)=FVOLASM(NDAYA,NLC)+DZC(K)*DXYP(L)*
c     &                  HP(L)*(( SAL(L,K)/CSERT(K,NS,1) )-1.0)
	      FVOLASM(NDAYA,NLC)=FVOLASM(NDAYA,NLC)+TSCDA*DZC(K)*
     &      DXYP(L)*GVCSCLP(L)*HP(L)*(1.0-( CSERT(K,NSID,1)/SAL(L,K) ))
            SAL(L,K)=TSCDA*CSERT(K,NSID,1)+(1.-TSCDA)*SAL(L,K)
          ENDIF
          ENDIF
        ENDIF
        IF(IWASM.EQ.1) WRITE(1,1111)N,NLC,ICDA(NLC),JCDA(NLC),NSID,
     &                 CSERT(K,NSID,1),SAL(L,K)
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          SAL(L,K)=TSCDA*SAL(LDATA,K)+(1.-TSCDA)*SAL(L,K)
	  ENDIF
	  SALASM=SALASM+HP(L)*DXYP(L)*(SAL(L,K)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDIF
C
        IF(ISCDA(2).GT.0)THEN
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=TEM(L,K)
        NSID=NCSERA(NLC,2)
        IF(IWASM.EQ.1) WRITE(1,1112)N,NLC,ICDA(NLC),JCDA(NLC),NSID,
     &                 CSERT(K,NSID,2),TEM(L,K)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,2).GT.0)THEN
	      FTEMASM(NDAYA,NLC)=FTEMASM(NDAYA,NLC)+TSCDA*DZC(K)*
     &       DXYP(L)*GVCSCLP(L)*HP(L)*(CSERT(K,NSID,2)-TEM(L,K))
            TEM(L,K)=TSCDA*CSERT(K,NSID,2)+(1.-TSCDA)*TEM(L,K)
          ENDIF
          ENDIF
        ENDIF
        IF(IWASM.EQ.1) WRITE(1,1112)N,NLC,ICDA(NLC),JCDA(NLC),NSID,
     &                 CSERT(K,NSID,2),TEM(L,K)
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          TEM(L,K)=TSCDA*TEM(LDATA,K)+(1.-TSCDA)*TEM(L,K)
	  ENDIF
	  TEMASM=TEMASM+HP(L)*DXYP(L)*(TEM(L,K)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDIF
C
        IF(ISCDA(3).GT.0)THEN
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=DYE(L,K)
        NSID=NCSERA(NLC,3)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NS.GT.0)THEN
          IF(CSERT(K,NSID,3).GT.0)THEN
            DYE(L,K)=TSCDA*CSERT(K,NSID,3)+(1.-TSCDA)*DYE(L,K)
          ENDIF
          ENDIF
        ENDIF
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          DYE(L,K)=TSCDA*DYE(LDATA,K)+(1.-TSCDA)*DYE(L,K)
	  ENDIF
	  DYEASM=DYEASM+HP(L)*DXYP(L)*(DYE(L,K)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDIF
C
        IF(ISCDA(4).GT.0)THEN
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=SFL(L,K)
        NSID=NCSERA(NLC,4)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,4).GT.0)THEN
            SFL(L,K)=TSCDA*CSERT(K,NSID,4)+(1.-TSCDA)*SFL(L,K)
          ENDIF
          ENDIF
        ENDIF
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          SFL(L,K)=TSCDA*SFL(LDATA,K)+(1.-TSCDA)*SFL(L,K)
	  ENDIF
	  SFLASM=SFLASM+HP(L)*DXYP(L)*(SFL(L,K)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDIF
C
        IF(ISCDA(5).GT.0)THEN
        DO NT=1,NTOX
        M=MSVTOX(NT)
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=TOX(L,K,NT)
        NSID=NCSERA(NLC,M)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,M).GT.0)THEN
            TOX(L,K,NT)=TSCDA*CSERT(K,NSID,M)+(1.-TSCDA)*TOX(L,K,NT)
          ENDIF
          ENDIF
        ENDIF
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          TOX(L,K,NT)=TSCDA*TOX(LDATA,K,NT)+(1.-TSCDA)*TOX(L,K,NT)
	  ENDIF
	  TOXASM(NT)=TOXASM(NT)
     &            +HP(L)*DXYP(L)*(TOX(L,K,NT)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDDO
        ENDIF
C
        IF(ISCDA(6).GT.0)THEN
        DO NS=1,NSED
        M=MSVSED(NS)
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=SED(L,K,NS)
        NSID=NCSERA(NLC,M)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,M).GT.0)THEN
            SED(L,K,NS)=TSCDA*CSERT(K,NSID,M)+(1.-TSCDA)*SED(L,K,NS)
          ENDIF
          ENDIF
        ENDIF
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          SED(L,K,NS)=TSCDA*SED(LDATA,K,NS)+(1.-TSCDA)*SED(L,K,NS)
	  ENDIF
	  SEDASM(NS)=SEDASM(NS)
     &            +HP(L)*DXYP(L)*(SED(L,K,NS)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDDO
        ENDIF
C
 6222 FORMAT(' TC,SNEW,SASSM,SOLD='4F10.2)
C
        IF(ISCDA(7).GT.0)THEN
        DO NX=1,NSND
        M=MSVSND(NX)
        DO K=1,KC
        DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
	  CONASMOLD=SND(L,K,NX)
        NSID=NCSERA(NLC,M)
        IF(ITPCDA(NLC).EQ.0)THEN
          IF(NSID.GT.0)THEN
          IF(CSERT(K,NSID,M).GT.0)THEN
            SND(L,K,NX)=TSCDA*CSERT(K,NSID,M)+(1.-TSCDA)*SND(L,K,NX)
          ENDIF
          ENDIF
        ENDIF
        IF(ITPCDA(NLC).EQ.1)THEN
          LDATA=LIJ(ICCDA(NLC),JCCDA(NLC))
          SND(L,K,NX)=TSCDA*SND(LDATA,K,NX)+(1.-TSCDA)*SND(L,K,NX)
	  ENDIF
	  SNDASM(NX)=SNDASM(NX)
     &            +HP(L)*DXYP(L)*(SND(L,K,NX)-CONASMOLD)*DZC(K)
        ENDDO
        ENDDO
        ENDDO
        ENDIF
C
        IF(IWASM.EQ.1)THEN
          CLOSE(1)
	  ENDIF
C
      IF(IS2TIM.GE.1) THEN
        IF(ISBAL.GE.1)THEN
          SALOUT=SALOUT-SALASM
          DYEOUT=DYEOUT-DYEASM
          DO NT=1,NTOX
            TOXOUT2T(NT)=TOXOUT2T(NT)-TOXASM(NT)
          ENDDO
          DO NS=1,NSED
            SEDOUT2T(NS)=SEDOUT2T(NS)-SEDASM(NS)
          ENDDO
          DO NS=1,NSND
            SNDOUT2T(NS)=SNDOUT2T(NS)-SNDASM(NS)
          ENDDO
	  ENDIF
	ENDIF
C
      ENDIF
C
 1111 FORMAT(' SAL '5I5,2F10.3)	
 1112 FORMAT(' TEM '5I5,2F10.3)
 1212 FORMAT(' N,NDAYA = ',2I12)	
C
C**********************************************************************C
C
C **  SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
C
C     IF(ISTRAN(2).GE.1) CALL CALHEAT(ISTL)
C
C**********************************************************************C
C
C **  DYE DECAY CALCULATION
C
C     CDYETMP=1./(1.+DELTD2*RKDYE)
C     CDYETMP=(1.-DELTD2*RKDYE)/(1.+DELTD2*RKDYE)
C     DO K=1,KC
C     DO L=2,LA
C     DYE(L,K)=CDYETMP*DYE(L,K)
C     ENDDO
C     ENDDO
C
C**********************************************************************C
C
C **  BOTTOM AND INTERNAL SEDIMENT SOURCE-SINK CALCULATION
C
C     IF(ISTRAN(4).GE.1) CALL CALSED(ISTL,1.0)
C     IF(ISTRAN(4).EQ.5) CALL CALSED3(ISTL,1.0,SED)
C
C**********************************************************************C
C
      RETURN
      END
