C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE HDMT
C
C **  SUBROUTINE HDMT EXECUTES THE FULL HYDRODYNAMIC AND MASS TRANSPORT
C **  TIME INTERGATION 
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
C
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      ENDIF
C
      FOURDPI=4./PI
C
C**********************************************************************C
C
C **  INITIALIZE COURNT NUMBER DIAGNOSTICS
C
      DO K=1,KC
      DO L=2,LA
       CFLUUU(L,K)=0.
       CFLVVV(L,K)=0.
       CFLWWW(L,K)=0.
       CFLCAC(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
      ILOGC=0
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
C **  CALCULATE VELOCITY GRADIENTS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)      
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      HMC(L)=0.25*(HMP(L)+HMP(L-1)+HMP(LS)+HMP(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     &         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      U1V(L)=0.25*(H1P(LS)*(U1(LSE,1)+U1(LS,1))
     &          +H1P(L)*(U1(L+1,1)+U1(L,1)))*H1VI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     &         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      V1U(L)=0.25*(H1P(L-1)*(V1(LNW,1)+V1(L-1,1))
     &          +H1P(L)*(V1(LN,1)+V1(L,1)))*H1UI(L)           
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
      IF(ISWAVE.EQ.1) CALL WAVEBL
      IF(ISWAVE.EQ.2) CALL WAVESXY
C
C**********************************************************************C
C
C **  FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      ISTL=3
	IS2TL=0
      CALL CALTBXY(ISTL,IS2TL)
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
C
      IF(ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
C
C----------------------------------------------------------------------C
C
      N=-1
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX1(L)=(AVCON1*H1UI(L)+STBX(L)*SQRT(V1U(L)*V1U(L)
     &        +U1(L,1)*U1(L,1)))*U1(L,1)
       TBY1(L)=(AVCON1*H1VI(L)+STBY(L)*SQRT(U1V(L)*U1V(L)
     &        +V1(L,1)*V1(L,1)))*V1(L,1)
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
C
C      IF(ISVEG.GE.1.AND.KC.GT.1)
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C       TBX1(L)=TBX1(L)+0.5*FXVEG(L,1)*SQRT(V1U(L)*V1U(L)
C     &        +U1(L,1)*U1(L,1))*U1(L,1)
C       TBY1(L)=TBY1(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*U1V(L)
C     &        +V1(L,1)*V1(L,1))*V1(L,1)
C       ENDDO
C      ENDDO
C      ENDIF
C
C**********************************************************************C
C
C **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      CALL CALTBXY(ISTL,IS2TL)
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE STRESSES
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     &       +U(L,1)*U(L,1)))*U(L,1)
       TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     &       +V(L,1)*V(L,1)))*V(L,1)
       ENDDO
      ENDDO
C
C      IF(ISVEG.GE.1.AND.KC.GT.1)
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
C     &        +U(L,1)*U(L,1))*U(L,1)
C       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
C     &        +V(L,1)*V(L,1))*V(L,1)
C       ENDDO
C      ENDDO
C      ENDIF
C
      N=0
      CALL CALTSXY
C
C **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
C
      DO L=2,LA
        HDFUFX(L)=1.
        HDFUFY(L)=1.
        HDFUF(L)=1.
      ENDDO
C
C **  SET CORNER CORRECTIONS
C
       DO L=2,LA
         WCOREST(L)=1.
	   WCORWST(L)=1.
	   WCORNTH(L)=1.
	   WCORSTH(L)=1.
	 ENDDO


C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.EQ.0)THEN
C  
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ1(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX1(L))**2
     &                        +(TVAR3N(L)+TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2
     &                         +(TVAR3S(L)+TSY1(L))**2)
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
C
      ENDIF
C
C     ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.GE.1)THEN
C  
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX1(L))**2
     &                        +(TVAR3N(L)+TBY1(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U1(L+1,1)+U1(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V1(LN,1)+V1(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ1(L,0 )=CTURB2*SQRT(TAUB2)       
         QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2
     &                         +(TVAR3S(L)+TSY1(L))**2)
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)       
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
C
      ENDIF
C
C     ENDIF
C
C **   SET GRAIN STRESS
C
       DO L=2,LA
	   TAUBSED(L)=QQ(L,0 )/CTURB2
	   TAUBSND(L)=QQ(L,0 )/CTURB2
       ENDDO
C
C**********************************************************************C
C
C **   SET SWITCHES FOR THREE TIME LEVEL STEP
C
       ISTL=3
       IS2TL=0
       DELT=DT2
       DELTD2=DT
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
C
C**********************************************************************C
C**********************************************************************C
C 
C **  BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT
C **  CALCULATION
C
C **  SET CYCLE COUNTER AND CALL TIMER
C
      NTIMER=0
      N=0
      TIMEDAY=TCON*TBEGIN/86400.
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      CALL TIMELOG(N,TIMEDAY) 
C     CALL DTIME (TARRAY)
C     WRITE(9,200)N, TARRAY(1),TARRAY(2)
C     CALL FLUSH(9)
  200 FORMAT('  N=',I10,5X,'USER TIME=',E12.4,5X,'SYSTEM TIME=',E12.4)
      NTIMER=1
C
C----------------------------------------------------------------------C
C
      DO 1000 N=1,NTS
C
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
C
C      NLOGTMP=2*NLOGTMP
      IF(ILOGC.EQ.NTSMMT)THEN
        CLOSE(8,STATUS='DELETE')
        OPEN(8,FILE='EFDCLOG.OUT',STATUS='UNKNOWN')
        IF(ISDRY.GT.0)THEN
          OPEN(1,FILE='DRYWET.LOG',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        ENDIF
        IF(ISCFL.EQ.1)THEN
          OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        ENDIF
        ILOGC=0
      ENDIF
      ILOGC=ILOGC+1
C     WRITE(8,6626)ILOGC
 6626 FORMAT(' ILOGC = ',I10)
C
C
      IF(N.LE.NLTS) SNLT=0.
      IF(N.GT.NLTS.AND.N.LE.NTTS)THEN
         NTMP1=N-NLTS
         NTMP2=NTTS-NLTS+1
         SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
      ENDIF
      IF(N.GT.NTTS) SNLT=1.
C
      IF(N.LE.NTSVB)THEN
       GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
      ELSE
       GP=GPO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE VOLUME, MASS, MOMENTUM, AND ENERGY BALANCE
C
C
      IF(IS2TIM.EQ.0) THEN
C
      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
         CALL CALBAL1
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV1
          ELSE
           CALL CBALOD1
         ENDIF
       ENDIF
C
C  ** INITIALIZE SEDIMENT BUDGET CALCULATION   (DLK 10/15)
C
      IF(NCTBC.NE.NTSTBC.AND.ISSBAL.GE.1)THEN
         CALL BUDGET1
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0)THEN
C           CALL BUDGEV1
C          ELSE
C           CALL BUDGOD1
C         ENDIF
       ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  REENTER HERE FOR TWO TIME LEVEL CORRECTION
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN 
	    IF(ISTOPT(0).EQ.0)CALL CALAVBOLD (ISTL)
	    IF(ISTOPT(0).GE.1)CALL CALAVB (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALAVB2 (ISTL)
      ENDIF
      IF(ISCRAY.EQ.0)THEN
        TAVB=TAVB+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TAVB=TAVB+T2TMP-T1TMP
        WTAVB=WTAVB+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
      IF(ISTL.EQ.3)THEN
        IF(ISWAVE.EQ.1) CALL WAVEBL
        IF(ISWAVE.EQ.2) CALL WAVESXY
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C NOTES ON VARIOUS VERSIONS OF CALEXP
C
C  CALEXP   -  PRODUCTION VERSION WITH HORIZONTAL MOMENTUM SOURCE
C              AND 3D IMPLICIT VEGETATION DRAG
C                OPERATES WITH CALPUV2
C                OPERATES WITH CALPUV5
C                OPERATES WITH CALPUV9
C                OPERATES WITH CALPUVA
C                OPERATES WITH CALPUVB
C                OPERATES WITH CALPUVC
C                OPERATES WITH CALPUVD
C                OPERATES WITH CALPUVE
C  CALEXP1  -  OLD VERSION OF CALEXP WITH HORIZONTAL DOMAIN 
C              DECOMPOSITION FOR SINGLE LAYER, NOT NEEDED FOR LINKING
C  CALEXP2  -  ????????????????????????????
C                OPERATES WITH CALPUV6
C                OPERATES WITH CALPUV7
C                OPERATES WITH CALPUV8
C  CALEXP3  -  SMOLARKIEWCZ-MARGOLIN SCHEME
C                OPERATES WITH CALPUV3
C  CALEXP9  -  PROTOTYPE OF COSMIC ADVECTION SCHEME BUILT FROM CALEXP
C              TO BE REPLACE OF MODIFIED BY COSEXP
C                OPERATES WITH CALPUV2
C                OPERATES WITH CALPUV5
C                OPERATES WITH CALPUV9
C                OPERATES WITH CALPUVA
C                OPERATES WITH CALPUVB
C                OPERATES WITH CALPUVC
C                OPERATES WITH CALPUVD
C                OPERATES WITH CALPUVE
C  CALEXP1D -  TWO TIME LEVEL 1D VERSION CALLED FROM HDMT1D
C  CALEXP2T -  TWO TIME LEVEL 3D VERSION CALLED FROM HDMT2T
C
CX     IF(KC.EQ.1.AND.NDM.GE.2)THEN
CX        CALL CALEXP1 (ISTL)
CX       ELSE    
        IF(ISCDMA.LE.2) CALL CALEXP (ISTL)
c        IF(ISCDMA.EQ.3) CALL CALEXP3 (ISTL)
c        IF(ISCDMA.EQ.4) CALL CALEXP3 (ISTL)
        IF(ISCDMA.EQ.5) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.6) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.9) CALL CALEXP9 (ISTL)
CX      ENDIF
C
      IF(ISCRAY.EQ.0)THEN
        TCEXP=TCEXP+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCEXP=TCEXP+T2TMP-T1TMP
        WTCEXP=WTCEXP+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,
C **  VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
C
      CALL CALCSER (ISTL)
      CALL CALVEGSER (ISTL)
      CALL CALQVS (ISTL)
      PSERT(0)=0.
      IF(NPSER.GE.1) CALL CALPSER (ISTL)
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     SETUP CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,1)
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING (S&M SOLUTION ONLY)
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.GE.3.AND.ISCDMA.LE.8)THEN
      IF(ISTL.EQ.3)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
C
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
      ENDDO
C
      DO L=2,LA
        FXE(L)=FXE(L)+DT*SUB(L)*DYU(L)*(TSX(L)-TSX1(L))
        FYE(L)=FYE(L)+DT*SVB(L)*DXV(L)*(TSY(L)-TSY1(L))
      ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C NOTES ON VARIOUS VERSIONS OF CALPUV
C
C  CALPUV    -  RIGID LID SOLVER  
C  CALPUVA   -  EXPERIMENTAL
C  CALPUVB   -  EXPERIMENTAL
C  CALPUVC   -  EXPERIMENTAL
C  CALPUVD   -  EXPERIMENTAL
C  CALPUVE   -  EVERGLADES DRYING AND WETTING VERSION
C  CALPUVF   -  EXPERIMENTAL
C  CALPUV2   -  OLD NONDRYING VERSION
C  CALPUV3   -  SMOLARKIEWCZ-MARGOLIN SCHEME
C  CALPUV5   -  PRODUCTION VERSION FOR DRYING AND WETTING
C  CALPUV6   -  EXPERIMENTAL DRYING AND WETTING
C  CALPUV7   -  EXPERIMENTAL DRYING AND WETTING
C  CALPUV8   -  EXPERIMENTAL DRYING AND WETTING
C  CALPUV9   -  PRODUCTION VERSION
C  CALPUV1D  -  TWO TIME LEVEL 1D VERSION CALLED FROM HDMT1D
C  CALPUV2T  -  TWO TIME LEVEL 3D VERSION CALLED FROM HDMT2T
C
CJH      IF(ISDRY.EQ.-1)THEN
CJH         CALL CALPUV(ISTL)
CJH         GOTO 5555
CJH      ENDIF
C
CJH      IF(ISEVER.GE.1)THEN
CJH        CALL CALPUVE(ISTL)
CJH        GOTO 5555
CJH      ENDIF
C
C NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATE WITH CALEXP OR CALEXP9
C
CJH      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV2(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.2) CALL CALPUV5(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV2(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.4) CALL CALPUV5(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.9)THEN
           IF(ISCHAN.EQ.0.AND.ISDRY.EQ.0) CALL CALPUV9(ISTL)
           IF(ISCHAN.GE.1.OR.ISDRY.GE.1) CALL CALPUV9C(ISTL)
CJH        ENDIF
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.99) CALL CALPUV9(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.10) CALL CALPUVA(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.11) CALL CALPUVB(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.12) CALL CALPUVC(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.13) CALL CALPUVD(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.14) CALL CALPUVF(ISTL)
CJH        IF(ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV5(ISTL)
CJH        IF(ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV5(ISTL)
CJH        IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV5(ISTL)
CJH        IF(ISDRY.EQ.99) CALL CALPUV5(ISTL)
CJH        GOTO 5555
CJH      ENDIF
C
C NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATES WITH CALEXP3
C
CJH      IF(ISCDMA.GT.2.AND.ISCDMA.LT.5)THEN 
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV3(ISTL)
CJH        IF(ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV3(ISTL)
CJH        GOTO 5555
CJH      ENDIF
C
C NOTE THE FOLLOWING BLOCK OF OPTIONS OPERATES WITH CALEXP2
C
CJH      IF(ISCDMA.GE.5.AND.ISCDMA.LE.8)THEN
CJH        IF(ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV6(ISTL)
CJH        IF(ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV6(ISTL)
CJH        IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL)
CJH        IF(ISDRY.EQ.99) CALL CALPUV7(ISTL)
CX      IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL) !7 MOVED TO 8
C       IF(ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV8(ISTL)
CJH      ENDIF
C
 5555 CONTINUE
C
      IF(ISCRAY.EQ.0)THEN
        TPUV=TPUV+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TPUV=TPUV+T2TMP-T1TMP
        WTPUV=WTPUV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  WRITE DIAGNOSTICS
C
C----------------------------------------------------------------------C
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      IF(ISLOG.GE.1)THEN
      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN
C     CALL FLUSH(8)
      ENDIF
C
C  17 FORMAT('  N,ITER,AVMA,AVMI,ABMA,ABMI',2I5,4(1X,F8.4))
   17 FORMAT('  N,ITER,RSQ,CFMAX,AVMA,AVMI,ABMA,ABMI',
     &        I7,I5,2E12.4,4(1X,F8.4))
C
      ERRMAX=MAX(ERRMAX,ERR)
      ERRMIN=MIN(ERRMIN,ERR)
      ITRMAX=MAX(ITRMAX,ITER)
      IRRMIN=MIN(ITRMIN,ITER)
C
C**********************************************************************C
C
C **  ADVANCE INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO K=1,KC
      DO L=2,LA
      UHDY2(L,K)=UHDY1(L,K)
      UHDY1(L,K)=UHDY(L,K)
      VHDX2(L,K)=VHDX1(L,K)
      VHDX1(L,K)=VHDX(L,K)
      U2(L,K)=U1(L,K)
      V2(L,K)=V1(L,K)
      U1(L,K)=U(L,K)
      V1(L,K)=V(L,K)
      W2(L,K)=W1(L,K)
      W1(L,K)=W(L,K)
      ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9)THEN
      IF(ISTL.EQ.3)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       ENDDO
      ENDDO
C
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
      ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
C
C----------------------------------------------------------------------C
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        CALL CALUVW (ISTL,IS2TL)
       ELSE
        DO ND=1,NDM
         LF=2+(ND-1)*LDM
         LL=LF+LDM-1
         DO L=LF,LL
          UHDY(L,1)=UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDX(L,1)=VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
         ENDDO
        ENDDO
        CALL CALUVW (ISTL,IS2TL)
      ENDIF
      IF(ISCRAY.EQ.0)THEN
        TUVW=TUVW+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TUVW=TUVW+T2TMP-T1TMP
        WTUVW=WTUVW+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     DIAGNOSTIC CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,2)
C
C**********************************************************************C
C
C **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS 
C **  AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      CALL CALCONC (ISTL,IS2TL)
C
C----------------------------------------------------------------------C
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      ENDDO
      ENDDO
C
      DO NS=1,NSED
       DO K=1,KB
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
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
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
C----------------------------------------------------------------------C
C
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
C
      IF(ISMMC.EQ.1)THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SAL(L,K).GT.SALMAX)THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SAL(L,K).LT.SALMIN)THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(DYE(L,K).GT.SALMAX)THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(DYE(L,K).LT.SALMIN)THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SFL(L,K).GT.SALMAX)THEN
       SALMAX=SFL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SFL(L,K).LT.SALMIN)THEN
       SALMIN=SFL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN
C
      ENDIF
C
C
      IF(ISMMC.EQ.2)THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(TEM(L,K).GT.SALMAX)THEN
       SALMAX=TEM(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(TEM(L,K).LT.SALMIN)THEN
       SALMIN=TEM(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN
C
      ENDIF
C
 6001 FORMAT('  N=',I10)
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
 6006 FORMAT('  SFLMAX=',F14.4,5X,'I,J,K=',(3I10))
 6007 FORMAT('  SFLMIN=',F14.4,5X,'I,J,K=',(3I10))
 6008 FORMAT('  TEMMAX=',F14.4,5X,'I,J,K=',(3I10))
 6009 FORMAT('  TEMMIN=',F14.4,5X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
C **  CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT 
C **  CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOULBE TIME
C **  STEP TRANSPORT FIELD
C
C----------------------------------------------------------------------C
C
      ITMP=0
      IF(ISTRAN(4).GE.1) ITMP=1
      IF(ISTRAN(8).GE.1) ITMP=1
      IF(ISICM.GE.1) ITMP=1
C
      IF(ITMP.EQ.1)THEN
      NTMP=MOD(N,2)
      IF(NTMP.EQ.0.AND.ISTL.EQ.3)THEN
C
C **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        UHDY2E(L)=0.
        VHDX2E(L)=0.
        HWQ(L)=0.25*(H2P(L)+2.*H1P(L)+HP(L))
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        UHDYWQ(L,K)=0.5*(UHDY1(L,K)+UHDY2(L,K))
        VHDXWQ(L,K)=0.5*(VHDX1(L,K)+VHDX2(L,K))
        UHDY2E(L)=UHDY2E(L)+UHDYWQ(L,K)*DZC(K)
        VHDX2E(L)=VHDX2E(L)+VHDXWQ(L,K)*DZC(K)
        UWQ(L,K)=0.5*(U1(L,K)+U2(L,K))
        VWQ(L,K)=0.5*(V1(L,K)+V2(L,K))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3E(L)=UHDY2E(L+1)
       TVAR3N(L)=VHDX2E(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
        TVAR2E(L,K)=UHDYWQ(L+1   ,K)
        TVAR2N(L,K)=VHDXWQ(LNC(L),K)
        ENDDO
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
        WWQ(L,K)=SWB(L)*(WWQ(L,K-1)
     &       -DZC(K)*(TVAR2E(L,K)-UHDYWQ(L,K)-TVAR3E(L)+UHDY2E(L)
     &       +TVAR2N(L,K)-VHDXWQ(L,K)-TVAR3N(L)+VHDX2E(L))*DXYIP(L))
     &        +SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
        ENDDO
       ENDDO       
      ENDDO       
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        HPPTMP=H2WQ(L)+DT2*DXYIP(L)*(QSUME(L)
     &      -(TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)))
        HWQ(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HWQ(L)
       ENDDO
      ENDDO
C 
C     ADD CHANNEL INTERACTIONS
C

      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        IF(MDCHTYP(NMD).EQ.1)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     &    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
        ENDIF            
        IF(MDCHTYP(NMD).EQ.2)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     &    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF            
        IF(MDCHTYP(NMD).EQ.3)THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
     &    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     &    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     &    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        ENDIF
        ENDDO
      ENDIF
C
C     END ADD CHANNEL INTERACTIONS
C
C      IF(ISTRAN(6).GE.1) CALL CALWQC(2)
      IF(ISTRAN(8).GE.1) CALL WQ3D
      IF(ISTRAN(4).GE.1) CALL CALSFT(2)
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        H2WQ(L)=HWQ(L)
       ENDDO
      ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING
C **  AN EQUATION OF STATE 
C
      IF(ISTL.EQ.3)THEN
       DO K=1,KC
       DO L=2,LA
       B1(L,K)=B(L,K)
       ENDDO
       ENDDO
      ENDIF
C
      IF(BSC.GT.1.E-6)THEN
        CALL CALBUOY
	ELSE
        DO K=1,KC
          DO L=2,LA
            B(L,K)=0.
          ENDDO
        ENDDO
	ENDIF
C
      IF(NCTBC.NE.NTSTBC.AND.ISBAL.GE.1)THEN
         CALL CALBAL4
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV4
          ELSE
           CALL CBALOD4
         ENDIF
      ENDIF
C
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)      
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     &         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     &         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES 
C **  AT TIME LEVEL (N)
C
      IF(ISTL.NE.2.AND.ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
C **  INTENSITIES 
C
C----------------------------------------------------------------------C
C
C      IF(ISTL.EQ.2)THEN
C
C      DO K=1,KC
C       DO L=2,LA
C        QQ(L,K)=SQRT(QQ(L,K)*QQ1(L,K))
C       ENDDO
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
      IF(ISCDMA.EQ.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
      ENDDO
C
      ELSE
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ1(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ1(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       ENDDO
      ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)
C
      CALL CALTBXY(ISTL,IS2TL)
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     &       +U(L,1)*U(L,1)))*U(L,1)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     &       +V(L,1)*V(L,1)))*V(L,1)
       ENDDO
      ENDDO
C
C      IF(ISVEG.GE.1.AND.KC.GT.1)
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
C     &        +U(L,1)*U(L,1))*U(L,1)
C       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
C     &        +V(L,1)*V(L,1))*V(L,1)
C       ENDDO
C      ENDDO
C      ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.EQ.0)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
C
      ENDIF  
C
C     ENDIF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
C     IF(KC.GT.1.OR.ISTRAN(4).GE.1)THEN
C
      IF(ISWAVE.GE.1)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     &                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)       
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     &                         +(TVAR3S(L)+TSY(L))**2)
       ENDDO
      ENDDO
C
      ENDIF 
C
C     ENDIF
C
C**********************************************************************C
C
C **  CALCULATE TURBULENT INTENSITY SQUARED 
C
      IF(ISCRAY.EQ.0)THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
      IF(KC.GT.1)THEN
        IF(ISQQ.EQ.1)THEN 
	    IF(ISTOPT(0).EQ.0)CALL CALQQ1OLD (ISTL)
	    IF(ISTOPT(0).GE.1)CALL CALQQ1 (ISTL)
        ENDIF
        IF(ISQQ.EQ.2) CALL CALQQ2 (ISTL)
      ENDIF
      IF(ISCRAY.EQ.0)THEN
        TQQQ=TQQQ+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQQQ=TQQQ+T2TMP-T1TMP
        WTQQQ=WTQQQ+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE MEAN MASS TRANSPORT FIELD
C
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.GE.1)THEN
          NTMP=MOD(N,2)
          IF(ISTL.EQ.3.AND.NTMP.EQ.0) CALL CALMMT
        ENDIF
      ENDIF
C
C      IF(ISSSMMT.NE.2) CALL CALMMT
C
C**********************************************************************C
C
C **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
C **  IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
C
C----------------------------------------------------------------------C
C
      IF(NCTBC.EQ.NTSTBC)THEN
       NCTBC=0
       ISTL=2
       DELT=DT
       DELTD2=0.5*DT
       DZDDELT=DZ/DELT
       ROLD=0.5
       RNEW=0.5
       GOTO 500
      ELSE
       NCTBC=NCTBC+1
       ISTL=3
       DELT=DT2
       DELTD2=DT
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO TIME SERIES FILES
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON  
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      IF(ISTMSR.GE.1)THEN
C        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCTMSR.GE.NWTMSR)THEN
            CALL TMSR
            ICALLTP=1
            NCTMSR=1  
           ELSE
            NCTMSR=NCTMSR+1
          ENDIF
C        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO DUMP FILES
C
      IF(ISDUMP.GE.1)THEN
        IF(TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP)THEN
          IF(NCDUMP.GE.NSDUMP)THEN
            CALL DUMP
            ICALLTP=1
            NCDUMP=1  
           ELSE
            NCDUMP=NCDUMP+1
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C ** WRITE BOTTOM VELOCITIES TO FILE FOR AJ MINE SEDIMENT BED MODEL
C     (DLK - ADDED 3/25/97)
      IF(IHYDOUT.GE.1)THEN
        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN
          IF(NCHYDOUT.EQ.NWHYDOUT)THEN
            CALL HYDOUT
            NCHYDOUT=1  
           ELSE
            NCHYDOUT=NCHYDOUT+1
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  OUTPUT ZERO DIMENSION VOLUME BALANCE
C 
C----------------------------------------------------------------------C
C
C      IF(ISDRY.GE.1)THEN
C        OPEN(1,FILE='ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
C        DO LS=1,LORMAX
C        IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN
C           WTM=VOLSEL(LS+1)-VOLZERD
C           WTMP=VOLZERD-VOLSEL(LS)
C           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
C           WTM=WTM/DELVOL
C           WTMP=WTMP/DELVOL
C           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
C           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
C        ENDIF
C        ENDDO
C        IF(ISDYNSTP.EQ.0)THEN
C          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
C        ELSE
C          TIME=TIMESEC/TCTMSR
C        ENDIF
C        WRITE(1,5304)TIME,SELZERD,ASFZERD,VOLZERD,VETZERD
C        CLOSE(1)
C      ENDIF
      ICALLTP=0       
C
 5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C                         
C**********************************************************************C
C
C **  WRITE VERTICAL SCALAR FIELD PROFILES 
C
      IF(ISVSFP.EQ.1)THEN
        IF(N.GE.NBVSFP.AND.N.LE.NSVSFP)THEN
          CALL VSFP
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE MEAN MASS TRANSPORT FIELD
C
      IF(ISSSMMT.NE.2)THEN
        IF(ISICM.EQ.0) CALL CALMMT
      ENDIF
C
C      IF(ISSSMMT.NE.2) CALL CALMMT
C
C**********************************************************************C
C
C **  ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
C
      IF(ISPD.EQ.1)THEN
        IF(N.GE.NPDRT) CALL DRIFTER
      ENDIF
      IF(ISLRPD.GE.1)THEN
        IF(ISCRAY.EQ.0)THEN
          T1TMP=SECNDS(0.0)
         ELSE
          T1TMP=SECOND( )
          CALL TIMEF(WT1TMP)
        ENDIF
        IF(ISLRPD.LE.2)THEN
          IF(N.GE.NLRPDRT(1)) CALL LAGRES
        ENDIF
        IF(ISLRPD.GE.3)THEN
          IF(N.GE.NLRPDRT(1)) CALL GLMRES
        ENDIF
        IF(ISCRAY.EQ.0)THEN
          TLRPD=TLRPD+SECNDS(T1TMP)
         ELSE
          T2TMP=SECOND( )
          CALL TIMEF(WT2TMP)
          TLRPD=TLRPD+T2TMP-T1TMP
          WTLRPD=WTLRPD+(WT2TMP-WT1TMP)*0.001
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES
C
      IF(ISBAL.GE.1)THEN
         CALL CALBAL5
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0)THEN
           CALL CBALEV5
          ELSE
           CALL CBALOD5
         ENDIF
       ENDIF
C
C   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
C
       IF(ISSBAL.GE.1)THEN
       CALL BUDGET5
       ENDIF
C       NTMP=MOD(N,2)
C       IF(NTMP.EQ.0)THEN
C         CALL BUDGEV5
C        ELSE
C         CALL BUDGOD5
C       ENDIF
C
C**********************************************************************C
C
C **  PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
C
CHARM      IF(ISHTA.EQ.1) CALL CALHTA
C
C**********************************************************************C
C
C **  CALCULATE DISPERSION COEFFICIENTS
C
C     IF(N.GE.NDISP)THEN
      IF(N.GE.NDISP.AND.NCTBC.EQ.1)THEN
       IF(ISDISP.EQ.2) CALL CALDISP2
       IF(ISDISP.EQ.3) CALL CALDISP3
      ENDIF
C
C**********************************************************************C
C
C **  PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
C
      IF(ISLSHA.EQ.1.AND.N.EQ.NCLSHA)THEN
       CALL LSQHARM 
       NCLSHA=NCLSHA+(NTSPTC/24)
      ENDIF
C
C**********************************************************************C
C
C **  PRINT INTERMEDIATE RESULTS 
C
C----------------------------------------------------------------------C
C
      IF(NPRINT .EQ. NTSPP)THEN
       NPRINT=1
       CALL OUTPUT1
      ELSE
       NPRINT=NPRINT+1
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING GRAPHICS FILES 
C
C----------------------------------------------------------------------C
C
      IF(N.GE.NCPPH.AND.ISPPH.GE.1)THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      ENDIF
C
C----------------------------------------------------------------------C
C
      IPLTTMP=0
      IF(ISVPH.EQ.1.OR.ISVPH.EQ.2)IPLTTMP=1
      IF(N.GE.NCVPH.AND.IPLTTMP.EQ.1)THEN
       CALL VELPLTH
       NCVPH=NCVPH+(NTSPTC/NPVPH)
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(N.GE.NCVPV.AND.ISVPV.GE.1)THEN
       CALL VELPLTV
       NCVPV=NCVPV+(NTSPTC/NPVPV)
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        TVAR1S(L,K)=TOX(L,K,1)
       ENDDO
      ENDDO
C
      IPLTTMP=0
      IF(ISSPH(1).EQ.1.OR.ISSPH(1).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(1).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(2).EQ.1.OR.ISSPH(2).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(2).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(3).EQ.1.OR.ISSPH(3).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(3).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(4).EQ.1.OR.ISSPH(4).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(4).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(5).EQ.1.OR.ISSPH(5).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(5).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(6).EQ.1.OR.ISSPH(6).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(6).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
      ENDIF
C
      IPLTTMP=0
      IF(ISSPH(7).EQ.1.OR.ISSPH(7).EQ.2)IPLTTMP=1
      IF(N.GE.NCSPH(7).AND.IPLTTMP.EQ.1)THEN
       IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
      ENDIF
C----------------------------------------------------------------------C
C
      DO ITMP=1,7
      IF(N.GE.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN
       CALL SALPLTV(ITMP)
       NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))
      ENDIF
      ENDDO
C----------------------------------------------------------------------C
C
C **  WRITE EFDC EXPLORER FORMAT OUTPUT
C
      IF(ISSPH(8).EQ.1.OR.ISBEXP.EQ.1)THEN
      IF(N.GE.NCSPH(8))THEN
       CALL EE_LINKAGE(0)
       NCSPH(8)=NCSPH(8)+(NTSPTC/NPSPH(8))
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING 3D HDF GRAPHICS FILES 
C
C----------------------------------------------------------------------C
C
      IF(N.EQ.NC3DO.AND.IS3DO.EQ.1)THEN
       CALL OUT3D
       NC3DO=NC3DO+(NTSPTC/NP3DO)
      ENDIF
C
C**********************************************************************C
C
C **  WRITE RESTART FILE EVERY ISRESTO M2 TIDAL CYCLES
C
      IF(ISRESTO.GE.1)THEN
        NRESTO=ISRESTO*NTSPTC
        ISSREST=MOD(N,NRESTO)
        IF(ISSREST.EQ.0)THEN
          CALL RESTOUT(0)
          IF(ISTRAN(8).GE.1)THEN
            IF(IWQRST.EQ.1) CALL WWQRST
            IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  RECORD TIME
C
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
C **  SUPPORTED ON OTHER SYSTEMS.
C
      IF(NTIMER.EQ.NTSPTC)THEN
      CALL TIMELOG(N,TIMEDAY)   
C     CALL DTIME (TARRAY)
C     WRITE(9,200)N, TARRAY(1),TARRAY(2)
C     CALL FLUSH(9)
      NTIMER=1
      ELSE
      NTIMER=NTIMER+1
      ENDIF
C
C**********************************************************************C
C
      IF(ISHOW.EQ.1) CALL SHOWVAL1
      IF(ISHOW.EQ.2) CALL SHOWVAL2
      IF(ISHOW.EQ.3) CALL SHOWVAL3
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  TIME LOOP COMPLETED
C
      IF(ISCRAY.EQ.0)THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      ENDIF
C
C**********************************************************************C
C**********************************************************************C
C
C **  CALCULATE VECTOR POTENTIAL AND VECTOR POTENTIAL TRANSPORTS
C **  USING RESULTS OF THE HARMONIC ANALYSIS
C
C     IF(ISVPTHA.NE.1) GOTO 2000
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)      
C     VPZ(L,K)=TCVP*SUB(L)*SUB(LS)*SVB(L)*SVB(L-1)*HMC(L)*
C    &         ((AMSU(L,K)+AMSU(LS,K))*(AMCV(L,K)+AMCV(L-1,K))
C    &         -(AMCU(L,K)+AMCU(LS,K))*(AMSV(L,K)+AMSV(L-1,K)))
C     ENDDO
C     ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)      
C     VPX(L,K)=TCVP*((AMSV(L,K)+AMSV(L,K+1))*(AMCW(L,K)+AMCW(LS,K))
C    &              -(AMCV(L,K)+AMCV(L,K+1))*(AMSW(L,K)+AMSW(LS,K)))
C     VPY(L,K)=TCVP*((AMSW(L,K)+AMSW(L-1,K))*(AMCU(L,K)+AMCU(L,K+1))
C    &              -(AMCW(L,K)+AMCW(L-1,K))*(AMSU(L,K)+AMSU(L,K+1)))
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     LN=LNC(L)      
C     UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L)-DZI*(VPY(L,K)-VPY(L,K-1))
C     VVPT(L,K)=DZI*(VPX(L,K)-VPX(L,K-1))-(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
C     ENDDO
C     ENDDO
C
C     DO K=1,KS
C     DO L=2,LA
C     LS=LSC(L)
C     LN=LNC(L)
C     WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)
C    &         -(VPX(LN,K)-VPX(L,K))/DYP(L)
C     ENDDO
C     ENDDO
C
C----------------------------------------------------------------------C
C
 2000 CONTINUE
C
C**********************************************************************C
C
C **  PRINT FINAL RESULTS
C
      CALL OUTPUT2
C
C**********************************************************************C
C
C **  WRITE RESTART FILE
C
      IF(ISRESTO.EQ.-1.OR.ISRESTO.EQ.-11)THEN
        CALL RESTOUT(0)
        IF(ISTRAN(8).GE.1)THEN
          IF(IWQRST.EQ.1) CALL WWQRST
          IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
        ENDIF
      ENDIF
      IF(ISRESTO.EQ.-2)THEN
        CALL RESTMOD
      ENDIF
C
C**********************************************************************C
C
C **  COMPLETE LEAST SQUARES HARMONIC ANALYSIS
C
      LSLSHA=1
      IF(ISLSHA.EQ.1) CALL LSQHARM 
C
C**********************************************************************C
C
C **  OUTPUT COURANT NUMBER DIAGNOSTICS
C
      OPEN(1,FILE='CFLMAX.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='CFLMAX.OUT')
C
      DO L=2,LA
       WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
       WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
       WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
       WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
      ENDDO
C
      CLOSE(1)
C
 1991 FORMAT(2I5,12F7.2)
 1992 FORMAT(10X,12F7.2)
C 
C**********************************************************************C
C
      RETURN
      END
