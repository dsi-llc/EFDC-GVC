C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE TMSRHOUS
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 02/19/2002        john hamrick       02/19/2002       john hamrick
C  changed tox bed output
C----------------------------------------------------------------------C
C
C **  SUBROUTINE TMSR WRITES TIME SERIES FILES FOR SURFACE ELEVATON,
C **  VELOCITY, CONCENTRATION, AND VOLUME SOURCES AT SPECIFIED  
C **  (I,J) POINTS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C     DIMENSION ATMP(KCM),BTMP(KCM)
C
      DIMENSION SEDSND(KCM),TXWF(KCM),TXWC(KCM),TXWP(KCM),TXBT(KBM),
     &          SEDSNDB(KBM),TXBF(KBM),TXBC(KBM),TXBP(KBM),PORH(KBM)
      DIMENSION QCHANUIJ(LCM),QCHANVIJ(LCM)
C
      CHARACTER(80) :: TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7,
     &         TITLE11,TITLE12,TITLE13,TITLE14,TITLE15,TITLE16,TITLE17,
     &         TITLE18,TITLE19,TFNTXWT,TFNTXWF,TFNTXWC,TFNTXWP,
     &         TFNTXBT,TFNTXBF,TFNTXBC,TFNTXBP 
      CHARACTER(10) :: CTUNIT
      CHARACTER(1) ::  CZTT(0:9)
      CHARACTER(1) ::  CCHTMF,CCHTMS
      CHARACTER(2) ::  CNTOX(25),CNSND(6),CNSBL(6)
      CHARACTER(5) ::  OUTLOC
C
C**********************************************************************C
C
      OUTLOC(1:2)='&&'
      IF(JSTMSR.NE.1) GOTO 300
C
C----------------------------------------------------------------------C
C
      IF(MLTMSR.GT.MLTMSRM)THEN
        WRITE (6,600)
        STOP
      ENDIF
      IF(MLTMSR.GT.99)THEN
        WRITE (6,601)
        STOP
      ENDIF
C
  600 FORMAT(' NUMBER OF TIME SER LOC, MLTMSR, EXCEEDS DIM, MLTMSRM')
  601 FORMAT(' NUMBER OF TIME SERIES LOCATIONS EXCEED 99')
C
            TAUW1=0.
            TAUW2=0.
c
      CZTT(0)='0'
      CZTT(1)='1'
      CZTT(2)='2'
      CZTT(3)='3'
      CZTT(4)='4'
      CZTT(5)='5'
      CZTT(6)='6'
      CZTT(7)='7'
      CZTT(8)='8'
      CZTT(9)='9'
C
      DO MLTM=1,MLTMSR
      MSDIG=MOD(MLTM,10)
      MTMP=MLTM-MSDIG
      MFDIG=MTMP/10
      CCHTMF=CZTT(MFDIG)
      CCHTMS=CZTT(MSDIG)
      CNTMSR(MLTM)= CCHTMF // CCHTMS
      ENDDO
C
      IF(TCTMSR.EQ.1.) CTUNIT='SECONDS'
      IF(TCTMSR.EQ.60.) CTUNIT='MINUTES'
      IF(TCTMSR.EQ.3600.) CTUNIT='HOURS'
      IF(TCTMSR.EQ.86400.) CTUNIT='DAYS'
C
       CNTOX( 1)= '01'
       CNTOX( 2)= '02'
       CNTOX( 3)= '03'
       CNTOX( 4)= '04'
       CNTOX( 5)= '05'
       CNTOX( 6)= '06'
       CNTOX( 7)= '07'
       CNTOX( 8)= '08'
       CNTOX( 9)= '09'
       CNTOX(10)= '10'
       CNTOX(11)= '11'
       CNTOX(12)= '12'
       CNTOX(13)= '13'
       CNTOX(14)= '14'
       CNTOX(15)= '15'
       CNTOX(16)= '16'
       CNTOX(17)= '17'
       CNTOX(18)= '18'
       CNTOX(19)= '19'
       CNTOX(20)= '20'
       CNTOX(21)= '21'
       CNTOX(22)= '22'
       CNTOX(23)= '23'
       CNTOX(24)= '24'
       CNTOX(25)= '25'
C
       CNSND( 1)= '01'
       CNSND( 2)= '02'
       CNSND( 3)= '03'
       CNSND( 4)= '04'
       CNSND( 5)= '05'
       CNSND( 6)= '06'
C
       CNSBL( 1)= '01'
       CNSBL( 2)= '02'
       CNSBL( 3)= '03'
       CNSBL( 4)= '04'
       CNSBL( 5)= '05'
       CNSBL( 6)= '06'
C

C **  WRITE HEADINGS
C
      TITLE1=' SALINITY (PSU) TIME SERIES, K=1,KC'
      TITLE2=' TEMPERATURE (DEG C) TIME SERIES, K=1,KC'
      TITLE3=' DYE CONC (KG/M**3) TIME SERIES, K=1,KC'
      TITLE4=' SED CONC (MG/LITER) TIME SERIES, K=1,KC'
      TITLE5=' TOXIC CONC (M/TOT VOL - UG/LITER) 1-4 BED,5-8 WC'
      TITLE6=' VISCOSITY (CM**2/S) TIME SERIES, K=1,KS'
      TITLE7=' DIFFUSIVITY (CM**2/S) TIME SERIES, K=1,KS'
      TITLE11=' SURFACE ELEVATION & DEPTH (METERS) TIME SERIES'
      TITLE12=' EXT MODE E,N VEL (CM/S) TBX TBY TB TBCG TBNG (CM/S)**2'
      TITLE13=' EXT MODE U,V TRANSPORT (M**3/S) TIME SERIES'
      TITLE14=' INT MODE EAST VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE15=' INT MODE NORTH VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE16=' EXT MODE VOLUME S/S (M**3/S) TIME SERIES'
      TITLE17=' INT MODE VOL S/S (M**3/S) TIME SERIES, K=1,KC'
      TITLE18=' SED BED LOAD QSX QSY (GM/S) CQSX CQSY (MG/L) '
      TITLE19=' BED TOP  KBT HBED(KBT) HBED(KBT-1) VOIDR FRAC SED/SND'
      TFNTXWT=' TOTAL TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
      TFNTXWF=' FREE DIS TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'  
      TFNTXWC=' DOC COMP TOXIC CONC WATER COL (M/TOT VOL), UG/LITER' 
      TFNTXWP=' TOT PART TOXIC CONC WATER COL (M/M), UG/GM' 
      TFNTXBT=' TOTAL TOXIC CONC SED BED (M/TOT VOL), UG/LITER' 
      TFNTXBF=' FREE DIS TOXIC CONC SED BED (M/PORE VOL), UG/LITER' 
      TFNTXBC=' DOC COMP TOXIC CONC SED BED (M/PORE VOL), UG/LITER' 
      TFNTXBP=' TOT PART TOXIC CONC SED BED (M/M), UG/GM' 
C
      IF(ISTMSR.EQ.2)THEN
      DO MLTM=1,MLTMSR
        IF(MTMSRC(MLTM).EQ.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            FNSAL(MLTM)='SALTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            FNTEM(MLTM)='TEMTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            FNDYE(MLTM)='DYETS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            FNSFL(MLTM)='SFLTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            FNSED(MLTM)='SEDTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
C          IF(ISTRAN(7).GE.1)THEN
C            FNSND(MLTM)='SNDTS' // CNTMSR(MLTM) // '.OUT'
C          ENDIF
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
            FNSND(MLTM,NX)='SND' // CNSND(NX) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNSBL(MLTM,NX)='SBL' // CNSBL(NX) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            ENDDO
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            FNDOX(MLTM)='DOXTS' // CNTMSR(MLTM) // '.OUT'
            FNTOC(MLTM)='TOCTS' // CNTMSR(MLTM) // '.OUT'
            FNNHX(MLTM)='NHXTS' // CNTMSR(MLTM) // '.OUT'
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
            FNTXWT(MLTM,NT)='TXWT' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXWF(MLTM,NT)='TXWF' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXWC(MLTM,NT)='TXWC' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXWP(MLTM,NT)='TXWP' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXBT(MLTM,NT)='TXBT' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXBF(MLTM,NT)='TXBF' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXBC(MLTM,NT)='TXBC' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            FNTXBP(MLTM,NT)='TXBP' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          FNAVV(MLTM)='AVVTS' // CNTMSR(MLTM) // '.OUT'
          FNAVB(MLTM)='AVBTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          FNSEL(MLTM)='SELTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF 
        IF(MTMSRUE(MLTM).EQ.1)THEN
          FNUVE(MLTM)='UVETS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
          FNUVT(MLTM)='UVTTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRU(MLTM).EQ.1)THEN
          FNU3D(MLTM)='U3DTS' // CNTMSR(MLTM) // '.OUT'
          FNV3D(MLTM)='V3DTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          FNQQE(MLTM)='QQETS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          FNQ3D(MLTM)='Q3DTS' // CNTMSR(MLTM) // '.OUT'
        ENDIF
      ENDDO
      JSTMSR=0
      ENDIF
C
      IF(JSTMSR.EQ.0) GOTO 300
C
      DO MLTM=1,MLTMSR
        WRITE(OUTLOC(3:5),107) MLTM
        IF(OUTLOC(3:3).EQ.' ') OUTLOC(3:3)='0'
        IF(OUTLOC(4:4).EQ.' ') OUTLOC(4:4)='0'
        IF(MTMSRC(MLTM).EQ.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            FNSAL(MLTM)='SALTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            CLOSE(11,STATUS='DELETE')
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            WRITE (11,100) TITLE1
            WRITE (11,101) CLTMSR(MLTM)
            WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (11,102) CTUNIT
            CLOSE(11)
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            FNTEM(MLTM)='TEMTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            CLOSE(21,STATUS='DELETE')
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            WRITE (21,100) TITLE2
            WRITE (21,101) CLTMSR(MLTM)
            WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (21,102) CTUNIT
            CLOSE(21)
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            FNDYE(MLTM)='DYETS' // CNTMSR(MLTM) // '.OUT'
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            FNDYE(MLTM)='SFLTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            JHFNSED=305
           IF(MLTM.EQ.1) THEN
            FNSED(MLTM)='SEDTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(JHFNSED,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            CLOSE(JHFNSED,STATUS='DELETE')
            OPEN(JHFNSED,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            WRITE (JHFNSED,100) TITLE4
            WRITE (JHFNSED,102) CTUNIT
           ENDIF
            WRITE (JHFNSED,104) OUTLOC,CLTMSR(MLTM)
            WRITE (JHFNSED,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah        CLOSE(41)
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
            JHFNSND=400+NX
           IF(MLTM.EQ.1) THEN
            FNSND(MLTM,NX)='SND'// CNSND(NX) // 'TS' // 
     &                        CNTMSR(MLTM) // '.OUT'
            OPEN(JHFNSND,FILE=FNSND(MLTM,NX),STATUS='UNKNOWN')
            CLOSE(JHFNSND,STATUS='DELETE')
            OPEN(JHFNSND,FILE=FNSND(MLTM,NX),STATUS='UNKNOWN')
            WRITE (JHFNSND,100) TITLE4
            WRITE (JHFNSND,102) CTUNIT
           ENDIF
            WRITE (JHFNSND,104) OUTLOC,CLTMSR(MLTM)
            WRITE (JHFNSND,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah        CLOSE(41)
            ENDDO
            DO NX=1,NSND
            JHFNSBL=400+NSND+NX
           IF(MLTM.EQ.1) THEN
            FNSBL(MLTM,NX)='SBL'// CNSBL(NX) // 'TS' // 
     &                        CNTMSR(MLTM) // '.OUT'
            OPEN(JHFNSBL,FILE=FNSBL(MLTM,NX),STATUS='UNKNOWN')
            CLOSE(JHFNSBL,STATUS='DELETE')
            OPEN(JHFNSBL,FILE=FNSBL(MLTM,NX),STATUS='UNKNOWN')
            WRITE (JHFNSBL,100) TITLE18
            WRITE (JHFNSBL,102) CTUNIT
           ENDIF
            WRITE (JHFNSBL,104) OUTLOC,CLTMSR(MLTM)
            WRITE (JHFNSBL,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah        CLOSE(41)
            ENDDO
C            FNSND(MLTM)='SNDTS' // CNTMSR(MLTM) // '.OUT'
C            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
C            CLOSE(41,STATUS='DELETE')
C            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
C            WRITE (41,100) TITLE4
C            WRITE (41,101) CLTMSR(MLTM)
C            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
C            WRITE (41,102) CTUNIT
C            CLOSE(41)
          ENDIF
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
            JHFNBED=600
           IF(MLTM.EQ.1) THEN
            FNBED(MLTM)='BEDTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(JHFNBED,FILE=FNBED(MLTM),STATUS='UNKNOWN')
            CLOSE(JHFNBED,STATUS='DELETE')
            OPEN(JHFNBED,FILE=FNBED(MLTM),STATUS='UNKNOWN')
            WRITE (JHFNBED,100) TITLE19
            WRITE (JHFNBED,102) CTUNIT
           ENDIF
            WRITE (JHFNBED,104) OUTLOC,CLTMSR(MLTM)
            WRITE (JHFNBED,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah        CLOSE(41)
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            FNDOX(MLTM)='DOXTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(41,FILE=FNDOX(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNDOX(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
            FNTOC(MLTM)='TOCTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(42,FILE=FNTOC(MLTM),STATUS='UNKNOWN')
            CLOSE(42,STATUS='DELETE')
            OPEN(42,FILE=FNTOC(MLTM),STATUS='UNKNOWN')
            WRITE (42,100) TITLE4
            WRITE (42,101) CLTMSR(MLTM)
            WRITE (42,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (42,102) CTUNIT
            CLOSE(42)
            FNNHX(MLTM)='NHXTS' // CNTMSR(MLTM) // '.OUT'
            OPEN(43,FILE=FNNHX(MLTM),STATUS='UNKNOWN')
            CLOSE(43,STATUS='DELETE')
            OPEN(43,FILE=FNNHX(MLTM),STATUS='UNKNOWN')
            WRITE (43,100) TITLE4
            WRITE (43,101) CLTMSR(MLTM)
            WRITE (43,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (43,102) CTUNIT
            CLOSE(43)
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
              JHFNTOX=310+NT
             IF(MLTM.EQ.1) THEN
              FNTOX(MLTM,NT)='TOX' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTOX,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTOX,STATUS='DELETE')
              OPEN(JHFNTOX,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTOX,100) TITLE5
              WRITE (JHFNTOX,102) CTUNIT
             ENDIF
              WRITE (JHFNTOX,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTOX,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
              JHFNTXWT=310+NTOX+NT
             IF(MLTM.EQ.1) THEN
              FNTXWT(MLTM,NT)='TXWT' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXWT,FILE=FNTXWT(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXWT,STATUS='DELETE')
              OPEN(JHFNTXWT,FILE=FNTXWT(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXWT,100) TFNTXWT
              WRITE (JHFNTXWT,102) CTUNIT
             ENDIF
              WRITE (JHFNTXWT,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXWT,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
              JHFNTXWF=310+NTOX*2+NT
             IF(MLTM.EQ.1) THEN
              FNTXWF(MLTM,NT)='TXWF' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXWF,FILE=FNTXWF(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXWF,STATUS='DELETE')
              OPEN(JHFNTXWF,FILE=FNTXWF(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXWF,100) TFNTXWF
              WRITE (JHFNTXWF,102) CTUNIT
             ENDIF
              WRITE (JHFNTXWF,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXWF,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
              JHFNTXWC=310+NTOX*3+NT
             IF(MLTM.EQ.1) THEN
              FNTXWC(MLTM,NT)='TXWC' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXWC,FILE=FNTXWC(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXWC,STATUS='DELETE')
              OPEN(JHFNTXWC,FILE=FNTXWC(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXWC,100) TFNTXWC
              WRITE (JHFNTXWC,102) CTUNIT
             ENDIF
              WRITE (JHFNTXWC,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXWC,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
              JHFNTXWP=310+NTOX*4+NT
             IF(MLTM.EQ.1) THEN
              FNTXWP(MLTM,NT)='TXWP' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXWP,FILE=FNTXWP(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXWP,STATUS='DELETE')
              OPEN(JHFNTXWP,FILE=FNTXWP(MLTM,NT),STATUS='UNKNOWN')
              WRITE(JHFNTXWP,100) TFNTXWP
              WRITE (JHFNTXWP,102) CTUNIT
             ENDIF
              WRITE(JHFNTXWP,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXWP,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
             JHFNTXBT=310+NTOX*5+NT
             IF(MLTM.EQ.1) THEN
              FNTXBT(MLTM,NT)='TXBT' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXBT,FILE=FNTXBT(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXBT,STATUS='DELETE')
              OPEN(JHFNTXBT,FILE=FNTXBT(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXBT,100) TFNTXBT
              WRITE (JHFNTXBT,102) CTUNIT
             ENDIF
              WRITE (JHFNTXBT,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXBT,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
             JHFNTXBF=310+NTOX*6+NT
             IF(MLTM.EQ.1) THEN
              FNTXBF(MLTM,NT)='TXBF' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXBF,FILE=FNTXBF(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXBF,STATUS='DELETE')
              OPEN(JHFNTXBF,FILE=FNTXBF(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXBF,100) TFNTXBF
              WRITE (JHFNTXBF,102) CTUNIT
             ENDIF
              WRITE (JHFNTXBF,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXBF,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
             JHFNTXBC=310+NTOX*7+NT
             IF(MLTM.EQ.1) THEN
              FNTXBC(MLTM,NT)='TXBC' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXBC,FILE=FNTXBC(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXBC,STATUS='DELETE')
              OPEN(JHFNTXBC,FILE=FNTXBC(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXBC,100) TFNTXBC
              WRITE (JHFNTXBC,102) CTUNIT
             ENDIF
              WRITE (JHFNTXBC,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXBC,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
              JHFNTXBP=310+NTOX*8+NT
             IF(MLTM.EQ.1) THEN
              FNTXBP(MLTM,NT)='TXBP' // CNTOX(NT) // 'TS' // 
     &                           CNTMSR(MLTM) // '.OUT'
              OPEN(JHFNTXBP,FILE=FNTXBP(MLTM,NT),STATUS='UNKNOWN')
              CLOSE(JHFNTXBP,STATUS='DELETE')
              OPEN(JHFNTXBP,FILE=FNTXBP(MLTM,NT),STATUS='UNKNOWN')
              WRITE (JHFNTXBP,100) TFNTXBP
              WRITE (JHFNTXBP,102) CTUNIT
             ENDIF
              WRITE (JHFNTXBP,104) OUTLOC,CLTMSR(MLTM)
              WRITE (JHFNTXBP,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah          CLOSE(51)
            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          FNAVV(MLTM)='AVVTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          WRITE (61,100) TITLE6
          WRITE (61,101) CLTMSR(MLTM)
          WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (61,102) CTUNIT
          CLOSE(61)
          FNAVB(MLTM)='AVBTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE7
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          JHFNSEL=610
         IF(MLTM.EQ.1) THEN
          FNSEL(MLTM)='SELTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(JHFNSEL,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          CLOSE(JHFNSEL,STATUS='DELETE')
          OPEN(JHFNSEL,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          WRITE (JHFNSEL,100) TITLE11
          WRITE (JHFNSEL,102) CTUNIT
         ENDIF 
          WRITE (JHFNSEL,104) OUTLOC,CLTMSR(MLTM)
          WRITE (JHFNSEL,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah      CLOSE(11)
        ENDIF
        IF(MTMSRUE(MLTM).EQ.1)THEN
          JHFNUVE=611
         IF(MLTM.EQ.1) THEN
          FNUVE(MLTM)='UVETS' // CNTMSR(MLTM) // '.OUT'
          OPEN(JHFNUVE,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          CLOSE(JHFNUVE,STATUS='DELETE')
          OPEN(JHFNUVE,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          WRITE (JHFNUVE,100) TITLE12
          WRITE (JHFNUVE,101) CLTMSR(MLTM)
         ENDIF
          WRITE (JHFNUVE,104)OUTLOC, CTUNIT
          WRITE (JHFNUVE,105)OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah      CLOSE(21)
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
          JHFNUVT=612
         IF(MLTM.EQ.1) THEN
          FNUVT(MLTM)='UVTTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(JHFNUVT,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          CLOSE(JHFNUVT,STATUS='DELETE')
          OPEN(JHFNUVT,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          WRITE (JHFNUVT,100) TITLE13
          WRITE (JHFNUVT,102) CTUNIT
         ENDIF
          WRITE (JHFNUVT,104) OUTLOC, CLTMSR(MLTM)
          WRITE (JHFNUVT,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah      CLOSE(31)
        ENDIF
        IF(MTMSRU(MLTM).EQ.1)THEN
          FNU3D(MLTM)='U3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          CLOSE(41,STATUS='DELETE')
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          WRITE (41,100) TITLE14
          WRITE (41,101) CLTMSR(MLTM)
          WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (41,102) CTUNIT
          CLOSE(41)
          FNV3D(MLTM)='V3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          CLOSE(51,STATUS='DELETE')
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          WRITE (51,100) TITLE15
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          CLOSE(51)
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          JHFNQQE=614
         IF(MLTM.EQ.1) THEN
          FNQQE(MLTM)='QQETS' // CNTMSR(MLTM) // '.OUT'
          OPEN(JHFNQQE,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          CLOSE(JHFNQQE,STATUS='DELETE')
          OPEN(JHFNQQE,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          WRITE (JHFNQQE,100) TITLE16
          WRITE (JHFNQQE,102) CTUNIT
         ENDIF
          WRITE (JHFNQQE,104) OUTLOC,CLTMSR(MLTM)
          WRITE (JHFNQQE,105) OUTLOC,ILTMSR(MLTM),JLTMSR(MLTM)
cjah      CLOSE(61)
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          FNQ3D(MLTM)='Q3DTS' // CNTMSR(MLTM) // '.OUT'
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE17
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
        ENDIF
      ENDDO
C
      JSTMSR=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      ELSE
        TIME=TIMESEC/TCTMSR
      ENDIF
C
      FOURDPI=4./PI
C
C **  STEP CURRENT TIME INTERVALS FOR WRITE SCENARIOS
C
      DO NTSSS=1,NTSSTSP
       DO MTSSS=1,MTSSTSP(NTSSS)
        IF(TIME.GE.TSSTRT(MTSSS,NTSSS))THEN
        IF(TIME.LE.TSSTOP(MTSSS,NTSSS))THEN
          MTSCUR(NTSSS)=MTSSS
        ENDIF
        ENDIF
       ENDDO
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO L=2,LA
          QCHANUIJ(L)=0.0
          QCHANVIJ(L)=0.0
        ENDDO
        DO NMD=1,MDCHH
          LMDCHHT=LMDCHH(NMD)
          LMDCHUT=LMDCHU(NMD)
          LMDCHVT=LMDCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            QCHANUIJ(LMDCHHT)=QCHANUIJ(LMDCHHT)+QCHANU(NMD)
            QCHANUIJ(LMDCHUT)=QCHANUIJ(LMDCHUT)-QCHANU(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            QCHANVIJ(LMDCHHT)=QCHANVIJ(LMDCHHT)+QCHANV(NMD)
            QCHANVIJ(LMDCHVT)=QCHANVIJ(LMDCHVT)-QCHANV(NMD)
          ENDIF
          IF(MDCHTYP(NMD).EQ.3)THEN
            QCHANUIJ(LMDCHHT)=QCHANUIJ(LMDCHHT)+QCHANU(NMD)
            QCHANUIJ(LMDCHUT)=QCHANUIJ(LMDCHUT)-QCHANU(NMD)
            QCHANVIJ(LMDCHHT)=QCHANVIJ(LMDCHHT)+QCHANV(NMD)
            QCHANVIJ(LMDCHVT)=QCHANVIJ(LMDCHVT)-QCHANV(NMD)
          ENDIF
        ENDDO
      ELSE
        DO L=2,LA
          QCHANUIJ(L)=0.0
          QCHANVIJ(L)=0.0
        ENDDO
      ENDIF
C
      DO MLTM=1,MLTMSR
       WRITE(OUTLOC(3:5),107) MLTM
       IF(OUTLOC(3:3).EQ.' ') OUTLOC(3:3)='0'
       IF(OUTLOC(4:4).EQ.' ') OUTLOC(4:4)='0'
       NTSSS=NTSSSS(MLTM)
       MTSCC=MTSCUR(NTSSS)
       IF(TIME.GE.TSSTRT(MTSCC,NTSSS))THEN
       IF(TIME.LE.TSSTOP(MTSCC,NTSSS))THEN
        I=ILTMSR(MLTM)
        J=JLTMSR(MLTM)
        L=LIJ(I,J)
        LN=LNC(L)
        IF(MTMSRC(MLTM).EQ.1)THEN
          IF(ISTRAN(1).GE.1)THEN
            OPEN(11,FILE=FNSAL(MLTM),POSITION='APPEND')
            WRITE(11,201)TIME,(SAL(L,K),K=1,KC)
            CLOSE(11)
          ENDIF
          IF(ISTRAN(2).GE.1)THEN
            OPEN(21,FILE=FNTEM(MLTM),POSITION='APPEND')
            WRITE(21,201)TIME,(TEM(L,K),K=1,KC)
            CLOSE(21)
          ENDIF
          IF(ISTRAN(3).GE.1)THEN
            OPEN(31,FILE=FNDYE(MLTM),POSITION='APPEND')
            WRITE(31,201)TIME,(DYE(L,K),K=1,KC)
            CLOSE(31)
          ENDIF
          IF(ISTRAN(4).GE.1)THEN
            OPEN(31,FILE=FNSFL(MLTM),POSITION='APPEND')
            WRITE(31,201)TIME,(SFL(L,K),K=1,KC)
            CLOSE(31)
          ENDIF
          IF(ISTRAN(6).GE.1)THEN
            JHFNSED=305
cjah        OPEN(41,FILE=FNSED(MLTM),POSITION='APPEND')
            IF(ISNDAL.EQ.2)THEN
                SEDBTMP=SEDBT(L,KBT(L))+SEDBT(L,KBT(L)-1)
            ELSE
                SEDBTMP=SEDBT(L,KBT(L))
            ENDIF
            WRITE(JHFNSED,301)OUTLOC,TIME,SEDBTMP,(SEDT(L,K),K=1,KC)
cjah        CLOSE(41)
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
            JHFNSND=400+NX
cjah        OPEN(JHFNSND,FILE=FNSND(MLTM,NX),POSITION='APPEND')
            IF(ISNDAL.EQ.2)THEN
                SNDBTMP=SNDB(L,KBT(L),NX)+SNDB(L,KBT(L)-1,NX)
            ELSE
                SNDBTMP=SNDB(L,KBT(L),NX)
            ENDIF
            WRITE(JHFNSND,301)OUTLOC,TIME,SNDBTMP,(SND(L,K,NX),K=1,KC),
     &                                 SNDEQSAV(L,NX)
cjah        CLOSE(41)
            ENDDO
            DO NX=1,NSND
            JHFNSBL=400+NSND+NX
cjah        OPEN(JHFNSBL,FILE=FNSBL(MLTM,NX),POSITION='APPEND')
c            CQBEDLOADX=0.
c            CQBEDLOADY=0.
            IF(UHDYE(L).NE.0.0)CQBEDLOADX(L,NX)=QSBDLDX(L,NX)/UHDYE(L)
            IF(VHDXE(L).NE.0.0)CQBEDLOADY(L,NX)=QSBDLDY(L,NX)/VHDXE(L)
            WRITE(JHFNSBL,301)OUTLOC,TIME,QSBDLDX(L,NX),QSBDLDY(L,NX),
     &               CQBEDLOADX(L,NX),CQBEDLOADY(L,NX),SNDFBL(L,NX)
cjah        CLOSE(41)
            ENDDO
C            OPEN(41,FILE=FNSND(MLTM,NX),POSITION='APPEND')
C            WRITE(41,201)TIME,SNDBT(L,KBT(L)),(SNDT(L,K),K=1,KC)
C            CLOSE(41)
          ENDIF
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
            JHFNBED=600
cjah        OPEN(41,FILE=FNBED(MLTM),POSITION='APPEND')
            KTMP=KBT(L)
            KTMP1=KBT(L)-1
            KTMP1=MAX(KTMP1,1)
            NSXD=NSED+NSND
c            WRITE(41,221)TIME,KTMP,HBED(L,KTMP),HBED(L,KTMP1),
            WRITE(JHFNBED,321)OUTLOC,TIME,KTMP,HBED(L,KTMP),
     &             HBED(L,KTMP1),
     &             HBEDA(L),VDRBED(L,KTMP),(VFRBED(L,KTMP,NX),NX=1,NSXD)
cjah        CLOSE(41)
          ENDIF
          IF(ISTRAN(8).GE.1)THEN
            OPEN(41,FILE=FNDOX(MLTM),POSITION='APPEND')
            WRITE(41,201)TIME,(WQV(L,K,19),K=1,KC)
            CLOSE(41)
            OPEN(42,FILE=FNTOC(MLTM),POSITION='APPEND')
            WRITE(42,201)TIME,(WQV(L,K,6),K=1,KC)
            CLOSE(42)
            OPEN(43,FILE=FNNHX(MLTM),POSITION='APPEND')
            WRITE(43,201)TIME,(WQV(L,K,14),K=1,KC)
            CLOSE(43)
          ENDIF
          IF(ISTRAN(5).GE.1)THEN
            DO NT=1,NTOX
C
cjah          OPEN(51,FILE=FNTOX(MLTM,NT),POSITION='APPEND')
cjah          OPEN(52,FILE=FNTXWT(MLTM,NT),POSITION='APPEND')
cjah          OPEN(53,FILE=FNTXWF(MLTM,NT),POSITION='APPEND')
cjah          OPEN(54,FILE=FNTXWC(MLTM,NT),POSITION='APPEND')
cjah          OPEN(55,FILE=FNTXWP(MLTM,NT),POSITION='APPEND')
cjah          OPEN(56,FILE=FNTXBT(MLTM,NT),POSITION='APPEND')
cjah          OPEN(57,FILE=FNTXBF(MLTM,NT),POSITION='APPEND')
cjah          OPEN(58,FILE=FNTXBC(MLTM,NT),POSITION='APPEND')
cjah          OPEN(59,FILE=FNTXBP(MLTM,NT),POSITION='APPEND')
              JHFNTOX=310+NT
              JHFNTXWT=310+NTOX+NT
              JHFNTXWF=310+NTOX*2+NT
              JHFNTXWC=310+NTOX*3+NT
              JHFNTXWP=310+NTOX*4+NT
              JHFNTXBT=310+NTOX*5+NT
              JHFNTXBF=310+NTOX*6+NT
              JHFNTXBC=310+NTOX*7+NT
              JHFNTXBP=310+NTOX*8+NT
C
              NDOC=NSED+NSND+1
              DO K=1,KC
                SEDSND(K)=SEDT(L,K)+SNDT(L,K)
                TXWF(K)=TOXFDFW(L,K,NT)*TOX(L,K,NT)
                TXWC(K)=TOXCDFW(L,K,NT)*TOX(L,K,NT)
                TXWP(K)=TOXPFTW(L,K,NT)*TOX(L,K,NT)
              ENDDO
              DO K=1,KB
                SEDSNDB(K)=SEDBT(L,K)+SNDBT(L,K)
                TXBF(K)=0.
                TXBC(K)=0.
                TXBP(K)=0.
                TXBT(K)=0.
              ENDDO
              DO K=1,KBT(L)
                PORH(K)=1.0/PORBED(L,K)
                TXBF(K)=TOXFDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
                TXBC(K)=TOXCDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
                TXBP(K)=TOXPFTB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
                TXBT(K)=TOXB(L,K,NT)/HBED(L,K)
              ENDDO
C
              K=KBT(L)
              WRITE(JHFNTOX,301)OUTLOC,TIME,TXBT(K),TXBF(K),TXBC(K),
     &           TXBP(K),TOX(L,1,NT),TXWF(1),TXWC(1),TXWP(1)
C                WRITE(51,201)TIME,TOXFDFB(L,K,NT),TOXCDFB(L,K,NT),
C     &          TOXPFTB(L,K,NT),TOXFDFW(L,1,NT),TOXCDFW(L,1,NT),
C     &          TOXPFTW(L,1,NT)
C
              DO K=1,KC
                IF(SEDSND(K).GT.0.0)TXWP(K)=1000.*TXWP(K)/SEDSND(K)
              ENDDO
C
              DO K=1,KBT(L)
                TXBP(K)=TXBP(K)*HBED(L,K)
                TXBF(K)=TXBF(K)*PORH(K)
                TXBC(K)=TXBC(K)*PORH(K)
                IF(SEDSNDB(K).GT.0.0)TXBP(K)=1000.*TXBP(K)/SEDSNDB(K)
              ENDDO
              WRITE(JHFNTXWT,301)OUTLOC,TIME,(TOX(L,K,NT),K=1,KC)
              WRITE(JHFNTXWF,301)OUTLOC,TIME,(TXWF(K),K=1,KC)
              WRITE(JHFNTXWC,301)OUTLOC,TIME,(TXWC(K),K=1,KC)
              WRITE(JHFNTXWP,301)OUTLOC,TIME,(TXWP(K),K=1,KC)
              WRITE(JHFNTXBT,301)OUTLOC,TIME,(TXBT(K),K=1,KB)
              WRITE(JHFNTXBF,301)OUTLOC,TIME,(TXBF(K),K=1,KB)
              WRITE(JHFNTXBC,301)OUTLOC,TIME,(TXBC(K),K=1,KB)
              WRITE(JHFNTXBP,301)OUTLOC,TIME,(TXBP(K),K=1,KB)
C
cjah          CLOSE(51)
cjah          CLOSE(52)
cjah          CLOSE(53)
cjah          CLOSE(54)
cjah          CLOSE(55)
cjah          CLOSE(56)
cjah          CLOSE(57)
cjah          CLOSE(58)
cjah          CLOSE(59)
C
            ENDDO
c            DO NT=1,NTOX
c            OPEN(51,FILE=FNTOXPF(MLTM,NT),POSITION='APPEND')
c            WRITE(51,201)TIME,TOXPFTB(L,KBT(L),NT),
c     &                  (TOXPFTW(L,K,NT),K=1,KC)
c            CLOSE(51)
c            ENDDO
c            DO NT=1,NTOX
c            OPEN(51,FILE=FNTOXFD(MLTM,NT),POSITION='APPEND')
c            WRITE(51,201)TIME,TOXPFTB(L,KBT(L),NT),
c     &                  (TOXPFTW(L,K,NT),K=1,KC)
c            CLOSE(51)
c            ENDDO
          ENDIF
        ENDIF
        IF(MTMSRA(MLTM).EQ.1)THEN
          OPEN(61,FILE=FNAVV(MLTM),POSITION='APPEND')
          OPEN(71,FILE=FNAVB(MLTM),POSITION='APPEND')
           DO K=1,KS
           ATMP(K)=10000.*AV(L,K)*HP(L)
           ENDDO
          WRITE(61,201)TIME,(ATMP(K),K=1,KS)
           DO K=1,KS
           ATMP(K)=10000.*AB(L,K)*HP(L)
           ENDDO
          WRITE(71,201)TIME,(ATMP(K),K=1,KS)
          CLOSE(61)
          CLOSE(71)
        ENDIF
        IF(MTMSRP(MLTM).EQ.1)THEN
          JHFNSEL=610
cjah      OPEN(11,FILE=FNSEL(MLTM),POSITION='APPEND')
          PPTMP=HP(L)+BELV(L)
          TMPVAL=VDWASTE(L)/DXYP(L)
C          HHTMP=PPTMP-BELV(L)
C          GWELTMP=AGWELV(L)-BELAGW(L)
          WRITE(JHFNSEL,301)OUTLOC,TIME,PPTMP,HP(L),BELV(L),HBEDA(L),
     &          ZELBEDA(L),TMPVAL,VDWASTE(L)
cjah      CLOSE(11)
        ENDIF 
        IF(MTMSRUE(MLTM).EQ.1)THEN
          JHFNUVE=611
cjah      OPEN(21,FILE=FNUVE(MLTM),POSITION='APPEND')
          UTMP1=50.*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=50.*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          IF(SPB(L).EQ.0)THEN
            UTMP1=2.*UTMP1
            VTMP1=2.*VTMP1
          ENDIF
          UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
          UTMP1=5000.*(TBX(L+1)+TBX(L))
          VTMP1=5000.*(TBY(LN)+TBY(L))
          TBEAST=CUE(L)*UTMP1+CVE(L)*VTMP1
          TBNORT=CUN(L)*UTMP1+CVN(L)*VTMP1
C          TAUBDYN=10000.*QQ(L,0)/CTURB2
          TAUBDYN=10000.*TAUB(L)
c          UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
c          VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
          TAUB2=TAUBDYN*TAUBDYN
          IF(ISWAVE.GT.0)THEN
            CURANG=ATAN2(VTMP,UTMP)
            TAUW1=10000.*QQWV1(L)
            TAUW2=10000.*QQWV2(L)
            TAUB2=TAUB2+0.5*(TAUW2*TAUW2)
     &           +FOURDPI*TAUBDYN*TAUW2*COS(CURANG-WACCWE(L))
          END IF
          TAUB2=MAX(TAUB2,0.)
          TAUTOT=SQRT(TAUB2)
          VELMAG=UTMP*UTMP+VTMP*VTMP
          TMPDRAG=0.0
          TAUBSEDDYN=10000.*TAUBSED(L)
          TAUBSNDDYN=10000.*TAUBSND(L)
          IF(VELMAG.GT.0.0)TMPDRAG=TAUTOT/VELMAG
C          WRITE(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUB,TAUW1,
C     &                 TAUW2,TAUTOT,TMPDRAG
          WRITE(JHFNUVE,301)OUTLOC,TIME,UTMP,VTMP,TBEAST,TBNORT,TAUBDYN,
     &                 TAUBSEDDYN,TAUBSNDDYN
cjah      CLOSE(21)
        ENDIF
        IF(MTMSRUT(MLTM).EQ.1)THEN
cjah      OPEN(31,FILE=FNUVT(MLTM),POSITION='APPEND')
          JHFNUVT=612
          QBEDLOADX=0.
          QBEDLOADY=0.
          CQBEDLOADXT=0.
          CQBEDLOADYT=0.
          DO NX=1,NSND
            QBEDLOADX=QBEDLOADX+QSBDLDX(L,NX)
            QBEDLOADY=QBEDLOADY+QSBDLDY(L,NX)
            CQBEDLOADXT=CQBEDLOADXT+CQBEDLOADX(L,NX)
            CQBEDLOADYT=CQBEDLOADYT+CQBEDLOADY(L,NX)
          ENDDO
c          IF(UHDYE(L).NE.0.0)CQBEDLOADX=QBEDLOADX/UHDYE(L)
c          IF(VHDXE(L).NE.0.0)CQBEDLOADY=QBEDLOADY/VHDXE(L)
          WRITE(JHFNUVT,301)OUTLOC,TIME,UHDYE(L),VHDXE(L),QBEDLOADX,
     &                 QBEDLOADY,
     &                 CQBEDLOADXT,CQBEDLOADYT
cjah      CLOSE(31)
        ENDIF
        IF(MTMSRU(MLTM).EQ.1)THEN
          OPEN(41,FILE=FNU3D(MLTM),POSITION='APPEND')
          OPEN(51,FILE=FNV3D(MLTM),POSITION='APPEND')
          RUVTMP=50.
          IF(SPB(L).EQ.0) RUVTMP=100.
          DO K=1,KC
           UTMP1=RUVTMP*(U(L+1,K)+U(L,K))
           VTMP1=RUVTMP*(V(LN,K)+V(L,K))
           ATMP(K)=CUE(L)*UTMP1+CVE(L)*VTMP1
           BTMP(K)=CUN(L)*UTMP1+CVN(L)*VTMP1
          ENDDO
          WRITE(41,201)TIME,(ATMP(K),K=1,KC)
          WRITE(51,201)TIME,(BTMP(K),K=1,KC)
          CLOSE(41)
          CLOSE(51)
        ENDIF
        IF(MTMSRQE(MLTM).EQ.1)THEN
          JHFNQQE=614
cjah      OPEN(61,FILE=FNQQE(MLTM),POSITION='APPEND')
          QRNT=DXYP(L)*RAINT(L)
          WRITE(JHFNQQE,301)OUTLOC,TIME,QSUME(L),QRNT,EVAPSW(L),
     &      EVAPGW(L),
     &      RIFTR(L), QCHANUIJ(L),QCHANVIJ(L),QDWASTE(L),VDWASTE(L)
cjah      CLOSE(61)
        ENDIF
        IF(MTMSRQ(MLTM).EQ.1)THEN
          OPEN(71,FILE=FNQ3D(MLTM),POSITION='APPEND')
          WRITE (71,201)TIME,(QSUM(L,K),K=1,KC)
          CLOSE(71)
        ENDIF
       ENDIF
       ENDIF
      ENDDO  
C
C**********************************************************************C
C
  100 FORMAT(A80)
  101 FORMAT('  AT LOCATION  ',A20)
  104 FORMAT(A5,1x,'  AT LOCATION  ',A20)
  102 FORMAT('  TIME IN FIRST COLUMN HAS UNITS OF ',A10)
  103 FORMAT('  CELL I,J = ',2I5)
  105 FORMAT(A5,1x,'  CELL I,J = ',2I5)
  107   FORMAT(I3)
  201 FORMAT(F12.5,12E12.4)
  221 FORMAT(F12.5,I5,14E12.4)
  301 FORMAT(a5,1x,F12.5,12E12.4)
  321 FORMAT(a5,1x,F12.5,I5,14E12.4)
C
C**********************************************************************C
C
      RETURN
      END
