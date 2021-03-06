C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE INITBIN
C
C**********************************************************************C
C
C M.R. MORTON    23 JUL 1998
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
C
C**********************************************************************C
C
C INITIALIZES BINARY FILE FOR EFDC OUTPUT.  PLACES CONTROL
C PARAMETERS FOR POST-PROCESSOR IN HEADER SECTION OF BINARY
C FILES WQWCAVG.BIN, WQWCMIN.BIN, AND WQWCMAX.BIN.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'

      PARAMETER(MXPARM=30)
      REAL TEND
      DIMENSION XLON(LCM), YLAT(LCM)
      INTEGER NPARM, NCELLS
      LOGICAL FEXIST
      CHARACTER*20 WQNAME(MXPARM)
      CHARACTER*10 WQUNITS(MXPARM)
      CHARACTER*3  WQCODE(MXPARM)

C
C THE FOLLOWING PARAMETERS ARE SPECIFIED IN EFDC.INP AND WQ3DWC.INP:
C KC      = NUMBER OF VERTICAL LAYERS
C IWQTSDT = NUMBER OF TIME STEPS PER DATA DUMP
C DT      = TIME STEP OF EFDC MODEL IN SECONDS
C LA      = NUMBER OF ACTIVE CELLS + 1 IN MODEL
C TBEGAN  = BEGINNING TIME OF RUN IN DAYS
C
C THE PARAMETER NPARM MUST BE CHANGED IF THE OUTPUT DATA
C IS CHANGED IN SUBROUTINE WWQTSBIN:
C NPARM   = NUMBER OF PARAMETERS WRITTEN TO BINARY FILE
C
C NREC    = NUMBER OF RECORDS WRITTEN TO BINARY FILE (ONE RECORD
C           IS A COMPLETE DATA DUMP FOR TIME INTERVAL IWQTSDT)
C
      IF(IDNOTRVA.EQ.0)THEN
        NPARM = 23
      ELSE
        NPARM = 24
      ENDIF
      NCELLS = LA-1
      NREC = 0
      TEND = TBEGIN
      MAXRECL = 32
      IF(NPARM .GE. 8)THEN
        MAXRECL = NPARM*4
      ENDIF
C
C THE FOLLOWING WATER QUALITY NAMES, UNITS, AND 3-CHARACTER CODES
C SHOULD BE MODIFIED TO MATCH THE PARAMETERS WRITTEN TO THE BINARY
C FILE IN SUBROUTINE WWQTSBIN.  THE CHARACTER STRINGS MUST BE
C EXACTLY THE LENGTH SPECIFIED BELOW IN ORDER FOR THE POST-PROCESSOR
C TO WORK CORRECTLY.
C
C BE SURE WQNAME STRINGS ARE EXACTLY 20-CHARACTERS LONG:
C------------------'         1         2'
C------------------'12345678901234567890'
      WQNAME( 1) = 'SALINITY            '
      WQNAME( 2) = 'CHLOROPHYLL-A       '
      WQNAME( 3) = 'DISSOLVED_OXYGEN    '
      WQNAME( 4) = 'TOT_ORG_CARBON      '
      WQNAME( 5) = 'DISS_ORG_CARBON     '
      WQNAME( 6) = 'PART_ORG_CARBON     '
      WQNAME( 7) = 'TOTAL_NITROGEN      '
      WQNAME( 8) = 'DISS_ORG_NITROGEN   '
      WQNAME( 9) = 'AMMONIA_NITROGEN    '
      WQNAME(10) = 'NITRATE_NITROGEN    '
      WQNAME(11) = 'PART_ORG_NITROGEN   '
      WQNAME(12) = 'TOTAL_PHOSPHORUS    '
      WQNAME(13) = 'DISS_ORG_PHOSPHORUS '
      WQNAME(14) = 'PART_ORG_PHOSPHORUS '
      WQNAME(15) = 'DISS_ORTHOPHOSPHATE '
      WQNAME(16) = 'SILICA_AVAIL_DISS   '
      WQNAME(17) = 'FECAL_COLIFORM_BAC  '
      WQNAME(18) = 'CHEM_OXYGEN_DEMAND  '
      WQNAME(19) = 'ALGAE_DIATOMS       '
      WQNAME(20) = 'GREEN_ALGAE         '
      WQNAME(21) = 'TOT_SUSPENDED_SOLIDS'
      WQNAME(22) = 'TOT_LIGHT_EXTINCTION'
      WQNAME(23) = 'WATER_TEMPERATURE   '
      WQNAME(24) = 'MACROALGAE          '

C      WQNAME( 7) = 'DISS_ORTHOPHOSPHATE '
C      WQNAME(13) = 'TOTAL_SILICA        '
C      WQNAME(14) = 'SILICA_UNAVAILABLE  '
C      WQNAME(15) = 'SILICA_AVAILABLE    '
C      WQNAME(17) = 'BOD5                '
C      WQNAME(22) = 'DAILY_MAX_DO        '
C      WQNAME(23) = 'DAILY_MIN_DO        '
C      WQNAME(28) = 'DAILY_MAX_CHL-A     '
C      WQNAME(29) = 'DAILY_MIN_CHL-A     '
C
C BE SURE WQUNITS STRINGS ARE EXACTLY 10-CHARACTERS LONG:
C-------------------'         1'
C-------------------'1234567890'
      WQUNITS( 1) = 'G/L       '
      WQUNITS( 2) = 'UG/L      '
      WQUNITS( 3) = 'MG/L      '
      WQUNITS( 4) = 'MG/L      '
      WQUNITS( 5) = 'MG/L      '
      WQUNITS( 6) = 'MG/L      '
      WQUNITS( 7) = 'MG/L      '
      WQUNITS( 8) = 'MG/L      '
      WQUNITS( 9) = 'MG/L      '
      WQUNITS(10) = 'MG/L      '
      WQUNITS(11) = 'MG/L      '
      WQUNITS(12) = 'MG/L      '
      WQUNITS(13) = 'MG/L      '
      WQUNITS(14) = 'MG/L      '
      WQUNITS(15) = 'MG/L      '
      WQUNITS(16) = 'MG/L      '
      WQUNITS(17) = 'MPN/100ML '
      WQUNITS(18) = 'MG/L      '
      WQUNITS(19) = 'UG/L      '
      WQUNITS(20) = 'UG/L      '
      WQUNITS(21) = 'MG/L      '
      WQUNITS(22) = 'PER_METER '
      WQUNITS(23) = 'DEGC      '
      WQUNITS(24) = 'UG/L      '

C      WQUNITS(22) = 'MG/L      '
C      WQUNITS(23) = 'MG/L      '
C      WQUNITS(24) = 'G/L       '
C      WQUNITS(25) = 'MG/L      '
C      WQUNITS(26) = 'MG/L      '
C      WQUNITS(27) = 'MG/L      '
C      WQUNITS(28) = 'UG/L      '
C      WQUNITS(29) = 'UG/L      '
C      WQUNITS(30) = 'G/SQ.M    '
C
C BE SURE WQCODE STRINGS ARE EXACTLY 3-CHARACTERS LONG:
C
C------------------'123'
      WQCODE( 1) = 'SAL'
      WQCODE( 2) = 'CHL'
      WQCODE( 3) = 'DOO'
      WQCODE( 4) = 'TOC'
      WQCODE( 5) = 'DOC'
      WQCODE( 6) = 'POC'
      WQCODE( 7) = 'TNN'
      WQCODE( 8) = 'DON'
      WQCODE( 9) = 'NH4'
      WQCODE(10) = 'NO3'
      WQCODE(11) = 'PON'
      WQCODE(12) = 'TPP'
      WQCODE(13) = 'DOP'
      WQCODE(14) = 'POP'
      WQCODE(15) = 'P4D'
      WQCODE(16) = 'SAD'
      WQCODE(17) = 'FCB'
      WQCODE(18) = 'COD'
      WQCODE(19) = 'DIA'
      WQCODE(20) = 'GRN'
      WQCODE(21) = 'TSS'
      WQCODE(22) = 'KET'
      WQCODE(23) = 'TEM'
      WQCODE(24) = 'MAC'

C      WQCODE( 7) = 'P4D'
C      WQCODE(13) = 'TSI'
C      WQCODE(14) = 'SUU'
C      WQCODE(15) = 'SAA'
C      WQCODE(17) = 'BD5'
C      WQCODE(22) = 'MXO'
C      WQCODE(23) = 'MNO'
C      WQCODE(28) = 'MXC'
C      WQCODE(29) = 'MNC'
C
C---------------------------------------------------------
C
C IF WQWCAVG.BIN ALREADY EXISTS, OPEN FOR APPENDING HERE.
C
      IF(ISWQAVG .EQ. 2)THEN
        INQUIRE(FILE='WQWCAVG.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCAVG.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
          WRITE(0,*) 'OLD FILE WQWCAVG.BIN FOUND...OPENING FOR APPEND'
          READ(2, REC=1) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM,
     +      NCELLS, KC
          NR1 = 1 + NPARM*3 + NCELLS*4 + (NCELLS*KC+1)*NREC + 1
          CLOSE(2)
        ELSE
          ISWQAVG=1
        ENDIF
      ENDIF
C---------------------------------------------------------
C WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
C SECTION OF THE WQWCAVG.BIN BINARY FILE:
C
C
C---------------------------------------------------------
C
C IF WQWCAVG.BIN ALREADY EXISTS, DELETE IT HERE.
C
      IF(ISWQAVG .EQ. 1)THEN
        TBEGAN = TBEGIN
        INQUIRE(FILE='WQWCAVG.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCAVG.BIN')
          CLOSE(UNIT=2, STATUS='DELETE')
          WRITE(0,*) 'OLD FILE WQWCAVG.BIN DELETED...'
        ENDIF

        OPEN(UNIT=2, FILE='WQWCAVG.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
C---------------------------------------------------------
C WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
C SECTION OF THE WQWCAVG.BIN BINARY FILE:
C
        WRITE(2) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM, NCELLS, KC
        DO I=1,NPARM
          WRITE(2) WQNAME(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQUNITS(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQCODE(I)
        ENDDO
C
C WRITE CELL I,J MAPPING REFERENCE TO HEADER SECTION OF BINARY FILE:
C
        DO L=2,LA
          WRITE(2) IL(L)
        ENDDO
        DO L=2,LA
          WRITE(2) JL(L)
        ENDDO
C
C **  READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
C **  CURVILINEAR PORTION OF THE GRID FROM FILE LXLY.INP:
C
        OPEN(1,FILE='LXLY.INP',STATUS='UNKNOWN')
C
        DO NS=1,4
          READ(1,1111)
        ENDDO
 1111   FORMAT(80X)
C
        DO LL=1,LVC
          READ(1,*) I,J,XUTME,YUTMN
          L=LIJ(I,J)
          XLON(L)=XUTME
          YLAT(L)=YUTMN
        ENDDO
        CLOSE(1)
C
C WRITE XLON AND YLAT OF CELL CENTERS TO HEADER SECTION OF
C BINARY OUTPUT FILE:
C
        DO L=2,LA
          WRITE(2) XLON(L)
        ENDDO
        DO L=2,LA
          WRITE(2) YLAT(L)
        ENDDO

        INQUIRE(UNIT=2, NEXTREC=NR1)
        CLOSE(2)
      ENDIF
C
C---------------------------------------------------------
C
C IF WQWCMIN.BIN ALREADY EXISTS, OPEN FOR APPENDING HERE.
C
      IF(ISWQMIN .EQ. 2)THEN
        INQUIRE(FILE='WQWCMIN.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCMIN.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
          WRITE(0,*) 'OLD FILE WQWCMIN.BIN FOUND...OPENING FOR APPEND'
          READ(2, REC=1) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM,
     +      NCELLS, KC
          NR2 = 1 + NPARM*3 + NCELLS*4 + (NCELLS*KC+1)*NREC + 1
          CLOSE(2)
        ELSE
          ISWQMIN=1
        ENDIF
      ENDIF
C
C---------------------------------------------------------
C
C IF WQWCMIN.BIN ALREADY EXISTS, DELETE IT HERE.
C
      IF(ISWQMIN .EQ. 1)THEN
        TBEGAN = TBEGIN
        INQUIRE(FILE='WQWCMIN.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCMIN.BIN')
          CLOSE(UNIT=2, STATUS='DELETE')
          WRITE(0,*) 'OLD FILE WQWCMIN.BIN DELETED...'
        ENDIF

        OPEN(UNIT=2, FILE='WQWCMIN.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
C--------------------------------------------------------------------
C WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
C SECTION OF THE WQWCMIN.BIN BINARY FILE:
C
        WRITE(2) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM, NCELLS, KC
        DO I=1,NPARM
          WRITE(2) WQNAME(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQUNITS(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQCODE(I)
        ENDDO
C
C WRITE CELL I,J MAPPING REFERENCE TO HEADER SECTION OF BINARY FILE:
C
        DO L=2,LA
          WRITE(2) IL(L)
        ENDDO
        DO L=2,LA
          WRITE(2) JL(L)
        ENDDO
C
C WRITE XLON AND YLAT OF CELL CENTERS TO HEADER SECTION OF
C BINARY OUTPUT FILE:
C
        DO L=2,LA
          WRITE(2) XLON(L)
        ENDDO
        DO L=2,LA
          WRITE(2) YLAT(L)
        ENDDO

        INQUIRE(UNIT=2, NEXTREC=NR2)
        CLOSE(2)
      ENDIF
C
C---------------------------------------------------------
C
C IF WQWCMAX.BIN ALREADY EXISTS, OPEN FOR APPENDING HERE.
C
      IF(ISWQMAX .EQ. 2)THEN
          INQUIRE(FILE='WQWCMAX.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCMAX.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
          WRITE(0,*) 'OLD FILE WQWCMAX.BIN FOUND...OPENING FOR APPEND'
          READ(2, REC=1) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM,
     +      NCELLS, KC
          NR3 = 1 + NPARM*3 + NCELLS*4 + (NCELLS*KC+1)*NREC + 1
          CLOSE(2)
        ELSE
          ISWQMAX=1
        ENDIF
      ENDIF
C
C---------------------------------------------------------
C
C IF WQWCMAX.BIN ALREADY EXISTS, DELETE IT HERE.
C
      IF(ISWQMAX .EQ. 1)THEN
        TBEGAN = TBEGIN
        INQUIRE(FILE='WQWCMAX.BIN', EXIST=FEXIST)
        IF(FEXIST)THEN
          OPEN(UNIT=2, FILE='WQWCMAX.BIN')
          CLOSE(UNIT=2, STATUS='DELETE')
          WRITE(0,*) 'OLD FILE WQWCMAX.BIN DELETED...'
        ENDIF

        OPEN(UNIT=2, FILE='WQWCMAX.BIN', ACCESS='DIRECT',
     +     FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=MAXRECL)
C---------------------------------------------------------
C WRITE CONTROL PARAMETERS FOR POST-PROCESSOR TO HEADER
C SECTION OF THE WQWCMAX.BIN BINARY FILE:
C
        WRITE(2) NREC, TBEGAN, TEND, DT, IWQTSDT, NPARM, NCELLS, KC
        DO I=1,NPARM
          WRITE(2) WQNAME(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQUNITS(I)
        ENDDO
        DO I=1,NPARM
          WRITE(2) WQCODE(I)
        ENDDO
C
C WRITE CELL I,J MAPPING REFERENCE TO HEADER SECTION OF BINARY FILE:
C
        DO L=2,LA
          WRITE(2) IL(L)
        ENDDO
        DO L=2,LA
          WRITE(2) JL(L)
        ENDDO
C
C WRITE XLON AND YLAT OF CELL CENTERS TO HEADER SECTION OF
C BINARY OUTPUT FILE:
C
        DO L=2,LA
          WRITE(2) XLON(L)
        ENDDO
        DO L=2,LA
          WRITE(2) YLAT(L)
        ENDDO

        INQUIRE(UNIT=2, NEXTREC=NR3)
        CLOSE(2)
      ENDIF

      RETURN
      END
