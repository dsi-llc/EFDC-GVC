C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE TMSRBIN
C
C**********************************************************************C
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
C**********************************************************************C
C
C M. MORTON  07 JUN 1999
C INITIALIZES THE BINARY HYDRODYNAMIC ARRAYS TO 0.0
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C SUM MOST RECENT VALUES INTO SUMMATION ARRAYS:
C
      DO LL=2,LA
        PPTMP = GI*P(LL)
        HHTMP = PPTMP - BELV(LL)
        SELSUM(LL) = SELSUM(LL) + PPTMP
        DEPSUM(LL) = DEPSUM(LL) + HHTMP
        LN=LNC(LL)
        UTMP1 = 50.*(UHDYE(LL+1) + UHDYE(LL))/(DYP(LL)*HP(LL))
        VTMP1 = 50.*(VHDXE(LN)   + VHDXE(LL))/(DXP(LL)*HP(LL))
        IF(SPB(LL) .EQ. 0)THEN
          UTMP1 = 2.*UTMP1
          VTMP1 = 2.*VTMP1
        ENDIF
        UTMP = CUE(LL)*UTMP1 + CVE(LL)*VTMP1
        VTMP = CUN(LL)*UTMP1 + CVN(LL)*VTMP1
        VXXSUM(LL) = VXXSUM(LL) + UTMP
        VYYSUM(LL) = VYYSUM(LL) + VTMP
C        QXXSUM(LL) = QXXSUM(LL) + UHDYE(LL)
C        QYYSUM(LL) = QYYSUM(LL) + VHDXE(LL)
        UTMP = MAX( UHDYE(LL), UHDYE(LL+1) )
        VTMP = MAX( VHDXE(LL), VHDXE(LN) )
        QXXSUM(LL) = QXXSUM(LL) + UTMP
        QYYSUM(LL) = QYYSUM(LL) + VTMP
      ENDDO
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN
        TIMTMP=TIMTMP/TCTMSR  
      ELSE
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF
      TIMEHYD = TIMEHYD + TIMTMP
      NHYCNT = NHYCNT + 1
C
C CHECK IF TIME TO DUMP RESULTS TO BINARY FILE:
C
      IF(N .GE. NBTMSR .AND. N .LE. NSTMSR)THEN
        IF(NCTMSR .EQ. NWTMSR)THEN
          NREC0 = NREC0+1
          TIMTMP = TIMEHYD / NHYCNT
          OPEN(UNIT=2, FILE='HYDTS.BIN',ACCESS='DIRECT',
     +       FORM='UNFORMATTED',STATUS='UNKNOWN', RECL=MAXRECL0)

          READ(2, REC=1) NDUM, XDUM, XDUM,
     +      XDT, IXDT, NPARM, NCELLS, NLAYERS
          NDUM=NDUM
          XDUM=XDUM
          WRITE(2, REC=1) NREC0, TBEGIN, TIMTMP,
     +      XDT, IXDT, NPARM, NCELLS, NLAYERS

          WRITE(2, REC=NR0) TIMTMP
          DO LL=2,LA
            SELSUM(LL) = SELSUM(LL)  / FLOAT(NHYCNT)
            DEPSUM(LL) = DEPSUM(LL)  / FLOAT(NHYCNT)
            VXXSUM(LL) = VXXSUM(LL)  / FLOAT(NHYCNT)
            VYYSUM(LL) = VYYSUM(LL)  / FLOAT(NHYCNT)
            QXXSUM(LL) = QXXSUM(LL)  / FLOAT(NHYCNT)
            QYYSUM(LL) = QYYSUM(LL)  / FLOAT(NHYCNT)

            WRITE(2) SELSUM(LL), DEPSUM(LL), VXXSUM(LL), VYYSUM(LL),
     *        QXXSUM(LL), QYYSUM(LL), BELV(LL), ZBR(LL)
          ENDDO
          INQUIRE(UNIT=2, NEXTREC=NR0)
          CLOSE(2)

          CALL HYDZERO
        ENDIF
      ENDIF
C
      RETURN
      END
