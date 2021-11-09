C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WSMTSBIN
C
C**********************************************************************C
C
C **  LAST MODIFIED BY M.R. MORTON ON 02 JUNE 1999
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
C WRITE SEDIMENT TIME-SERIES OUTPUT TO BINARY FILE.
C AVERAGES BENTHIC FLUX RATES OVER ISMTSDT TIME STEPS (E.G., DAILY AVG).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'

      IF(ISSDBIN .GT. 0)THEN
        IF( MOD(ITNWQ,ISMTSDT) .EQ. 0 )THEN
          NREC4 = NREC4+1
          TIMTMP = TIMEBF / FLOAT(NBFCNT)
          OPEN(UNIT=2, FILE='WQSDTS.BIN',ACCESS='DIRECT',
     +       FORM='UNFORMATTED',STATUS='UNKNOWN', RECL=MAXRECL4)

          READ(2, REC=1) NDUM, XDUM, XDUM,
     +      XDT, IXDT, NPARM, NCELLS, NLAYERS
          NDUM=NDUM
          XDUM=XDUM
          WRITE(2, REC=1) NREC4, TBEGIN, TIMTMP,
     +      XDT, IXDT, NPARM, NCELLS, NLAYERS
          WRITE(2, REC=NR6) TIMTMP
          DO LL=2,LA
            BFO2SUM(LL)  = BFO2SUM(LL)  / FLOAT(NBFCNT)
            BFNH4SUM(LL) = BFNH4SUM(LL) / FLOAT(NBFCNT)
            BFNO3SUM(LL) = BFNO3SUM(LL) / FLOAT(NBFCNT)
            BFPO4SUM(LL) = BFPO4SUM(LL) / FLOAT(NBFCNT)
            BFSADSUM(LL) = BFSADSUM(LL) / FLOAT(NBFCNT)
            BFCODSUM(LL) = BFCODSUM(LL) / FLOAT(NBFCNT)
            BFSMTSUM(LL) = BFSMTSUM(LL) / FLOAT(NBFCNT)
            BFBSTSUM(LL) = BFBSTSUM(LL) / FLOAT(NBFCNT)

            WRITE(2) BFO2SUM(LL), BFNH4SUM(LL), BFNO3SUM(LL),
     +        BFPO4SUM(LL), BFSADSUM(LL), BFCODSUM(LL), BFSMTSUM(LL),
     +        BFBSTSUM(LL)
          ENDDO
          INQUIRE(UNIT=2, NEXTREC=NR6)
          CLOSE(2)
          CALL WQZERO4
        ENDIF
      ENDIF
C
      RETURN
      END
