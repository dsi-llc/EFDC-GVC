C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSMRST
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
C
C**********************************************************************C
C
C READ ICS FROM RESTART FILE FROM INSMRST.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      LOGICAL FEXIST
      INTEGER NN
      REAL XTIME
C
C CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE
C THE ASCII FILE INSTEAD.
C
      INQUIRE(FILE='WQSDRST.BIN', EXIST=FEXIST)
      IF(.NOT. FEXIST)THEN
        OPEN(1,FILE='WQSDRST.INP',STATUS='UNKNOWN')
        READ(1,999)
        READ(1,999)
C
        DO M=2,LA
          READ(1,*) L,(SMPON(L,NW),NW=1,NSMG),
     *      (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
     *      SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
     *      SMBST(L),SMT(L)
        ENDDO
C
        CLOSE(1)
      ELSE
C        OPEN(UNIT=1, FILE='WQSDRST.BIN', ACCESS='TRANSPARENT',
        OPEN(UNIT=1, FILE='WQSDRST.BIN', 
     +     FORM='UNFORMATTED', STATUS='UNKNOWN')
        READ(1) NN, XTIME
        XTIME=XTIME
        WRITE(0,911) NN, XTIME
911     FORMAT(' READING BINARY WQSDRST.BIN FILE ...    NN, TIME = ',
     +      I7, F11.5)
        DO M=2,LA
          READ(1) L
          READ(1) (SMPON(L,NW),NW=1,NSMG),
     *      (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
     *      SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
     *      SMBST(L),SMT(L)
        ENDDO
        CLOSE(1)
      ENDIF
C
   90 FORMAT(I5, 18E12.4)
  999 FORMAT(1X)
C
      RETURN
      END
