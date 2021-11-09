C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQRST
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
C READ ICS FROM RESTART FILE FROM INWQRST.
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
      LK=(LA-1)*KC
C
      INQUIRE(FILE='WQWCRST.BIN', EXIST=FEXIST)
      IF(.NOT. FEXIST)THEN
        OPEN(1,FILE='WQWCRST.INP',STATUS='UNKNOWN')
        READ(1,999)
        READ(1,999)
C J.S.
        NWQV0=NWQV
        IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
        DO M=1,LK
C       WRITE(8,*) 'M=  ',M       ! DEBUG, JI, 9/5/97
          READ(1,* ) L,K,(WQV(L,K,NW),NW=1,NWQV0)
        ENDDO
C
        CLOSE(1)
      ELSE
C        OPEN(UNIT=1, FILE='WQWCRST.BIN', ACCESS='TRANSPARENT',
        OPEN(UNIT=1, FILE='WQWCRST.BIN', 
     +     FORM='UNFORMATTED', STATUS='UNKNOWN')
        READ(1) NN, XTIME
        XTIME=XTIME
        WRITE(0,911) NN, XTIME
911     FORMAT(' READING BINARY WQWCRST.BIN FILE ...    NN, TIME = ',
     +     I7, F11.5)
        NWQV0=NWQV
        IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
        DO M=1,LK
          READ(1) L, K
          DO NW=1,NWQV0
            READ(1) WQV(L, K, NW)
          ENDDO
        ENDDO
        CLOSE(1)
      ENDIF
C
C INITIALIZE MACROALGAE SO BIOMASS ONLY EXISTS IN BOTTOM LAYER:
      IF(IDNOTRVA.GT.0)THEN
C      DO L=2,LA
C       WQV(L,1,22)=SMAC(LL)*WQV(L,K,22)
C      ENDDO
        IF(KC.GT.1)THEN
          DO K=2,KC
           DO L=2,LA
            WQV(L,K,22)=0.
           ENDDO
          ENDDO
        ENDIF
      ENDIF
C
   90 FORMAT(2I5, 21E12.4)
  999 FORMAT(1X)
C
      RETURN
      END
