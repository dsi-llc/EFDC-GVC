C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQRST
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
C WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT IWQORST.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C      REAL XTEM
C      LOGICAL FEXIST
C
C WRITE ASCII RESTART FILE:
C
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
      WRITE(1,101) N,TIME
      WRITE(1,102)
C
C J.S.
      NWQV0=NWQV
      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
      DO L=2,LA
        DO K=1,KC
          WRITE(1,90) L,K,(WQV(L,K,NW),NW=1,NWQV0)
        ENDDO
      ENDDO
C
      CLOSE(1)
C
C ALSO WRITE BINARY RESTART FILE:
C
C      INQUIRE(FILE='WQWCRST.BIN', EXIST=FEXIST)
C      IF(FEXIST)THEN
C        OPEN(UNIT=1, FILE='WQWCRST.BIN', ACCESS='TRANSPARENT',
C     +    FORM='UNFORMATTED', STATUS='UNKNOWN')
C        CLOSE(UNIT=1, STATUS='DELETE')
C      ENDIF
C      OPEN(UNIT=1, FILE='WQWCRST.BIN', ACCESS='TRANSPARENT',
C     +   FORM='UNFORMATTED', STATUS='UNKNOWN')
C      WRITE(1) N, TIME
C      NWQV0=NWQV
C      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
C      DO L=2,LA
C        DO K=1,KC
C          WRITE(1) L, K
C          DO NW=1,NWQV0
C            WRITE(1) WQV(L,K,NW)
C          ENDDO
C        ENDDO
C      ENDDO
C      CLOSE(1)
C
   90 FORMAT(2I5, 1P, 22E12.4)
  101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',I10,F12.5)
  102 FORMAT('C   L    K  BC          BD          BG          ',
     &       'RPOC        LPOC        DOC         ',
     &       'RPOP        LPOP        DOP         PTO4        ',
     &       'RPON        LPON        DON         AMN         ',
     &       'NIT         SU          SA          COD         ',
     &       'DO          TAM         FCB        MALG')
C
      RETURN
      END
