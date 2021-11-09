C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WSMRST
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
C WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT ISMORST.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C      LOGICAL FEXIST
C
C WRITE ASCII RESTART FILE:
C
      OPEN(1,FILE='WQSDRST.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WQSDRST.OUT',STATUS='UNKNOWN')
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
      WRITE(1,101) N,TIME
      WRITE(1,888)
C
      DO L=2,LA
        WRITE(1,90) L,(SMPON(L,NW),NW=1,NSMG),
     *    (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
     *    SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
     *    SMBST(L),SMT(L)
      ENDDO
C
      CLOSE(1)
C
C ALSO WRITE BINARY RESTART FILE:
C
C      INQUIRE(FILE='WQSDRST.BIN', EXIST=FEXIST)
C      IF(FEXIST)THEN
C        OPEN(UNIT=1, FILE='WQSDRST.BIN', ACCESS='TRANSPARENT',
C     +    FORM='UNFORMATTED', STATUS='UNKNOWN')
C        CLOSE(UNIT=1, STATUS='DELETE')
C      ENDIF
C      OPEN(UNIT=1, FILE='WQSDRST.BIN', ACCESS='TRANSPARENT',
C     +   FORM='UNFORMATTED', STATUS='UNKNOWN')
C      WRITE(1) N, TIME
C      DO L=2,LA
C        WRITE(1) L
C        WRITE(1) (SMPON(L,NW),NW=1,NSMG),
C     *    (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
C     *    SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
C     *    SMBST(L),SMT(L)
C      ENDDO
C      CLOSE(1)
C
C   90 FORMAT(I5, 18E12.4)
   90 FORMAT(I5, 1P, 18E12.4)
  101 FORMAT('CC  SM RESTART FILE TIME STEP, TIME = ',I10,F13.5)
  888 FORMAT('    L',
     *'       GPON1       GPON2       GPON3       GPOP1       GPOP2',
     *'       GPOP3       GPOC1       GPOC2       GPOC3       G1NH4',
     *'       G2NH4       G2NO3       G2PO4       G2H2S        GPSI',
     *'        G2SI        GBST          GT')
C
      RETURN
      END
