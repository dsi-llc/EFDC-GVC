C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQNC
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
C WRITE INFORMATION OF NEGATIVE WQ STATE VARIABLES (UNIT IWQONC).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER*5 WQVN(NWQVM)
C
      IF(IWQNC.EQ.1)THEN
C
      OPEN(1,FILE=NCOFN,STATUS='UNKNOWN',POSITION='APPEND')
C
      DATA WQVN/
     * 'BC   ','BD   ','BG   ','RPOC ','LPOC ','DOC  ','RPOP ','LPOP ',
     * 'DOP  ','PO4T ','RPON ','LPON ','DON  ','NH4  ','NO3  ','SU   ',
     * 'SA   ','COD  ','O2   ','TAM  ','FCB  ','MALG '/
C
      DO L=2,LA
        DO K=1,KC
          DO NW=1,NWQV
            IF(WQV(L,K,NW).LT.0.0) WRITE(1,90) WQVN(NW),
     *        ITNWQ,L,IL(L),JL(L),K,WQV(L,K,NW)
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(1)
C
      ENDIF
C
   90 FORMAT(A5, I8, 4I5, E11.3)
C
      RETURN
      END
