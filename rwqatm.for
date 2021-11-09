C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQATM
C
C**********************************************************************C
C
C MRM **  ADDED BY MIKE MORTON  8 JUNE 1998
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
C ** COMPUTES WET ATMOSPHERIC DEPOSITION USING CONSTANT CONCENTRATIONS
C ** FOR THE 21 STATE VARIABLES MULTIPLIED BY THE RAINFALL FLOW RATE
C ** ENTERING EACH GRID CELL.  COMPUTED LOADS ARE IN G/DAY.
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C  CV2 = CONVERSION TO GET UNITS OF G/DAY
C  WQATM(NW) HAS UNITS OF MG/L
C  RAINT(L) HAS UNITS OF M/SEC
C  DXYP(L) HAS UNITS OF M2
C  WQATML(L,KC,NW) HAS UNITS OF G/DAY
      CV2=86400.0
      DO NW=1,NWQV
        DO L=2,LA
         WQATML(L,KC,NW)=WQATM(NW)*RAINT(L)*DXYP(L)*CV2
        ENDDO
      ENDDO
C
      IF(ITNWQ.EQ.0)THEN
C
       OPEN(1,FILE='WQATM.DIA',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WQATM.DIA',STATUS='UNKNOWN')
C
       IF(ISDYNSTP.EQ.0)THEN
         TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
       ELSE
         TIME=TIMESEC/86400.
       ENDIF
       WRITE(1,112) N,TIME
C
       DO L=2,LA
         WRITE(1,110) IL(L),JL(L),(WQATML(L,KC,NW),NW=1,NWQV)
       ENDDO
C
      CLOSE(1)
C
      ENDIF
C
  110 FORMAT(1X,2I4,2X,1P,7E11.3,/,15X,7E11.3,/,15X,7E11.3)
  112 FORMAT('# WET ATMOSPHERIC DEPOSITION DIAGNOSTIC FILE',/,
     &   ' N, TIME = ', I10, F12.5/)
C
C**********************************************************************C
C
      RETURN
      END
