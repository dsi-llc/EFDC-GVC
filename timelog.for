C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY PAUL M CRAIG ON 21 JULY 2003
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
      SUBROUTINE TIMELOG(N,TIMEDAY)
C
      CHARACTER*8 MRMDATE,MRMTIME*10

C WRITE OUT MODEL TIME STEP AND SUN/PC SYSTEM CLOCK TIME TO TIME.LOG FILE:
CJH      CALL DATE(MRMDATE)                                              
CJH      CALL TIME(MRMTIME)                                              
      ! *** WRITE OUT MODEL TIME STEP AND SYSTEM CLOCK TIME TO TIME.LOG
      CALL DATE_AND_TIME(MRMDATE,MRMTIME)
      WRITE(9,100)N,TIMEDAY,MRMDATE,MRMTIME

  100 FORMAT(' ','N =',I10,5X,'TIMEDAY =',F12.4,5X,'DATE = ',A8,
     &            5X,'TIME = ',A10)
      RETURN
      END
