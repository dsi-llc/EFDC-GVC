C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
CTT   SUBROUTINE RWQSUN
C
C**********************************************************************C
C
C READ IN TEMPORALLY VARYING PARAMETERS FOR DAILY SOLAR RADIATION (WQI0)
C AND FRACTIONAL DAYLENGTH (WQFD) (UNIT INWQSUN).
C
C**********************************************************************C
C
CTT     INCLUDE 'EFDC.PAR'
CTT      INCLUDE 'EFDC.CMN'
C
CTT      CHARACTER TITLE(3)*79, SUNCONT*3
C
CTT      OPEN(1,FILE=SUNFN,STATUS='UNKNOWN')
CTT      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
CTT      IF(IWQTSUN.EQ.0)THEN
CTT        READ(1,50) (TITLE(M),M=1,3)
CTT        WRITE(2,999)
CTT        WRITE(2,50) (TITLE(M),M=1,3)
CTT      ENDIF
C
CTT      WRITE(2,60)'* IO & FD AT            ', IWQTSUN,
CTT     *  ' TH DAY FROM MODEL START'
C
CTT      READ(1,999)
CTT      READ(1,50) TITLE(1)
CTT      WRITE(2,50) TITLE(1)
CTT      READ(1,53) WQI0,WQFD
CTT      WRITE(2,53) WQI0,WQFD
C
CTT      IF(IWQTSUN.EQ.0)THEN
CTT        WQI1 = WQI0
CTT        WQI2 = WQI0
CTT      ENDIF
C
CTT      READ(1,52) IWQTSUN, SUNCONT
CTT      WRITE(2,52) IWQTSUN, SUNCONT
CTT      IF(SUNCONT.EQ.'END')THEN
CTT        CLOSE(1)
CTT        IWQSUN = 0
CTT      ENDIF
C
CTT      CLOSE(2)
C
CTT  999 FORMAT(1X)
CTT   50 FORMAT(A79)
CTT   52 FORMAT(I7, 1X, A3)
CTT   53 FORMAT(10F8.3)
CTT   60 FORMAT(/, A24, I5, A24)
C
CTT      RETURN
CTT      END
C
C**********************************************************************C
C
      SUBROUTINE RWQSUN
C
C**********************************************************************C
C
C **  NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
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
C **  READS AND INTERPOLATES DAILY AVERAGE SOLAR RADIATION AND
C **  DAYLIGHT FRACTION
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(ITNWQ.GT.0) GOTO 1000
C
C**********************************************************************C
C
C **  READ IN DAILY AVERAGE SOLAR SW RAD SERIES FROM FILE 'SUNDAY.INP'
C
C----------------------------------------------------------------------C
C
      OPEN(1,FILE='SUNDAY.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,7
      READ(1,1)
      ENDDO
C
      M=0
      ISPAR=1
C      MCSUNDAY=1 TEMP USE ISPAR FOR MCSUNDAY
      READ(1,*,IOSTAT=ISO)NSUNDAY,TCSUNDAY,
     &                   TASUNDAY,RMULADJ,ADDADJ
      IF(ISO.GT.0) GOTO 900
      DO M=1,NSUNDAY
        READ(1,*,IOSTAT=ISO)TSSRD(M),SOLSRD(M),SOLFRD(M)
        IF(ISO.GT.0) GOTO 900
        TSSRD(M)=TCSUNDAY*( TSSRD(M)+TASUNDAY )
        SOLSRD(M)=RMULADJ*(SOLSRD(M)+ADDADJ) * PARADJ
      ENDDO
C
      CLOSE(1)
C
      GOTO 901
C
  900 CONTINUE
      WRITE(6,601)M
      STOP
C
  901 CONTINUE
C
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR FILE SUNDAY.INP ')
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  DAILY AVERAGE SOLAR SW RADIATION INTERPOLTATION FOR WATER QUALITY
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      ELSE
        TIME=TIMESEC/86400.
      ENDIF
C
       M1=ISPAR
C      TEMP USE ISPAR FOR MCSUNDAY
C
  100  CONTINUE
       M2=M1+1
       IF(TIME.GT.TSSRD(M2))THEN
         M1=M2
         GOTO 100
        ELSE
         ISPAR=M1
C      TEMP USE ISPAR FOR MCSUNDAY
       ENDIF
C
       TDIFF=TSSRD(M2)-TSSRD(M1)
       WTM1=(TSSRD(M2)-TIME)/TDIFF
       WTM2=(TIME-TSSRD(M1))/TDIFF
       SOLSRDT=WTM1*SOLSRD(M1)+WTM2*SOLSRD(M2)
       SOLFRDT=WTM1*SOLFRD(M1)+WTM2*SOLFRD(M2)
C
C**********************************************************************C
C
      RETURN
      END
