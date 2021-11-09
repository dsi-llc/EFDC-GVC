C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE HYDOUT
C
C **  ADDED BY DON KINGERY, CH2M-HILL, CH2M-HILL ON 25 MARCH 1997
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
C **  SUBROUTINE HYDOUT WRITES TIME SERIES FILES FOR BOTTOM VELOCITY
C **  OVER THE MODEL DOMAIN FOR INPUT TO THE SEDIMENT BED MODEL
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4
      CHARACTER*14 TUNITS
C
C**********************************************************************C
C
      IF(JSHYDOUT.NE.1) GOTO 300
C
C----------------------------------------------------------------------C
C
C  WRITE HEADING
C
      TITLE1='BOTTOM VELOCITY OUTPUT FOR SEDIMENT BED MODEL'
      TITLE2='   DLAT          DLONG'
      TITLE3='        UVEL     VVEL     SPD     DIR '
      TITLE4='       (CM/S)   (CM/S)   (CM/S)  (DEG)'
C
      OPEN(11,FILE='SEDBED.HYD',STATUS='UNKNOWN')
      CLOSE(11,STATUS='DELETE')
      OPEN(11,FILE='SEDBED.HYD',STATUS='UNKNOWN')
      WRITE (11,100) TITLE1
      WRITE (11,110) TITLE2,TITLE3
      WRITE (11,120) TITLE4
  100 FORMAT(A80)
  110 FORMAT(1X,A22,A38)
  120 FORMAT(23X,A38)
      CLOSE(11)
      JSHYDOUT=0

  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      ELSE
        TIME=TIMESEC/TCTMSR
      ENDIF
      IF(TCTMSR.EQ.1.0)THEN
        TUNITS=' SECONDS'
       ELSE IF(TCTMSR.EQ.60.)THEN
        TUNITS=' MINUTES'
       ELSE IF(TCTMSR.EQ.3600.)THEN
        TUNITS=' HOURS'
       ELSE IF(TCTMSR.EQ.86400.)THEN
        TUNITS=' DAYS'
       ELSE
        TUNITS=' UNKNOWN UNITS'
      ENDIF 
C
      OPEN(11, FILE='SEDBED.HYD',POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(11,200)TIME,TUNITS
  200 FORMAT(1X, 'TIME:',F12.3,A14)
      DO L=2,LC-1
      LN=LNC(L)
      RUVTMP=50.
      IF(SPB(L).EQ.0) RUVTMP=100.
        UTMP1=RUVTMP*(U(L+1,1)+U(L,1))
        VTMP1=RUVTMP*(V(LN,1)+V(L,1))
        UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
        VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
        SPD=SQRT(UTMP**2+VTMP**2)
        IF(UTMP.EQ.0.0.AND.VTMP.EQ.0.0)THEN
        DIR=0
        GOTO 210
      ENDIF
      DIR=ATAN2(UTMP,VTMP)
      DIR=57.2958*DIR
      IF(DIR.LT.0.)DIR=DIR+360.
  210   WRITE(11,201)DLAT(L),DLON(L),UTMP,VTMP,SPD,DIR
      ENDDO
      CLOSE(11)
  201 FORMAT(E12.7,3X,E12.7,3X,3(F6.1,3X),F4.0)

      RETURN
      END
