C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SURFPLTE
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
C **  SUBROUTINE SURFPLT WRITES FILES TO CONTOUR FREE SURFACE 
C **  ELEVATION
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      integer IVER
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
C **  OUTPUT EFDC EXPLORER FORMAT
C
c      IF(IPPHXY.EQ.3)THEN
C
      LINES=LA-1
      IF(JSPPH.EQ.1)THEN
        OPEN(92,FILE='SURFCON.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
c######################################################################
c     HQI Change for EFDC Explorer compatibility
c    .        FORM='UNFORMATTED')
     .        FORM='BINARY')
c######################################################################
        CLOSE(92,STATUS='DELETE')
        OPEN(92,FILE='SURFCON.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
c######################################################################
c     HQI Change for EFDC Explorer compatibility
c    .        FORM='UNFORMATTED')
     .        FORM='BINARY')
        IVER=101
c       IVER=102
c######################################################################

        WRITE(92)IVER,IC,JC,LINES
        JSPPH=0 
C        CLOSE(92) 
      ENDIF
C  
C      OPEN(92,FILE='SURFCON.OUT',POSITION='APPEND'
C     &             ,FORM='UNFORMATTED')
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON 
      ENDIF
C
      IF(ISDYNSTP.EQ.0)THEN  
        DELT=DT  
      ELSE  
        DELT=DTDYN  
      ENDIF  
      WRITE (92)N,TIME,DELT
C  
      DO L=2,LA  
        WRITE(92)HP(L)  
      ENDDO  
C  
C      CLOSE(92)  
C
c      ENDIF
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,9E14.5)
  201 FORMAT(9E14.5)
  250 FORMAT(12E12.4)
CMRM  200 FORMAT(2I5,1X,1P,8E13.5) 
CMRM  250 FORMAT(1P,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
