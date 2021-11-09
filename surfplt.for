C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SURFPLT
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
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF(IPPHXY.LE.2)THEN
C
      IF(JSPPH.NE.1) GOTO 300
C
      OPEN(10,FILE='SURFCON.OUT')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='SURFCON.OUT')
      TITLE='INSTANTANEOUS SURFACE ELEVATION CONTOURS'
C
      LINES=LA-1
      LEVELS=1
      DBS=0.
C
      WRITE (10,99) TITLE
      WRITE (10,101)LINES,LEVELS
      WRITE (10,250)DBS
      CLOSE(10)
      JSPPH=0
C
  300 CONTINUE
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      OPEN(10,FILE='SURFCON.OUT',POSITION='APPEND')
      WRITE (10,100)N,TIME
C
      IF(IS1DCHAN.EQ.0)THEN
	  IF(IPPHXY.EQ.0)THEN
        DO L=2,LA
         SURFEL=BELV(L)+HP(L)
         WRITE(10,201)SURFEL,BELV(L),HP(L),
     &                 HBED(L,KBT(L)),HBEDA(L)
        ENDDO
	  ENDIF
	  IF(IPPHXY.EQ.1)THEN
        DO L=2,LA
         SURFEL=BELV(L)+HP(L)
         WRITE(10,200)IL(L),JL(L),SURFEL,BELV(L),HP(L),
     &                 HBED(L,KBT(L)),HBEDA(L)
        ENDDO
	  ENDIF
	  IF(IPPHXY.EQ.2)THEN
        DO L=2,LA
         SURFEL=BELV(L)+HP(L)
         WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L),HP(L),
     &                 HBED(L,KBT(L)),HBEDA(L)
        ENDDO
	  ENDIF
       ELSE
        DO L=2,LA
         SURFEL=GI*P(L)
         WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L),HP(L),
     &                FADXP(L),FADYP(L)
        ENDDO
      ENDIF
C
      CLOSE(10)
C
      ENDIF
C
C**********************************************************************C
C
C **  OUTPUT EFDC EXPLORER FORMAT
C
      IF(IPPHXY.EQ.3) CALL SURFPLTE
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
