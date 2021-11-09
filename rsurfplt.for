C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSURFPLT
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
      IF(JSRPPH.NE.1) GOTO 300
C
      OPEN(10,FILE='RSURFCN.OUT',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='RSURFCN.OUT',STATUS='UNKNOWN')
      TITLE='INSTANTANEOUS SURFACE ELEVATION CONTOURS'
C
      LINES=LA-1
      LEVELS=1
      DBS=0.
C
      WRITE (10,99) TITLE
      WRITE (10,100)LINES,LEVELS
      WRITE (10,250)DBS
      CLOSE(10)
      JSRPPH=0
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
      OPEN(10,FILE='RSURFCN.OUT',POSITION='APPEND',STATUS='UNKNOWN')
      WRITE (10,100)N,TIME
C
      DO L=2,LA
      SURFEL=HLPF(L)+BELV(L)
      WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL
      ENDDO
C
      CLOSE(10)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
CMRM  200 FORMAT(2I5,1X,1P,6E13.5) 
CMRM  250 FORMAT(1P,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
