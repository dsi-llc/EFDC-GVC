      SUBROUTINE SURFPLT  
C  
C CHANGE RECORD  
C **  SUBROUTINE SURFPLT WRITES FILES TO CONTOUR FREE SURFACE  
C **  ELEVATION  
C  
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'

      CHARACTER*80 TITLE  
C  
C *** EE BEGIN BLOCK  
C  
      INTEGER*4    VER  
C  
C *** EE END BLOCK  
C  
      IF(IPPHXY.LE.2)THEN  
        IF(JSPPH.NE.1) GOTO 300  
        OPEN(10,FILE='SURFCON.OUT')  
        CLOSE(10,STATUS='DELETE')  
        OPEN(10,FILE='SURFCON.OUT')  
        TITLE='INSTANTANEOUS SURFACE ELEVATION CONTOURS'  
        LINES=LA-1  
        LEVELS=1  
        DBS=0.  
        WRITE (10,99) TITLE  
        WRITE (10,101)LINES,LEVELS  
        WRITE (10,250)DBS  
        CLOSE(10)  
        JSPPH=0  
  300   CONTINUE  
        IF(ISDYNSTP.EQ.0)THEN  
          TIME=DT*FLOAT(N)+TCON*TBEGIN  
          TIME=TIME/TCON  
        ELSE  
          TIME=TIMESEC/TCON  
        ENDIF  
        OPEN(10,FILE='SURFCON.OUT',POSITION='APPEND')  
        WRITE (10,100)N,TIME  
        IF(IPPHXY.EQ.0)THEN  
          DO L=2,LA  
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,201)SURFEL,BELV(L),HP(L),  
     &            HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        IF(IPPHXY.EQ.1)THEN  
          DO L=2,LA  
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,200)IL(L),JL(L),SURFEL,BELV(L),HP(L),  
     &            HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        IF(IPPHXY.EQ.2)THEN  
          DO L=2,LA  
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L)
     &                 ,HP(L),HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        CLOSE(10)  
      ENDIF  
C  
   99 FORMAT(A80)  
  100 FORMAT(I10,F12.4)  
  101 FORMAT(2I10)  
  200 FORMAT(2I5,1X,9E14.5)  
  201 FORMAT(9E14.5)  
  250 FORMAT(12E12.4)  
      RETURN  
      END  

