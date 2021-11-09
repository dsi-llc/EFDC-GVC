C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE VELPLTHE
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY PAUL M. CRAIG ON 21 JULY 2003
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 05/02/2002        John Hamrick       05/01/2002       John Hamrick
C  added real flags RSSBCE(L),RSSBCW(L),RSSBCN(L),RSSBCS(L)
C  to modified  the outputed cell center velocity for cells have source/sinks
C----------------------------------------------------------------------C
C
C **  SUBROUTINE VELPLTH WRITES A HORIZONTAL INSTANTANEOUS VELOCITY
C **  VECTOR FILE 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      REAL DBS(10)
      DIMENSION UPROFILE(KCM),VPROFILE(KCM)
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7
C *** EE BEGIN BLOCK
      INTEGER*4    IVER
C *** EE END BLOCK
C
C***********************************************************************C
C *** EE BEGIN BLOCK
C *** OUTPUT EFDC EXPLORER FORMAT.  DO NOT CHANGE OUTPUTS!  
C ***                               MUST EXACTLY MATCH EFDC_EXPLORER INPUTS!
C
C      IF(IVPHXY.EQ.3)THEN
C
        IF(JSVPH.EQ.1)THEN  
          LINES=LA-1  
          OPEN(10,FILE='VELVECH.OUT',STATUS='UNKNOWN',
     .         ACCESS='SEQUENTIAL',FORM='BINARY')
          IVER=102  
          WRITE(10)IVER,IC,JC,KC,LINES
          CLOSE(10)  

C**********************************************************************C
C *** HQI BEGIN BLOCK
          OPEN(94,FILE='TAUVECH.OUT',STATUS='UNKNOWN',
     &            ACCESS='SEQUENTIAL',FORM='BINARY')
          CLOSE(94,STATUS='DELETE')
          OPEN(94,FILE='TAUVECH.OUT',STATUS='UNKNOWN',
     &            ACCESS='SEQUENTIAL',FORM='BINARY')
          WRITE(94)IVER,IC,JC,KC,LINES
C *** HQI END BLOCK
C**********************************************************************C

          JSVPH=0  
        ENDIF
C  
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)+TCON*TBEGIN
        ELSE
          TIME=TIMESEC
        ENDIF
        TIME=TIME/86400.  ! PMC OUTPUT IN DAYS
C
        IF(ISDYNSTP.EQ.0)THEN  
          DELT=DT  
        ELSE  
          DELT=DTDYN  
        ENDIF  
C  
        OPEN(10,FILE='VELVECH.OUT',POSITION='APPEND',STATUS='OLD',
     &       FORM='BINARY')

        WRITE (10)N,TIME,DELT  
C  
        IF(ISVPH.EQ.1)THEN  
          DO L=2,LA  
            LN=LNC(L)  
            DO K=1,KC  
              UTMPS=0.5*STCUV(L)*(RSSBCE(L)*U(L+1,K)+RSSBCW(L)*U(L,K))  ! m/s
              VTMPS=0.5*STCUV(L)*(RSSBCN(L)*V(LN ,K)+RSSBCS(L)*V(L,K))  
              UPROFILE(K)=CUE(L)*UTMPS+CVE(L)*VTMPS  
              VPROFILE(K)=CUN(L)*UTMPS+CVN(L)*VTMPS  
            ENDDO  
            WRITE(10)(UPROFILE(K),VPROFILE(K),W(L,K),K=1,KC)
          ENDDO  
        ENDIF  

C**********************************************************************C
C *** HQI BEGIN BLOCK
        WRITE (94)N,TIME,DELT
C
        IF(ISVPH.EQ.1)THEN
          DO L=2,LA
            LN=LNC(L)
            UTMP1=5000.*(TBX(L+1)+TBX(L))
            VTMP1=5000.*(TBY(LN)+TBY(L))
            TBEAST=CUE(L)*UTMP1+CVE(L)*VTMP1
            TBNORT=CUN(L)*UTMP1+CVN(L)*VTMP1
            TMPV1=10000.*TAUB(L)
            TMPV2=10000.*TAUBSED(L)
            TMPV3=10000.*TAUBSND(L)
            WRITE(94)TBEAST,TBNORT,TMPV1,TMPV2,TMPV3
          ENDDO
        ENDIF
C *** HQI END BLOCK
C**********************************************************************C
C  
        CLOSE(10)  
C
C      ENDIF
C *** EE END BLOCK
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,8E14.6)
  201 FORMAT(8E14.6)
  250 FORMAT(12E12.4)
CMRM  200 FORMAT(2I5,1X,1p,6E13.5) 
CMRM  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
