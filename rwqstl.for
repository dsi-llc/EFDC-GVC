C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQSTL(IWQTSTL)
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
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR SETTLING
C VELOCITIES OF ALGAE, RPOM, LPOM & PARTICULATE METAL (UNIT INWQSTL).
C ALSO SPATIALLY/TEMPORALLY VARYING REAERATION ADJUSTMENT FACTOR.
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER TITLE(3)*79, STLCONT*3
C
      OPEN(1,FILE=STLFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(IWQTSTL.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
      WRITE(2,60)'* SETTLING VELOCITY AT  ', IWQTSTL,
     *  ' TH DAY FROM MODEL START'

      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO I=1,IWQZ
        READ(1,*) MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),
     *    WQWSLP(I),WQWSS(I), WQWSM
        WRITE(2,51) MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),WQWSLP(I),
     *    WQWSS(I), WQWSM
      ENDDO

      READ(1,52) IWQTSTL, STLCONT
      WRITE(2,52) IWQTSTL, STLCONT
      IF(STLCONT.EQ.'END')THEN
        CLOSE(1)
        IWQSTL = 0
      ENDIF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I3, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
