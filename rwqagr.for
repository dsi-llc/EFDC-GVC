C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQAGR(IWQTAGR)
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
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR ALGAL
C GROWTH, RESP. & PREDATION RATES, AND BASE LIGHT EXTINCT. COEFF.
C (UNIT INWQAGR).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER TITLE(3)*79, AGRCONT*3
C
      OPEN(1,FILE=AGRFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(IWQTAGR.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
C
      WRITE(2,60)'* ALGAL KINETIC VALUE AT', IWQTAGR,
     *  ' TH DAY FROM MODEL START'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO I=1,IWQZ
C        READ(1,51) MM,WQPMC(I),WQPMD(I),WQPMG(I),WQBMRC(I),
C     *    WQBMRD(I),WQBMRG(I),WQPRRC(I),WQPRRD(I),WQPRRG(I),WQKEB(I),
C     *    WQSDCOEF(I)
C        WRITE(2,51) MM,WQPMC(I),WQPMD(I),WQPMG(I),WQBMRC(I),
C     *    WQBMRD(I),WQBMRG(I),WQPRRC(I),WQPRRD(I),WQPRRG(I),WQKEB(I),
C     *    WQSDCOEF(I)
        READ(1,*) MM, WQPMC(I),WQPMD(I),WQPMG(I),WQPMM(I),WQBMRC(I),
     *    WQBMRD(I),WQBMRG(I),WQBMRM(I),WQPRRC(I),WQPRRD(I),
     *    WQPRRG(I),WQPRRM(I),WQKEB(I)
        WRITE(2,51) MM, WQPMC(I),WQPMD(I),WQPMG(I),WQPMM(I),WQBMRC(I),
     *    WQBMRD(I),WQBMRG(I),WQBMRM(I),WQPRRC(I),WQPRRD(I),
     *    WQPRRG(I),WQPRRM(I),WQKEB(I)
      ENDDO
C
      READ(1,52) IWQTAGR, AGRCONT
      WRITE(2,52) IWQTAGR, AGRCONT
C
      IF(AGRCONT.EQ.'END')THEN
        CLOSE(1)
        IWQAGR = 0
      ENDIF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 14F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
