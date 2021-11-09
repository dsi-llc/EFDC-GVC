C
C**********************************************************************C
C
      SUBROUTINE RWQICI(IWQTICI)
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
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING ICS (UNIT INWQICI).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION XWQV(NWQVM)
      CHARACTER TITLE(3)*79, ICICONT*3
C
      OPEN(1,FILE=ICIFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(IWQTICI.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
C
      WRITE(2,60)'* INITIAL CONDITIONS AT ', IWQTICI,
     *  ' TH DAY FROM MODEL START'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO M=2,LA
        READ(1,84) I,J,K,(XWQV(NW),NW=1,NWQV)
        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
          PRINT*, 'I, J, K, LINE# = ', I,J,K,M-1
          STOP 'ERROR!! INVALID (I,J) IN FILE 1'
        ENDIF
        L=LIJ(I,J)
        DO NW=1,NWQV
          WQV(L,K,NW)=XWQV(NW)
        ENDDO
        WRITE(2,84) I,J,K,(WQV(L,K,NW),NW=1,NWQV)
      ENDDO
C
C: WQCHLX=1/WQCHLX
C
      DO L=2,LA
        DO K=1,KC
          WQCHL(L,K) = WQV(L,K,1)*WQCHLC + WQV(L,K,2)*WQCHLD
     *      + WQV(L,K,3)*WQCHLG
          IF(IWQSRP.EQ.1)THEN
            O2WQ = MAX(WQV(L,K,19), 0.0)
            WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQV(L,K,20) )
            WQTAMP(L,K) = WQV(L,K,20) - WQTAMD
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))
           ELSE IF(IWQSRP.EQ.2)THEN
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDT(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDT(L,K))
           ELSE
            WQPO4D(L,K) = WQV(L,K,10)
            WQSAD(L,K)  = WQV(L,K,17)
          ENDIF
        ENDDO
      ENDDO
C
      READ(1,52) IWQTICI, ICICONT
      WRITE(2,52) IWQTICI, ICICONT
      IF(ICICONT.EQ.'END')THEN
        CLOSE(1)
        IWQICI = 0
      ENDIF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 21E12.4)
C
      RETURN
      END
