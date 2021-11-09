C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQBEN(IWQTBEN)
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
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR BENTHIC
C FLUXES OF PO4D,NH4,NO3,SAD,COD,O2,TAM (UNIT INWQBEN).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION XBFPO4D(NSMZM),XBFNH4(NSMZM),XBFNO3(NSMZM)
      DIMENSION XBFCOD(NSMZM),XBFO2(NSMZM),XBFSAD(NSMZM)
      CHARACTER TITLE(3)*79, BENCONT*3
C
      OPEN(1,FILE=SUNFN,STATUS='UNKNOWN')
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(IWQTBEN.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
      WRITE(2,60)'* BENTHIC FLUXES AT     ', IWQTBEN,
     *  ' TH DAY FROM MODEL START'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
C
      DO I=1,IWQZ
        READ(1,51) MM,XBFPO4D(I),XBFNH4(I),XBFNO3(I),
     *    XBFSAD(I),XBFCOD(I),XBFO2(I)
        WRITE(2,51) MM,XBFPO4D(I),XBFNH4(I),XBFNO3(I),XBFSAD(I),
     *    XBFCOD(I),XBFO2(I)
      ENDDO
C
      DO L=2,LA
        IWQ = IWQZMAP(L,1)
        WQBFPO4D(L)=XBFPO4D(IWQ)
        WQBFNH4(L)=XBFNH4(IWQ)
        WQBFNO3(L)=XBFNO3(IWQ)
        WQBFSAD(L)=XBFSAD(IWQ)
        WQBFCOD(L)=XBFCOD(IWQ)
        WQBFO2(L)=XBFO2(IWQ)
      ENDDO
C
      READ(1,52) IWQTBEN, BENCONT
      WRITE(2,52) IWQTBEN, BENCONT
C
      IF(BENCONT.EQ.'END')THEN
        CLOSE(1)
        IWQBEN = 0
      ENDIF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
