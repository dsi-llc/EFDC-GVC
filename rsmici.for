C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSMICI(ISMTICI)
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
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING ICS (UNIT INSMICI).
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION XSMPON(NSMGM),XSMPOP(NSMGM),XSMPOC(NSMGM)
      CHARACTER TITLE(3)*79, ICICONT*3
C
      OPEN(1,FILE='WQSDICI.INP',STATUS='OLD')
C
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(ISMTICI.EQ.0)THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      ENDIF
C
      WRITE(2,60)'* INITIAL CONDITIONS AT ', ISMTICI,
     *  ' TH DAY FROM MODEL START'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
C
      DO M=2,LA
CQUESTION        READ(INSMRST,90) I,J,(XSMPON(NW),NW=1,NSMG),
        READ(1,90) I,J,(XSMPON(NW),NW=1,NSMG),
     *    (XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),XSM1NH4,
     *    XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        WRITE(2,90) I,J,(XSMPON(NW),NW=1,NSMG),
     *    (XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),XSM1NH4,
     *    XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
          PRINT*, 'I, J, LINE# = ', I,J,M-1
          STOP 'ERROR!! INVALID (I,J) IN FILE WQSDICI.INP'
        ENDIF
        L=LIJ(I,J)
        DO MM=1,NSMG
          SMPON(L,MM)=XSMPON(MM)
          SMPOP(L,MM)=XSMPOP(MM)
          SMPOC(L,MM)=XSMPOC(MM)
        ENDDO
        SM1NH4(L)=XSM1NH4
        SM2NH4(L)=XSM2NH4
        SM2NO3(L)=XSM2NO3
        SM2PO4(L)=XSM2PO4
        SM2H2S(L)=XSM2H2S
        SMPSI(L) =XSMPSI
        SM2SI(L) =XSM2SI
        SMBST(L) =XSMBST
        SMT(L)   =XSMT
      ENDDO
C
      READ(1,52) ISMTICI, ICICONT
      WRITE(2,52) ISMTICI, ICICONT
C
      IF(ICICONT.EQ.'END')THEN
        CLOSE(1)
        ISMICI = 0
      ENDIF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 22F8.4)
   90 FORMAT(2I5, 22E16.4)
C
      RETURN
      END
