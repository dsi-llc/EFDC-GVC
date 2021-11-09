C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQZERO3
C
C**********************************************************************C
C
C M. MORTON  29 APR 1999
C INITIALIZES THE LIMITATION AND D.O. COMPONENT ANALYSIS ARRAYS:
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DO LL=2,LA
        DO K=1,KC
C
C ZERO THE DIURNAL DO VARIABLES:
C
          XLIMIC(LL,K) = 0.0
          XLIMID(LL,K) = 0.0
          XLIMIG(LL,K) = 0.0
          XLIMIM(LL,K) = 0.0
          XLIMVM(LL,K) = 0.0
          XLIMDM(LL,K) = 0.0
          XLIMNC(LL,K) = 0.0
          XLIMND(LL,K) = 0.0
          XLIMNG(LL,K) = 0.0
          XLIMNM(LL,K) = 0.0
          XLIMPC(LL,K) = 0.0
          XLIMPD(LL,K) = 0.0
          XLIMPG(LL,K) = 0.0
          XLIMPM(LL,K) = 0.0
          XLIMTC(LL,K) = 0.0
          XLIMTD(LL,K) = 0.0
          XLIMTG(LL,K) = 0.0
          XLIMTM(LL,K) = 0.0
          XDOSAT(LL,K) = 0.0
          XDODEF(LL,K) = 0.0
          XDOPSL(LL,K) = 0.0
          XDOSOD(LL,K) = 0.0
          XDOKAR(LL,K) = 0.0
          XDODOC(LL,K) = 0.0
          XDONIT(LL,K) = 0.0
          XDOCOD(LL,K) = 0.0
          XDOPPB(LL,K) = 0.0
          XDORRB(LL,K) = 0.0
          XDOPPM(LL,K) = 0.0
          XDORRM(LL,K) = 0.0
          XDOTRN(LL,K) = 0.0
          XDOALL(LL,K) = 0.0
          XDODZ(LL,K)  = 0.0
        ENDDO
      ENDDO
      TIMESUM3 = 0.0
      NLIM = 0
C
      RETURN
      END
