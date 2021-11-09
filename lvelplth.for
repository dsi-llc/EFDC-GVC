C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LVELPLTH
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
C **  SUBROUTINE LVELPLTH WRITES A HORIZONTAL LAGRANGIAN MEAN VELOCITY
C **  VECTOR FILE AND A HORIZONTAL AVERAGED LAGRANGIAN MEAN VELOCITY 
C **  FILE 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      REAL DBS(10)
      CHARACTER*80 TITLE1,TITLE2
C
C**********************************************************************C
C
C **  WRITE HEADINGS
C
      NRLRPD=0
      NBLRPD=0
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      IF(ISRED(L).EQ.1)NRLRPD=NRLRPD+1
      IF(ISRED(L).EQ.0)NBLRPD=NBLRPD+1
      ENDDO
C
      TITLE1='HORIZONTAL LAGRANGIAN MEAN VELOCITY '
      TITLE2='HORIZ AVERAGED LAGRANGIAN MEAN VEL '
      IF(IPLRPD.EQ.1) LINES=NLRPD
      IF(IPLRPD.EQ.2) LINES=NRLRPD
      IF(IPLRPD.EQ.3) LINES=NBLRPD
      LINES=NLRPD
      LEVELS=2
      DBS(1)=0.
      DBS(2)=99.
C
      OPEN(10,FILE='LMVVECH.OUT',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='LMVVECH.OUT',STATUS='UNKNOWN')
      WRITE (10,99) TITLE1
      WRITE (10,101)LINES,LEVELS
      WRITE (10,250)(DBS(L),L=1,LEVELS)
C
      IF(IPLRPD.EQ.1)THEN
      DO M=1,MLRPDRT
      WRITE (10,100)NLRPDRT(M)
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
      VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
      VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
      VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
      UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
      VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
      VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
      VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
      WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDDO
      ENDDO
      ENDIF
C
      IF(IPLRPD.EQ.2)THEN
      DO M=1,MLRPDRT
      WRITE (10,100)NLRPDRT(M)
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      IF(ISRED(L).EQ.1)THEN
       UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDIF
      ENDDO
      ENDDO
      ENDIF
C
      IF(IPLRPD.EQ.3)THEN
      DO M=1,MLRPDRT
      WRITE (10,100)NLRPDRT(M)
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      IF(ISRED(L).EQ.0)THEN
       UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDIF
      ENDDO
      ENDDO
      ENDIF
C
      CLOSE(10)
C
C**********************************************************************C
C
      OPEN(10,FILE='ALMVVCH.OUT',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='ALMVVCH.OUT',STATUS='UNKNOWN')
      WRITE (10,99) TITLE2
      WRITE (10,101)LINES,LEVELS
      WRITE (10,250)(DBS(L),L=1,LEVELS)
C
      IF(IPLRPD.EQ.1)THEN
      M=MLRAVG
      WRITE (10,100)NTS
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
      VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
      VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
      VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
      UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
      VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
      VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
      VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
      WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDDO
      ENDIF
C
      IF(IPLRPD.EQ.2)THEN
      M=MLRAVG
      WRITE (10,100)NTS
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      IF(ISRED(L).EQ.1)THEN
       UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDIF
      ENDDO
      ENDIF
C
      IF(IPLRPD.EQ.3)THEN
      M=MLRAVG
      WRITE (10,100)NTS
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
      IF(ISRED(L).EQ.0)THEN
       UTMP=100.*STCUV(L)*XLRPD(NLR,KC,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,KC,M)
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=100.*STCUV(L)*XLRPD(NLR,1,M)
       VTMP=100.*STCUV(L)*YLRPD(NLR,1,M)
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)I,J,DLON(L),DLAT(L),VELEKC,VELNKC,
     &              VELEKB,VELNKB
      ENDIF
      ENDDO
      ENDIF
C
      CLOSE(10)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10)
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
