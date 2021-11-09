C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALAREA (IAREA)
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
C**********************************************************************C
C
C ** SUBROUTINE CALAREA CALCULATES 1D CROSS SECTION AREA, WETTED
C ** PERIMETER AND DA/DH FORM DEPTH
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
CDIAG      IF(IAREA.LT.3)THEN
CDIAG        OPEN(1,FILE='CALAREA.DIA')
CDIAG        CLOSE(1,STATUS='DELETE')
CDIAG      ENDIF
C
      IF(IAREA.EQ.1) GOTO 100
      IF(IAREA.EQ.2) GOTO 200
      IF(IAREA.EQ.3) GOTO 300
      IF(IAREA.EQ.4) GOTO 400
C
C**********************************************************************C
C
C **  AREA, WETTED PERIMETER AND DA/DH INTERPOLATION FROM DEPTH
C **  FOR IAREA=1 (INITIALIZATION CALL)
C
  100 CONTINUE
C      
      DO L=2,LA
C
      M1=1
  101 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(H1P(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 101
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-H1P(L))/TDIFF
      WTM2=(H1P(L)-EHXYS(M1,L))/TDIFF
C
      IF(LCT(L).EQ.6)THEN
        FADYP1(L)=WTM1*AREADY(M1,L)+WTM2*AREADY(M2,L)
        WPDYP1(L)=WTM1*WPERDY(M1,L)+WTM2*WPERDY(M2,L)
        SRFYP1(L)=WTM1*SURFDY(M1,L)+WTM2*SURFDY(M2,L)
        DADH1(L)=(AREADY(M2,L)-AREADY(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
      ENDIF
C
      IF(LCT(L).EQ.7)THEN
        FADXP1(L)=WTM1*AREADX(M1,L)+WTM2*AREADX(M2,L)
        WPDXP1(L)=WTM1*WPERDX(M1,L)+WTM2*WPERDX(M2,L)
        SRFXP1(L)=WTM1*SURFDX(M1,L)+WTM2*SURFDX(M2,L)
        DADH1(L)=(AREADX(M2,L)-AREADX(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
      ENDIF
C
      ENDDO
C
      DO L=2,LA
C
      IF(LCT(L).EQ.6)THEN
        FADYU1(L)=0.5*(FADYP1(L)+FADYP1(L-1))
        WPDYU1(L)=0.5*(WPDYP1(L)+WPDYP1(L-1))
        SRFYU1(L)=0.5*(SRFYP1(L)+SRFYP1(L-1))
      ENDIF
C
      IF(LCT(L).EQ.7)THEN
        FADXV1(L)=0.5*(FADXP1(L)+FADXP1(LSC(L)))
        WPDXV1(L)=0.5*(WPDXP1(L)+WPDXP1(LSC(L)))
        SRFXV1(L)=0.5*(SRFXP1(L)+SRFXP1(LSC(L)))
      ENDIF
C
      ENDDO
C
      DO NJX=1,NJUNX
       L=LJUNX(NJX)
       IF(JUNTPX(NJX).EQ.1)THEN 
         FADXV1(LNC(L))=FADXP1(LNC(L))
         WPDXV1(LNC(L))=WPDXP1(LNC(L))
         SRFXV1(LNC(L))=SRFXP1(LNC(L))
       ENDIF
       IF(JUNTPX(NJX).EQ.2)THEN
         FADXV1(L)=FADXP1(LSC(L)) 
         WPDXV1(L)=WPDXP1(LSC(L)) 
         SRFXV1(L)=SRFXP1(LSC(L)) 
       ENDIF
       IF(JUNTPX(NJX).EQ.3)THEN
         FADXV1(LNC(L))=FADXP1(LNC(L))
         WPDXV1(LNC(L))=WPDXP1(LNC(L))
         SRFXV1(LNC(L))=SRFXP1(LNC(L))
         FADXV1(L)=FADXP1(LSC(L)) 
         WPDXV1(L)=WPDXP1(LSC(L)) 
         SRFXV1(L)=SRFXP1(LSC(L)) 
       ENDIF 
      ENDDO
C
      DO NJY=1,NJUNY
       L=LJUNY(NJY)
       IF(JUNTPY(NJY).EQ.1)THEN
         FADYU1(L+1)=FADYP1(L+1)
         WPDYU1(L+1)=WPDYP1(L+1)
         FADYU1(L+1)=FADYP1(L+1)
       ENDIF 
       IF(JUNTPY(NJY).EQ.2)THEN
         FADYU1(L)=FADYP1(L-1) 
         WPDYU1(L)=WPDYP1(L-1) 
         FADYU1(L)=FADYP1(L-1) 
       ENDIF 
       IF(JUNTPY(NJY).EQ.3)THEN
         FADYU1(L+1)=FADYP1(L+1)
         WPDYU1(L+1)=WPDYP1(L+1)
         FADYU1(L+1)=FADYP1(L+1)
         FADYU1(L)=FADYP1(L-1) 
         WPDYU1(L)=WPDYP1(L-1) 
         FADYU1(L)=FADYP1(L-1) 
       ENDIF 
      ENDDO
C
      GOTO 500
C                         
C**********************************************************************C
C
C **  AREA, WETTED PERIMETER AND DA/DH INTERPOLATION FROM DEPTH
C **  IAREA=2 UPDATE CALL
C 
  200 CONTINUE
C     
      DO L=2,LA
C
      M1=1
  201 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(HP(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 201
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-HP(L))/TDIFF
      WTM2=(HP(L)-EHXYS(M1,L))/TDIFF
C
      IF(LCT(L).EQ.6)THEN
        FADYP(L)=WTM1*AREADY(M1,L)+WTM2*AREADY(M2,L)
        WPDYP(L)=WTM1*WPERDY(M1,L)+WTM2*WPERDY(M2,L)
        SRFYP(L)=WTM1*SURFDY(M1,L)+WTM2*SURFDY(M2,L)
        DADH(L)=(AREADY(M2,L)-AREADY(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
      ENDIF
C
      IF(LCT(L).EQ.7)THEN
        FADXP(L)=WTM1*AREADX(M1,L)+WTM2*AREADX(M2,L)
        WPDXP(L)=WTM1*WPERDX(M1,L)+WTM2*WPERDX(M2,L)
        SRFXP(L)=WTM1*SURFDX(M1,L)+WTM2*SURFDX(M2,L)
        DADH(L)=(AREADX(M2,L)-AREADX(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
      ENDIF
C
      ENDDO
C
      DO L=2,LA
C
      IF(LCT(L).EQ.6)THEN
        FADYU(L)=0.5*(FADYP(L)+FADYP(L-1))
        WPDYU(L)=0.5*(WPDYP(L)+WPDYP(L-1))
        SRFYU(L)=0.5*(SRFYP(L)+SRFYP(L-1))
      ENDIF
C
      IF(LCT(L).EQ.7)THEN
        FADXV(L)=0.5*(FADXP(L)+FADXP(LSC(L)))
        WPDXV(L)=0.5*(WPDXP(L)+WPDXP(LSC(L)))
        SRFXV(L)=0.5*(SRFXP(L)+SRFXP(LSC(L)))
      ENDIF
C
      ENDDO
C
      DO NJX=1,NJUNX
       L=LJUNX(NJX)
       IF(JUNTPX(NJX).EQ.1)THEN 
         FADXV(LNC(L))=FADXP(LNC(L))
         WPDXV(LNC(L))=WPDXP(LNC(L))
         SRFXV(LNC(L))=SRFXP(LNC(L))
       ENDIF
       IF(JUNTPX(NJX).EQ.2)THEN
         FADXV(L)=FADXP(LSC(L)) 
         WPDXV(L)=WPDXP(LSC(L)) 
         SRFXV(L)=SRFXP(LSC(L)) 
       ENDIF
       IF(JUNTPX(NJX).EQ.3)THEN
         FADXV(LNC(L))=FADXP(LNC(L))
         WPDXV(LNC(L))=WPDXP(LNC(L))
         SRFXV(LNC(L))=SRFXP(LNC(L))
         FADXV(L)=FADXP(LSC(L)) 
         WPDXV(L)=WPDXP(LSC(L)) 
         SRFXV(L)=SRFXP(LSC(L)) 
       ENDIF 
      ENDDO
C
      DO NJY=1,NJUNY
       L=LJUNY(NJY)
       IF(JUNTPY(NJY).EQ.1)THEN
         FADYU(L+1)=FADYP(L+1)
         WPDYU(L+1)=WPDYP(L+1)
         FADYU(L+1)=FADYP(L+1)
       ENDIF 
       IF(JUNTPY(NJY).EQ.2)THEN
         FADYU(L)=FADYP(L-1) 
         WPDYU(L)=WPDYP(L-1) 
         FADYU(L)=FADYP(L-1) 
       ENDIF 
       IF(JUNTPY(NJY).EQ.3)THEN
         FADYU(L+1)=FADYP(L+1)
         WPDYU(L+1)=WPDYP(L+1)
         FADYU(L+1)=FADYP(L+1)
         FADYU(L)=FADYP(L-1) 
         WPDYU(L)=WPDYP(L-1) 
         FADYU(L)=FADYP(L-1) 
       ENDIF 
      ENDDO
C
      GOTO 500
C                         
C**********************************************************************C
C
C **  DEPTH, WETTED PERIMETER AND DA/DH INTERPOLATION FROM AREA
C **  IAREA=3 (INVERSE INTERPOLATION CALL)
C
  300 CONTINUE
C
      OPEN(1,FILE='CALAREA.DIA',POSITION='APPEND')
C
CDIAG      WRITE(1,1001)N
C
      DO L=2,LA
C
      IF(LCT(L).EQ.6)THEN
        AREATMP=FADYP(L)
        M1=1
  301   CONTINUE
        M2=M1+1
       IF(M2.GT.NXYSDAT(L)) GOTO 500
        IF(AREATMP.GT.AREADY(M2,L))THEN
          M1=M2
          GOTO 301
        ENDIF      
        TDIFF=AREADY(M2,L)-AREADY(M1,L)
        WTM1=(AREADY(M2,L)-AREATMP)/TDIFF
        WTM2=(AREATMP-AREADY(M1,L))/TDIFF
C        FADYP(L)=WTM1*AREADY(M1,L)+WTM2*AREADY(M2,L)
        WPDYP(L)=WTM1*WPERDY(M1,L)+WTM2*WPERDY(M2,L)
        HP(L)=WTM1*EHXYS(M1,L)+WTM2*EHXYS(M2,L)
        DADH(L)=(AREADY(M2,L)-AREADY(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
        SURF=HP(L)+BELV(L)
CDIAG        WRITE(1,1003)IL(L),JL(L),M1,M2,FADYP(L),AREADY(M1,L),
CDIAG     &              AREADY(M2,L),HP(L),EHXYS(M1,L),EHXYS(M2,L),SURF              
      ENDIF
C
      ENDDO
C
CDIAG      WRITE(1,1002)N
C
      DO L=2,LA
C
      IF(LCT(L).EQ.7)THEN
        AREATMP=FADXP(L)
        M1=1
  302   CONTINUE
        M2=M1+1
        IF(M2.GT.NXYSDAT(L)) GOTO 600
        IF(AREATMP.GT.AREADX(M2,L))THEN
          M1=M2
          GOTO 302
        ENDIF      
        TDIFF=AREADX(M2,L)-AREADX(M1,L)
        WTM1=(AREADX(M2,L)-AREATMP)/TDIFF
        WTM2=(AREATMP-AREADX(M1,L))/TDIFF
C        FADXP(L)=WTM1*AREADX(M1,L)+WTM2*AREADX(M2,L)
        WPDXP(L)=WTM1*WPERDX(M1,L)+WTM2*WPERDX(M2,L)
        HP(L)=WTM1*EHXYS(M1,L)+WTM2*EHXYS(M2,L)
        DADH(L)=(AREADX(M2,L)-AREADX(M1,L))/(EHXYS(M2,L)-EHXYS(M1,L))
        SURF=HP(L)+BELV(L)
CDIAG        WRITE(1,1003)IL(L),JL(L),M1,M2,FADXP(L),AREADX(M1,L),
CDIAG     &             AREADX(M2,L),HP(L),EHXYS(M1,L),EHXYS(M2,L),SURF              
      ENDIF
C
      ENDDO
C
      DO L=2,LA
C
      IF(LCT(L).EQ.6)THEN
        FADYU(L)=0.5*(FADYP(L)+FADYP(L-1))
        WPDYU(L)=0.5*(WPDYP(L)+WPDYP(L-1))
        SRFYU(L)=0.5*(SRFYP(L)+SRFYP(L-1))
CDIAG        WRITE(1,1004)IL(L),JL(L),FADYU(L),WPDYU(L),SRFYU(L)
      ENDIF
C
      IF(LCT(L).EQ.7)THEN
        FADXV(L)=0.5*(FADXP(L)+FADXP(LSC(L)))
        WPDXV(L)=0.5*(WPDXP(L)+WPDXP(LSC(L)))
        SRFXV(L)=0.5*(SRFXP(L)+SRFXP(LSC(L)))
      ENDIF
C
      ENDDO
C
CDIAG      CLOSE(1)
C
      DO NJX=1,NJUNX
       L=LJUNX(NJX)
       IF(JUNTPX(NJX).EQ.1)THEN 
         FADXV(LNC(L))=FADXP(LNC(L))
         WPDXV(LNC(L))=WPDXP(LNC(L))
         SRFXV(LNC(L))=SRFXP(LNC(L))
       ENDIF
       IF(JUNTPX(NJX).EQ.2)THEN
         FADXV(L)=FADXP(LSC(L)) 
         WPDXV(L)=WPDXP(LSC(L)) 
         SRFXV(L)=SRFXP(LSC(L)) 
       ENDIF
       IF(JUNTPX(NJX).EQ.3)THEN
         FADXV(LNC(L))=FADXP(LNC(L))
         WPDXV(LNC(L))=WPDXP(LNC(L))
         SRFXV(LNC(L))=SRFXP(LNC(L))
         FADXV(L)=FADXP(LSC(L)) 
         WPDXV(L)=WPDXP(LSC(L)) 
         SRFXV(L)=SRFXP(LSC(L)) 
       ENDIF 
      ENDDO
C
      DO NJY=1,NJUNY
       L=LJUNY(NJY)
       IF(JUNTPY(NJY).EQ.1)THEN
         FADYU(L+1)=FADYP(L+1)
         WPDYU(L+1)=WPDYP(L+1)
         FADYU(L+1)=FADYP(L+1)
       ENDIF 
       IF(JUNTPY(NJY).EQ.2)THEN
         FADYU(L)=FADYP(L-1) 
         WPDYU(L)=WPDYP(L-1) 
         FADYU(L)=FADYP(L-1) 
       ENDIF 
       IF(JUNTPY(NJY).EQ.3)THEN
         FADYU(L+1)=FADYP(L+1)
         WPDYU(L+1)=WPDYP(L+1)
         FADYU(L+1)=FADYP(L+1)
         FADYU(L)=FADYP(L-1) 
         WPDYU(L)=WPDYP(L-1) 
         FADYU(L)=FADYP(L-1) 
       ENDIF 
      ENDDO
C
      GOTO 500
C                         
C**********************************************************************C
C
C **  AREA, WETTED PERIMETER AND DA/DH INTERPOLATION FROM DEPTH
C **  FOR OPEN BOUNDARY CELLS, IAREA=4 UPDATE CALL
C 
  400 CONTINUE
C
      DO LL=1,NPBW
      L=LPBW(LL)
      M1=1
  401 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(HP(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 401
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-HP(L))/TDIFF
      WTM2=(HP(L)-EHXYS(M1,L))/TDIFF
      FADYP(L)=WTM1*AREADY(M1,L)+WTM2*AREADY(M2,L)
      ENDDO
C
      DO LL=1,NPBE
      L=LPBE(LL)
      M1=1
  402 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(HP(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 402
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-HP(L))/TDIFF
      WTM2=(HP(L)-EHXYS(M1,L))/TDIFF
      FADYP(L)=WTM1*AREADY(M1,L)+WTM2*AREADY(M2,L)
      ENDDO
C
      DO LL=1,NPBS
      L=LPBS(LL)
      M1=1
  403 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(HP(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 403
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-HP(L))/TDIFF
      WTM2=(HP(L)-EHXYS(M1,L))/TDIFF
      FADXP(L)=WTM1*AREADX(M1,L)+WTM2*AREADX(M2,L)
      ENDDO
C
      DO LL=1,NPBN
      L=LPBN(LL)
      M1=1
  404 CONTINUE
      M2=M1+1
      IF(M2.GT.NXYSDAT(L)) GOTO 600
      IF(HP(L).GT.EHXYS(M2,L))THEN
        M1=M2
        GOTO 404
      ENDIF      
      TDIFF=EHXYS(M2,L)-EHXYS(M1,L)
      WTM1=(EHXYS(M2,L)-HP(L))/TDIFF
      WTM2=(HP(L)-EHXYS(M1,L))/TDIFF
      FADXP(L)=WTM1*AREADX(M1,L)+WTM2*AREADX(M2,L)
      ENDDO
C
C**********************************************************************C
C
  500 GOTO 700
C
  600 CONTINUE
C
      WRITE(6,601)N,IL(L),JL(L),IAREA
      WRITE(8,601)N,L,IAREA
      STOP
C
  700 CONTINUE
C
  601 FORMAT(' 1D CHAN ERROR N,I,J,IA = ',I10,3I5)
 1001 FORMAT(I10,' X CHAN L,M1,M2,A,A1,A2,H,H1,H2,SEL')
 1002 FORMAT(I10,' Y CHAN L,M1,M2,A,A1,A2,H,H1,H2,SEL')
 1003 FORMAT(4I5,7F12.4)
 1004 FORMAT(2I5,7F12.4)
C
      RETURN
      END
