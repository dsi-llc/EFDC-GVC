C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE VELPLTV
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
C **  SUBROUTINE VELPLTV WRITES A FIL FOR VERTICAL PLANE CONTOURING
C **  OF VELOCITY NORMAL TO AN ARBITARY SEQUENCE OF (I,J) POINTS AND
C **  AND VERTICAL PLANE TANGENTIAL-VERTICAL VELOCITY VECTORS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      REAL VELN (KCM,100)
      REAL VELT (KCM,100)
      REAL WZCORD (KCM,100)
      CHARACTER*80 TITLE1,TITLE2
C
C**********************************************************************C
C
      IF(JSVPV.NE.1) GOTO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      LEVELS=KC
      TITLE1='INSTANTANEOUS NORMAL VELOCITY CONTOURS'
      TITLE2='INSTANTANEOUS TANGENTIAL VELOCITY VECTORS'
C
      IF(ISECVPV.GE.1)THEN
        OPEN(11,FILE='VELCNV1.OUT')
        OPEN(21,FILE='VELVCV1.OUT')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        OPEN(11,FILE='VELCNV1.OUT')
        OPEN(21,FILE='VELVCV1.OUT')
      ENDIF 
      IF(ISECVPV.GE.2)THEN
        OPEN(12,FILE='VELCNV2.OUT')
        OPEN(22,FILE='VELVCV2.OUT')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        OPEN(12,FILE='VELCNV2.OUT')
        OPEN(22,FILE='VELVCV2.OUT')
      ENDIF
      IF(ISECVPV.GE.3)THEN
        OPEN(13,FILE='VELCNV3.OUT')
        OPEN(23,FILE='VELVCV3.OUT')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        OPEN(13,FILE='VELCNV3.OUT')
        OPEN(23,FILE='VELVCV3.OUT')
      ENDIF
      IF(ISECVPV.GE.4)THEN
        OPEN(14,FILE='VELCNV4.OUT')
        OPEN(24,FILE='VELVCV4.OUT')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        OPEN(14,FILE='VELCNV4.OUT')
        OPEN(24,FILE='VELVCV4.OUT')
      ENDIF 
      IF(ISECVPV.GE.5)THEN
        OPEN(15,FILE='VELCNV5.OUT')
        OPEN(25,FILE='VELVCV5.OUT')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        OPEN(15,FILE='VELCNV5.OUT')
        OPEN(25,FILE='VELVCV5.OUT')
      ENDIF
      IF(ISECVPV.GE.6)THEN
        OPEN(16,FILE='VELCNV6.OUT')
        OPEN(26,FILE='VELVCV6.OUT')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        OPEN(16,FILE='VELCNV6.OUT')
        OPEN(26,FILE='VELVCV6.OUT')
      ENDIF
      IF(ISECVPV.GE.7)THEN
        OPEN(17,FILE='VELCNV7.OUT')
        OPEN(27,FILE='VELVCV7.OUT')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        OPEN(17,FILE='VELCNV7.OUT')
        OPEN(27,FILE='VELVCV7.OUT')
      ENDIF 
      IF(ISECVPV.GE.8)THEN
        OPEN(18,FILE='VELCNV8.OUT')
        OPEN(28,FILE='VELVCV8.OUT')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        OPEN(18,FILE='VELCNV8.OUT')
        OPEN(28,FILE='VELVCV8.OUT')
      ENDIF
      IF(ISECVPV.GE.9)THEN
        OPEN(19,FILE='VELCNV9.OUT')
        OPEN(29,FILE='VELVCV9.OUT')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        OPEN(19,FILE='VELCNV9.OUT')
        OPEN(29,FILE='VELVCV9.OUT')
      ENDIF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LINES=NIJVPV(IS)
      WRITE (LUN1,99)TITLE1,CVTITLE(LUN1)
      WRITE (LUN2,99)TITLE2,CVTITLE(LUN2)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      CLOSE(LUN1)
      CLOSE(LUN2)
      ENDDO
C
      JSVPV=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      IF(ISECVPV.GE.1)THEN
        OPEN(11,FILE='VELCNV1.OUT',POSITION='APPEND')
        OPEN(21,FILE='VELVCV1.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.2)THEN
        OPEN(12,FILE='VELCNV2.OUT',POSITION='APPEND')
        OPEN(22,FILE='VELVCV2.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.3)THEN
        OPEN(13,FILE='VELCNV3.OUT',POSITION='APPEND')
        OPEN(23,FILE='VELVCV3.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.4)THEN
        OPEN(14,FILE='VELCNV4.OUT',POSITION='APPEND')
        OPEN(24,FILE='VELVCV4.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.5)THEN
        OPEN(15,FILE='VELCNV5.OUT',POSITION='APPEND')
        OPEN(25,FILE='VELVCV5.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.6)THEN
        OPEN(16,FILE='VELCNV6.OUT',POSITION='APPEND')
        OPEN(26,FILE='VELVCV6.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.7)THEN
        OPEN(17,FILE='VELCNV7.OUT',POSITION='APPEND')
        OPEN(27,FILE='VELVCV7.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.8)THEN
        OPEN(18,FILE='VELCNV8.OUT',POSITION='APPEND')
        OPEN(28,FILE='VELVCV8.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.9)THEN
        OPEN(19,FILE='VELCNV9.OUT',POSITION='APPEND')
        OPEN(29,FILE='VELVCV9.OUT',POSITION='APPEND')
      ENDIF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      WRITE (LUN1,100)N,TIME
      WRITE (LUN2,100)N,TIME
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       LN=LNC(L)
       LS=LSC(L)
        DO K=1,KC
        VELN(K,NN)=50.*((U(L+1,K)+U(L,K))*COSC+(V(LN,K)
     &                                         +V(L,K))*SINC)
        VELT(K,NN)=-50.*((U(L+1,K)+U(L,K))*SINC-(V(LN,K)
     &                                          +V(L,K))*COSC)
        WZCORD(K,NN)=50.*(W(L,K)+W(L,K-1))+GI*ZZ(K)*(DTI*(P(L)-P1(L))
     &         +50.*(U(L+1,K)*(P(L+1)-P(L))*DXIU(L+1)
     &              +U(L,K)*(P(L)-P(L-1))*DXIU(L)
     &              +V(LN,K)*(P(LN)-P(L))*DYIV(LN)
     &              +V(L,K)*(P(L)-P(LS))*DYIV(L)))
     &         +50.*(1.-ZZ(K))*(U(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     &                         +U(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     &                         +V(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     &                         +V(L,K)*(BELV(L)-BELV(LS))*DYIV(L))
C    &         -50.*(1.-ZZ(K))*(U(L+1,K)*(HMP(L+1)-HMP(L))*DXIU(L+1)
C    &                         +U(L,K)*(HMP(L)-HMP(L-1))*DXIU(L)
C    &                         +V(LN,K)*(HMP(LN)-HMP(L))*DYIV(LN)
C    &                         +V(L,K)*(HMP(L)-HMP(LS))*DYIV(L))
C
        ENDDO
       ENDDO
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       ZETA=P(L)*GI-SBPLTV(1)*(HMP(L)+BELV(L))
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(VELN(K,NN),K=1,KC)
       WRITE(LUN2,250)(VELT(K,NN),K=1,KC)
       WRITE(LUN2,250)(WZCORD(K,NN),K=1,KC)
       ENDDO
      CLOSE(LUN1)
      CLOSE(LUN2)
      ENDDO
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10,F12.4)
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
