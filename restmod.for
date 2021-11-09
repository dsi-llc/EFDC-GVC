C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTMOD
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
C **  SUBROUTINE RESTOUT WRITES A RESTART FILE
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      DIMENSION LIJMOD(100)
C
C**********************************************************************C
C
      OPEN(99,FILE='RESTART.OUT',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE='RESTART.OUT',STATUS='UNKNOWN')
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      WRITE(99,909)N,TIME
C
      OPEN(1,FILE='RESTMOD.INP',STATUS='UNKNOWN')
      READ(1,*)NIJMOD
      DO NNIJ=1,NIJMOD
      READ(1,*)ITMP,JTMP
      LIJMOD(NNIJ)=LIJ(ITMP,JTMP)
      ENDDO
      CLOSE(1)
C
      DO L=2,LA
      LSMOD=1
      DO NNIJ=1,NIJMOD
      IF(L.EQ.LIJMOD(NNIJ)) LSMOD=0
      ENDDO
      IF(LSMOD.EQ.1)THEN
C       WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L)
        WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
        WRITE(99,907)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
        WRITE(99,907)(U(L,K),K=1,KC)
        WRITE(99,907)(U1(L,K),K=1,KC)
        WRITE(99,907)(V(L,K),K=1,KC)
        WRITE(99,907)(V1(L,K),K=1,KC)
        WRITE(99,907)(QQ(L,K),K=0,KC)
        WRITE(99,907)(QQ1(L,K),K=0,KC)
        WRITE(99,907)(QQL(L,K),K=0,KC)
        WRITE(99,907)(QQL1(L,K),K=0,KC)
        WRITE(99,907)(DML(L,K),K=0,KC)
        IF(ISCO(1).EQ.1)THEN
          WRITE(99,907)(SAL(L,K),K=1,KC)
          WRITE(99,907)(SAL1(L,K),K=1,KC)
        ENDIF
        IF(ISCO(2).EQ.1)THEN
          WRITE(99,907)(TEM(L,K),K=1,KC)
          WRITE(99,907)(TEM1(L,K),K=1,KC)
        ENDIF
        IF(ISCO(3).EQ.1)THEN
          WRITE(99,907)(DYE(L,K),K=1,KC)
          WRITE(99,907)(DYE1(L,K),K=1,KC)
        ENDIF
        IF(ISCO(4).EQ.1)THEN
          WRITE(99,907)SEDB(L,1,1),(SED(L,K,1),K=1,KC)
          WRITE(99,907)SEDB1(L,1,1),(SED1(L,K,1),K=1,KC)
        ENDIF
        IF(ISCO(5).EQ.1)THEN
          WRITE(99,907)(SFL(L,K),K=1,KC)
          WRITE(99,907)(SFL2(L,K),K=1,KC)
        ENDIF
      ENDIF
      ENDDO
C
      DO M=1,5
      IF(ISCO(M).EQ.1)THEN
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDIF
      ENDDO
C
C**********************************************************************C
C
      CLOSE(99)
C
C     IF(ISRESTO.EQ.-2) CALL OUT3D
C
C**********************************************************************C
C
  906 FORMAT(5E17.8)
  907 FORMAT(13E17.8)
CMRM      906 FORMAT(1P,5E17.8)
CMRM      907 FORMAT(1P,13E17.8)
  908 FORMAT(12I10)
  909 FORMAT(I20,4X,F12.4)
C
C**********************************************************************C
C
      RETURN
      END
