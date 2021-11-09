C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CBALOD3 
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
C **  SUBROUTINES CBALOD CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C     DIMENSION CONT(LCM,KCM)
C
C**********************************************************************C
C
C **  ACCUMULATE INTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      VOLOUTO=VOLOUTO-QSUME(L)
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NQSIJ
      L=LQS(LL)
C
      PPEOUTO=PPEOUTO-QSS(K,LL)*G*( 0.5*(BELV(L)+BELV(L-1))
     &       +0.125*(HP(L)+H2P(L)+HP(L-1)+H2P(L-1))*(Z(K)+Z(K-1)) )
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(1).GE.1)THEN
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SAL(L,K)
       ENDDO 
      ENDDO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
       DO K=1,KC
       SALOUTO=SALOUTO
     &           -MAX(QSS(K,NS),0.)*CQS(K,NS,1)
     &           -MIN(QSS(K,NS),0.)*SAL(L,K) 
     &           -MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,1)
     &           -MIN(QSERT(K,NQSTMP),0.)*SAL(L,K)
       ENDDO
      ENDDO
C
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
       SALOUTO=SALOUTO+QCTLT(K,NCTL)*CONT(LU,K)
     &         -RQWD*QCTLT(K,NCTL)*CONT(LU,K)
       ENDDO
      ENDDO
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       SALOUTO=SALOUTO+
     &   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         SALOUTO=SALOUTO-
     &     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
       ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(3).GE.1)THEN
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=DYE(L,K)
       ENDDO 
      ENDDO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
       DO K=1,KC
       DYEOUTO=DYEOUTO
     &           -MAX(QSS(K,NS),0.)*CQS(K,NS,3)
     &           -MIN(QSS(K,NS),0.)*DYE(L,K) 
     &           -MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,3)
     &           -MIN(QSERT(K,NQSTMP),0.)*DYE(L,K)
       ENDDO
      ENDDO
C
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF(ID.EQ.0.AND.JD.EQ.0)THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      ENDIF
       DO K=1,KC
      DYEOUTO=DYEOUTO+QCTLT(K,NCTL)*CONT(LU,K)
     &        -RQWD*QCTLT(K,NCTL)*CONT(LU,K)
       ENDDO
      ENDDO
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       DYEOUTO=DYEOUTO+
     &   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         DYEOUTO=DYEOUTO-
     &     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
       ENDIF
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      RETURN
      END
