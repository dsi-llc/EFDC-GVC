C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BUDGET3 
C
C **  ADDED BY DON KINGERY, CH2M-HILL ON 15 OCTOBER 1996
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
C **  SUBROUTINES BUDGETN CALCULATE SEDIMENT BUDGET (TOTAL SEDIMENTS)
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
      VOLCONT=0.
      VOLMAST=0.
C
      DO L=2,LA
      VOLMIN=VOLMIN+QSUME(L)
      VOLCONT=VOLCONT+QSUME(L)
      ENDDO

      IF(ISTRAN(1).GE.1)THEN
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SAL1(L,K)
       ENDDO 
      ENDDO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
       DO K=1,KC
       SMASSIN=SMASSIN
     &           +MAX(QSS(K,NS),0.)*CQS(K,NS,1)
     &           +MIN(QSS(K,NS),0.)*SAL1(L,K) 
     &           +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,1)
     &           +MIN(QSERT(K,NQSTMP),0.)*SAL1(L,K)
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
       SMASSIN=SMASSIN+QCTLT(K,NCTL)*CONT(LU,K)
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
       SMASSIN=SMASSIN+
     &   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         SMASSIN=SMASSIN-
     &     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
       ENDIF
      ENDDO
C
      ENDIF

C
C    ACCUMULATE INTERNAL SOURCES AND SINKS FOR SED   (DLK 10/15)
C
      IF(ISTRAN(6).GE.1)THEN
C
      DO NN=1,NSED
      M=MSVSED(NN)
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SED(L,K,NN)
       ENDDO 
      ENDDO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       SEDIN=SEDIN
     &           +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     &           +MIN(QSS(K,NS),0.)*SED1(L,K,NN) 
     &           +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     &           +MIN(QSERT(K,NQSTMP),0.)*SED1(L,K,NN)
       IF(NN.EQ.1) VOLMAST=VOLMAST+QSERT(K,NQSTMP)+QSS(K,NS)
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
      SEDIN=SEDIN+QCTLT(K,NCTL)*CONT(LU,K)
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
       SEDIN=SEDIN+
     &   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         SEDIN=SEDIN-
     &     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
       ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C    ACCUMULATE INTERNAL SOURCES AND SINKS FOR SND   (DLK 10/15)
C
      IF(ISTRAN(7).GE.1)THEN
C
      DO NN=1,NSND
      M=MSVSND(NN)
       DO K=1,KC
        DO L=2,LC
        CONT(L,K)=SND(L,K,NN)
        ENDDO 
       ENDDO
C
       DO NS=1,NQSIJ
       L=LQS(NS)
       NQSTMP=NQSERQ(NS)
       NCSTMP=NCSERQ(NS,M)
        DO K=1,KC
        SEDIN=SEDIN
     &           +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     &           +MIN(QSS(K,NS),0.)*SND1(L,K,NN) 
     &           +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     &           +MIN(QSERT(K,NQSTMP),0.)*SND1(L,K,NN)
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
       SEDIN=SEDIN+QCTLT(K,NCTL)*CONT(LU,K)
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
        SEDIN=SEDIN+
     &    ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
        IF(LD.NE.1.OR.LD.NE.LC)THEN
          SEDIN=SEDIN-
     &     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
        ENDIF
       ENDDO
       ENDDO
C
      ENDIF
C
C     WRITE(6,600)VOLCONT,VOLMAST
  600 FORMAT(' VOLCON,VOLMAS = ',2E14.6)

      RETURN
      END
