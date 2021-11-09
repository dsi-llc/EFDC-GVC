C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T3(IBALSTDT)
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
C 05/01/2002        john hamrick       05/01/2002       john hamrick
C  subroutine added for 2 time-level balances including sed,snd,tox
C 06/06/2002        john hamrick       06/06/2002       john hamrick
C  modified snd mass balance with respect to bed load outflow
C 06/07/2002        john hamrick       06/07/2002       john hamrick
C  added qdwaste to water mass balance
C----------------------------------------------------------------------C
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      END IF
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
        VOLOUT=VOLOUT-DELT*(QSUME(L)+QGW(L)-QDWASTE(L))
      ENDDO
C
      DO L=2,LA
        WVOLOUT=WVOLOUT-DELT*(QSUME(L)-QDWASTE(L))
      ENDDO
C
      IF(IBALSTDT.EQ.1)THEN
        DO L=2,LA
          WVOLOUT=WVOLOUT-DTSED*QMORPH(L)
C          WVOLOUT=WVOLOUT-DTSED*(QWTRBEDA(L)+QWTRBEDA1(L)
C    &                          +QWTRBED(L,KBT(L)))
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NQSIJ
      L=LQS(LL)
C
      PPEOUT=PPEOUT-DELT*QSS(K,LL)*G*( 0.5*(BELV(L)+BELV(L-1))
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
       SALOUT=SALOUT
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,1)
     &           -DELT*MIN(QSS(K,NS),0.)*SAL(L,K) 
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,1)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*SAL(L,K)
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
       DO K=1,KC
       SALOUT=SALOUT+DELT*QCTLT(K,NCTL)*CONT(LU,K)
       ENDDO
      ENDIF
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
       SALOUT=SALOUT+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         SALOUT=SALOUT-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))
     &      +QWRSERT(NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
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
      NCSTMP=NCSERQ(NS,3)
       DO K=1,KC
       DYEOUT=DYEOUT
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,3)
     &           -DELT*MIN(QSS(K,NS),0.)*DYE(L,K) 
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,3)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*DYE(L,K)
       DYEOUT2T=DYEOUT2T
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,3)
     &           -DELT*MIN(QSS(K,NS),0.)*DYE(L,K) 
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,3)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*DYE(L,K)
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
       DO K=1,KC
      DYEOUT=DYEOUT+DELT*QCTLT(K,NCTL)*CONT(LU,K)
      DYEOUT2T=DYEOUT2T+DELT*QCTLT(K,NCTL)*CONT(LU,K)
       ENDDO
      ENDIF
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
       DYEOUT=DYEOUT+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       DYEOUT2T=DYEOUT2T+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         DYEOUT=DYEOUT-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))
     &      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
         DYEOUT2T=DYEOUT2T-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))
     &      +QWRSERT(NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
       ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(5).GE.1)THEN
C
      DO NT=1,NTOX
      M=MSVTOX(NT)
cSO   WRITE(8,*)'nt m ',NT,M
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=TOX(L,K,NT)
       ENDDO 
      ENDDO
C
C  TOXOUT2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
c  TO WATER COLUMN VOLUME SOURCES AND SINKS
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       TOXOUT2T(NT)=TOXOUT2T(NT)
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     &           -DELT*MIN(QSS(K,NS),0.)*TOX(L,K,NT) 
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*TOX(L,K,NT)
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
       DO K=1,KC
       TOXOUT2T(NT)=TOXOUT2T(NT)+DELT*QCTLT(K,NCTL)*TOX(LU,K,NT)
       ENDDO
      ENDIF
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
       TOXOUT2T(NT)=TOXOUT2T(NT)+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
         TOXOUT2T(NT)=TOXOUT2T(NT)-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))
     &      +QWRSERT(NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
       ENDIF
      ENDDO
C
C  TOXBLB2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
c  DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
C
      IF(IBALSTDT.EQ.1)THEN
        IF(NSBDLDBC.GT.0) THEN 
          TOXBLB2T(NT)=TOXBLB2T(NT)+DTSED*TOXBLB(NT)
        ENDIF
C
C  TOXFLUXW2T(NT) IS WATER COLUMN SIDE TOXIC FLUX DUE TO SUSPENDED LOAD
C    (POSITIVE INTO WATER COLUMN)
C  TOXFLUXB2T(NT) IS BED SIDE TOXIC FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)
C  TADFLUX2T(NT) IS PORE WATER ADVECTION+DIFFUSION FLUX (POSITIVE INTO WATER COLUMN)
C  TOXFBL2T(NT) IS NET TOXIC FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT
C    (SHOULD EQUAL TOXBLB2T(NT)
C
        DO L=2,LA
c        TOXFTMP=TOXF(L,0,NT)+TOXFB(L,KBT(L),NT)
c        TOXFLUXW2T(NT)=TOXFLUXW2T(NT)+DELT*DXYP(L)*TOXFTMP
c        TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DELT*DXYP(L)*(TOXFTMP
c     &     +TOXFBL(L,NT))
c        TADFLUX2T(NT)=TADFLUX2T(NT)+DELT*DXYP(L)*TADFLUX(L,NT)
          TOXFLUXW2T(NT)=TOXFLUXW2T(NT)+DTSED*DXYP(L)*TOXF(L,0,NT)
          TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DTSED*DXYP(L)*TOXFB(L,KBT(L),NT)
          TADFLUX2T(NT)=TADFLUX2T(NT)+DTSED*DXYP(L)*TADFLUX(L,NT)
C        TOXFBL2T(NT)=TOXFBL2T(NT)+DELT*DXYP(L)*TOXFBL(L,NT)
        ENDDO
C
        TOXFBL2T(NT)=TOXFBL2T(NT)+DTSED*TOXFBLT(NT)
C
      ENDIF
C
      ENDDO
c
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(6).GE.1)THEN
C
      DO NSX=1,NSED
      M=MSVSED(NSX)
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SED(L,K,NSX)
       ENDDO 
      ENDDO
C
C SEDOUT2T(NSX) IS IS NET COHESIVE MASS GOING OUT OF DOMAIN DUE
C   TO WATER COLUMN VOLUME SOURCES AND SINKS
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       SEDOUT2T(NSX)=SEDOUT2T(NSX)
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     &           -DELT*MIN(QSS(K,NS),0.)*SED(L,K,NSX)
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*SED(L,K,NSX)
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
       DO K=1,KC
       SEDOUT2T(NSX)=SEDOUT2T(NSX)+DELT*QCTLT(K,NCTL)*CONT(LU,K)
       ENDDO
      ENDIF
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
       SEDOUT2T(NSX)=SEDOUT2T(NSX)+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
       SEDOUT2T(NSX)=SEDOUT2T(NSX)-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))
     &      +QWRSERT(NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
       ENDIF
      ENDDO
C
C SEDFLUX2T(NSX) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED
C   TO WATER COLUMN
C
      IF(IBALSTDT.EQ.1)THEN
      DO L=2,LA
        SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(L)*SEDF(L,0,NSX)
      ENDDO
      ENDIF
C
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(7).GE.1)THEN
C
      DO NSX=1,NSND
      M=MSVSND(NSX)
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SND(L,K,NSX)
       ENDDO 
      ENDDO
C
C  SNDOUT2T(NSX) IS NET NONCOHESIVE MASS GOING OUT OF DOMAIN DUE
c  TO WATER COLUMN VOLUME SOURCES AND SINKS
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       SNDOUT2T(NSX)=SNDOUT2T(NSX)
     &           -DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     &           -DELT*MIN(QSS(K,NS),0.)*SND(L,K,NSX)
     &           -DELT*MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     &           -DELT*MIN(QSERT(K,NQSTMP),0.)*SND(L,K,NSX)
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
        DO K=1,KC
          SNDOUT2T(NSX)=SNDOUT2T(NSX)+DELT*QCTLT(K,NCTL)*CONT(LU,K)
        ENDDO
      ENDIF
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
       SNDOUT2T(NSX)=SNDOUT2T(NSX)+
     &   DELT*( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC)THEN
       SNDOUT2T(NSX)=SNDOUT2T(NSX)-
     &     DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))
     &      +QWRSERT(NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
       ENDIF
      ENDDO
C
C  SBLOUT2T(NSX) IS NET NONCOHESIVE SEDIMENT MASS GOING OUT OF DOMAIN DUE
c  DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
C
      IF(IBALSTDT.EQ.1)THEN
      IF(NSBDLDBC.GT.0) THEN 
        DO NSB=1,NSBDLDBC
          LUTMP=LSBLBCU(NSB)
          LDTMP=LSBLBCD(NSB)
          IF(LDTMP.EQ.0) THEN
              SBLOUT2T(NSX)=SBLOUT2T(NSX)+
     &        DTSED*QSBDLDOT(LUTMP,NSX)
          ENDIF
        ENDDO
      ENDIF
      ENDIF
C
C  SNDFLUX2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX DUE TO SUSPENDED LOAD 
C    (POSITIVE INTO WATER COLUMN)
C  SNDFBL2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX FROM BED ASSOCIATED WITH 
C    BED LOAD TRANSPORT (SHOULD EQUAL SBLOUT2T(NSX))
C
c      TMPVAL=SNDFBL2T(NSX)
      IF(IBALSTDT.EQ.1)THEN
      DO L=2,LA
        SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(L)*(SNDF(L,0,NSX)
     &                -SNDFBL(L,NSX))
        SNDFBL2T(NSX)=SNDFBL2T(NSX)+DTSED*DXYP(L)*SNDFBL(L,NSX)
      ENDDO
      ENDIF
c      TMPVAL=SNDFBL2T(NSX)-TMPVAL
c       WRITE(8,800)N,NSX,SNDFBL2T(NSX),TMPVAL
C
      ENDDO
C
      ENDIF
C
  800 FORMAT('N,NS,SNDFBL2T,DEL',2I5,2E14.5)
C
C----------------------------------------------------------------------C
C
      RETURN
      END
