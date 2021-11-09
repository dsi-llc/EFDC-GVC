C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T3B(IBALSTDT)
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
      IF(IBALSTDT.EQ.1)THEN
        DO L=2,LA
          WVOLOUT=WVOLOUT-DTSED*QMORPH(L)
          BVOLOUT=BVOLOUT+DTSED*QMORPH(L)
          VOLMORPH2T=VOLMORPH2T+DTSED*QMORPH(L)
C          WVOLOUT=WVOLOUT-DTSED*(QWTRBEDA(L)+QWTRBEDA1(L)
C    &                          +QWTRBED(L,KBT(L)))
        ENDDO
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
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT)
     &                    +DTSED*DXYP(LBANK)*TOXFBEBKB(LBANK,NT)
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT)
     &                    +DTSED*DXYP(LCHAN)*TOXFBECHB(LCHAN,NT)
          ENDDO
        ENDIF
c
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
C SEDFLUX2T(NSX) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED
C   TO WATER COLUMN
C
      IF(IBALSTDT.EQ.1)THEN
c
        DO L=2,LA
          SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(L)*SEDF(L,0,NSX)
        ENDDO
c
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
             SEDFLUX2T(NSX)=SEDFLUX2T(NSX)
     &                    +DTSED*DXYP(LBANK)*SEDFBEBKB(LBANK,NSX)
            SEDFLUX2T(NSX)=SEDFLUX2T(NSX)
     &                    +DTSED*DXYP(LCHAN)*SEDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
c
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
c
        DO L=2,LA
          SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(L)*(SNDF(L,0,NSX)
     &                -SNDFBL(L,NSX))
          SNDFBL2T(NSX)=SNDFBL2T(NSX)+DTSED*DXYP(L)*SNDFBL(L,NSX)
        ENDDO
c
        IF(ISBKERO.GE.1)THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
             SNDFLUX2T(NSX)=SNDFLUX2T(NSX)
     &                    +DTSED*DXYP(LBANK)*SNDFBEBKB(LBANK,NSX)
            SNDFLUX2T(NSX)=SNDFLUX2T(NSX)
     &                    +DTSED*DXYP(LCHAN)*SNDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
c
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
