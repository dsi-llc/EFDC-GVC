C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BANKEROSED
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE BANKEROSED CALCULATES SEDIMENT TRANSPORT DUE TO BANK
C **  EROSION.  TRANSPORT IS FROM BANK BED CELL TO CHANNEL BED AND 
C **  WATER COLUMN CELLS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      COMMON/SSEDTOX1/ CTMPDRY(LCM),CSHIELDS50(LCM),
     &                 USTAR(LCM),UCELLCTR(LCM),VCELLCTR(LCM),
     &                 QWBDTOP(LCM),QSBDTOP(LCM),ZETATOP(LCM),
     &                 ZBEDGT(LCM),QSBDLDP(LCM),RBPSBL(LCM),
     &                 TOXBBALO(LCM),TOXBBALN(LCM),TOXWBALO(LCM),
     &                 ROUSE(LCM),ZEQ(LCM),ZEQI(LCM),ZEQD(LCM),
     &                 ZEQDI(LCM),SNDEQ(LCM),SNDEQB(LCM),
     &                 TOXWBALN(LCM),FACSUSL(LCM),FACBEDL(LCM),
     &                 PEXP(LCM,NSNM),PHID(LCM,NSNM),
     &                 SEDFPA(LCM,NSCM),SNDFPA(LCM,NSNM)
C
      COMMON/SSEDTOX1A/ CBEDTOTAL(LCM),QCELLCTR(LCM),HGDH(LCM),
     &                  FRACCOH(LCM,KBM),FRACNON(LCM,KBM),USTARSED(LCM),
     &                  USTARSND(LCM),QWATPA(LCM),QSSDPA(LCM)
C
      COMMON/SSEDTOX2/ CDECAYW(LCM,KCM),CDECAYB(LCM,KBM),STRSE(LCM,KBM),
     &                 HYDCN(LCM,KBM),COEFK(LCM,KBM),COEFSK(LCM,KBM),
     &                 PRESE(LCM,KBM),PRESH(LCM,KBM),PREST(LCM,KBM),
     &                 STRST(LCM,KBM),DZBTR1(LCM,KBM),SEDDIA50(LCM,KBM),
     &                 SEDBALL(LCM,KBM),ZBEDC(LCM,KBM),SGSM1(LCM,KBM),
     &                 DSTRSE(LCM,KBM),DZBTR(LCM,KBM),STRSEM(LCM,KBM),
     &                 SEDDIA90(LCM,KBM),SEDGEOSTD(LCM,KBM),
     &                 SEDDIAGS(LCM,KBM),ZOTOTAL(LCM),ZOGRAIN(LCM)
C
      COMMON/SSEDTOX3/ ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),CUPP(LCM,KBM+1),
     &                 RRHS(LCM,KBM+1),TOXTMP(LCM,KBM+1),
     &                 GAMTMP(LCM,KBM+1),ACOEF(LCM,0:KBM),
     &                 QCOEF(LCM,0:KBM),ZBEDG(LCM,0:KBM)
C
      COMMON/SSEDTOX4/ WSETA(LCM,0:KSM,NSTM),SEDS(LCM,KCM,NSCM),
     &                 SNDS(LCM,KCM,NSNM),TOXS(LCM,KCM,NTXM),
     &                 SEDBS(LCM,KBM,NSCM),SNDBS(LCM,KBM,NSNM),
     &                 TOXBS(LCM,KBM,NTXM)
C
      COMMON/SSEDTOX5/ NSP2(NTXM),DERRB(KBM),STRESSS(0:KSM)
C
      COMMON/SSEDTOX6/ DELT,DELTI,DSEDGMM,FOURDPI,SEDMDGM,S2TL,S3TL,
     &                 CORDT,DIASED,GPDIASED,BEDEX
C
      COMMON/SSEDTOX7/ ISTL,IS2TL,ISUD
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C   
 7777 FORMAT(2I5,10E12.4)
C
C**********************************************************************C
C
C  INITIALIZE BANK EROSION VARIABLES
C
      DO L=1,LC
        QSBDTOPBEBKB(L)=0.
        QSBDTOPBECHB(L)=0.
        QSBDTOPBECHW(L)=0.
        QWBDTOPBEBKB(L)=0.
        QWBDTOPBECHB(L)=0.
        QWBDTOPBECHW(L)=0.
      ENDDO
C
      DO NS=1,NSED
      DO L=1,LC
        SEDFBEBKB(L,NS)=0.
        SEDFBECHB(L,NS)=0.
        SEDFBECHW(L,NS)=0.
      ENDDO
      ENDDO
C
      DO NS=1,NSND
      DO L=1,LC
        SNDFBEBKB(L,NS)=0.
        SNDFBECHB(L,NS)=0.
        SNDFBECHW(L,NS)=0.
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO L=1,LC
        TOXFBEBKB(L,NT)=0.
        TOXFBECHB(L,NT)=0.
        TOXFBECHW(L,NT)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C  LOAD SEDIMENT FLUXES
C
      IF(ISTRAN(6).GT.0)THEN
        DO NS=1,NSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            BKEROBKB=VFRBED(LBANK,K,NS)*FBESER(NP)
     &               *BESERT(NBESERN(NP))
            BKEROCHW=FWCBESERT(NBESERN(NP))
            BKEROCHB=1.-FWCBESERT(NBESERN(NP))
C**********************************************************************C
C HQI Change to restrict bank erosion to cells with available sediment
C RM, 01/08/07
C            SEDFBEBKB(NP,NS)=BKEROBKB*DXYIP(LBANK)
C            SEDFBECHB(NP,NS)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
C            SEDFBECHW(NP,NS)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
            IF((DELT*BKEROBKB*DXYIP(LBANK)).GT.SEDB(LBANK,KBT(LBANK),NS)
     &        ) THEN
              SEDFBEBKB(NP,NS)=0.
              SEDFBECHB(NP,NS)=0.
              SEDFBECHW(NP,NS)=0.
            ELSE              
              SEDFBEBKB(NP,NS)=BKEROBKB*DXYIP(LBANK)
              SEDFBECHB(NP,NS)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
              SEDFBECHW(NP,NS)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
            ENDIF
C End HQI Change
C**********************************************************************C
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GT.0)THEN
        DO NX=1,NSND
        NS=NSED+NX
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            BKEROBKB=VFRBED(LBANK,K,NS)*FBESER(NP)
     &               *BESERT(NBESERN(NP))
            BKEROCHW=FWCBESERT(NBESERN(NP))
            BKEROCHB=1.-FWCBESERT(NBESERN(NP))
C**********************************************************************C
C HQI Change to restrict bank erosion to cells with available sediment
C RM, 01/08/07
c     SNDFBEBKB(NP,NX)=BKEROBKB*DXYIP(LBANK)
c     SNDFBECHB(NP,NX)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
c     SNDFBECHW(NP,NX)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
            IF((DELT*BKEROBKB*DXYIP(LBANK)).GT.SNDB(LBANK,KBT(LBANK),NX)
     &        ) THEN
              SNDFBEBKB(NP,NX)=0.
              SNDFBECHB(NP,NX)=0.
              SNDFBECHW(NP,NX)=0.
            ELSE              
              SNDFBEBKB(NP,NX)=BKEROBKB*DXYIP(LBANK)
              SNDFBECHB(NP,NX)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
              SNDFBECHW(NP,NX)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
            ENDIF
C End HQI Change
C**********************************************************************C
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C  UPDATE BED AND WATER COLUMN SEDIMENT CONCENTRATION
C
      IF(ISTRAN(6).GT.0)THEN
        DO NS=1,NSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            WVEL=DELT*HPI(LCHAN)*DZIC(1)
            SEDB(LBANK,KBT(LBANK),NS)=SEDB(LBANK,KBT(LBANK),NS)
     &                               -DELT*SEDFBEBKB(NP,NS)
            SEDB(LCHAN,KBT(LCHAN),NS)=SEDB(LCHAN,KBT(LCHAN),NS)
     &                               -DELT*SEDFBECHB(NP,NS)
            SED(LCHAN,1,NS)=SED(LCHAN,1,NS)+WVEL*SEDFBECHW(NP,NS)
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GT.0)THEN
        DO NS=1,NSND
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            WVEL=DELT*HPI(LCHAN)*DZIC(1)
            SNDB(LBANK,KBT(LBANK),NS)=SNDB(LBANK,KBT(LBANK),NS)
     &                               -DELT*SNDFBEBKB(NP,NS)
            SNDB(LCHAN,KBT(LCHAN),NS)=SNDB(LCHAN,KBT(LCHAN),NS)
     &                               -DELT*SNDFBECHB(NP,NS)
            SND(LCHAN,1,NS)=SND(LCHAN,1,NS)+WVEL*SNDFBECHW(NP,NS)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FROM BANK
C
C  COHESIVE
C
      IF(ISTRAN(6).GT.0)THEN
      IF(IBMECH.EQ.1.AND.SEDVRDT.LT.0.0)THEN
C
        DO NS=1,NSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP)
     &      +0.001*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
            QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP)
     &      +0.001*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
        ENDDO
        ENDDO
C
      ELSE
C
        DO NS=1,NSED
          DSEDGMM=1./(1.E6*SSG(NS))
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP)
     &        +DSEDGMM*SEDFBEBKB(NP,NS)
            QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP)
     &        +DSEDGMM*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)
          ENDDO
        ENDDO
C
      ENDIF
      ENDIF
C
C  NONCOHESIVE
C
      IF(ISTRAN(7).GT.0)THEN
      IF(IBMECH.EQ.1.AND.SEDVRDT.LT.0.0)THEN
C
        DO NS=1,NSND
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP)
     &      +0.001*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
            QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP)
     &      +0.001*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
          ENDDO
        ENDDO
C
      ELSE
C
        DO NS=1,NSND
          DSEDGMM=1./(1.E6*SSG(NS+NSED))
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            K=KBT(LBANK)
            QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP)
     &        +DSEDGMM*SNDFBEBKB(NP,NS)
            QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP)
     &        +DSEDGMM*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)
          ENDDO
        ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FOR CHANNEL BED AND 
C  WATER COLUMN
C
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        BKEROCHW=FWCBESERT(NBESERN(NP))
        BKEROCHB=1.-FWCBESERT(NBESERN(NP))
        TMPVAL=DXYP(LBANK)*DXYIP(LCHAN)
        QSBDTOPBECHB(NP)=-TMPVAL*BKEROCHB*QSBDTOPBEBKB(NP)
        QWBDTOPBECHB(NP)=-TMPVAL*BKEROCHB*QWBDTOPBEBKB(NP)
        QSBDTOPBECHW(NP)=TMPVAL*BKEROCHW*QSBDTOPBEBKB(NP)
        QWBDTOPBECHW(NP)=TMPVAL*BKEROCHW*QWBDTOPBEBKB(NP)
      ENDDO
C
C**********************************************************************C
C
      RETURN
      END
