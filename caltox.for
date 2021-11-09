C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTOX
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
C **  SUBROUTINE CALSND CALCULATES NONCOHESIVER SEDIMENT SETTLING,
C **  DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
czzdiff      COMMON/PMC/CCSHEAR(LCM),DSTAR(NSTM),BDLDFACTOR(NSTM)  
czzdiff      REAL*4 TMPVAL,TMPVAL1,TMPVAL2  
czzdiff      REAL*4 DSTAR  
czzdiff      REAL*4 CCUSTAR(LCM)  
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
      DIMENSION TOXFPA(LCM)
C
C**********************************************************************C
C
C **  UPDATE FRACTION OF PARTICULATE ORGANIC CARBON IN BED
C
      IVAL=0
      DO NT=1,NTOX
        IF(ISTOC(NT).GE.2)IVAL=1
      ENDDO
C
      IF(IVAL.EQ.1.AND.ISTPOCB.EQ.4)THEN
c        IF(NSND.LE.3) CALL SETFPOCB3NC(1)
c         IF(NSND.EQ.4) CALL SETFPOCB4NC(1)
        CALL SETFPOCB(1)
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN WATER COLUMN
C
C **  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
C **  TOXPARW(NS,NT) = PARTITION COEFFICIENT IN WATER COLUMN
C **  TOXPFTW(L,K,NT) = TOTAL PARTICULATE FRACTION IN WATER COLUMN
C                       USED AS TEMPORARY VARIBLE IN THIS AND 
C                       FOLLOWING CODE BLOCK
C 
C      IF(TOXPARW(1,1).LT.1.E-6) RETURN
      IF(ISTRAN(5).GE.1)THEN
C
        DO NT=1,NTOX
          NSP2(NT)=NSED+NSND
          IF(ISTOC(NT).EQ.1) NSP2(NT)=NSP2(NT)+2
          IF(ISTOC(NT).EQ.2) NSP2(NT)=NSP2(NT)+1
        END DO
C
        DO NT=1,NTOX
          DO NS=1,NSP2(NT)
            DO K=1,KC
              DO L=2,LA
                TOXPFW(L,K,NS,NT)=0.
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          IF(ISTRAN(6).GE.1)THEN
C partition to cohesive
            DO NS=1,NSED
              IF(ITXPARW(NS,NT).EQ.0)THEN
                DO K=1,KC
                  DO L=2,LA
                   TOXPFW(L,K,NS,NT)=SED(L,K,NS)*STFPOCW(L,K,NS)
     &                              *TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
              IF(ITXPARW(NS,NT).EQ.1)THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KC
                  DO L=2,LA
                    TMPVAL=1.
                    IF(SED(L,K,NS).GT.0.) TMPVAL=SED(L,K,NS)**TMPEXP
                    TOXPFW(L,K,NS,NT)=TMPVAL*SED(L,K,NS)
     &                               *STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
C partition to noncohesive
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
              NS=NX+NSED
              IF(ITXPARW(NS,NT).EQ.0)THEN
                DO K=1,KC
                  DO L=2,LA
                    TOXPFW(L,K,NS,NT)=SND(L,K,NX)*STFPOCW(L,K,NS)
     &                               *TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
              IF(ITXPARW(NS,NT).EQ.1)THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KC
                  DO L=2,LA
                    TMPVAL=1.
                    IF(SND(L,K,NX).GT.0.) TMPVAL=SND(L,K,NX)**TMPEXP
                    TOXPFW(L,K,NS,NT)=TMPVAL*SND(L,K,NX)
     &                               *STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
C partition (complex to dissolved organic carbon)
cjmh040203          IF(ISTOC(NT).EQ.2)THEN
          IF(ISTOC(NT).EQ.1.OR.ISTOC(NT).EQ.2)THEN
            NS=1+NSED+NSND
            IF(ITXPARWC(1,NT).EQ.0)THEN
              DO K=1,KC
                DO L=2,LA
                  TOXPFW(L,K,NS,NT)=STDOCW(L,K)*TOXPARWC(1,NT)
                ENDDO
              ENDDO
            ENDIF
            IF(ITXPARWC(1,NT).EQ.1)THEN
              TMPEXP=CONPARWC(1,NT)
              DO K=1,KC
                DO L=2,LA
                  TMPVAL=1.
                  IF(STDOCW(L,K).GT.0.) TMPVAL=STDOCW(L,K)**TMPEXP
                  TOXPFW(L,K,NS,NT)=TMPVAL*STDOCW(L,K)*TOXPARWC(1,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C partition to particulate organic carbon
cjmh040203          IF(ISTOC(NT).EQ.3)THEN
          IF(ISTOC(NT).EQ.1)THEN
              NS=2+NSED+NSND
            IF(ITXPARWC(2,NT).EQ.0)THEN
              DO K=1,KC
                DO L=2,LA
                  TOXPFW(L,K,NS,NT)=STPOCW(L,K)*TOXPARWC(2,NT)
                ENDDO
              ENDDO
            ENDIF
            IF(ITXPARW(NS,NT).EQ.1)THEN
              TMPEXP=CONPARW(NS,NT)
              DO K=1,KC
                DO L=2,LA
                  TMPVAL=1.
                  IF(STPOCW(L,K).GT.0.) TMPVAL=STPOCW(L,K)**TMPEXP
                  TOXPFW(L,K,NS,NT)=TMPVAL*STPOCW(L,K)*TOXPARWC(2,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C !!  TOXPFTW IS TEMPORARILY USED TO STORE TOTAL SORBED
C
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXPFTW(L,K,NT)=0.
            ENDDO
          ENDDO
          DO NS=1,NSP2(NT)
            DO K=1,KC
              DO L=2,LA
                TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)+TOXPFW(L,K,NS,NT)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          DO NS=1,NSP2(NT)
            DO K=1,KC
              DO L=2,LA
                TOXPFW(L,K,NS,NT)=TOXPFW(L,K,NS,NT)/(1.+TOXPFTW(L,K,NT))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        NFD=NSED+NSND+1
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXFDFW(L,K,NT)=1./(1.+TOXPFTW(L,K,NT))
              TOXCDFW(L,K,NT)=TOXPFW(L,K,NFD,NT)
            ENDDO
          ENDDO
        ENDDO
C
C      OPEN(1,FILE='PARTWATER.OUT')
C      CLOSE(1,STATUS='DELETE')
C      OPEN(1,FILE='PARTWATER.OUT')
C      DO L=2,LA
C         WRITE(1,1907)IL(L),JL(L),(TOXPFW(L,1,NS,1),NS=1,NFD),
C     &                TOXPFTW(L,1,1)
C      ENDDO
C      CLOSE(1)
C
C      NVAL=NSED+NSND
C      OPEN(1,FILE='CARBWATER.OUT')
C      CLOSE(1,STATUS='DELETE')
C      OPEN(1,FILE='CARBWATER.OUT')
C      DO L=2,LA
C         WRITE(1,1907)IL(L),JL(L),STDOCW(L,1),
C     &                (STFPOCW(L,1,NS),NS=1,NVAL)
C      ENDDO
C      CLOSE(1)
C
      ENDIF
C
 1907 FORMAT(2I6,10E13.4)
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN SEDIMENT BED
C
C **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED
C **  TOXPARB(NS,NT) = PARTITION COEFFICIENT IN SEDIMENT BED
C **  TOXPFTB(L,NT) = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
C                       USED AS TEMPORARY VARIBLE IN THIS AND 
C                       FOLLOWING CODE BLOCK
C 
      IF(ISTRAN(5).GE.1)THEN
C
C
        DO NT=1,NTOX
          NSP2(NT)=NSED+NSND
          IF(ISTOC(NT).EQ.1) NSP2(NT)=NSP2(NT)+2
          IF(ISTOC(NT).EQ.2) NSP2(NT)=NSP2(NT)+1
        END DO
C
        DO NT=1,NTOX
          DO NS=1,NSP2(NT)
            DO K=1,KB
              DO L=2,LA
                TOXPFB(L,K,NS,NT)=0.
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
c partition to cohesives
          IF(ISTRAN(6).GE.1)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*STFPOCB(L,K,NS)
     &                             *TOXPARB(NS,NT)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
c partition to noncohesives
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
              NS=NX+NSED
              DO K=1,KB
                DO L=2,LA
                  TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*STFPOCB(L,K,NS)
     &                             *TOXPARB(NS,NT)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
c partition (complex) to DOC
          IF(ISTOC(NT).EQ.1.OR.ISTOC(NT).EQ.2)THEN
            NS=1+NSED+NSND
            DO K=1,KB
              DO L=2,LA
C0516                TOXPFB(L,K,NS,NT)=HBED(L,K)
                TOXPFB(L,K,NS,NT)=PORBED(L,K)*HBED(L,K)
     &                           *STDOCB(L,K)*TOXPARBC(1,NT)
              ENDDO
            ENDDO
          ENDIF
c partition to POC
          IF(ISTOC(NT).EQ.1)THEN
            NS=2+NSED+NSND
            DO K=1,KB
              DO L=2,LA
                TOXPFB(L,K,NS,NT)=HBED(L,K)
     &                           *STPOCB(L,K)*TOXPARBC(2,NT)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
C !!  TOXPFTB IS TEMPORARILY USED TO STORE TOTAL SORBED
C
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXPFTB(L,K,NT)=0.
            ENDDO
          ENDDO
          DO NS=1,NSP2(NT)
            DO K=1,KB
              DO L=2,LA
                TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          DO NS=1,NSP2(NT)
            DO K=1,KB
              DO L=2,LA
c                IF(SEDBALL(L,K).GT.0.0)THEN
c                  HBEDTMP=HBED(L,K)+1.E-9
                IF(HBED(L,K).GT.0.0)THEN
                  TOXPFB(L,K,NS,NT)=TOXPFB(L,K,NS,NT)/
     &                   (PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
C                write(8,8888)N,L,K,NS,TOXPFB(L,K,NS,NT)
                ELSE
                  TOXPFB(L,K,NS,NT)=0.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
        NFD=NSED+NSND+1
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              IF(HBED(L,K).GT.0.0)THEN
                TOXFDFB(L,K,NT)=PORBED(L,K)*HBED(L,K)
     &                       /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
                TOXCDFB(L,K,NT)=TOXPFB(L,K,NFD,NT)
C     &                       /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
              ELSE
                TOXFDFB(L,K,NT)=0.
                TOXCDFB(L,K,NT)=0.
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C      OPEN(1,FILE='PARTBED.OUT')
C      CLOSE(1,STATUS='DELETE')
C      OPEN(1,FILE='PARTBED.OUT')
C        DO L=2,LA
C         K=KBT(L)
C         WRITE(1,1907)IL(L),JL(L),(TOXPFB(L,K,NS,1),NS=1,NFD),
C     &                TOXPFTB(L,K,1)
C      ENDDO
C      CLOSE(1)
C
C      NVAL=NSED+NSND
C      OPEN(1,FILE='CARBBED.OUT')
C      CLOSE(1,STATUS='DELETE')
C      OPEN(1,FILE='CARBBED.OUT')
C      DO L=2,LA
C         K=KBT(L)
C         WRITE(1,1907)IL(L),JL(L),STDOCB(L,K),
C     &                            (STFPOCB(L,K,NS),NS=1,NVAL)
C      ENDDO
C      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE PARTICULATE TOXIC CONTAMINANT SETTLING 
C **  AND BED EXCHANGE FLUXES
C
C **  TOXF(L,KC,NT) = TOXIC CONTAMINANT SETTLING AND BED EXCHANGE
C                        FLUX.  USED AS TEMPORARY VARIABLE IN THIS 
C                        AND FOLLOWING CODE BLOCK
C
      IF(ISTRAN(5).GE.1)THEN
C
        DO NT=1,NTOX
          DO K=0,KC
            DO L=2,LA
              TOXF(L,K,NT)=0.
            ENDDO
          ENDDO
        ENDDO
C
        IF(KC.GE.2)THEN
          DO NT=1,NTOX
            IF(ISTRAN(6).GE.1)THEN
c particule cohesive flux
              DO NS=1,NSED
                IF(ITXPARW(NS,NT).EQ.0)THEN
                  DO K=1,KS
                    DO L=2,LA
                      TOXF(L,K,NT)=TOXF(L,K,NT)+SEDF(L,K,NS)
     &                            *STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                    ENDDO
                  ENDDO
                ENDIF
                IF(ITXPARW(NS,NT).EQ.1)THEN
                  TMPEXP=CONPARW(NS,NT)
                  DO K=1,KS
                    DO L=2,LA
                      TMPVAL=1.
                      IF(SED(L,K+1,NS).GT.0.) TMPVAL=
     &                                        SED(L,K+1,NS)**TMPEXP
                      TOXF(L,K,NT)=TOXF(L,K,NT)+TMPVAL*SEDF(L,K,NS)
     &                            *STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
            IF(ISTRAN(7).GE.1)THEN
c particule noncohesive flux
              DO NX=1,NSND
                NS=NX+NSED
                IF(ITXPARW(NS,NT).EQ.0)THEN
                  DO K=1,KS
                    DO L=2,LA
                      TOXF(L,K,NT)=TOXF(L,K,NT)+SNDF(L,K,NX)
     &                            *STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                    ENDDO
                  ENDDO
                ENDIF
                IF(ITXPARW(NS,NT).EQ.1)THEN
                  TMPEXP=CONPARW(NS,NT)
                  DO K=1,KS
                    DO L=2,LA
                      TMPVAL=1.
                      IF(SND(L,K+1,NX).GT.0.) TMPVAL=
     &                                        SND(L,K+1,NX)**TMPEXP
                      TOXF(L,K,NT)=TOXF(L,K,NT)+TMPVAL*SNDF(L,K,NX)
     &                            *STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C
        IF(KC.GE.2)THEN
          DO NT=1,NTOX
            DO K=1,KS
              DO L=2,LA
                TOXF(L,K,NT)=TOXF(L,K,NT)/(1.+TOXPFTW(L,K+1,NT))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        DO NT=1,NTOX
          IF(ISTRAN(6).GE.1)THEN
c particule cohesive flux
            DO NS=1,NSED
              IF(ITXPARW(NS,NT).EQ.0)THEN
                DO L=2,LA
                  TOXF(L,0,NT)=TOXF(L,0,NT)+MIN(SEDF(L,0,NS),0.)
     &                        *STFPOCW(L,1,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDIF
              IF(ITXPARW(NS,NT).EQ.1)THEN
                TMPEXP=CONPARW(NS,NT)
                DO L=2,LA
                  TMPVAL=1.
                  IF(SED(L,1,NS).GT.0.) TMPVAL=SED(L,1,NS)**TMPEXP
                  TOXF(L,0,NT)=TOXF(L,0,NT)+TMPVAL*MIN(SEDF(L,0,NS),0.)
     &                        *STFPOCW(L,1,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
c particule noncohesive flux
            DO NX=1,NSND
              NS=NX+NSED
              IF(ITXPARW(NS,NT).EQ.0)THEN
                DO L=2,LA
                  SNDFEFF=SNDF(L,0,NX)-SNDFBL(L,NX)
                  TOXF(L,0,NT)=TOXF(L,0,NT)+
     &            MIN(SNDFEFF,0.)*STFPOCW(L,1,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDIF
              IF(ITXPARW(NS,NT).EQ.1)THEN
                TMPEXP=CONPARW(NS,NT)
                DO L=2,LA
                  TMPVAL=1.
                  IF(SND(L,1,NX).GT.0.) TMPVAL=SND(L,1,NX)**TMPEXP
                  SNDFEFF=SNDF(L,0,NX)-SNDFBL(L,NX)
                  TOXF(L,0,NT)=TOXF(L,0,NT)+
     &            TMPVAL*MIN(SNDFEFF,0.)*STFPOCW(L,1,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDIF
           ENDDO
          ENDIF
        ENDDO
C
        DO NT=1,NTOX
          DO L=2,LA
            TOXF(L,0,NT)=TOXF(L,0,NT)/(1.+TOXPFTW(L,1,NT))
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          DO L=2,LA
            TOXFB(L,KBT(L),NT)=0.
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          IF(ISTRAN(6).GE.1)THEN
c particule cohesive flux
            DO NS=1,NSED
              DO L=2,LA
                TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)+
     &          MAX(SEDF(L,0,NS),0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(NS,NT)
              ENDDO
            ENDDO
          ENDIF
          IF(ISTRAN(7).GE.1)THEN
c particule noncohesive flux
            DO NX=1,NSND
              NS=NX+NSED
              DO L=2,LA
                SNDFEFF=SNDF(L,0,NX)-SNDFBL(L,NX)
                TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)
     &              +MAX(SNDFEFF,0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(NS,NT)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
C---START BANK EROSION-------------------------
C
        IF(ISBKERO.GE.1)THEN
C
          DO NT=1,NTOX
            DO L=2,LA
              TOXFBEBKB(L,NT)=0.
              TOXFBECHB(L,NT)=0.
              TOXFBECHW(L,NT)=0.
            ENDDO
          ENDDO
C
          DO NT=1,NTOX
            IF(ISTRAN(6).GE.1)THEN
c particule cohesive flux
              DO NS=1,NSED
                DO NP=1,NBEPAIR
                  LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
                  LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
                  KBANK=KBT(LBANK)
                  TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)
     &              +SEDFBEBKB(NP,NS)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                  TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)
     &              +SEDFBECHB(NP,NS)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                  TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)
     &              +SEDFBECHW(NP,NS)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                ENDDO
              ENDDO
            ENDIF
            IF(ISTRAN(7).GE.1)THEN
c particule noncohesive flux
              DO NX=1,NSND
                NS=NX+NSED
                DO NP=1,NBEPAIR
                  LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
                  LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
                  KBANK=KBT(LBANK)
                  TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)
     &              +SNDFBEBKB(NP,NX)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                  TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)
     &              +SNDFBECHB(NP,NX)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                  TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)
     &              +SNDFBECHW(NP,NX)*STFPOCB(LBANK,KBANK,NS)
     &                *TOXPARB(NS,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
C
        ENDIF
C
C---END BANK EROSION-------------------------
C
        DO NT=1,NTOX
          DO L=2,LA
c            IF(SEDBALL(L,KBT(L)).GT.0.0)THEN 
c              HBEDTMP=HBED(L,KBT(L))+1.E-9
            IF(HBED(L,KBT(L)).GT.0.0)THEN 
              TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)/
     &           (PORBED(L,KBT(L))*HBED(L,KBT(L))+TOXPFTB(L,KBT(L),NT))
            ELSE
              TOXFB(L,KBT(L),NT)=0.
            ENDIF
          ENDDO
        ENDDO
C
C---START BANK EROSION-------------------------
C
        IF(ISBKERO.GE.1)THEN
C
        DO NT=1,NTOX
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK=KBT(LBANK)
            IF(HBED(LBANK,KBANK).GT.0.0)THEN 
              TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)/
     &           (PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)
     &           +TOXPFTB(LBANK,KBANK,NT))
              TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)/
     &           (PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)
     &           +TOXPFTB(LBANK,KBANK,NT))
              TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)/
     &           (PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)
     &           +TOXPFTB(LBANK,KBANK,NT))
            ELSE
              TOXFBEBKB(NP,NT)=0.
              TOXFBECHB(NP,NT)=0.
              TOXFBECHW(NP,NT)=0.
            ENDIF
          ENDDO
        ENDDO
C
        ENDIF
C
C---END BANK EROSION-------------------------
C
C ** DIAGNOSTICS OF FLUX
C
        IF(ISDTXBUG.EQ.1)THEN
          IF(N.EQ.1)THEN
            OPEN(2,FILE='TOXFLX.DIA')
            CLOSE(2,STATUS='DELETE')
            OPEN(2,FILE='TOXFLX.DIA')
            DO L=2,LA
              WRITE(2,2222)IL(L),JL(L),HBED(L,KBT(L)),
     &        TOXPFTB(L,KBT(L),1),TOXFB(L,KBT(L),1),TOXF(L,0,1)
            ENDDO
            CLOSE(2)
          ENDIF
        ENDIF
C
        IF(ISDTXBUG.EQ.2)THEN
            OPEN(2,FILE='TOXFLX.DIA',ACCESS='APPEND')
            DO L=2,LA
              WRITE(2,2223)N,IL(L),JL(L),HBED(L,KBT(L)),
     &   SEDFBEBKB(L,1),SEDFBECHB(L,1),SEDFBECHW(L,1)
               WRITE(2,2224)
     &   TOXFBEBKB(L,1),TOXFBECHB(L,1),TOXFBECHW(L,1)  

            ENDDO
            CLOSE(2)
        ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN AND BED
C **  NOTING THAT TO THIS POINT TOXPFTW AND TOXPFTB HAVE BEEN USED
C **  TO TEMPORILY SORTED THE SORBED PORTION
C
      IF(ISTRAN(5).GE.1)THEN
C
        NFD=NSED+NSND+1
C
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXPFTW(L,K,NT)=( TOXPFTW(L,K,NT)/(1.+TOXPFTW(L,K,NT)) )
     &                        -TOXPFW(L,K,NFD,NT)
            ENDDO
          ENDDO
        ENDDO
C
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
c              IF(SEDBALL(L,K).GT.0.0)THEN
c                HBEDTMP=HBED(L,K)+1.E-9
              IF(HBED(L,K).GT.0.0)THEN
                TOXPFTB(L,K,NT)=( TOXPFTB(L,K,NT)
     &                      /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT)) )
     &                      -TOXPFB(L,K,NFD,NT)
              ELSE
                TOXPFTB(L,K,NT)=0.
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
      ENDIF
C
C
C**********************************************************************C
C
C FIXED FOR BED LOAD JMH 5/22/02
c
C **  DETERIME TOXIC FLUX FROM BED LOAD SORBED MATERIAL
C
      DO NT=1,NTOX
        TOXBLB(NT)=0.0
        DO L=2,LA
          TOXFBL(L,NT)=0.0
        ENDDO
      ENDDO
C 
      DO NT=1,NTOX
        DO NX=1,NSND
          NS=NX+NSED
          DO L=2,LA
            IF(LMASKDRY(L))THEN
              LE=L+1
              LW=l-1
              LN=LNC(L)
              LS=LSC(L)
              TMPTOXC=0.0
              TMPTOXE=0.0
              TMPTOXW=0.0
              TMPTOXN=0.0
              TMPTOXS=0.0
C
              K=KBT(L)
              IF(SNDB(L,K,NX).GT.0.0) THEN
                TMPTOXC=TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)*HBED(L,K)
     &                 /SNDB(L,K,NX)
              ENDIF
              IF(SUB(LE).GT.0.5) THEN
                K=KBT(LE)
                IF(SNDB(LE,K,NX).GT.0.0) THEN
                  TMPTOXE=TOXPFB(LE,K,NS,NT)*TOXB(LE,K,NT)*HBED(LE,K)
     &                 /SNDB(LE,K,NX)
                ENDIF
              ENDIF
              IF(SUB(L).GT.0.5) THEN
                K=KBT(LW)
                IF(SNDB(LW,K,NX).GT.0.0) THEN
                  TMPTOXW=TOXPFB(LW,K,NS,NT)*TOXB(LW,K,NT)*HBED(LW,K)
     &                 /SNDB(LW,K,NX)
                ENDIF
              ENDIF
              IF(SVB(LN).GT.0.5) THEN
                K=KBT(LN)
                IF(SNDB(LN,K,NX).GT.0.0) THEN
                  TMPTOXN=TOXPFB(LN,K,NS,NT)*TOXB(LN,K,NT)*HBED(LN,K)
     &                 /SNDB(LN,K,NX)
                ENDIF
              ENDIF
              IF(SVB(L).GT.0.5) THEN
                K=KBT(LS)
                IF(SNDB(LS,K,NX).GT.0.0) THEN
                  TMPTOXS=TOXPFB(LS,K,NS,NT)*TOXB(LS,K,NT)*HBED(LS,K)
     &                 /SNDB(LS,K,NX)
                ENDIF
              ENDIF
C
C ABOVE FROM HOUSATONIC, CHECK FOR EQUIVALENCE OF BELOW
C
C	        K=KBT(L)
C	        IF(SNDB(L,K,NX).GT.0.0) THEN
C	          TMPTOXC=TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)*HBED(L,K)
C     &                 /SNDB(L,K,NX)
C	        ENDIF
C			K=KBT(LE)
C	        IF(SUB(LE).GT.0.5.AND.SNDB(LE,K,NX).GT.0.0) THEN
C	          TMPTOXE=TOXPFB(LE,K,NS,NT)*TOXB(LE,K,NT)*HBED(LE,K)
C     &                 /SNDB(LE,K,NX)
C	        ENDIF
C	        K=KBT(LW)
C	        IF(SUB(L).GT.0.5.AND.SNDB(LW,K,NX).GT.0.0) THEN
C	          TMPTOXW=TOXPFB(LW,K,NS,NT)*TOXB(LW,K,NT)*HBED(LW,K)
C     &                 /SNDB(LW,K,NX)
C	        ENDIF
C	        K=KBT(LN)
C	        IF(SVB(LN).GT.0.5.AND.SNDB(LN,K,NX).GT.0.0) THEN
C	          TMPTOXN=TOXPFB(LN,K,NS,NT)*TOXB(LN,K,NT)*HBED(LN,K)
C     &                 /SNDB(LN,K,NX)
C	        ENDIF
C	        K=KBT(LS)
C	        IF(SVB(L).GT.0.5.AND.SNDB(LS,K,NX).GT.0.0) THEN
C	          TMPTOXS=TOXPFB(LS,K,NS,NT)*TOXB(LS,K,NT)*HBED(LS,K)
C     &                 /SNDB(LS,K,NX)
C	        ENDIF
C
              TOXFBL(L,NT)=TOXFBL(L,NT)+DXYIP(L)*(
     &        SUB(L+1)*MAX(QSBDLDX(L+1,NX),0.)*TMPTOXC
     &       +SUB(L+1)*MIN(QSBDLDX(L+1,NX),0.)*TMPTOXE
     &       -SUB(L  )*MAX(QSBDLDX(L  ,NX),0.)*TMPTOXW
     &       -SUB(L  )*MIN(QSBDLDX(L  ,NX),0.)*TMPTOXC
     &       +SVB(LN )*MAX(QSBDLDY(LN ,NX),0.)*TMPTOXC
     &       +SVB(LN )*MIN(QSBDLDY(LN ,NX),0.)*TMPTOXN
     &       -SVB(L  )*MAX(QSBDLDY(L  ,NX),0.)*TMPTOXS
     &       -SVB(L  )*MIN(QSBDLDY(L  ,NX),0.)*TMPTOXC )

              IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
                ATOXPFBLX(L,NX,NT)=ATOXPFBLX(L,NX,NT)
     &                   +DELT*SUB(L  )*MAX(QSBDLDX(L  ,NX),0.)*TMPTOXW
     &                   +DELT*SUB(L  )*MIN(QSBDLDX(L  ,NX),0.)*TMPTOXC
                  ATOXPFBLY(L,NX,NT)=ATOXPFBLY(L,NX,NT)
     &                   +DELT*SVB(L  )*MAX(QSBDLDY(L  ,NX),0.)*TMPTOXS
     &                   +DELT*SVB(L  )*MIN(QSBDLDY(L  ,NX),0.)*TMPTOXC
              ENDIF

            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
 8822 FORMAT(3I5,E14.5)
C
      IF(IS2TIM.GE.1) THEN
      IF(ISBAL.GE.1)THEN
      IF(NSBDLDBC.GT.0) THEN 
        DO NSB=1,NSBDLDBC
          LUTMP=LSBLBCU(NSB)
          LDTMP=LSBLBCD(NSB)
          IF(LDTMP.EQ.0) THEN
            DO NT=1,NTOX
              DO NX=1,NSND
                NS=NX+NSED
                TMPTOXC=0.0
                 K=KBT(LUTMP)
                IF(SNDB(LUTMP,K,NX).GT.0.0) THEN
                  TMPTOXC=TOXPFB(LUTMP,K,NS,NT)*TOXB(LUTMP,K,NT)
     &                     *HBED(LUTMP,K)/SNDB(LUTMP,K,NX)
                ENDIF
                  TOXBLB(NT)=TOXBLB(NT)+QSBDLDOT(LUTMP,NX)*TMPTOXC
                TOXFBL(LUTMP,NT)=TOXFBL(LUTMP,NT)+DXYIP(LUTMP)* 
     &                       QSBDLDOT(LUTMP,NX)*TMPTOXC
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DO NSB=1,NSBDLDBC
          LUTMP=LSBLBCU(NSB)
          LDTMP=LSBLBCD(NSB)
          IF(LDTMP.NE.0) THEN
            DO NT=1,NTOX
              DO NX=1,NSND
                NS=NX+NSED
                TMPTOXC=0.0
                 K=KBT(LUTMP)
                IF(SNDB(LUTMP,K,NX).GT.0.0) THEN
                  TMPTOXC=TOXPFB(LUTMP,K,NS,NT)*TOXB(LUTMP,K,NT)
     &                     *HBED(LUTMP,K)/SNDB(LUTMP,K,NX)
                ENDIF
                TOXFBL(LUTMP,NT)=TOXFBL(LUTMP,NT)+DXYIP(LUTMP)* 
     &                       QSBDLDOT(LUTMP,NX)*TMPTOXC
                TOXFBL(LDTMP,NT)=TOXFBL(LDTMP,NT)-DXYIP(LDTMP)* 
     &                       QSBDLDOT(LUTMP,NX)*TMPTOXC
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      ENDIF
      ENDIF
C
      DO NT=1,NTOX
        TOXFBLT(NT)=0.
        DO L=2,LA
          TOXFBLT(NT)=TOXFBLT(NT)+DXYP(L)*TOXFBL(L,NT)
        ENDDO
      ENDDO
C
 8862 FORMAT('N,NX,SNDFBLTOT,QSBLLDXDY   ='2I5,2E14.5)
 8899 FORMAT('N,TOXFBLT(NT),TOXBLB(NT)='I5,2E14.5)
C
C END FIXED FOR BED LOAD JMH 5/22/02
C
C**********************************************************************C
C
C **  ADJUST TOXIC FLUXES ACROSS WATER COLUMN - BED INTERFACE TO 
C **  INCLUDE WATER ENTRAINMENT AND EXPULSION ASSOCIATED WITH 
C **  DEPOSITION AND RESUSPENSION
C 
      DO NT=1,NTOX
        DO L=2,LA
          TMPVAL=( MIN(QSBDTOP(L),0.0)+MIN(QWBDTOP(L),0.0) )
          TOXF(L,0,NT)=TOXF(L,0,NT)+TMPVAL*(1.-TOXPFTW(L,1,NT))
        ENDDO
      ENDDO
C
      DO NT=1,NTOX
        DO L=2,LA
          K=KBT(L)
          TMPVAL=( MAX(QSBDTOP(L),0.0)+MAX(QWBDTOP(L),0.0) )/HBED(L,K)
          TOXFB(L,K,NT)=TOXFB(L,K,NT)+TMPVAL*(1.-TOXPFTB(L,K,NT))
        ENDDO
      ENDDO
C
C
C---START BANK EROSION-------------------------
C
        IF(ISBKERO.GE.1)THEN
C
        DO NT=1,NTOX
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK=KBT(LBANK)
            TMPVAL=(1.-TOXPFTB(LBANK,KBANK,NT))/HBED(LBANK,KBANK)
            TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)
     &           +TMPVAL*(QSBDTOPBEBKB(NP)+QWBDTOPBEBKB(NP))
            TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)
     &           +TMPVAL*(QSBDTOPBECHB(NP)+QWBDTOPBECHB(NP))
            TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)
     &           +TMPVAL*(QSBDTOPBECHW(NP)+QWBDTOPBECHW(NP))
          ENDDO
        ENDDO
C
        IF(ISDTXBUG.EQ.2)THEN
            OPEN(2,FILE='TOXFLX.DIA',ACCESS='APPEND')
            DO L=2,LA
              WRITE(2,2224)
     &   TOXFBEBKB(L,1),TOXFBECHB(L,1),TOXFBECHW(L,1)
            ENDDO
            CLOSE(2)
        ENDIF
C
        ENDIF
C
C---END BANK EROSION-------------------------
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.1)THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(1)
            AA11=WVEL-TOXF(L,0,NT)
            AA12=-TOXFB(L,KBT(L),NT)
            AA21=TOXF(L,0,NT)
            AA22=DELTI+TOXFB(L,KBT(L),NT)
            BB11=WVEL*TOX(L,1,NT)
cjmh216            BB22=DELTI*TOXB1(L,KBT(L),NT)
C            BB22=DELTI*TOXB(L,KBT(L),NT)
C FIXED FOR BED LOAD JMH 5/22/02
            BB22=DELTI*TOXB(L,KBT(L),NT)-TOXFBL(L,NT)
C END FIXED FOR BED LOAD JMH 5/22/02
            DETI=1./(AA11*AA22-AA12*AA21)
            TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
            TOXB1(L,KBT(L),NT)=TOXB(L,KBT(L),NT)
            TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
            TOXF(L,0,NT)=TOXF(L,0,NT)*TOX(L,1,NT)
            TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)
            TOXTOTMP=HP(L)*TOX(L,1,NT)+TOXB(L,KBT(L),NT)
c         WRITE(8,676)N,L,TOX(L,1,NT),TOXB(L,KBT(L),NT),TOXTOTMP,
c     $      TOX1(L,1,NT),TOXB1(L,KBT(L),NT),TOXF(L,0,NT),
c     &      TOXFB(L,KBT(L),NT)
C adjust wc and bed toxic consistent with flux
             TOX(L,1,NT)=(BB11+TOXF(L,0,NT)+TOXFB(L,KBT(L),NT))/WVEL
             TOXB(L,KBT(L),NT)=DELT*(BB22-TOXF(L,0,NT)
     &                                  -TOXFB(L,KBT(L),NT))
c         WRITE(8,677)N,L,TOX(L,1,NT),TOXB(L,KBT(L),NT)
C end adjust wc and bed toxic consistent with flux
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
             DO NS=1,NSED+NSND+2
             TMPFLUX=DELT*DXYP(L)*TOXF(L,0,NT)*TOXPFW(L,1,NS,NT)
     &         +DELT*DXYP(L)*TOXFB(L,KBT(L),NT)*TOXPFB(L,KBT(L),NS,NT)
             ATOXPFBW(L,NS,NT)=ATOXPFBW(L,NS,NT)+TMPFLUX
             ATOXPFBWP(L,NS,NT)=ATOXPFBWP(L,NS,NT)+MAX(TMPFLUX,0.0)
             ATOXPFBWN(L,NS,NT)=ATOXPFBWN(L,NS,NT)+MIN(TMPFLUX,0.0)
             ENDDO
             TMPFLUX=DELT*DXYP(L)*TOXF(L,0,NT)*TOXFDFW(L,1,NT)
     &         +DELT*DXYP(L)*TOXFB(L,KBT(L),NT)*TOXFDFB(L,KBT(L),NT)
             ATOXFDFBW(L,NT)=ATOXFDFBW(L,NT)+TMPFLUX
             ATOXFDFBWP(L,NT)=ATOXFDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWN(L,NT)=ATOXFDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             TMPFLUX=DELT*DXYP(L)*TOXF(L,0,NT)*TOXCDFW(L,1,NT)
     &         +DELT*DXYP(L)*TOXFB(L,KBT(L),NT)*TOXCDFB(L,KBT(L),NT)    
             ATOXCDFBW(L,NT)=ATOXCDFBW(L,NT)+TMPFLUX
             ATOXCDFBWP(L,NT)=ATOXCDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWN(L,NT)=ATOXCDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ENDIF
          ENDDO
C
C  676 FORMAT('N,I,J,TW,TB,TB1,TF,TFB=',3I5,5E13.4)
  676 FORMAT('N,L,T,TB,TT,T1.TB1,F,FB=',2I5,8E13.4)
  677 FORMAT('N,L,T,TB               =',2I5,8E13.4)
C
C----------------------------------------------------------------------C
C
          IF(ISBKERO.GE.1)THEN
            DO NP=1,NBEPAIR
              LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
              LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
              KBANK=KBT(LBANK)
              KCHAN=KBT(LCHAN)
              TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)
     &                           *TOXB(LBANK,KBANK,NT)
              TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)
     &                           *TOXB(LBANK,KBANK,NT)
              TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)
     &                           *TOXB(LBANK,KBANK,NT)
                TOXB(LBANK,KBANK,NT)=TOXB(LBANK,KBANK,NT)
     &          -DELT*TOXFBEBKB(NP,NT)
                TOXB(LCHAN,KCHAN,NT)=TOXB(LCHAN,KCHAN,NT)
     &          -DELT*TOXFBECHB(NP,NT)
                TOX(LCHAN,1,NT)=TOX(LCHAN,1,NT)
     &          +DELT*DZIC(1)*HPI(LCHAN)*TOXFBECHW(NP,NT)
c
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
               ATOXFBEBKB(LBANK,NT)=ATOXFBEBKB(LBANK,NT)
     &                           -DELT*DXYP(LBANK)*TOXFBEBKB(NP,NT)
                 ATOXFBECHB(LCHAN,NT)=ATOXFBECHB(LCHAN,NT)
     &                           -DELT*DXYP(LCHAN)*TOXFBECHB(NP,NT)
               ATOXFBECHW(LCHAN,NT)=ATOXFBECHW(LCHAN,NT)
     &                           +DELT*DXYP(LCHAN)*TOXFBECHW(NP,NT)
             ENDIF               
C
            ENDDO    
          ENDIF
C
          IF(ISDTXBUG.EQ.2)THEN
            OPEN(2,FILE='TOXFLX.DIA',ACCESS='APPEND')
            DO L=2,LA
              WRITE(2,2224)TOXB(L,KBT(L),1),TOX(L,1,1)
            ENDDO
            CLOSE(2)
          ENDIF
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.2)THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
          K=2
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(K)
            CLEFT=WVEL-TOXF(L,K-1,NT)
            CRIGHT=WVEL*TOX(L,K,NT)
            TOX(L,K,NT)=CRIGHT/CLEFT
          ENDDO
C
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(1)
            AA11=WVEL-TOXF(L,0,NT)
            AA12=-TOXFB(L,KBT(L),NT)
            AA21=TOXF(L,0,NT)
            AA22=DELTI+TOXFB(L,KBT(L),NT)
            BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,KC,NT)
C            BB22=DELTI*TOXB(L,KBT(L),NT)
C FIXED FOR BED LOAD JMH 5/22/02
            BB22=DELTI*TOXB(L,KBT(L),NT)-TOXFBL(L,NT)
C END FIXED FOR BED LOAD JMH 5/22/02
            DETI=1./(AA11*AA22-AA12*AA21)
            TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
            TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     &                      +S2TL*TOXB1(L,KBT(L),NT)
            TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
            TOXF(L,0,NT)=TOXF(L,0,NT)*TOX(L,1,NT)
            TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)
C adjust wc and bed toxic consistent with flux
             TOX(L,1,NT)=(BB11+TOXF(L,0,NT)+TOXFB(L,KBT(L),NT))/WVEL
             TOXB(L,KBT(L),NT)=DELT*(BB22-TOXF(L,0,NT)
     &                                  -TOXFB(L,KBT(L),NT))
C end adjust wc and bed toxic consistent with flux
          ENDDO
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.3)THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
          K=3
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(K)
            CLEFT=WVEL-TOXF(L,K-1,NT)
            CRIGHT=WVEL*TOX(L,K,NT)
            TOX(L,K,NT)=CRIGHT/CLEFT
          ENDDO
C
          K=2
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(K)
            CLEFT=WVEL-TOXF(L,K-1,NT)
            CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
            TOX(L,K,NT)=CRIGHT/CLEFT
          ENDDO
C
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(1)
            AA11=WVEL-TOXF(L,0,NT)
            AA12=-TOXFB(L,KBT(L),NT)
            AA21=TOXF(L,0,NT)
            AA22=DELTI+TOXFB(L,KBT(L),NT)
            BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,KC-1,NT)
C            BB22=DELTI*TOXB(L,KBT(L),NT)
C FIXED FOR BED LOAD JMH 5/22/02
            BB22=DELTI*TOXB(L,KBT(L),NT)-TOXFBL(L,NT)
C END FIXED FOR BED LOAD JMH 5/22/02
            DETI=1./(AA11*AA22-AA12*AA21)
            TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
            TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     &                      +S2TL*TOXB1(L,KBT(L),NT)
            TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
            TOXF(L,0,NT)=TOXF(L,0,NT)*TOX(L,1,NT)
            TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)
C adjust wc and bed toxic consistent with flux
             TOX(L,1,NT)=(BB11+TOXF(L,0,NT)+TOXFB(L,KBT(L),NT))/WVEL
             TOXB(L,KBT(L),NT)=DELT*(BB22-TOXF(L,0,NT)
     &                                  -TOXFB(L,KBT(L),NT))
C end adjust wc and bed toxic consistent with flux
          ENDDO
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC.GE.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(KC.GE.3)K1P1=2
      IF(ISTRAN(5).GE.1.AND.KC.GE.3)THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
          K=KC
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(K)
            CLEFT=WVEL-TOXF(L,K-1,NT)
            CRIGHT=WVEL*TOX(L,K,NT)
            TOX(L,K,NT)=CRIGHT/CLEFT
          ENDDO
C
          DO K=KS,2,-1
            DO L=2,LA
              WVEL=DELTI*HP(L)*DZC(K)
              CLEFT=WVEL-TOXF(L,K-1,NT)
              CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
              TOX(L,K,NT)=CRIGHT/CLEFT
            ENDDO
          ENDDO
C
          DO L=2,LA
            WVEL=DELTI*HP(L)*DZC(1)
            AA11=WVEL-TOXF(L,0,NT)
            AA12=-TOXFB(L,KBT(L),NT)
            AA21=TOXF(L,0,NT)
            AA22=DELTI+TOXFB(L,KBT(L),NT)
            BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,K1P1,NT)
C            BB22=DELTI*TOXB(L,KBT(L),NT)
C FIXED FOR BED LOAD JMH 5/22/02
            BB22=DELTI*TOXB(L,KBT(L),NT)-TOXFBL(L,NT)
C END FIXED FOR BED LOAD JMH 5/22/02
            DETI=1./(AA11*AA22-AA12*AA21)
            TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
            TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     &                      +S2TL*TOXB1(L,KBT(L),NT)
cjmh216            TOXBTMP=DETI*(AA11*BB22-AA21*BB11)
cjmh216            IF(TOXBTMP.LT.0.0)THEN
cjmh216              TOXBTMP=TOXB1(L,KBT(L),NT)
cjmh216              TOX(L,1,NT)=TOXS(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)/WVEL
cjmh216            ENDIF
cjmh216            TOXB(L,KBT(L),NT)=TOXBTMP
            TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
            TOXF(L,0,NT)=TOXF(L,0,NT)*TOX(L,1,NT)
            TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)
c         WRITE(8,676)N,L,TOX(L,1,NT),TOXB(L,KBT(L),NT),TOXTOTMP,
c     $      TOX1(L,1,NT),TOXB1(L,KBT(L),NT),TOXF(L,0,NT),
c     &      TOXFB(L,KBT(L),NT)
C adjust wc and bed toxic consistent with flux
             TOX(L,1,NT)=(BB11+TOXF(L,0,NT)+TOXFB(L,KBT(L),NT))/WVEL
             TOXB(L,KBT(L),NT)=DELT*(BB22-TOXF(L,0,NT)
     &                                  -TOXFB(L,KBT(L),NT))
c         WRITE(8,677)N,L,TOX(L,1,NT),TOXB(L,KBT(L),NT)
C end adjust wc and bed toxic consistent with flux
          ENDDO
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C      
C**********************************************************************C
C
C **  ADD PARENT TO ACTIVE LAYER TOXIC TRANSPORT
C
      IF(ISNDAL.GE.2)THEN
C
      IF(ISTRAN(5).GE.1)THEN
      DO NT=1,NTOX
C
      DO L=2,LA
        TOXFPA(L)=0.0
      ENDDO
C
      DO NS=1,NSED
        DO L=2,LA
          KTOPTP=KBT(L)
          KTOPM1=KBT(L)-1
          FTPOS=0.0
          FTNEG=0.0
          IF(SEDFPA(L,NS).GT.0.0.AND.SEDB(L,KTOPM1,NS).GT.0.0)
     &      FTPOS=SEDFPA(L,NS)*TOXPFB(L,KTOPM1,NS,NT)/SEDB(L,KTOPM1,NS)
          IF(SEDFPA(L,NS).LT.0.0.AND.SEDB(L,KTOPTP,NS).GT.0.0)
     &      FTNEG=SEDFPA(L,NS)*TOXPFB(L,KTOPTP,NS,NT)/SEDB(L,KTOPTP,NS)
          TOXFPA(L)=TOXFPA(L)+FTPOS*TOXB(L,KTOPM1,NT)
     &                       +FTNEG*TOXB(L,KTOPTP,NT)
        ENDDO
      ENDDO
C
      DO NX=1,NSND
        NS=NX+NSED
        DO L=2,LA
          KTOPTP=KBT(L)
          KTOPM1=KBT(L)-1
          FTPOS=0.0
          FTNEG=0.0
          IF(SNDFPA(L,NX).GT.0.0.AND.SNDB(L,KTOPM1,NX).GT.0.0)
     &      FTPOS=SNDFPA(L,NX)*TOXPFB(L,KTOPM1,NS,NT)/SNDB(L,KTOPM1,NX)
          IF(SNDFPA(L,NX).LT.0.0.AND.SNDB(L,KTOPTP,NX).GT.0.0)
     &      FTNEG=SNDFPA(L,NX)*TOXPFB(L,KTOPTP,NS,NT)/SNDB(L,KTOPTP,NX)
          TOXFPA(L)=TOXFPA(L)+FTPOS*TOXB(L,KTOPM1,NT)
     &                       +FTNEG*TOXB(L,KTOPTP,NT)
        ENDDO
      ENDDO
C
      DO L=2,LA
        KTOPTP=KBT(L)
        KTOPM1=KBT(L)-1
        FTPOS=0.0
        FTNEG=0.0
        IF(QWATPA(L).GT.0.0)
     &    FTPOS=QSSDPA(L)*(1.-TOXPFTB(L,KTOPM1,NT))
        IF(QWATPA(L).LT.0.0)
     &      FTNEG=QSSDPA(L)*(1.-TOXPFTB(L,KTOPTP,NT))
        TOXFPA(L)=TOXFPA(L)+FTPOS*TOXB(L,KTOPM1,NT)
     &                       +FTNEG*TOXB(L,KTOPTP,NT)
      ENDDO
C
      DO L=2,LA
        KTOPTP=KBT(L)
        KTOPM1=KBT(L)-1
        TOXB(L,KTOPTP,NT)=TOXB(L,KTOPTP,NT)+DELT*TOXFPA(L)
        TOXB(L,KTOPM1,NT)=TOXB(L,KTOPM1,NT)-DELT*TOXFPA(L)
      ENDDO
C
      ENDDO
      ENDIF
      ENDIF
C      
C**********************************************************************C
C
 8888 FORMAT(4I5,7E14.5)
 2222 FORMAT(2I5,7E14.5)
 2223 FORMAT(3I5,7E14.5)
 2224 FORMAT(29X,7E14.5)
C      
C**********************************************************************C
C
      RETURN
      END
