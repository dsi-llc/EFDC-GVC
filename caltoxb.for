C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTOXB
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
c
      DIMENSION DIFTOXBW(LCM),PARTDIF(LCM,KBM)
C
      DIMENSION TOXBSMB(LCM,KBM),TOXSMB(LCM)
C
C**********************************************************************C
C
C
C **  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
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
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*STFPOCB(L,K,NS)
     &                           *TOXPARB(NS,NT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NX+NSED
            DO K=1,KB
              DO L=2,LA
                TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*STFPOCB(L,K,NS)
     &                           *TOXPARB(NS,NT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(ISTOC(NT).EQ.1.OR.ISTOC(NT).EQ.2)THEN
          NS=1+NSED+NSND
          DO K=1,KB
            DO L=2,LA
C0516              TOXPFB(L,K,NS,NT)=HBED(L,K)
               TOXPFB(L,K,NS,NT)=PORBED(L,K)*HBED(L,K)
     &                           *STDOCB(L,K)*TOXPARBC(1,NT)
            ENDDO
          ENDDO
        ENDIF
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
                IF(HBED(L,K).GT.0.0)THEN
                  TOXPFB(L,K,NS,NT)=TOXPFB(L,K,NS,NT)/
     &                   (PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
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
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
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
C**********************************************************************C
C
C **  ADVECT AND DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER 
C **  COLUMN LAYER
C
      IF(ISTRAN(5).GE.1)THEN
        DO NT=1,NTOX
C
C
C **  ADD PARTICLE MIXING AND SCALE PARTICAL DIFFUSION FOR SOLUTION
C
          DO K=1,KB
	      DO L=2,LA
	        PARTDIF(L,K)=0.0
            ENDDO
	    ENDDO
C
          IF(ISPMXZ(NT).GE.1)THEN
C
            CALL PARTMIX(NT)
            DO K=1,KB
              DO L=2,LA
                PARTDIF(L,K)=PARTMIXZ(L,K)
              ENDDO
            ENDDO
C
          ELSE
C
            DO L=2,LA
              DEPINBED=0.
              PARTDIF(L,KBT(L))=0.0
              IF(KBT(L).GT.2)THEN
                DO K=KBT(L),2,-1
                  KM=K-1
                  DEPINBED=DEPINBED+HBED(L,K)
                  PARTDIF(L,KM)=0.0
                  IF(DEPINBED.LT.DPDIFTOX(NT))THEN
                    PARTDIF(L,KM)=2.*PDIFTOX(NT)/
     &              (2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
                  ENDIF
                ENDDO
              ENDIF
              IF(KBT(L).EQ.2)THEN
                K=2
                KM=K-1
                DEPINBED=DEPINBED+HBED(L,K)
                PARTDIF(L,KM)=0.0
                IF(DEPINBED.LT.DPDIFTOX(NT))THEN
                  PARTDIF(L,KM)=2.*PDIFTOX(NT)/
     &            (2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
                ENDIF
              ENDIF
            ENDDO
C
          ENDIF
C 
C **  SETUP IMPLICIT TRI-DIAGONAL SYSTEM FOR PORE WATER ADVECTION AND DIFFUSION
C
          IF(IADTOXDP.EQ.1)THEN
C
C           CALL DOUBLE PRECISION SOLVER
            CALL ADTOXDP(NT,DELTI,PARTDIF)
C
C
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
             DO L=2,LA
             TMPBOTW=TOXFDFW(L,1,NT)+TOXCDFW(L,1,NT)
             TMPBOTB=TOXFDFB(L,KBT(L),NT)+TOXCDFB(L,KBT(L),NT)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXFDFW(L,1,NT)/TMPBOTW
             IF(TMPBOTB.GT.0.0) FACTORB=TOXFDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB
             ATOXFDFBW(L,NT)=ATOXFDFBW(L,NT)+TMPFLUX
             ATOXFDFBWP(L,NT)=ATOXFDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWN(L,NT)=ATOXFDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXFDFBWAD(L,NT)=ATOXFDFBWAD(L,NT)+TMPFLUX
             ATOXFDFBWADP(L,NT)=ATOXFDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWADN(L,NT)=ATOXFDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXCDFW(L,1,NT)/TMPBOTW 
             IF(TMPBOTB.GT.0.0) FACTORB=TOXCDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB   
             ATOXCDFBW(L,NT)=ATOXCDFBW(L,NT)+TMPFLUX
             ATOXCDFBWP(L,NT)=ATOXCDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWN(L,NT)=ATOXCDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXCDFBWAD(L,NT)=ATOXCDFBWAD(L,NT)+TMPFLUX
             ATOXCDFBWADP(L,NT)=ATOXCDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWADN(L,NT)=ATOXCDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             ENDDO   
             ENDIF
C
          ELSE
C
C **  ELIMINATE BED TO WATER COLUMN DIFFUSION FOR DRY CELLS
C
            DO L=2,LA
              DIFTOXBW(L)=0.0
            ENDDO
C
            DO L=2,LA
              IF(LMASKDRY(L)) DIFTOXBW(L)=DIFTOXS(NT)
            ENDDO
C
            DO L=2,LA
C
              DIFBWFAC=2./HBED(L,KBT(L))
              IF(ISDIFBW(NT).EQ.1)DIFBWFAC=1.0
              TOXBBALO(L)=0.
              KBTP1=KBT(L)+1
              KBTM1=KBT(L)-1
              ALOW(L,1)=0.
              CUPP(L,KBTP1)=0.
C
              DO K=1,KBTM1
                CUPP(L,K)=MIN(QWTRBED(L,K),0.)
     &          -(DIFTOX(NT)+PARTDIF(L,K))*(PORBED(L,K)+PORBED(L,K+1))/
     &                  (HBED(L,K)+HBED(L,K+1))
              ENDDO
              CUPP(L,KBT(L))=MIN(QWTRBED(L,KBT(L)),0.)
     &          -DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
C
              DO K=2,KBT(L)
                ALOW(L,K)=-MAX(QWTRBED(L,K-1),0.)
     &          -(DIFTOX(NT)+PARTDIF(L,K-1))*(PORBED(L,K-1)+PORBED(L,K))
     &                  /(HBED(L,K-1)+HBED(L,K))
              ENDDO
              ALOW(L,KBTP1)=-MAX(QWTRBED(L,KBT(L)),0.)
     &                 -DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
C
              DO K=1,KBT(L)
                BMNN(L,K)=DELTI*HBED(L,K)*PORBED(L,K)/
     &                  (1.-TOXPFTB(L,K,NT))
              ENDDO
              BMNN(L,KBTP1)=DELTI*DZC(1)*HP(L)/(1.-TOXPFTW(L,1,NT))
C
              BMNN(L,1)=BMNN(L,1)+MAX(QWTRBED(L,1),0.)
     &        +(DIFTOX(NT)+PARTDIF(L,1))*(PORBED(L,2)+PORBED(L,1))/
     &                  (HBED(L,2)+HBED(L,1))
              DO K=2,KBTM1
                BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED(L,K),0.)
     &         +(DIFTOX(NT)+PARTDIF(L,K))*(PORBED(L,K+1)+PORBED(L,K))/
     &                  (HBED(L,K+1)+HBED(L,K))
     &                         -MIN(QWTRBED(L,K-1),0.)
     &         +(DIFTOX(NT)+PARTDIF(L,K-1))*(PORBED(L,K-1)+PORBED(L,K))/
     &                  (HBED(L,K-1)+HBED(L,K))
              ENDDO
              K=KBT(L)
              BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED(L,K),0.)
     &           +DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
     &                         -MIN(QWTRBED(L,K-1),0.)
     &           +(DIFTOX(NT)+PARTDIF(L,K-1))*(PORBED(L,K-1)
     &                  +PORBED(L,K))/(HBED(L,K-1)+HBED(L,K))
              K=KBTP1
              BMNN(L,K)=BMNN(L,K)-MIN(QWTRBED(L,K-1),0.)
     &          +DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
C
              DO K=1,KBT(L)
                RRHS(L,K)=DELTI*TOXB(L,K,NT)
                TOXBBALO(L)=TOXBBALO(L)+TOXB(L,K,NT)
              ENDDO
              RRHS(L,1)=RRHS(L,1)+MAX(QWTRBED(L,0),0.)*CONGW(L,NT+4)
              RRHS(L,KBTP1)=DELTI*DZC(1)*HP(L)*TOX(L,1,NT)
              TOXWBALO(L)=DZC(1)*HP(L)*TOX(L,1,NT)
C
            ENDDO
C
C **  TRI-DIAGONAL SOLVER
C
            DO L=2,LA
              KBTP1=KBT(L)+1
              BETTMP=BMNN(L,1)
              TOXTMP(L,1)=RRHS(L,1)/BETTMP
              DO KK=2,KBTP1
                GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
                BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
                TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/
     &                     BETTMP
              ENDDO
              DO KK=KBT(L),1,-1
                TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
              ENDDO
            ENDDO
C
C **  CONVERT SCALED SOLUTION VARIABLES AND CALCULATE FINAL MASS
C
            DO L=2,LA
              TOXBBALN(L)=0.0
              KBTP1=KBT(L)+1
              DO K=1,KBT(L)
                TOXBSMB(L,K)=TOXB(L,K,NT)
                TOXB(L,K,NT)=HBED(L,K)*PORBED(L,K)*TOXTMP(L,K)/
     &                     (1.-TOXPFTB(L,K,NT))
                TOXBBALN(L)=TOXBBALN(L)+TOXB(L,K,NT)
              ENDDO
              TOXSMB(L)=TOX(L,1,NT)
              TOX(L,1,NT)=TOXTMP(L,KBTP1)/(1.-TOXPFTW(L,1,NT))
              TOXWBALN(L)=DZC(1)*HP(L)*TOX(L,1,NT)
            ENDDO
C
            IF(IADTOXCOR.EQ.0)THEN
C
              DO L=2,LA
                TADFLUX(L,NT)=DELTI*(DZC(1)*HP(L)*TOX(L,1,NT)
     &                       -TOXWBALO(L))
              ENDDO

C
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
             DO L=2,LA
             TMPBOTW=TOXFDFW(L,1,NT)+TOXCDFW(L,1,NT)
             TMPBOTB=TOXFDFB(L,KBT(L),NT)+TOXCDFB(L,KBT(L),NT)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXFDFW(L,1,NT)/TMPBOTW
             IF(TMPBOTB.GT.0.0) FACTORB=TOXFDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB
             ATOXFDFBW(L,NT)=ATOXFDFBW(L,NT)+TMPFLUX
             ATOXFDFBWP(L,NT)=ATOXFDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWN(L,NT)=ATOXFDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXFDFBWAD(L,NT)=ATOXFDFBWAD(L,NT)+TMPFLUX
             ATOXFDFBWADP(L,NT)=ATOXFDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWADN(L,NT)=ATOXFDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXCDFW(L,1,NT)/TMPBOTW 
             IF(TMPBOTB.GT.0.0) FACTORB=TOXCDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB   
             ATOXCDFBW(L,NT)=ATOXCDFBW(L,NT)+TMPFLUX
             ATOXCDFBWP(L,NT)=ATOXCDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWN(L,NT)=ATOXCDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXCDFBWAD(L,NT)=ATOXCDFBWAD(L,NT)+TMPFLUX
             ATOXCDFBWADP(L,NT)=ATOXCDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWADN(L,NT)=ATOXCDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             ENDDO   
             ENDIF
C
            ENDIF
C
C **  CORRECT MASS ERROR AND DETERMINE NET FLUX FROM BED TO WATER COLUMN
C
            IF(IADTOXCOR.EQ.1)THEN
C
              DO L=2,LA
                ERRT=TOXBBALN(L)+TOXWBALN(L)-TOXBBALO(L)-TOXWBALO(L)
                ERRB=ERRT*TOXBBALN(L)/(TOXBBALN(L)+TOXWBALN(L))
                ERRW=ERRT-ERRB
                TOX(L,1,NT)=TOX(L,1,NT)-ERRW/(DZC(1)*HP(L))
                DO K=1,KBT(L)
                  DERRB(K)=TOXB(L,K,NT)/TOXBBALN(L)
                ENDDO
                DO K=1,KBT(L)
                  TOXB(L,K,NT)=TOXB(L,K,NT)-DERRB(K)*ERRB
                ENDDO
C               TADFLUX(L,NT)=DELTI*(TOX(L,1,NT)-TOXWBALO(L))*HP(L)
                TADFLUX(L,NT)=DELTI*(DZC(1)*HP(L)*TOX(L,1,NT)
     &                       -TOXWBALO(L))
C               WRITE(8,8999)N,L,ERRT,ERRB,ERRW,TADFLUX(L,NT)
                ENDDO

C
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
             DO L=2,LA
             TMPBOTW=TOXFDFW(L,1,NT)+TOXCDFW(L,1,NT)
             TMPBOTB=TOXFDFB(L,KBT(L),NT)+TOXCDFB(L,KBT(L),NT)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXFDFW(L,1,NT)/TMPBOTW
             IF(TMPBOTB.GT.0.0) FACTORB=TOXFDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB
             ATOXFDFBW(L,NT)=ATOXFDFBW(L,NT)+TMPFLUX
             ATOXFDFBWP(L,NT)=ATOXFDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWN(L,NT)=ATOXFDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXFDFBWAD(L,NT)=ATOXFDFBWAD(L,NT)+TMPFLUX
             ATOXFDFBWADP(L,NT)=ATOXFDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWADN(L,NT)=ATOXFDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXCDFW(L,1,NT)/TMPBOTW 
             IF(TMPBOTB.GT.0.0) FACTORB=TOXCDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB   
             ATOXCDFBW(L,NT)=ATOXCDFBW(L,NT)+TMPFLUX
             ATOXCDFBWP(L,NT)=ATOXCDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWN(L,NT)=ATOXCDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXCDFBWAD(L,NT)=ATOXCDFBWAD(L,NT)+TMPFLUX
             ATOXCDFBWADP(L,NT)=ATOXCDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWADN(L,NT)=ATOXCDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             ENDDO   
             ENDIF
C
            ENDIF
C
            IF(IADTOXCOR.EQ.2)THEN
C
              DO L=2,LA
                ERRT=TOXBBALN(L)+TOXWBALN(L)-TOXBBALO(L)-TOXWBALO(L)
                ERRB=ERRT*TOXBBALN(L)/(TOXBBALN(L)+TOXWBALN(L))
                ERRW=ERRT-ERRB
                TOX(L,1,NT)=TOX(L,1,NT)-ERRW/(DZC(1)*HP(L))
                SUMDERRB=0.
                DO K=1,KBT(L)
                  DERRB(K)=ABS(TOXB(L,K,NT)-TOXBSMB(L,K))/HBED(L,K)
                  SUMDERRB=SUMDERRB+DERRB(K)
                ENDDO
                IF(SUMDERRB.GT.0.0)THEN
                  DO K=1,KBT(L)
                    DERRB(K)=DERRB(K)/SUMDERRB
                  ENDDO
                ENDIF
                DO K=1,KBT(L)
                  TOXB(L,K,NT)=TOXB(L,K,NT)-DERRB(K)*ERRB
                ENDDO
C               TADFLUX(L,NT)=DELTI*(TOX(L,1,NT)-TOXWBALO(L))*HP(L)
                TADFLUX(L,NT)=DELTI*(DZC(1)*HP(L)*TOX(L,1,NT)
     &                       -TOXWBALO(L))
                ENDDO

C
             IF(ISTOXALL.EQ.1.AND.NT.EQ.1)THEN
             DO L=2,LA
             TMPBOTW=TOXFDFW(L,1,NT)+TOXCDFW(L,1,NT)
             TMPBOTB=TOXFDFB(L,KBT(L),NT)+TOXCDFB(L,KBT(L),NT)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXFDFW(L,1,NT)/TMPBOTW
             IF(TMPBOTB.GT.0.0) FACTORB=TOXFDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB
             ATOXFDFBW(L,NT)=ATOXFDFBW(L,NT)+TMPFLUX
             ATOXFDFBWP(L,NT)=ATOXFDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWN(L,NT)=ATOXFDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXFDFBWAD(L,NT)=ATOXFDFBWAD(L,NT)+TMPFLUX
             ATOXFDFBWADP(L,NT)=ATOXFDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXFDFBWADN(L,NT)=ATOXFDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             FACTORW=0.0 
             FACTORB=0.0 
             IF(TMPBOTW.GT.0.0) FACTORW=TOXCDFW(L,1,NT)/TMPBOTW 
             IF(TMPBOTB.GT.0.0) FACTORB=TOXCDFB(L,KBT(L),NT)/TMPBOTB 
             TMPFLUX=DELT*DXYP(L)*MIN(TADFLUX(L,NT),0.0)*FACTORW
     &         +DELT*DXYP(L)*MAX(TADFLUX(L,NT),0.0)*FACTORB   
             ATOXCDFBW(L,NT)=ATOXCDFBW(L,NT)+TMPFLUX
             ATOXCDFBWP(L,NT)=ATOXCDFBWP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWN(L,NT)=ATOXCDFBWN(L,NT)+MIN(TMPFLUX,0.0)
             ATOXCDFBWAD(L,NT)=ATOXCDFBWAD(L,NT)+TMPFLUX
             ATOXCDFBWADP(L,NT)=ATOXCDFBWADP(L,NT)+MAX(TMPFLUX,0.0)
             ATOXCDFBWADN(L,NT)=ATOXCDFBWADN(L,NT)+MIN(TMPFLUX,0.0)
             ENDDO   
             ENDIF
C
            ENDIF
C
          ENDIF
C
        ENDDO
      ENDIF
C
 8999 FORMAT(' TAD  ',2I10,5E14.5,2F10.5)
C
C**********************************************************************C
C
      RETURN
      END
