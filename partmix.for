C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE PARTMIX(NT)
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
      DIMENSION PARTMIXAVG(KBM)
C
C**********************************************************************C
C
      DO K=1,KB
        DO L=2,LA
          PARTMIXZ(L,K)=0.0
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C  OLD PARTICLE MIXING IS NOW OPTOIN ISPMXZ=2
C
      IF(ISPMXZ(NT).EQ.2)THEN
 
      DO L=2,LA
        DEPINBED=0.
        LZ=LPMXZ(L)
        IF(KBT(L).GT.2)THEN
          DO K=KBT(L),2,-1
            KM=K-1
            DEPINBED=DEPINBED+HBED(L,K)
            DO ND=1,NPMXPTS-1
              NDP=ND+1
              IF(DEPINBED.GE.PMXDEPTH(ND,LZ).AND.
     &           DEPINBED.LT.PMXDEPTH(NDP,LZ))THEN
                WT1=DEPINBED-PMXDEPTH(ND,LZ)
                WT2=PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL=PMXDEPTH(ND+1,LZ)-PMXDEPTH(ND,LZ)
                PARTMIXZ(L,KM)=(WT2*PMXCOEF(ND,LZ)
     &                        +WT1*PMXCOEF(NDP,LZ))/TMPVAL
              ENDIF 
            ENDDO
          ENDDO
        ELSE
          K=KBT(L)
          KM=K-1
          DEPINBED=DEPINBED+HBED(L,K)
          DO ND=1,NPMXPTS-1
            NDP=ND+1
            IF(DEPINBED.GE.PMXDEPTH(ND,LZ).AND.
     &         DEPINBED.LT.PMXDEPTH(NDP,LZ))THEN
              WT1=DEPINBED-PMXDEPTH(ND,LZ)
              WT2=PMXDEPTH(NDP,LZ)-DEPINBED
              TMPVAL=PMXDEPTH(ND+1,LZ)-PMXDEPTH(ND,LZ)
              PARTMIXZ(L,KM)=(WT2*PMXCOEF(ND,LZ)
     &                        +WT1*PMXCOEF(NDP,LZ))/TMPVAL
            ENDIF 
          ENDDO
        ENDIF
      ENDDO
C
      DO L=2,LA
        DO K=KBT(L),2,-1
          KM=K-1
          PARTMIXZ(L,KM)=2.*PARTMIXZ(L,KM)/
     &      (2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
        ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C  NEW PARTICLE MIXING IS NOW OPTOIN ISPMXZ=1
C
      IF(ISPMXZ(NT).EQ.1)THEN
C
      DO L=2,LA
      LZ=LPMXZ(L)
C
        DO K=1,KB
          PARTMIXAVG(K)=0.0
          PARTMIXZ(L,K)=0.0
        ENDDO 
C
        DEPINBED=0.
        DO K=KBT(L),1,-1
          DELHBED=HBED(L,K)/10.
          DO KK=1,10
            DEPINBED=DEPINBED+DELHBED
            DO ND=1,NPMXPTS-1
              NDP=ND+1
              IF(DEPINBED.GE.PMXDEPTH(ND,LZ).AND.
     &          DEPINBED.LT.PMXDEPTH(NDP,LZ))THEN
                WT1=DEPINBED-PMXDEPTH(ND,LZ)
                WT2=PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL=PMXDEPTH(ND+1,LZ)-PMXDEPTH(ND,LZ)
                PARTMIXAVG(K)=PARTMIXAVG(K)
     &             +(WT2*PMXCOEF(ND,LZ)+WT1*PMXCOEF(NDP,LZ))/TMPVAL
              ENDIF 
            ENDDO
          ENDDO
          PARTMIXAVG(K)=PARTMIXAVG(K)/10.
        ENDDO
C
        DO K=1,KBT(L)-1
          TERM1=HBED(L,K)/PARTMIXAVG(K)
          TERM2=HBED(L,K+1)/PARTMIXAVG(K+1)
          TERM3=HBED(L,K)+HBED(L,K+1)
          PARTMIXZ(L,K)=TERM3/(TERM1+TERM2)
        ENDDO 
C
      ENDDO
C
      DO L=2,LA
        DO K=KBT(L),2,-1
          KM=K-1
          PARTMIXZ(L,KM)=2.*PARTMIXZ(L,KM)/
     &      (2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
        ENDDO
      ENDDO
C

      ENDIF
C
C**********************************************************************C
C
      RETURN
      END
C
