C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE ADTOXDP(NT,DELTI,PARTDIF)
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
C **  SUBROUTINE ADTOXDP SOLVES TOXIC PORE WATER ADVECTION AND 
C **  DIFFUSION IN DOUBLE PRECISION
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C      REAL*4
C       DELTI                   PASSED AS ARGUMENT    REQUIRES LOCAL DP VERSION
C      BETTMP                  LOCAL                         
C       DIFBWFAC                LOCAL
C
C      REAL*4
C      DZ(KCM)                 IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
C
C      REAL*4
C      DIFTOX(NTXM)            IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
C      DIFTOXS(NTXM)           IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
C
C      REAL*4 
C      HP(LCM)                 IN GLOBAL COMMON      REQUIRES LOCAL DP VERSION
C      TOXBBALO(LCM)           LOCAL
C      DIFTOXBW(LCM)           LOCAL
C      TOXBBALN(LCM)           LOCAL
C      TOXWBALO(LCM)           LOCAL
C      TOXWBALN(LCM)           LOCAL
C
C       REAL*4
C      HBED(LCM,KBM)           IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C      QWTRBED(LCM,0:KBM)      IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C       PARTDIF(LCM,KBM)        PASSED AS ARGUMENT   REQUIRES LOCAL DP VERSION
C       PORBED(LCM,KBM)         IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C
C      REAL*4
C      ALOW(LCM,KBM+1)         LOCAL  
C      BMNN(LCM,KBM+1)         LOCAL  
C      CUPP(LCM,KBM+1)         LOCAL  
C      RRHS(LCM,KBM+1)         LOCAL  
C      GAMTMP(LCM,KBM+1)       LOCAL  
C      TOXTMP(LCM,KBM+1)       LOCAL
C
C      REAL*4
C      TADFLUX(LCM,NTXM)       IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP  
C      CONGW(LCM,NSTVM)        IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C                                                FOR CONGW(L,NT+4)
C
C      REAL*4
C      TOX(LCM,KCM,NTXM)       IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP
C      TOXPFTW(LCM,KCM,NTXM)   IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C
C      REAL*4
C      TOXB(LCM,KBM,NT)        IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION  RETURN SP
C      TOXPFTB(LCM,KKBM,NTXM)  IN GLOBAL COMMON     REQUIRES LOCAL DP VERSION
C
C**********************************************************************C
C
C  **  DECLARE PASSED ARGUMENTS
C
      REAL*4 DELTI
C
      REAL*4 PARTDIF
      DIMENSION PARTDIF(LCM,KBM)
C
C  **  DECLARE LOCAL DOUBLE PRECISION VARIABLES
C
      REAL*8 DELTI_DP,BETTMP,DIFBWFAC
C
      REAL*8 DZC_DP
      DIMENSION DZC_DP(KCM)
C
      REAL*8 DIFTOX_DP,DIFTOXS_DP
      DIMENSION DIFTOX_DP(NTXM),DIFTOXS_DP(NTXM)
C
      REAL*8 DIFTOXBW,TOXBBALO,TOXBBALN,TOXWBALO,TOXWBALN,HP_DP
      DIMENSION DIFTOXBW(LCM),TOXBBALO(LCM),TOXBBALN(LCM),
     &                        TOXWBALO(LCM),TOXWBALN(LCM),HP_DP(LCM)
C
      REAL*8 PARTDIF_DP,HBED_DP,PORBED_DP
      DIMENSION PARTDIF_DP(LCM,KBM),HBED_DP(LCM,KBM),PORBED_DP(LCM,KBM)
C
      REAL*8 QWTRBED_DP
      DIMENSION QWTRBED_DP(LCM,0:KBM)
C
      REAL*8 ALOW,BMNN,CUPP,RRHS,GAMTMP,TOXTMP
      DIMENSION ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),
     &       CUPP(LCM,KBM+1),RRHS(LCM,KBM+1),GAMTMP(LCM,KBM+1),
     &       TOXTMP(LCM,KBM+1)
C
      REAL*8 TADFLUX_DP
      DIMENSION TADFLUX_DP(LCM,NTXM)
C        
      REAL*8 CONGW_DP    
      DIMENSION CONGW_DP(LCM,NSTVM)     
C
      REAL*8 TOX_DP,TOXPFTW_DP
      DIMENSION TOX_DP(LCM,KCM,NTXM),TOXPFTW_DP(LCM,KCM,NTXM)
C
      REAL*8 TOXB_DP,TOXPFTB_DP
      DIMENSION TOXB_DP(LCM,KBM,NTXM),TOXPFTB_DP(LCM,KBM,NTXM)
C
C**********************************************************************C
C
C  **  LOAD DOUBLE PRECISION VARIABLS
C
       DELTI_DP=DBLE(DELTI)
       DZC_DP(1)=DBLE(DZC(1))
       DIFTOX_DP(NT)=DBLE(DIFTOX(NT))
       DIFTOXS_DP(NT)=DBLE(DIFTOXS(NT))
C
       DO L=2,LA
         HP_DP(L)=DBLE(HP(L))
         QWTRBED_DP(L,0)=DBLE(QWTRBED(L,0))
         TOX_DP(L,1,NT)=DBLE(TOX(L,1,NT))
         TOXPFTW_DP(L,1,NT)=DBLE(TOXPFTW(L,1,NT))
C         TADFLUX_DP(L,NT)=DBLE(TADFLUX(L,NT))
         CONGW_DP(L,NT+4)=DBLE(CONGW(L,NT+4))
       ENDDO
C
       DO K=1,KB
       DO L=2,LA
         HBED_DP(L,K)=DBLE(HBED(L,K))
         QWTRBED_DP(L,K)=DBLE(QWTRBED(L,K))
         PORBED_DP(L,K)=DBLE(PORBED(L,K))
         PARTDIF_DP(L,K)=DBLE(PARTDIF(L,K))
         TOXB_DP(L,K,NT)=DBLE(TOXB(L,K,NT))
         TOXPFTB_DP(L,K,NT)=DBLE(TOXPFTB(L,K,NT))
       ENDDO
       ENDDO
C
C**********************************************************************C
C
          DO L=2,LA
            DIFTOXBW(L)=0.0
          ENDDO
C
          DO L=2,LA
            IF(LMASKDRY(L)) DIFTOXBW(L)=DIFTOXS_DP(NT)
          ENDDO
C
            DO L=2,LA
C 
              DIFBWFAC=2./HBED_DP(L,KBT(L))
              IF(ISDIFBW(NT).EQ.1)DIFBWFAC=1.0
              TOXBBALO(L)=0.
              KBTP1=KBT(L)+1
              KBTM1=KBT(L)-1
              ALOW(L,1)=0.
              CUPP(L,KBTP1)=0.
C
              DO K=1,KBTM1
C###############################################################################
C HQI Change to implement particle mixing as a rate in units of length/time 
C instead of a mixing coefficient in units of length**2/time
C RM, 09/28/05
c                CUPP(L,K)=MIN(QWTRBED_DP(L,K),0.)
c     &          -(DIFTOX_DP(NT)+PARTDIF_DP(L,K))*(PORBED_DP(L,K)
c     &                  +PORBED_DP(L,K+1))/
c     &                  (HBED_DP(L,K)+HBED_DP(L,K+1))
                CUPP(L,K)=MIN(QWTRBED_DP(L,K),0.)
     &          -((DIFTOX_DP(NT)/(HBED_DP(L,K)+HBED_DP(L,K+1)))+
     &          PARTDIF_DP(L,K))*(PORBED_DP(L,K)+PORBED_DP(L,K+1))
C###############################################################################
              ENDDO
              CUPP(L,KBT(L))=MIN(QWTRBED_DP(L,KBT(L)),0.)
     &          -DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
C
              DO K=2,KBT(L)
C###############################################################################
C HQI Change to implement particle mixing as a rate in units of length/time 
C instead of a mixing coefficient in units of length**2/time
C RM, 09/28/05
c                ALOW(L,K)=-MAX(QWTRBED_DP(L,K-1),0.)
c     &          -(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
c     &                  +PORBED_DP(L,K))
c     &                  /(HBED_DP(L,K-1)+HBED_DP(L,K))
                ALOW(L,K)=-MAX(QWTRBED_DP(L,K-1),0.)
     &          -((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K)))+
     &          PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
C###############################################################################
              ENDDO
              ALOW(L,KBTP1)=-MAX(QWTRBED_DP(L,KBT(L)),0.)
     &                 -DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
C
              DO K=1,KBT(L)
                BMNN(L,K)=DELTI_DP*HBED_DP(L,K)*PORBED_DP(L,K)/
     &                  (1.-TOXPFTB_DP(L,K,NT))
              ENDDO
              BMNN(L,KBTP1)=DELTI_DP*DZC_DP(1)*HP_DP(L)/
     &                  (1.-TOXPFTW_DP(L,1,NT))
C
C###############################################################################
C HQI Change to implement particle mixing as a rate in units of length/time 
C instead of a mixing coefficient in units of length**2/time
C RM, 09/28/05
c              BMNN(L,1)=BMNN(L,1)+MAX(QWTRBED_DP(L,1),0.)
c     &        +(DIFTOX_DP(NT)+PARTDIF_DP(L,1))*(PORBED_DP(L,2)
c     &                  +PORBED_DP(L,1))/
c     &                  (HBED_DP(L,2)+HBED_DP(L,1))
              BMNN(L,1)=BMNN(L,1)+MAX(QWTRBED_DP(L,1),0.)
     &        +((DIFTOX_DP(NT)/(HBED_DP(L,2)+HBED_DP(L,1)))+
     &        PARTDIF_DP(L,1))*(PORBED_DP(L,2)+PORBED_DP(L,1))
C###############################################################################
              DO K=2,KBTM1
C###############################################################################
C HQI Change to implement particle mixing as a rate in units of length/time 
C instead of a mixing coefficient in units of length**2/time
C RM, 09/28/05
c                BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
c     &         +(DIFTOX_DP(NT)+PARTDIF_DP(L,K))*(PORBED_DP(L,K+1)
c     &                  +PORBED_DP(L,K))/
c     &                  (HBED_DP(L,K+1)+HBED_DP(L,K))
c     &                         -MIN(QWTRBED_DP(L,K-1),0.)
c     &         +(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
c     &                  +PORBED_DP(L,K))/
c     &                  (HBED_DP(L,K-1)+HBED_DP(L,K))
                BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
     &         +((DIFTOX_DP(NT)/(HBED_DP(L,K+1)+HBED_DP(L,K)))+
     &         PARTDIF_DP(L,K))*(PORBED_DP(L,K+1)+PORBED_DP(L,K))
     &                         -MIN(QWTRBED_DP(L,K-1),0.)
     &         +((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K)))+
     &         PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
C###############################################################################
              ENDDO
              K=KBT(L)
C###############################################################################
C HQI Change to implement particle mixing as a rate in units of length/time 
C instead of a mixing coefficient in units of length**2/time
C RM, 09/28/05
c              BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
c     &           +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
c     &                         -MIN(QWTRBED_DP(L,K-1),0.)
c     &           +(DIFTOX_DP(NT)+PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)
c     &                  +PORBED_DP(L,K))/(HBED_DP(L,K-1)+HBED_DP(L,K))
              BMNN(L,K)=BMNN(L,K)+MAX(QWTRBED_DP(L,K),0.)
     &           +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
     &                         -MIN(QWTRBED_DP(L,K-1),0.)
     &           +((DIFTOX_DP(NT)/(HBED_DP(L,K-1)+HBED_DP(L,K)))
     &           +PARTDIF_DP(L,K-1))*(PORBED_DP(L,K-1)+PORBED_DP(L,K))
C###############################################################################
              K=KBTP1
              BMNN(L,K)=BMNN(L,K)-MIN(QWTRBED_DP(L,K-1),0.)
     &          +DIFBWFAC*DIFTOXBW(L)*PORBED_DP(L,KBT(L))
C
              DO K=1,KBT(L)
                RRHS(L,K)=DELTI_DP*TOXB_DP(L,K,NT)
                TOXBBALO(L)=TOXBBALO(L)+TOXB_DP(L,K,NT)
              ENDDO
              RRHS(L,1)=RRHS(L,1)+MAX(QWTRBED_DP(L,0),0.)
     &                  *CONGW_DP(L,NT+4)
              RRHS(L,KBTP1)=DELTI_DP*DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
              TOXWBALO(L)=DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
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
C                TOXBSMB(L,K)=TOXB_DP(L,K,NT)
                TOXB_DP(L,K,NT)=HBED_DP(L,K)*PORBED_DP(L,K)*TOXTMP(L,K)/
     &                     (1.-TOXPFTB_DP(L,K,NT))
                TOXBBALN(L)=TOXBBALN(L)+TOXB_DP(L,K,NT)
              ENDDO
C              TOXSMB(L)=TOX_DP(L,1,NT)
              TOX_DP(L,1,NT)=TOXTMP(L,KBTP1)/(1.-TOXPFTW_DP(L,1,NT))
              TOXWBALN(L)=DZC_DP(1)*HP_DP(L)*TOX_DP(L,1,NT)
            ENDDO
C
C **  CALCULATED ADVECTION-DIFFUSION FLUX FOR MASS BALANCE
C
            DO L=2,LA
              TADFLUX_DP(L,NT)=DELTI_DP*(DZC_DP(1)*HP_DP(L)
     &                       *TOX_DP(L,1,NT)-TOXWBALO(L))
            ENDDO
C
C**********************************************************************C
C
C **  RETURN UPDATED VARIABLES IN SINGLE PRECISION BACK TO GLOBAL 
C **  COMMON
C
C
       DO L=2,LA
         TOX(L,1,NT)=SNGL(TOX_DP(L,1,NT))
         TADFLUX(L,NT)=SNGL(TADFLUX_DP(L,NT))
       ENDDO
C
       DO K=1,KB
       DO L=2,LA
         TOXB(L,K,NT)=SNGL(TOXB_DP(L,K,NT))
       ENDDO
       ENDDO
C
C**********************************************************************C
C
      RETURN
      END
