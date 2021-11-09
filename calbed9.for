C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALBED9
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
C **  SUBROUTINE CALBED9 CALCULATES CALCULATES BED CONSOLIDATION 
C     WHERE A DIFFERENT TYPE OF CONSOLIDATION CAN BE USED FOR EACH
C     CELL
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
      DIMENSION SNDHYDCN(LCM,KBM)
C
C**********************************************************************C
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
        HBEDMIN=1.E-4
        IF(ISTRAN(7).GE.1)THEN
          HBEDMIN=MAX(HBEDMIN,SNDDMX)
        END IF
C
C----------------------------------------------------------------------C
C
C ** CONSTANT POROSITY BED
C
C        IF(IBMECH.EQ.0)THEN
C
C ** UPDATE TOP LAYER THICKNESS TO MAINTIAN CONSTANT POROSITY
C
          VOIDCON1=BEDPORC/(1.-BEDPORC)
C
          DO L=2,LA
            IF(LCONSOL(L).EQ.0)THEN
              K=KBT(L)
              HBEDTMP=(1.+VOIDCON1)*HBED(L,K)/(1.+VDRBED(L,K))
              TMPVALO=VDRBED(L,K)*HBED(L,K)/(1.+VDRBED(L,K))
              TMPVALN=VOIDCON1*HBEDTMP/(1.+VOIDCON1)
              QWBDTOP(L)=DELTI*(TMPVALO-TMPVALN)
              HBED(L,K)=HBEDTMP
              QWTRBED(L,K)=QWBDTOP(L)+QGW(L)/DXYP(L)
            ENDIF
          ENDDO
C
          DO K=0,KBT(L)-1
            DO L=2,LA
              IF(LCONSOL(L).EQ.0)THEN
                QWTRBED(L,K)=QGW(L)/DXYP(L)
              ENDIF
            ENDDO
          END DO
C
C ** ADD OR REMOVE LAYERS
C
C          CALL CALBLAY
C
C        ENDIF
C
C----------------------------------------------------------------------C
C
C ** SIMPLE CONSOLIDATING BED
C
C        IF(IBMECH.EQ.1)THEN
C
C ** DETERMINE TIME DIFFERENCE AND UPDATE VOID RATIO
C
C BEGIN JMH FIXED IBMECH.EQ.1  OPTION 12/30/02
C
C **  IF SEDVRDT.GT.0.0001 CONSOLIDATE TO SEDVRM (THE MINIMUM VOID RATIO)
C
          IF(SEDVRDT.GT.0.0001)THEN
		  TMPEXP=EXP(-DELT/SEDVRDT)
            DO K=1,KB
              DO L=2,LA
                IF(LCONSOL(L).EQ.1)THEN
                  IF(K.LE.KBT(L))THEN
                    VDRBED1(L,K)=VDRBED(L,K)
                    HBED1(L,K)=HBED(L,K)
                    VDRBED(L,K)=SEDVDRM+(VDRBED1(L,K)-SEDVDRM)*TMPEXP
                    TMPTOP=1.+VDRBED(L,K)
                    TMPBOT=1.+VDRBED1(L,K)
                    HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
	    ENDIF
C
C **  IF SEDVRDT.GT.0.0001 CONSOLIDATE TO SEDVRM INSTANTANEOUSLY
C
          IF(SEDVRDT.GE.0.0.AND.SEDVRDT.LE.0.0001)THEN
		  TMPEXP=0.0
            DO K=1,KB
              DO L=2,LA
                IF(LCONSOL(L).EQ.1)THEN
                  IF(K.LE.KBT(L))THEN
                    VDRBED1(L,K)=VDRBED(L,K)
                    HBED1(L,K)=HBED(L,K)
                    VDRBED(L,K)=SEDVDRM
                    TMPTOP=1.+VDRBED(L,K)
                    TMPBOT=1.+VDRBED1(L,K)
                    HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
	    ENDIF
C
C
C **  IF SEDVRDT.LT.0.0 MAINTAIN INITIAL VOID RATIO (SAVED IN VDRBED2)
C
C          IF(SEDVRDT.LT.0.0) THEN
C		  TMPEXP=1.0
C           DO K=1,KB
C              DO L=2,LA
C                IF(K.LE.KBT(L))THEN
C                  VDRBED1(L,K)=VDRBED(L,K)
C                  HBED1(L,K)=HBED(L,K)
C                  VDRBED(L,K)=VDRBED2(L,K)
C                  TMPTOP=1.+VDRBED(L,K)
C                  TMPBOT=1.+VDRBED1(L,K)
C                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
C                ENDIF
C              ENDDO
C            ENDDO
C	    ENDIF
C
          IF(SEDVRDT.LT.0.0) THEN
		  TMPEXP=1.0
C            DO K=1,KB
              DO L=2,LA
                IF(LCONSOL(L).EQ.1)THEN
                  K=KBT(L)
                  VDRBED1(L,K)=VDRBED(L,K)
                  HBED1(L,K)=HBED(L,K)
                  VDRBED(L,K)=VDRBED2(L,K)
c	            IF(VDRBED(L,K).NE.VDRBED1(L,K))THEN
                  TMPTOP=1.+VDRBED(L,K)
                  TMPBOT=1.+VDRBED1(L,K)
                  HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
c                 ENDIF
                ENDIF
              ENDDO
              DO L=2,LA
                IF(LCONSOL(L).EQ.1)THEN
                  K=KBT(L)-1
                  IF(K.GT.0)THEN
                    VDRBED1(L,K)=VDRBED(L,K)
                    HBED1(L,K)=HBED(L,K)
                    VDRBED(L,K)=VDRBED2(L,K)
c	              IF(VDRBED(L,K).NE.VDRBED1(L,K))THEN
                    TMPTOP=1.+VDRBED(L,K)
                    TMPBOT=1.+VDRBED1(L,K)
                    HBED(L,K)=TMPTOP*HBED1(L,K)/TMPBOT
c                   ENDIF
                  ENDIF
                ENDIF
              ENDDO
C            ENDDO
	    ENDIF
C
C END JMH FIXED IBMECH.EQ.1  OPTION 12/30/02
C
C ** UPDATE POROSITY
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).EQ.1)THEN
                IF(K.LE.KBT(L))THEN
                  PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                  PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** UPDATE PORE WATER FLOWS
C
          DO L=2,LA
            IF(LCONSOL(L).EQ.1)THEN
              QWTRBED(L,0)=QGW(L)/DXYP(L)
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).EQ.1)THEN
                IF(K.LE.KBT(L))THEN
                  TMPVAL=HBED(L,K)/(1.+VDRBED(L,K))
                  QWTRBED(L,K)=QWTRBED(L,K-1)
     &            -DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** ADD OR REMOVE LAYERS
C
C          CALL CALBLAY
C
C        ENDIF
C
C----------------------------------------------------------------------C
C
C ** SET FLAG FOR FINITE STRAIN CONSOLIDATION
C
        IFLAG=0
	  IF(ISTRAN(6).GE.1)IFLAG=IFLAG+1
	  IF(ISTRAN(7).GE.1)IFLAG=IFLAG+1
C
C----------------------------------------------------------------------C
C
C ** FINITE STRAIN CONSOLIDATING HOMOGENEOUS BED
C
C        IF(IBMECH.GE.2.AND.IFLAG.EQ.1)THEN
        IF(IFLAG.EQ.1)THEN
C
          WDENKGM3=1.E3
          WDENGMM3=1.E6
C
C ++  SET PHYSICAL VERTICAL COORDINATES OF THE BED
C
C     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
C     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDG(L,0)=ZELBEDA(L)
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
                ELSE  
                  ZBEDG(L,K)=HBED(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDGT(L)=ZBEDG(L,KBT(L))
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
                ELSE  
                  ZBEDC(L,K)=0.5*ZBEDG(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS
C     DZBTR = TRANSFORMED THICKNESS OF BED LAYER
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                  DZBTR1(L,K)=HBED1(L,K)/(1.+VDRBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED
C       EFFECTIVE STRESS USING FSTRSE
C     CALCULATE DERIVATIVE OF EFFECTIVE STRESS 
C       WITH RESPECT TO VOID RATIO, DSTRSE USING
C       FUNCTION FDSTRSE
C     CALCULATE HYDRAULIC CONDUCTIVITY DIVIED BY (1+VOID),
C       HYDCN =USING FUNCTION FHYDCN 
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  STRSE(L,K)=
     &            FSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)
                  DSTRSE(L,K)=
     &            FDSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)
                  HYDCN(L,K)=
     &            FHYDCN(VDRBED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                KBTM1=KBT(L)-1
                IF(K.LE.KBTM1)THEN
                  TMPVAL=( DZBTR(L,K)/HYDCN(L,K) )
     &            +( DZBTR(L,K+1)/HYDCN(L,K+1) )
                  COEFK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
c fix 0602                STTMP=DSTRSE(L,K)*HYDCN(L,K)
c fix 0602                STTMPP=DSTRSE(L,K+1)*HYDCN(L,K+1)
c fix 0602                TMPVAL=( DZBTR(L,K)/STTMP )
c fix 0602     &          +( DZBTR(L,K+1)/STTMPP )
c fix 0602                COEFSK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
                  DSTRESET=(DZBTR(L,K)*DSTRSE(L,K+1)
     &                   +DZBTR(L,K+1)*DSTRSE(L,K))
     &                   /(DZBTR(L,K)+DZBTR(L,K+1))
                  COEFSK(L,K)=DSTRESET*COEFK(L,K)
                ENDIF
                IF(K.EQ.KBT(L))THEN  
                  COEFK(L,K)=HYDCN(L,K)
                  COEFSK(L,K)=DSTRSE(L,K)*HYDCN(L,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE PRESSURE COMPONENTS
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                SGSM1(L,K)=0.
              ENDIF
            ENDDO
          ENDDO
C
          IF(ISTRAN(6).GT.0)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  IF(LCONSOL(L).GE.2)THEN
                    IF(K.LE.KBT(L))THEN
                      SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISTRAN(7).GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              DO K=1,KB
                DO L=2,LA
                  IF(LCONSOL(L).GE.2)THEN
                    IF(K.LE.KBT(L))THEN
                      SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  SGSM1(L,K)=SGSM1(L,K)-1.
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LT.KBT(L))THEN
                  ACOEF(L,K)=2.*COEFSK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
                  TMPVALK=DZBTR(L,K)*SGSM1(L,K)
                  TMPVALKP=DZBTR(L,K+1)*SGSM1(L,K+1)
                  QCOEF(L,K)=(TMPVALK+TMPVALKP)*COEFK(L,K)/
     &                     (DZBTR(L,K)+DZBTR(L,K+1))
                ELSE
                  ACOEF(L,K)=0.0
                  QCOEF(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              K=KBT(L)
              ACOEF(L,K)=2.*COEFSK(L,K)/DZBTR(L,K)
              QCOEF(L,K)=SGSM1(L,K)*COEFK(L,K)
              QCOEF(L,0)=QGW(L)*DXYIP(L)
              QWTRBED(L,0)=QGW(L)*DXYIP(L)
            ENDIF
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ALOW(L,1)=0.
              CUPP(L,KBT(L))=0.
              DO K=1,KBT(L)-1
                CUPP(L,K)=-DELT*ACOEF(L,K)/DZBTR(L,K)
              ENDDO
              DO K=2,KBT(L)
                ALOW(L,K)=-DELT*ACOEF(L,K-1)/DZBTR(L,K)
              ENDDO
              DO K=1,KBT(L)
                BMNN(L,K)=1.0-ALOW(L,K)-CUPP(L,K)
              ENDDO
              K=KBT(L)
              BMNN(L,K)=BMNN(L,K)+DELT*ACOEF(L,K)/DZBTR(L,K)
              DO K=1,KBT(L)
                RRHS(L,K)=VDRBED(L,K)
     &          +DELT*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)
                VDRBED1(L,K)=VDRBED(L,K)
                HBED1(L,K)=HBED(L,K)
              ENDDO
            ENDIF
          ENDDO
C
C         IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.2)THEN
                K=KBT(L)
C                RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*SEDVDRD/DZBTR(L,K)
                RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*VDRDEPO(1)
     &                   /DZBTR(L,K)
              ENDIF
            ENDDO
C         ENDIF
C
C         IF(IBMECH.EQ.3) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.3)THEN
                K=KBT(L)
                RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*(VDRBED(L,K)
     &                 +(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)
              ENDIF
            ENDDO
C         ENDIF
C
C          DO L=2,LA
C           DO K=1,KBT(L)
C            WRITE(8,*)K,ALOW(L,K),BMNN(L,K),CUPP(L,K),RRHS(L,K)
C           END DO
C          END DO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              BETTMP=BMNN(L,1)
              TOXTMP(L,1)=RRHS(L,1)/BETTMP
              DO KK=2,KBT(L)
                GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
                BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
                TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/
     &                     BETTMP
              ENDDO
              DO KK=KBT(L)-1,1,-1
                 TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
              ENDDO
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LT.KBT(L))THEN
                  QWTRBED(L,K)=-ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))
     &                      +QCOEF(L,K)
                ELSE
                  QWTRBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C         IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.2)THEN
                K=KBT(L)
                QWTRBED(L,K)=-ACOEF(L,K)*(SEDVDRD-TOXTMP(L,K))
     &                       +QCOEF(L,K)
              ENDIF
            ENDDO
C         ELSE
            DO L=2,LA
              IF(LCONSOL(L).EQ.3)THEN
                K=KBT(L)
                QWTRBED(L,K)=-ACOEF(L,K)*(VDRBED(L,K)
     &                +(STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K))+QCOEF(L,K)
              ENDIF
            ENDDO
C         END IF
C
C ++  CALCULATE VOID RATIOS  
C     VDRBED =  VOID RATIO OF BED LAYER 
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  VDRBED(L,K)=VDRBED1(L,K)
     &             -DELT*(QWTRBED(L,K)-QWTRBED(L,K-1))/DZBTR1(L,K)
                ELSE
                  VDRBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  UPDATE LAYER THICKNESS
C     HBED = BED LAYER THICKNESS
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  HBED(L,K)=HBED1(L,K)*(1.+VDRBED(L,K))
     &            /(1.+VDRBED1(L,K))
                ELSE
                  HBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              QWTRBED(L,0)=QGW(L)/DXYP(L)
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  QWTRBED(L,K)=QWTRBED(L,K-1)
     &             -DELTI*(HBED(L,K)-HBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
C     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
                ELSE  
                  ZBEDG(L,K)=HBED(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDGT(L)=ZBEDG(L,KBT(L))
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
                ELSE  
                  ZBEDC(L,K)=0.5*ZBEDG(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** UPDATE POROSITY
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                  PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
                ELSE
                  PORBED(L,K)=0.0
                  PORBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** ADD OR REMOVE LAYERS
C
C          CALL CALBLAY
C
        ENDIF
C
C----------------------------------------------------------------------C
C
C ** FINITE STRAIN CONSOLIDATING NON-HOMOGENEOUS BED
C
C       IF(IBMECH.GE.2.AND.IFLAG.EQ.2)THEN
        IF(IFLAG.EQ.2)THEN
C
          WDENKGM3=1.E3
          WDENGMM3=1.E6
C
C ++  SET PHYSICAL VERTICAL COORDINATES OF THE BED
C
C     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
C     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDG(L,0)=ZELBEDA(L)
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
                ELSE  
                  ZBEDG(L,K)=HBED(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDGT(L)=ZBEDG(L,KBT(L))
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
                ELSE  
                  ZBEDC(L,K)=0.5*ZBEDG(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS
C     DZBTR = TRANSFORMED THICKNESS OF BED LAYER
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                  DZBTR1(L,K)=HBED1(L,K)/(1.+VDRBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED
C       EFFECTIVE STRESS USING FSTRSE
C     CALCULATE DERIVATIVE OF EFFECTIVE STRESS 
C       WITH RESPECT TO VOID RATIO, DSTRSE USING
C       FUNCTION FDSTRSE
C     CALCULATE HYDRAULIC CONDUCTIVITY DIVIED BY (1+VOID),
C       HYDCN =USING FUNCTION FHYDCN 
C
C ++  NONCOHESIVE HYDRAULIC CONDUCTIVITY BASED ON R. R. RUMER, CHAP 3,
C ++  EQ 13 AND 14, IN 'FLOW THROUGH POROUS MEDIA' ED. R. J. M. DE WIEST, 
C ++  ACADEMIC PRESS, 1969
C
C     KH=2854*(n**2)*(d**2)/(1-n) = 2854*(e**2)*(d**2)/(1+e)
C     for KH in m/s and d in meters
C
 
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  STRSE(L,K)=
     &              FSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)
                  DSTRSE(L,K)=
     &              FDSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)
                  HYDCN(L,K)=
     &            FHYDCN(VDRBEDSED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)
	            TMPVAL=VDRBEDSND(L,K)/(1.+VDRBEDSND(L,K))
	            SNDHYDCN(L,K)=2854.0*TMPVAL*TMPVAL*(SEDDIA50(L,K)**2)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C 
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  COEFK(L,K)=(FRACCOH(L,K)/HYDCN(L,K))
     &                +(FRACNON(L,K)/SNDHYDCN(L,K))
                 ENDIF
               ENDIF
            ENDDO
          ENDDO
C 
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  HYDCN(L,K)=1./COEFK(L,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C 
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                COEFK(L,K)=0.
              ENDIF
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                KBTM1=KBT(L)-1
                IF(K.LE.KBTM1)THEN
                  TMPVAL=( DZBTR(L,K)/HYDCN(L,K) )
     &            +( DZBTR(L,K+1)/HYDCN(L,K+1) )
                  COEFK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
c fix 0602                STTMP=DSTRSE(L,K)*HYDCN(L,K)
c fix 0602                STTMPP=DSTRSE(L,K+1)*HYDCN(L,K+1)
c fix 0602                TMPVAL=( DZBTR(L,K)/STTMP )
c fix 0602     &          +( DZBTR(L,K+1)/STTMPP )
c fix 0602                COEFSK(L,K)=(DZBTR(L,K)+DZBTR(L,K+1))/TMPVAL
                  DSTRESET=(DZBTR(L,K)*DSTRSE(L,K+1)
     &                   +DZBTR(L,K+1)*DSTRSE(L,K))
     &                   /(DZBTR(L,K)+DZBTR(L,K+1))
                  COEFSK(L,K)=DSTRESET*COEFK(L,K)
                ENDIF
                IF(K.EQ.KBT(L))THEN  
                  COEFK(L,K)=HYDCN(L,K)
                  COEFSK(L,K)=DSTRSE(L,K)*HYDCN(L,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  CALCULATE PRESSURE COMPONENTS
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                SGSM1(L,K)=0.
              ENDIF
            ENDDO
          ENDDO
C
          IF(ISTRAN(6).GT.0)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  IF(LCONSOL(L).GE.2)THEN
                    IF(K.LE.KBT(L))THEN
                      SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISTRAN(7).GT.0)THEN
            DO NX=1,NSND
              NS=NSED+NX
              DO K=1,KB
                DO L=2,LA
                  IF(LCONSOL(L).GE.2)THEN
                    IF(K.LE.KBT(L))THEN
                      SGSM1(L,K)=SGSM1(L,K)+SSG(NS)*VFRBED(L,K,NS)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  SGSM1(L,K)=SGSM1(L,K)-1.
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LT.KBT(L))THEN
                  ACOEF(L,K)=2.*COEFSK(L,K)/(DZBTR(L,K)+DZBTR(L,K+1))
                  TMPVALK=DZBTR(L,K)*SGSM1(L,K)
                  TMPVALKP=DZBTR(L,K+1)*SGSM1(L,K+1)
                  QCOEF(L,K)=(TMPVALK+TMPVALKP)*COEFK(L,K)/
     &                     (DZBTR(L,K)+DZBTR(L,K+1))
                ELSE
                  ACOEF(L,K)=0.0
                  QCOEF(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              K=KBT(L)
              ACOEF(L,K)=2.*COEFSK(L,K)/DZBTR(L,K)
              QCOEF(L,K)=SGSM1(L,K)*COEFK(L,K)
              QCOEF(L,0)=QGW(L)*DXYIP(L)
              QWTRBED(L,0)=QGW(L)*DXYIP(L)
            ENDIF
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ALOW(L,1)=0.
              CUPP(L,KBT(L))=0.
              DO K=1,KBT(L)-1
                CUPP(L,K)=-DELT*ACOEF(L,K)/DZBTR(L,K)
              ENDDO
              DO K=2,KBT(L)
                ALOW(L,K)=-DELT*ACOEF(L,K-1)/DZBTR(L,K)
              ENDDO
              DO K=1,KBT(L)
                BMNN(L,K)=FRACCOH(L,K)-ALOW(L,K)-CUPP(L,K)
              ENDDO
              K=KBT(L)
              BMNN(L,K)=BMNN(L,K)+DELT*ACOEF(L,K)/DZBTR(L,K)
              DO K=1,KBT(L)
                RRHS(L,K)=FRACCOH(L,K)*VDRBEDSED(L,K)
     &          +DELT*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)
                VDRBED1(L,K)=VDRBED(L,K)
                HBED1(L,K)=HBED(L,K)
              ENDDO
            ENDIF
          ENDDO
C
C         IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.2)THEN
                K=KBT(L)
                RRHS(L,K)=RRHS(L,K)+DELT*ACOEF(L,K)*VDRDEPO(1)
     &          /DZBTR(L,K)
              ENDIF
            ENDDO
C         ENDIF
C
C         IF(IBMECH.EQ.3) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.3)THEN
                K=KBT(L)
                RRHS(L,K)=RRHS(L,K)
     &          +DELT*ACOEF(L,K)*(VDRBEDSED(L,K)
     &                   +(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)
              ENDIF
            ENDDO
C         ENDIF
C
C          DO L=2,LA
C           DO K=1,KBT(L)
C            WRITE(8,*)K,ALOW(L,K),BMNN(L,K),CUPP(L,K),RRHS(L,K)
C           END DO
C          END DO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              BETTMP=BMNN(L,1)
              TOXTMP(L,1)=RRHS(L,1)/BETTMP
              DO KK=2,KBT(L)
                GAMTMP(L,KK)=CUPP(L,KK-1)/BETTMP
                BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
                TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/
     &                     BETTMP
              ENDDO
              DO KK=KBT(L)-1,1,-1
                TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
              ENDDO
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LT.KBT(L))THEN
                  QWTRBED(L,K)=-ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))
     &                      +QCOEF(L,K)
                ELSE
                  QWTRBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C         IF(IBMECH.EQ.2) THEN
            DO L=2,LA
              IF(LCONSOL(L).EQ.2)THEN
                K=KBT(L)
                QWTRBED(L,K)=-ACOEF(L,K)*(VDRDEPO(1)-TOXTMP(L,K))
     &                     +QCOEF(L,K)
              ENDIF
            ENDDO
C         ELSE
            DO L=2,LA
              IF(LCONSOL(L).EQ.3)THEN
                K=KBT(L)
                QWTRBED(L,K)=-ACOEF(L,K)*(VDRBEDSED(L,K)
     &                +(STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K))+QCOEF(L,K)
              ENDIF
            ENDDO
C         END IF
C
C ++  CALCULATE VOID RATIOS  
C     VDRBED =  VOID RATIO OF BED LAYER 
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
	            IF(FRACCOH(L,K).GT.0.0)THEN
                    VDRBEDSED(L,K)=VDRBEDSED(L,K)
     &              -DELT*(QWTRBED(L,K)-QWTRBED(L,K-1))
     &              /(FRACCOH(L,K)*DZBTR1(L,K))
                  VDRBEDSED(L,K)=MAX(VDRBEDSED(L,K),VDRBEDSND(L,K))
                  ELSE
                    VDRBEDSED(L,K)=0.0
                  ENDIF
                ELSE
                  VDRBEDSED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ++  UPDATE VOID RATIO
C 
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  VDRBED(L,K)=FRACCOH(L,K)*VDRBEDSED(L,K)
     &                +FRACNON(L,K)*VDRBEDSND(L,K)
                ELSE
                  VDRBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C 
C ++  UPDATE LAYER THICKNESS
C     HBED = BED LAYER THICKNESS
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  HBED(L,K)=HBED1(L,K)*(1.+VDRBED(L,K))
     &            /(1.+VDRBED1(L,K))
                ELSE
                  HBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              QWTRBED(L,0)=QGW(L)/DXYP(L)
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  QWTRBED(L,K)=QWTRBED(L,K-1)
     &               -DELTI*(HBED(L,K)-HBED1(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C     ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
C     ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
                ELSE  
                  ZBEDG(L,K)=HBED(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO L=2,LA
            IF(LCONSOL(L).GE.2)THEN
              ZBEDGT(L)=ZBEDG(L,KBT(L))
            ENDIF
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
                ELSE  
                  ZBEDC(L,K)=0.5*ZBEDG(L,K)  
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** UPDATE POROSITY
C
          DO K=1,KB
            DO L=2,LA
              IF(LCONSOL(L).GE.2)THEN
                IF(K.LE.KBT(L))THEN
                  PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
                  PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
                ELSE
                  PORBED(L,K)=0.0
                  PORBED(L,K)=0.0
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
C ** ADD OR REMOVE LAYERS
C
C          CALL CALBLAY
C
        ENDIF
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C
C**********************************************************************C
C
      RETURN
      END
