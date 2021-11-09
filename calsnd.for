C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALSND
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
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
C**********************************************************************C
C
C **  SET MAXIMUM NONCOHESIVE SEDIMENT DIAMETER
C
      SNDDMX=0.
      IF(ISTRAN(7).GE.1)THEN
        DO NX=1,NSND
          NS=NSED+NX
          SNDDMX=MAX(SNDDMX,SEDDIA(NS))
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  SET NONCOHESIVE ARMORING PARAMETERS
C
      IF(ISTRAN(7).GE.1)THEN
C
        IF(ISNDAL.EQ.0)THEN
C
          DO K=1,KB
            DO L=2,LA
              SIGPHI(L,K)=0.
            ENDDO
          ENDDO
C
          DO NX=1,NSND
            DO L=2,LA
              PEXP(L,NX)=1.
              PHID(L,NX)=1.
            ENDDO
          ENDDO
C
        ENDIF
C
        IF(ISNDAL.GE.1)THEN
C
        ISGPFLAG=0
	  ISEHFLAG=0
        DO NX=1,NSND
          NS=NX+NSED
          IF(ISNDEQ(NS).EQ.1)ISGPFLAG=1
	    IF(ISBDLD(NS).GE.2)ISEHFLAG=1
	  ENDDO
C
C **  SET SIGPHI FOR GARCIA AND PARKER (1991) EQS 46 AND 47
C
        IF(ISGPFLAG.EQ.1)THEN
C
C **  GP - CONVERT SEDIMENT DIAMETERS IN M TO MM AND SET PHI SIZE
C
          DO NX=1,NSND
            NS=NSED+NX
            SEDPHI(NS)=-LOG(1000.*SEDDIA(NS))/LOG(2.)
          ENDDO
C
C **  GP - SET MEAN PHI FOR TOP LAYER OF BED
C
          DO L=2,LA
            TVAR3W(L)=0.
            TVAR3E(L)=0.
            SIGPHI(L,KBT(L))=0.
          ENDDO
C
          DO NX=1,NSND
           NS=NSED+NX
           DO L=2,LA
	      KTOP=KBT(L)
            TVAR3W(L)=TVAR3W(L)+SEDPHI(NS)*VFRBED(L,KTOP,NS)
            TVAR3E(L)=TVAR3E(L)+VFRBED(L,KTOP,NS)
           ENDDO
          ENDDO
C
          DO L=2,LA
           IF(TVAR3E(L).LE.0.) TVAR3E(L)=1. 
           TVAR3W(L)=TVAR3W(L)/TVAR3E(L)
          ENDDO
C
          DO NX=1,NSND
            NS=NSED+NX
              DO L=2,LA
	         KTOP=KBT(L)
               SIGPHI(L,KTOP)=SIGPHI(L,KTOP)+((SEDPHI(NS)-TVAR3W(L))**2)
     &         *VFRBED(L,KTOP,NS)/TVAR3E(L)
              ENDDO
          ENDDO
C
          DO L=2,LA
	      KTOP=KBT(L)
            SIGPHI(L,KTOP)=SQRT(SIGPHI(L,KTOP))
          ENDDO
C
        ENDIF
C        
C **  END CALCULATION OF SIGPHI FOR GARCIA AND PARKER (1991)
C
C **  SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN AND WU,WANG,&JI
C
	  IF(ISEHFLAG.EQ.1.OR.ISGPFLAG.EQ.1)THEN
C
          DO NX=1,NSND
            DO L=2,LA
              PEXP(L,NX)=0.0
              PHID(L,NX)=0.0
            ENDDO
          ENDDO
C
          DO L=2,LA
            SNDBT(L,KBT(L))=0.0
          ENDDO
C
          DO NX=1,NSND
            NS=NSED+NX
            DO L=2,LA
              K=KBT(L)
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
              DO NXX=1,NSND
                NSS=NSED+NXX
                PEXP(L,NX)=PEXP(L,NX)+SNDB(L,K,NXX)*SEDDIA(NS)
     &                    /(SEDDIA(NS)+SEDDIA(NSS))
                PHID(L,NX)=PHID(L,NX)+SNDB(L,K,NXX)*SEDDIA(NSS)
     &                    /(SEDDIA(NS)+SEDDIA(NSS))
              ENDDO
            ENDDO
          ENDDO
C
          DO NX=1,NSND
            DO L=2,LA
              K=KBT(L)
              IF(SNDBT(L,K).GT.0.0) THEN
                PEXP(L,NX)=PEXP(L,NX)/SNDBT(L,K)
                PHID(L,NX)=PHID(L,NX)/SNDBT(L,K)
              ELSE
                PEXP(L,NX)=1.0
                PHID(L,NX)=1.0
              END IF
            ENDDO
          ENDDO

C          DO L=2,LA
C	      WRITE(8,888)IL(L),JL(L),(PEXP(L,NX),NX=1,NSND),
C     &                  (PHID(L,NX),NX=1,NSND)
C	    ENDDO
C
        ENDIF
C
C **  END SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN
C

        ENDIF
C
      ENDIF
C
  888 FORMAT(2I5,6E12.4)
C
C**********************************************************************C
C
C **  D50 NOW SET IN SSEDTOX
C **  SET MEAN D50
C
C      IF(ISTRAN(7).GE.1)THEN
C
C        DO K=1,KB
C          DO L=2,LA
C            SEDDIA50(L,K)=0.
C            SNDBT(L,K)=0.
C          ENDDO
C        ENDDO
C
C        DO NX=1,NSND
C          NS=NSED+NX
C          DO K=1,KB
C            DO L=2,LA
C              SEDDIA50(L,K)=SEDDIA50(L,K)+SEDDIA(NS)*SNDB(L,K,NX)
C              SEDDIA50(L,K)=SEDDIA50(L,K)+SNDB(L,K,NX)*LOG(SEDDIA(NS))
C              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
C            ENDDO
C          ENDDO
C        ENDDO
C
C        DO K=1,KB
C          DO L=2,LA
C            IF(SNDBT(L,K).GT.0.)THEN
C              SEDDIA50(L,K)=SEDDIA50(L,K)/SNDBT(L,K)
C              SEDDIA50(L,K)=EXP(SEDDIA50(L,K)/SNDBT(L,K))
C            ELSE
C              SEDDIA50(L,K)=0.
C            ENDIF
C          ENDDO
C        ENDDO
C
C**********************************************************************C
C
C **   SET CRITICAL SHILED'S PARAMETER FOR D50
C
        DO L=2,LA
          CALL SETSHLD(DUM1,CSHIELDS50(L),SEDDIA50(L,KBT(L)),
     &                     SSG(NSED+1),DUM3,DUM4)
        ENDDO
C
CZZDIFF C **  CALCULATE CELL CENTER SHEAR STRESS  
CZZDIFF C  
CZZDIFF         DO L=2,LA  
CZZDIFF           TAUB=MAX(QQ(L,0),QQMIN)/CTURB2  
CZZDIFF           TAUB2=TAUB*TAUB+0.5*(QQWV2(L)*QQWV2(L))  
CZZDIFF           TAUB=SQRT(TAUB2)  
CZZDIFF           USTAR=SQRT(TAUB)  
CZZDIFF           CCSHEAR(L)=TAUB  
CZZDIFF           CCUSTAR(L)=USTAR  
CZZDIFF         ENDDO  
CZZDIFF C
C
C      ENDIF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.1)THEN
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM=1./(1.E6*SSG(NS))
          DIASED=SEDDIA(NS)
          DIASED3=3.*DIASED
          GPDIASED=G*(SSG(NS)-1.)*DIASED
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          K=0
C
          IF(ISNDVW.EQ.0)THEN
            DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
            ENDDO
          ENDIF
C
          IF(ISNDVW.GE.1)THEN
            DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)*
     &                     CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
            ENDDO
          ENDIF
C
          IF(IROUSE(NX).EQ.0)THEN
            DO L=2,LA
              IF(USTAR(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTAR(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ELSE
            DO L=2,LA
              IF(USTARSND(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTARSND(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ENDIF
C
            DO L=2,LA
              WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
            ENDDO
C
C----------------------------------------------------------------------C
C
          DO L=1,LC
            QSBDLDP(L)=0.
            QSBDLDX(L,NX)=0.
            QSBDLDY(L,NX)=0.
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            SNDFBL(L,NX)=0.
            RBPSBL(L)=1.0
          ENDDO
C
          CALL BEDLOAD(NX,NS)
C
C----------------------------------------------------------------------C
C
C **  SUSPENDED TRANSPORT HORIZONTAL LOOP
C
          DO L=2,LA
            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)
	    ENDDO
C
          DO L=2,LA
            ZEQ(L)=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUBSND(L),
     &        SEDDIA50(L,KBT(L)),HP(L),ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))
            ZEQMIN=0.5*DZC(1)
            ZEQ(L)=MIN(ZEQ(L),ZEQMIN)
            ZEQI(L)=1./ZEQ(L)
	    ENDDO
C
          DO L=2,LA
            SNDEQB(L)=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),
     &        TAUBSND(L),SEDDIA50(L,KBT(L)),SIGPHI(L,KBT(L)),ZEQ(L),
     &        VDRBED(L,KBT(L)),ISNDEQ(NS),ISNDAL)
	    ENDDO
C
C **  APPLIED LIMITOR TO GARCIA AND PARKER
C
          IF(ISNDEQ(NS).EQ.1)THEN
            IF(ISLTAUC(NS).EQ.1)THEN
		    CSHIELDS=TCSHIELDS(NS)
	        DO L=2,LA
	          IF(ISEDEFF.EQ.2)THEN
	            TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                  CSHIELDS=TMPVAL*CSHIELDS
	          ENDIF
	          SHIELDS=TAUBSND(L)/GPDIASED
	          IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	        ENDDO
	      ENDIF
	      IF(ISLTAUC(NS).EQ.2)THEN
	        DO L=2,LA
		      CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
	          IF(ISEDEFF.EQ.2)THEN
	            TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                  CSHIELDS=TMPVAL*CSHIELDS
	          ENDIF
	          SHIELDS=TAUBSND(L)/GPDIASED
	          IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
              ENDDO
	      ENDIF
	      IF(ISLTAUC(NS).EQ.3)THEN
	        DO L=2,LA
		      CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
	          IF(ISEDEFF.EQ.2)THEN
	            TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                  CSHIELDS=TMPVAL*CSHIELDS
	          ENDIF
	          SHIELDS=TAUBSND(L)/GPDIASED
	          IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
              ENDDO
	      ENDIF
          ENDIF
C
	    IF(ISEDEFF.EQ.1)THEN
            DO L=2,LA
              SNDEQB(L)=SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))
	      ENDDO
	    ENDIF
C

          DO L=2,LA
            IF(ROUSE(L).LT.0.999.OR.ROUSE(L).GT.1.001)THEN
              TOP=(ZEQ(L)**(ROUSE(L)-1.))-1.
              BOT=(1.-ROUSE(L))*(ZEQI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
              SNDEQSAV(L,NX)=SNDEQ(L)
             ELSE
              TOP=LOG(ZEQI(L))
              BOT=(ZEQI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
              SNDEQSAV(L,NX)=SNDEQ(L)
            ENDIF
	    ENDDO
C
          DO L=2,LA
            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)
            WESEMX=DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
C **  SET ROUSE PARAMETER AND EQUILIBRUIM CONCENTRATION
C            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTARSND(L),2)
C10/6            ZEQ=DIASED3*HPI(L)
C            D50TMP=SEDDIA50(L,KBT(L))
CXXX
C            ZEQ(L)=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUBSND(L),
C     &        SEDDIA50(L,KBT(L)),HP(L),ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))
C            ZEQMIN=0.5*DZC(1)
C            ZEQ(L)=MIN(ZEQ(L),ZEQMIN)
C            ZEQI(L)=1./ZEQ(L)
CXXX
            WSFAC=2.*(1.+ROUSE(L))/(2.+ROUSE(L)*(1.-ZEQ(L)))
C10/6            D50TMP=SEDDIA50(L,KBT(L))
C10/6            IF(ISNDAL.EQ.1) D50TMP=DIASED
CXXX
C            SIGP=SIGPHI(L,KBT(L))
C10/6            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
C10/6     &                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
C            SNDEQB(L)=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),
C     &        TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),ISNDEQ(NS),ISNDAL)
CXXX
C            IF(ROUSE(L).LT.0.999.OR.ROUSE(L).GT.1.001)THEN
C              TOP=(ZEQ(L)**(ROUSE(L)-1.))-1.
C              BOT=(1.-ROUSE(L))*(ZEQI(L)-1.)
C              SNDEQ(L)=SNDEQB(L)*TOP/BOT
C              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
C             ELSE
C              TOP=LOG(ZEQI(L))
C              BOT=(ZEQI(L)-1.)
C              SNDEQ(L)=SNDEQB(L)*TOP/BOT
C              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
C            ENDIF
C **  SET RESUSPENSION FLUX
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,0,NS)*SNDEQ(L)
C **  SET DEPOSITION VELOCITY 
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
c1104            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
            SNDBTMP=SNDB(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
     &             -DELT*SNDFBL(L,NX)
C **  ADDED BED LOAD FLUX TO SUSPENDED LOAD FLUX WITH ABOVE LINE
            IF(SNDBTMP.LT.0.0)THEN
c1104              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDF(L,0,NX)=DELTI*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            ENDIF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     &                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            SNDF(L,0,NX)=SNDF(L,0,NX)+SNDFBL(L,NX)
C
            IF(ISHOUSATONIC.EQ.1)THEN
            IF(IBMECH.EQ.1.AND.SEDVRDT.LT.0.0)THEN
              QSBDTOP(L)=QSBDTOP(L)+0.001*SNDF(L,0,NX)/SDENAVG(L,KBT(L))
              QWBDTOP(L)=QWBDTOP(L)
     &           +0.001*VDRBED(L,KBT(L))*SNDF(L,0,NX)/SDENAVG(L,KBT(L))
            ELSE
              QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)
              QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     &                  +VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
            ENDIF
            ENDIF
C
            IF(ISHOUSATONIC.EQ.0)THEN
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     &                  +VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
            ENDIF
          ENDDO
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C
C ++ FOLLOWING APPEARS IN ZZ CODE
C
C **  SUSPENDED TRANSPORT HORIZONTAL LOOP  
C  
CZZ           DO L=2,LA  
CZZ             SNDF(L,1,NX)=0.  
CZZ             PROBDEP=0.  
CZZ             WESE=0.  
CZZ             TAUBC=QQ(L,0)/CTURB2  
CZZ             TAUB=TAUBC  
CZZ             IF(TAUB.GT.0.0)THEN  
CZZ               IF(ISWAVE.GT.0)THEN  
CZZ                 UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12  
CZZ                 VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))  
CZZ                 CURANG=ATAN2(VTMP,UTMP)  
CZZ                 TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))+
CZZ      &              FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))  
CZZ                 TAUB2=MAX(TAUB2,0.)  
CZZ                 TAUB=SQRT(TAUB2)  
CZZ               ENDIF  
CZZ C **  SET MAXIMUM EROSION RATE  
CZZ               WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)  
CZZ               IF(KBT(L).NE.1)WESEMX=3.*WESEMX  
CZZ               WESEMX=MAX(WESEMX,0.)  
CZZ C **  SET ROUSE PARAMETER AND EQUILIBRUIM CONCENTRATION  
CZZ               USTAR=SQRT(TAUB)  
CZZ               ROUSE=WSETA(L,0,NS)/(VKC*USTAR)  
CZZ C           ZEQ=DIASED3*HPI(L)  
CZZ               D50TMP=SEDDIA50(L,KBT(L))  
CZZ               ZEQ=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUB,D50TMP,HP(L)
CZZ      &            ,ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))  
CZZ               ZEQMIN=0.5*DZC(1)  
CZZ               ZEQ=MIN(ZEQ,ZEQMIN)  
CZZ               ZEQI=1./ZEQ  
CZZ               WSFAC=2.*(1.+ROUSE)/(2.+ROUSE*(1.-ZEQ))  
CZZ C           D50TMP=SEDDIA50(L,KBT(L))  
CZZ C           IF(ISNDAL.EQ.1) D50TMP=DIASED  
CZZ               SIGP=SIGPHI(L,KBT(L))  
CZZ C           SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,  
CZZ C    &                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))  
CZZ               SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS)
CZZ      &            ,TAUB,D50TMP,SIGP,ZEQ,ISNDEQ(NS),ISNDAL)  
CZZ               IF(ROUSE.LT.0.999.OR.ROUSE.GT.1.001)THEN  
CZZ                 TOP=(ZEQ**(ROUSE-1.))-1.  
CZZ                 BOT=(1.-ROUSE)*(ZEQI-1.)  
CZZ                 SNDEQ=SNDEQB*TOP/BOT  
CZZ                 SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)  
CZZ               ELSE  
CZZ                 TOP=LOG(ZEQI)  
CZZ                 BOT=(ZEQI-1.)  
CZZ                 SNDEQ=SNDEQB*TOP/BOT  
CZZ                 SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)  
CZZ               ENDIF  
CZZ C **  SET RESUSPENSION FLUX  
CZZ               WESE=WSFAC*CTMPDRY(L)*WSETA(L,0,NS)*SNDEQ  
CZZ C **  SET DEPOSITION VELOCITY  
CZZ               WSETMP=WSFAC*WSETA(L,0,NS)  
CZZ               WVEL=DELT*HPI(L)*DZIC(1)  
CZZ               CLEFT=1.+WSETMP*WVEL  
CZZ               CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL  
CZZ               SND(L,1,NX)=CRIGHT/CLEFT  
CZZ               SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE  
CZZ             ELSE  
CZZ               SNDF(L,0,NX)=0.0  
CZZ             ENDIF  
CZZ             SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)-DELT*SNDFBL(L,NX)  
CZZ C **  ADDED BED LOAD FLUX TO SUSPENDED LOAD FLUX WITH ABOVE LINE  
CZZ             IF(SNDBTMP.LT.0.0)THEN  
CZZ               SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)-SNDFBL(L,NX)  
CZZ               SNDBTMP=0.0  
CZZ               SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL  
CZZ             ENDIF  
CZZ             SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)+S2TL*SNDB1(L,
CZZ      &          KBT(L),NX)  
CZZ             SNDB(L,KBT(L),NX)=SNDBTMP  
CZZ             SNDF(L,0,NX)=SNDF(L,0,NX)+SNDFBL(L,NX)  
CZZ             QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)  
CZZ             QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*(VDRBED(L,KBT(L))*MAX(SNDF(
CZZ      &          L,0,NX),0.)+SNDVDRD*MIN(SNDF(L,0,NX),0.))  
CZZ           ENDDO  
CZZ C  
CZZ         ENDDO  
CZZ       ENDIF
C
C ++ END APPEARS IN ZZ CODE
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.2)THEN
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM=1./(1.E6*SSG(NS))
          DIASED=SEDDIA(NS)
          DIASED3=3.*DIASED
          GPDIASED=G*(SSG(NS)-1.)*DIASED
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISNDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISNDVW.GE.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)*
     &                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(IROUSE(NX).EQ.0)THEN
            DO L=2,LA
              IF(USTAR(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTAR(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ELSE
            DO L=2,LA
              IF(USTARSND(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTARSND(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ENDIF
C
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
              ENDDO
            ENDDO
C
C----------------------------------------------------------------------C
C
          DO L=1,LC
            QSBDLDP(L)=0.
            QSBDLDX(L,NX)=0.
            QSBDLDY(L,NX)=0.
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            SNDFBL(L,NX)=0.
            RBPSBL(L)=1.0
          ENDDO
C
          CALL BEDLOAD(NX,NS)
C
C----------------------------------------------------------------------C
C
C **  HORIZONTAL LOOPS
C
          K=2
          DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
          ENDDO
C
          DO L=2,LA
            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)
	    ENDDO
C
          DO L=2,LA
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)
            WESEMX=DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
C **  SET ROUSE PARAMETER AND EQUILIBRUIM CONCENTRATION
c            USTARSND(L)=SQRT(TAUBSND(L))
C            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTARSND(L),2)
C10/6            ZEQ=DIASED3*HPI(L)
C            D50TMP=SEDDIA50(L,KBT(L))
            ZEQ(L)=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUBSND(L),
     &        SEDDIA50(L,KBT(L)),HP(L),ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))
            ZEQMIN=0.5*DZC(1)
            ZEQ(L)=MIN(ZEQ(L),ZEQMIN)
            ZEQD(L)=ZEQ(L)/DZC(1)
            ZEQDI(L)=1./ZEQD(L)
            WSFAC=1.
C10/6            D50TMP=SEDDIA50(L,KBT(L))
C10/6            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
C10/6            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
C10/6     &                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            SNDEQB(L)=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),
     &      TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),
     &      VDRBED(L,KBT(L)),ISNDEQ(NS),ISNDAL)
C
C **  APPLIED LIMITOR TO GARCIA AND PARKER
C
          IF(ISNDEQ(NS).EQ.1)THEN
            IF(ISLTAUC(NS).EQ.1)THEN
		    CSHIELDS=TCSHIELDS(NS)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.2)THEN
		    CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.3)THEN
		    CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
          ENDIF
C
	    IF(ISEDEFF.EQ.1)
     &      SNDEQB(L)=SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))
C
            IF(ROUSE(L).LT.0.999.OR.ROUSE(L).GT.1.001)THEN
              TOP=(ZEQD(L)**(ROUSE(L)-1.))-1.
              BOT=(1.-ROUSE(L))*(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ELSE
              TOP=LOG(ZEQDI(L))
              BOT=(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ENDIF
C **  SET RESUSPENSION FLUX
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,0,NS)*SNDEQ(L)
C **  SET DEPOSITION VELOCITY 
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
c1104            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
            SNDBTMP=SNDB(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
     &             -DELT*SNDFBL(L,NX)
C **  ADDED BED LOAD FLUX TO SUSPENDED LOAD FLUX WITH ABOVE LINE
            IF(SNDBTMP.LT.0.0)THEN
c1104              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDF(L,0,NX)=DELTI*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            ENDIF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     &                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            SNDF(L,0,NX)=SNDF(L,0,NX)+SNDFBL(L,NX)
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     &                  +VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.EQ.2 
C
          IF(ISTOPT(7).EQ.1)THEN
C
            DO L=2,LA
              CRNUM=1.+DELT*WSETA(L,1,NS)*HPI(L)*DZIC(KC)
              GRADSED=(SND(L,KC,NX)-SND(L,1,NX))/(DZC(KC)+DZC(1))
              SEDAVG=0.5*(SND(L,KC,NX)+SND(L,KC,NX)+1.E-16)
              WSETA(L,1,NS)=-CRNUM*DZC(KC)*WSETA(L,1,NS)*GRADSED/SEDAVG
            ENDDO
C
            DO L=2,LA
              AA11=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
              AA12=-MAX(WSETA(L,1,NS),0.)
              AA21=MIN(WSETA(L,1,NS),0.)
              AA22=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,1,NS),0.)
              BB11=DELTI*DZC(1)*HP(L)*SND(L,1,NX)
              BB22=DELTI*DZC(KC)*HP(L)*SND(L,KC,NX)
              DETI=1./(AA11*AA22-AA12*AA21)
              SND(L,1,NX)=DETI*( BB11*AA22-BB22*AA12 )
              SND(L,KC,NX)=DETI*( AA11*BB22-AA21*BB11 )
            ENDDO  
C
          ENDIF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=2
C
          DO L=2,LA
            SNDF(L,1,NX)=DELTI*DZC(KC)*HP(L)*(SND(L,KC,NX)
     &         -SNDS(L,KC,NX))
          ENDDO  
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.3)THEN
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM=1./(1.E6*SSG(NS))
          DIASED=SEDDIA(NS)
          DIASED3=3.*DIASED
          GPDIASED=G*(SSG(NS)-1.)*DIASED
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISNDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISNDVW.GE.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)*
     &                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(IROUSE(NX).EQ.0)THEN
            DO L=2,LA
              IF(USTAR(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTAR(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ELSE
            DO L=2,LA
              IF(USTARSND(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTARSND(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ENDIF
C
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
              ENDDO
            ENDDO
C
C----------------------------------------------------------------------C
C
          DO L=1,LC
            QSBDLDP(L)=0.
            QSBDLDX(L,NX)=0.
            QSBDLDY(L,NX)=0.
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            SNDFBL(L,NX)=0.
            RBPSBL(L)=1.0
          ENDDO
C
          CALL BEDLOAD(NX,NS)
C
C----------------------------------------------------------------------C
C
C **  HORIZONTAL LOOPS
C
          K=3
          DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
          ENDDO
C
          K=2
          DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
          ENDDO
C
          DO L=2,LA
            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)
	    ENDDO
C
          DO L=2,LA
c fix            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)
            WESEMX=DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
C **  SET ROUSE PARAMETER AND EQUILIBRUIM CONCENTRATION
c            TMPCGSQRT(TAUBSND(L))
C            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTARSND(L),2)
C10/6            ZEQ=DIASED3*HPI(L)
C            D50TMP=SEDDIA50(L,KBT(L))
            ZEQ(L)=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUBSND(L),
     &        SEDDIA50(L,KBT(L)),HP(L),ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))
            ZEQMIN=0.5*DZC(1)
            ZEQ(L)=MIN(ZEQ(L),ZEQMIN)
            ZEQD(L)=ZEQ(L)/DZC(1)
            ZEQDI(L)=1./ZEQD(L)
            WSFAC=1.
C10/6            D50TMP=SEDDIA50(L,KBT(L))
C10/6            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
C10/6            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
C10/6     &                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            SNDEQB(L)=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),
     &      TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),
     &      VDRBED(L,KBT(L)),ISNDEQ(NS),ISNDAL)
C
C **  APPLIED LIMITOR TO GARCIA AND PARKER
C
          IF(ISNDEQ(NS).EQ.1)THEN
            IF(ISLTAUC(NS).EQ.1)THEN
		    CSHIELDS=TCSHIELDS(NS)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.2)THEN
		    CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.3)THEN
		    CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
          ENDIF
C

	    IF(ISEDEFF.EQ.1)
     &      SNDEQB(L)=SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))
C
            IF(ROUSE(L).LT.0.999.OR.ROUSE(L).GT.1.001)THEN
              TOP=(ZEQD(L)**(ROUSE(L)-1.))-1.
              BOT=(1.-ROUSE(L))*(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ELSE
              TOP=LOG(ZEQDI(L))
              BOT=(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ENDIF
C **  SET RESUSPENSION FLUX
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,0,NS)*SNDEQ(L)
C **  SET DEPOSITION VELOCITY 
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
c1104            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
            SNDBTMP=SNDB(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
     &             -DELT*SNDFBL(L,NX)
C **  ADDED BED LOAD FLUX TO SUSPENDED LOAD FLUX WITH ABOVE LINE
            IF(SNDBTMP.LT.0.0)THEN
c1104              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDF(L,0,NX)=DELTI*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            ENDIF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     &                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            SNDF(L,0,NX)=SNDF(L,0,NX)+SNDFBL(L,NX)
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     &                  +VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.EQ.3
C
          IF(ISTOPT(7).EQ.1)THEN
C
            DO K=1,2
              DO L=2,LA
                CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
                GRADSED=(SND(L,K+1,NX)-SND(L,K,NX))/(DZC(K+1)+DZC(K))
                SEDAVG=0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)
                WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/
     &                         SEDAVG
              ENDDO
            ENDDO
C
C     TVAR1S=LOWER DIAGONAL
            DO L=2,LA
              TVAR1S(L,1)=0
            ENDDO
            DO K=2,KC
              DO L=2,LA
                TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
              ENDDO
            ENDDO
C     TVAR1N=UPPER DIAGONAL
            DO L=2,LA
              TVAR1N(L,KC)=0
            ENDDO
            DO K=1,KS
              DO L=2,LA
                TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
              ENDDO
            ENDDO
C     TVAR1W=MAIN DIAGONAL
            DO L=2,LA
              TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
              TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
            ENDDO
            DO K=2,KS
              DO L=2,LA
                TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     &                -MIN(WSETA(L,K,NS),0.)
              ENDDO
            ENDDO
C     TVAR1E=RIGHT HAND SIDE
            DO K=1,KC
              DO L=2,LA
                TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SND(L,K,NX)
              ENDDO
            ENDDO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
            DO L=2,LA
              TVAR3S(L)=TVAR1W(L,1)
            ENDDO
            DO L=2,LA
              TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
            ENDDO
            DO K=2,KC
              DO L=2,LA
                TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
                TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
                TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     &                   TVAR3S(L)
              ENDDO
            ENDDO
            DO K=KS,1,-1
              DO L=2,LA
                TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
              ENDDO
            ENDDO
            DO K=1,KC
              DO L=2,LA
                SND(L,K,NX)=TVAR2N(L,K)
              ENDDO
            ENDDO
C
          ENDIF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=3
C
          DO L=2,LA
            SNDF(L,KC-1,NX)=DELTI*DZC(KC)*HP(L)*(SND(L,KC,NX)
     &                -SNDS(L,KC,NX))
          ENDDO  
          DO L=2,LA
            SNDF(L,1,NX)=DELTI*DZC(KC-1)*HP(L)*(SND(L,KC-1,NX)
     &                 -SNDS(L,KC-1,NX))+SNDF(L,KC-1,NX)
          ENDDO  
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC.GT.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.GT.3)THEN
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM=1./(1.E6*SSG(NS))
          DIASED=SEDDIA(NS)
          DIASED3=3.*DIASED
          GPDIASED=G*(SSG(NS)-1.)*DIASED
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISNDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISNDVW.GE.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)*
     &                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(IROUSE(NX).EQ.0)THEN
            DO L=2,LA
              IF(USTAR(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTAR(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ELSE
            DO L=2,LA
              IF(USTARSND(L).GT.0.0) THEN
                ROUSE(L)=WSETA(L,0,NS)/(VKC*USTARSND(L))
              ELSE
                ROUSE(L)=250000.*WSETA(L,0,NS)
              END IF
	      ENDDO
          ENDIF
C
C
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
              ENDDO
            ENDDO
C
C----------------------------------------------------------------------C
C
          DO L=1,LC
            QSBDLDP(L)=0.
            QSBDLDX(L,NX)=0.
            QSBDLDY(L,NX)=0.
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            SNDFBL(L,NX)=0.
            RBPSBL(L)=1.0
          ENDDO
C
          CALL BEDLOAD(NX,NS)
C
C----------------------------------------------------------------------C
C
C **  HORIZONTAL LOOP
C
          K=KC
          DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
          ENDDO
C
          DO K=KS,2,-1
            DO L=2,LA
              WVEL=DELT*HPI(L)*DZIC(K)
              CLEFT=1.+WSETA(L,K-1,NS)*WVEL
              CRIGHT=MAX(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL
              SND(L,K,NX)=CRIGHT/CLEFT
              SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
            ENDDO
          ENDDO
C
          DO L=2,LA
            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)
	    ENDDO
C
          DO L=2,LA
c            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)
            WESEMX=DELTI*CTMPDRY(L)*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
C **  SET ROUSE PARAMETER AND EQUILIBRUIM CONCENTRATION
c            TMPCGSQRT(TAUBSND(L))
C            FACSUSL(L)=FSEDMODE(WSETA(L,0,NS),USTARSND(L),2)
C10/6            ZEQ=DIASED3*HPI(L)
C            D50TMP=SEDDIA50(L,KBT(L))
            ZEQ(L)=CSNDZEQ(DIASED,GPDIASED,TAUR(NS),TAUBSND(L),
     &        SEDDIA50(L,KBT(L)),HP(L),ISNDEQ(NS),SSG(NS),WSETA(L,0,NS))
            ZEQMIN=0.5*DZC(1)
            ZEQ(L)=MIN(ZEQ(L),ZEQMIN)
            ZEQD(L)=ZEQ(L)/DZC(1)
            ZEQDI(L)=1./ZEQD(L)
            WSFAC=1.
C10/6            D50TMP=SEDDIA50(L,KBT(L))
C10/6            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
C10/6            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
C10/6     &                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            SNDEQB(L)=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),
     &      TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),
     &      VDRBED(L,KBT(L)),ISNDEQ(NS),ISNDAL)
C
C **  APPLIED LIMITOR TO GARCIA AND PARKER
C
          IF(ISNDEQ(NS).EQ.1)THEN
            IF(ISLTAUC(NS).EQ.1)THEN
		    CSHIELDS=TCSHIELDS(NS)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.2)THEN
		    CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
	      IF(ISLTAUC(NS).EQ.3)THEN
		    CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
	        IF(ISEDEFF.EQ.2)THEN
	          TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
                CSHIELDS=TMPVAL*CSHIELDS
	        ENDIF
	        SHIELDS=TAUBSND(L)/GPDIASED
	        IF(SHIELDS.LT.CSHIELDS) SNDEQB(L)=0.
	      ENDIF
          ENDIF
C
	    IF(ISEDEFF.EQ.1)
     &        SNDEQB(L)=SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))
C
            IF(ROUSE(L).LT.0.999.OR.ROUSE(L).GT.1.001)THEN
              TOP=(ZEQD(L)**(ROUSE(L)-1.))-1.
              BOT=(1.-ROUSE(L))*(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ELSE
              TOP=LOG(ZEQDI(L))
              BOT=(ZEQDI(L)-1.)
              SNDEQ(L)=SNDEQB(L)*TOP/BOT
              SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)
            ENDIF
C **  SET RESUSPENSION FLUX
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,0,NS)*SNDEQ(L)
C **  SET DEPOSITION VELOCITY 
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
c1104            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
            SNDBTMP=SNDB(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
     &             -DELT*SNDFBL(L,NX)
C **  ADDED BED LOAD FLUX TO SUSPENDED LOAD FLUX WITH ABOVE LINE
            IF(SNDBTMP.LT.0.0)THEN
c1104              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDF(L,0,NX)=DELTI*SNDB(L,KBT(L),NX)-SNDFBL(L,NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            ENDIF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     &                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            SNDF(L,0,NX)=SNDF(L,0,NX)+SNDFBL(L,NX)
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SNDF(L,0,NX)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     &                  +VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.GT.3 
C
          IF(ISTOPT(7).EQ.1)THEN
C
            DO K=1,KS
              DO L=2,LA
                CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
                GRADSED=(SND(L,K+1,NX)-SND(L,K,NX))/(DZC(K+1)+DZC(K))
                SEDAVG=0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)
                WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/
     &                        SEDAVG
              ENDDO
            ENDDO
C
C     TVAR1S=LOWER DIAGONAL
            DO L=2,LA
              TVAR1S(L,1)=0
            ENDDO
            DO K=2,KC
              DO L=2,LA
                TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
              ENDDO
            ENDDO
C     TVAR1N=UPPER DIAGONAL
            DO L=2,LA
              TVAR1N(L,KC)=0
            ENDDO
            DO K=1,KS
              DO L=2,LA
                TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
              ENDDO
            ENDDO
C     TVAR1W=MAIN DIAGONAL
            DO L=2,LA
              TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
              TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
            ENDDO
            DO K=2,KS
              DO L=2,LA
                TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     &                -MIN(WSETA(L,K,NS),0.)
              ENDDO
            ENDDO
C     TVAR1E=RIGHT HAND SIDE
            DO K=1,KC
              DO L=2,LA
                TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SND(L,K,NX)
              ENDDO
            ENDDO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
            DO L=2,LA
              TVAR3S(L)=TVAR1W(L,1)
            ENDDO
            DO L=2,LA
             TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
            ENDDO
            DO K=2,KC
              DO L=2,LA
                TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
                TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
                TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     &                   TVAR3S(L)
              ENDDO
            ENDDO
            DO K=KS,1,-1
              DO L=2,LA
                TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
              ENDDO
            ENDDO
            DO K=1,KC
              DO L=2,LA
                SND(L,K,NX)=TVAR2N(L,K)
              ENDDO
            ENDDO
C
          ENDIF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC.GE.3
C
          DO L=2,LA
            SNDF(L,KS,NX)=DELTI*DZC(KC)*HP(L)*
     &                   (SND(L,KC,NX)-SNDS(L,KC,NX))
          ENDDO 
C 
          DO K=KS-1,1,-1
            DO L=2,LA
              SNDF(L,K,NX)=DELTI*DZC(K+1)*HP(L)*
     &                   (SND(L,K+1,NX)-SNDS(L,K+1,NX))+SNDF(L,K+1,NX)
            ENDDO  
          ENDDO  
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
        OPEN(1,FILE='NEGSEDSND.OUT',POSITION='APPEND')
C
        DO NS=1,NSND
          DO K=1,KC
            DO L=2,LA
              IF(SND(L,K,NS).LT.0.)THEN
                WRITE(1,107)TIME,NS,IL(L),JL(L),K,SND(L,K,NS)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        DO NS=1,NSND
          DO L=2,LA
            IF(SNDB(L,KBT(L),NS).LT.0.)THEN
              WRITE(1,108)TIME,NS,IL(L),JL(L),KBT(L),SNDB(L,KBT(L),NS),
     &             SNDF(L,0,NS)
            ENDIF
          ENDDO
        ENDDO
C
        CLOSE(1)
C
C **  ACCUMULATE NET POSTIVE AND NEGATIVE NONCOHESIVE SEDIMENT FLUXES
C
        DO NS=1,NSND
          DO L=2,LA
            SNDFDTAP(L,NS)=SNDFDTAP(L,NS)+DELT*MAX(SNDF(L,0,NS),0.0)
            SNDFDTAN(L,NS)=SNDFDTAN(L,NS)+DELT*MIN(SNDF(L,0,NS),0.0)
          ENDDO
        ENDDO
C
  107 FORMAT(' TIME,NS,I,J,K,NEGSND = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' TIME,NS,I,J,NEGSNDB,SNDF = ',F12.4,4I5,4E13.4)       
C
C**********************************************************************C
C
      RETURN
      END
