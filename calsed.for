C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALSED
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
C **  SUBROUTINE CALSED CALCULATES COHESIVER SEDIMENT SETTLING,
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
     &                 DSTRSE(LCM,KBM),DZBTR(LCM,KBM),STRSEM(LCM,KBM)
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
C **  COHESIVE SEDIMENT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.1)THEN
        DO NS=1,NSED
          DSEDGMM=1./(1.E6*SSG(NS))
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          K=0
C
          IF(ISEDVW.EQ.0)THEN
            DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.1)THEN
            DO L=2,LA
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),0.0,ISEDVW)
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.2)THEN
	      IF(ISWAVE.EQ.0)THEN
              DO L=2,LA
                STRESS=QQ(L,0)/CTURB2
                SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
              ENDDO
	      ELSE
              DO L=2,LA
                TAUBC=QQ(L,0)/CTURB2
                UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
                VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
                CURANG=ATAN2(VTMP,UTMP)
                TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
                TAUB2=MAX(TAUB2,0.)
                STRESS=SQRT(TAUB2)
                SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
	        ENDDO
	      ENDIF
          ENDIF
C
          IF(ISEDVW.GE.3.AND.ISEDVW.LE.4)THEN
	      IF(ISWAVE.EQ.0)THEN
              DO L=2,LA
                STRESS=0.5*QQ(L,0)/CTURB2
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
              ENDDO
	      ELSE
              DO L=2,LA
                TAUBC=QQ(L,0)/CTURB2
                UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
                VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
                CURANG=ATAN2(VTMP,UTMP)
                TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &             +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
                TAUB2=MAX(TAUB2,0.)
                STRESS=0.5*SQRT(TAUB2)+0.5*QQ(L,1)/CTURB2
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
              ENDDO
	      ENDIF
          ENDIF
C
          IF(ISEDVW.EQ.5)THEN
            DO L=2,LA
              UUSTARTMP=QCELLCTR(L)*SQRT(QQ(L,0)/CTURB2)
              STRESS=SQRT(HPI(L)*HPI(L)*UUSTARTMP)
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
            ENDDO
          ENDIF
C
          IF(ISEDVW.GE.1)THEN
            DO L=2,LA
              WSETA(L,K,NS)=MAX(WSETA(L,K,NS),WSEDO(NS))
            ENDDO
          ENDIF
C
            DO L=2,LA
              WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
            ENDDO
c
C----------------------------------------------------------------------C
C
C **  HORIZONTAL LOOP
C
          DO L=2,LA
            SEDF(L,1,NS)=0.
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            WESEMX=DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUBSED(L).GT.TAURB(L,KBT(L)))THEN
C **  MASS EROSION
              WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
              WESE=MIN(WESE,WESEMX)
            ELSE
              IF(TAUBSED(L).GT.TAURS(L,KBT(L)))THEN
C **  SURFACE EROSION
                WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
c1104                WESE=MIN(WESE,WESEMX)
                TAURTMP=TAURS(L,KBT(L))
C#######################################################################
C     HQI Change, 01/05/04, SO and RM
C     Change to implement HQI Sed-flume analysis based critical shear 
C     stress options
c                IF(IWRSP(1).GE.2) TAURTMP=TAUR(1)
                IF((IWRSP(1).GE.2).AND.(IWRSP(1).LT.99)) TAURTMP=TAUR(1)
cjh24dec2008
	          IF(IWRSP(1).ge.99)THEN
	            TAURTMP=TAUNS(L,KBT(L))
	          ENDIF
C#######################################################################
                TAUE=(TAUBSED(L)-TAURS(L,KBT(L)))/TAURTMP
                TAUE=MAX(TAUE,0.0)
	          IF(ISLTAUC(NS).EQ.1.AND.TAUBSED(L).LT.TAURTMP) TAUE=0.
	          TMPSEDHID=1.0
	          IF(ISTRAN(7).GE.1.AND.COSEDHID(1).NE.0.0) 
     &            TMPSEDHID=(FRACCOH(L,KBT(L)))**COSEDHID(1)
C#######################################################################
C     HQI change to implement spatially varying coefficient and exponent
C     in resuspension formulation based on Ed G analysis of Sed-flume data
C     RM 12/11/03
C     HAMRICK MODIFIED FOR GENERALITY
C
                IF(IWRSP(1).LT.99)THEN
                      WESE=TMPSEDHID*WESE*(TAUE**TEXP(NS))
                ELSE

cSO               IF(L.LE.265) then  ! Woods Pond
cSO                 WESE=TMPSEDHID*WESE*(TAUE**1.04)
cSO               ELSE               ! North of Woods Pond
cSO                WESE=TMPSEDHID*WESE*(TAUE**0.886)
cSO               ENDIF

cSO   HQI change to implement spatially varying coefficient and exponent
cSO   in resuspension formulation based re-analysis of Sedflume data with
cSO   a new partition of domain (see SSEDTOX.FOR, line 587)
cSO 05/12/04
c                 IF ( L.LE.2042 ) THEN
c                   WESE=TMPSEDHID*WESE*(TAUE**0.949) ! All of 5C
c                 ELSE                                !      + Woods Pond
c                   WESE=TMPSEDHID*WESE*(TAUE**1.31 ) ! North of (All of
c                 END IF                              !   5C + Woods Pond)

cSO 05/14/04
CJH                  IF ( L.LE.2042 ) THEN
CJH                    WESE=TMPSEDHID*WESE*(TAUE**0.949) ! All of 5C
CJH                  ELSE                                !      + Woods Pond
CJH                    WESE=TMPSEDHID*WESE*(TAUE**1.59 ) ! North of (All of
CJH                  END IF                              !   5C + Woods Pond)

                      WESE=TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )

                ENDIF
C#######################################################################
c1104
                WESE=MIN(WESE,WESEMX)
c1104
              ELSE
C **  NO EROSION 
                WESE=0.0
              ENDIF
            ENDIF
C **  SET PROBABILITY OF DEPOSITION 
C#######################################################################
C     HQI Change, 04/25/04
C     Change to compute probability of deposition using total bed shear 
C     rather than cohesive grain shear
c            IF(TAUBSED(L).LT.TAUD(NS)) 
c     &        PROBDEP=(TAUD(NS)-TAUBSED(L))/TAUD(NS)
C
            IF(IWRSP(1).LT.99)THEN
              IF(TAUB(L).LT.TAUD(NS)) 
     &           PROBDEP=(TAUD(NS)-TAUB(L))/TAUD(NS)
            ELSE
              IF(TAUB(L).LT.TAUDS(L)) 
     &           PROBDEP=(TAUDS(L)-TAUB(L))/TAUDS(L)
            ENDIF
C
C#######################################################################
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
cjmh216            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
            SEDBTMP=SEDB(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
C           IF(SEDBTMP.LT.0.0)THEN
C             SEDF(L,0,NS)=0.
C             SEDBTMP=SEDB1(L,KBT(L),NS)
C             SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
C           ENDIF
C           SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
C     &                        +S2TL*SEDB1(L,KBT(L),NS)
C           SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0)THEN
c1104              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
c1104
              SEDF(L,0,NS)=DELTI*SEDB(L,KBT(L),NS)
c1104
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)+(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            ENDIF
cjmh216            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
cjmh216     &                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB1(L,KBT(L),NS)=SEDB(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     &                  +VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.2)THEN
        DO NS=1,NSED
          DSEDGMM=1./(1.E6*SSG(NS))
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISEDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C 
          IF(ISEDVW.EQ.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),0.0,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.2)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
            ENDDO
            K=1
            DO L=2,LA
              LN=LNC(L)
              SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     &              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
            ENDDO
          ENDIF
C
          IF(ISEDVW.GE.3)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &            +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
            ENDDO
            K=1
            DO L=2,LA
              LN=LNC(L)
              STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     &              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 ) 
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
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
C **  HORIZONTAL LOOPS
C
          K=2
          DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
          ENDDO
C
          DO L=2,LA
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            WESEMX=DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUBSED(L).GT.TAURB(L,KBT(L)))THEN
C **  MASS EROSION
              WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
              WESE=MIN(WESE,WESEMX)
            ELSE
              IF(TAUBSED(L).GT.TAURS(L,KBT(L)))THEN
C **  SURFACE EROSION
                WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
c1104                WESE=MIN(WESE,WESEMX)
                TAURTMP=TAURS(L,KBT(L))
	          IF(IWRSP(1).GE.2)TAURTMP=TAUR(1)
		          IF(IWRSP(1).ge.99)TAURTMP=TAUNS(L,KBT(L))
                TAUE=(TAUBSED(L)-TAURS(L,KBT(L)))/TAURTMP
                TAUE=MAX(TAUE,0.0)
	          IF(ISLTAUC(NS).EQ.1.AND.TAUBSED(L).LT.TAURTMP) TAUE=0.
			  TMPSEDHID=1.0
	          IF(ISTRAN(7).GE.1.AND.COSEDHID(1).NE.0.0) 
     &            TMPSEDHID=(FRACCOH(L,KBT(L)))**COSEDHID(1)
                IF(IWRSP(1).LT.99)THEN
                  WESE=TMPSEDHID*WESE*(TAUE**TEXP(NS))
                ELSE
                  WESE=TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
                ENDIF
c1104
                WESE=MIN(WESE,WESEMX)
c1104
              ELSE
C **  NO EROSION 
                WESE=0.0
              ENDIF
            ENDIF
C **  SET PROBABILITY OF DEPOSITION 
            TAUBHYDRO=QQ(L,0)/CTURB2
            IF(ISPROBDEP(NS).EQ.0)THEN
              IF(TAUBSED(L).LT.TAUD(NS)) 
     &          PROBDEP=(TAUD(NS)-TAUBSED(L))/TAUD(NS)
            ENDIF
            IF(ISPROBDEP(NS).EQ.1)THEN
              IF(TAUBHYDRO.LT.TAUD(NS)) 
     &          PROBDEP=(TAUD(NS)-TAUBHYDRO)/TAUD(NS)
            ENDIF
            IF(ISPROBDEP(NS).EQ.2)THEN
	        IF(TAUBSED(L).LE.TAUD(NS))THEN
                PROBDEP=1.0
	        ELSE
                PROBDEP=FPROBDEP(TAUD(NS),TAUBSED(L))
              ENDIF
            ENDIF
            IF(ISPROBDEP(NS).EQ.3)THEN
	        IF(TAUBHYDRO.LE.TAUD(NS))THEN
                PROBDEP=1.0
	        ELSE
                PROBDEP=FPROBDEP(TAUD(NS),TAUBHYDRO)
              ENDIF
            ENDIF
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
cjmh216            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
            SEDBTMP=SEDB(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
C           IF(SEDBTMP.LT.0.0)THEN
C             SEDF(L,0,NS)=0.
C             SEDBTMP=SEDB1(L,KBT(L),NS)
C             SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
C           ENDIF
C           SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
C     &                        +S2TL*SEDB1(L,KBT(L),NS)
C           SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0)THEN
c1104              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
c1104
              SEDF(L,0,NS)=DELTI*SEDB(L,KBT(L),NS)
c1104
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)
     &                       +(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            ENDIF
cjmh216            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
cjmh216     &                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB1(L,KBT(L),NS)=SEDB(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     &                  +VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.EQ.2 
C
          IF(ISTOPT(6).EQ.1)THEN
C
            DO L=2,LA
              CRNUM=1.+DELT*WSETA(L,1,NS)*HPI(L)*DZIC(KC)
              GRADSED=(SED(L,KC,NS)-SED(L,1,NS))/(DZC(KC)+DZC(1))
              SEDAVG=0.5*(SED(L,KC,NS)+SED(L,1,NS)+1.E-16)
              WSETA(L,1,NS)=-CRNUM*DZC(KC)*WSETA(L,1,NS)*GRADSED/SEDAVG
            ENDDO
C
            DO L=2,LA
              AA11=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
              AA12=-MAX(WSETA(L,1,NS),0.)
              AA21=MIN(WSETA(L,1,NS),0.)
              AA22=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,1,NS),0.)
              BB11=DELTI*DZC(1)*HP(L)*SED(L,1,NS)
              BB22=DELTI*DZC(KC)*HP(L)*SED(L,KC,NS)
              DETI=1./(AA11*AA22-AA12*AA21)
              SED(L,1,NS)=DETI*( BB11*AA22-BB22*AA12 )
              SED(L,KC,NS)=DETI*( AA11*BB22-AA21*BB11 )
            ENDDO 
C
          ENDIF 
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=2
C
          DO L=2,LA
            SEDF(L,1,NS)=DELTI*DZC(KC)*HP(L)*(SED(L,KC,NS)
     &     -SEDS(L,KC,NS))
          ENDDO  
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.3)THEN
        DO NS=1,NSED
          DSEDGMM=1./(1.E6*SSG(NS))
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISEDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),0.0,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.2)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
            ENDDO
            DO K=1,KS
              DO L=2,LA
                LN=LNC(L)
                SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     &              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.GE.3)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
            ENDDO
            DO K=1,KS
              DO L=2,LA
                LN=LNC(L)
                STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     &              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
              ENDDO
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
C **  HORIZONTAL LOOPS
C
          K=3
          DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
          ENDDO
C
          K=2
          DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)-SEDF(L,K,NS)*WVEL
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
          ENDDO
C
          DO L=2,LA
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            WESEMX=DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            IF(TAUBSED(L).GT.TAURB(L,KBT(L)))THEN
C **  MASS EROSION
               WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
               WESE=MIN(WESE,WESEMX)
            ELSE
              IF(TAUBSED(L).GT.TAURS(L,KBT(L)))THEN
C **  SURFACE EROSION
                WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
c1104                WESE=MIN(WESE,WESEMX)
                TAURTMP=TAURS(L,KBT(L))
	          IF(IWRSP(1).GE.2)TAURTMP=TAUR(1)
	          IF(IWRSP(1).GE.99)TAURTMP=TAUNS(L,KBT(L))
                TAUE=(TAUBSED(L)-TAURS(L,KBT(L)))/TAURTMP
                TAUE=MAX(TAUE,0.0)
	          IF(ISLTAUC(NS).EQ.1.AND.TAUBSED(L).LT.TAURTMP) TAUE=0.
                TMPSEDHID=1.0
	          IF(ISTRAN(7).GE.1.AND.COSEDHID(1).NE.0.0) 
     &            TMPSEDHID=(FRACCOH(L,KBT(L)))**COSEDHID(1)
                IF(IWRSP(1).LT.99)THEN
                  WESE=TMPSEDHID*WESE*(TAUE**TEXP(NS))
                ELSE
                  WESE=TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
                ENDIF
c1104
                WESE=MIN(WESE,WESEMX)
c1104
              ELSE
C **  NO EROSION 
                WESE=0.0
              ENDIF
            ENDIF
C **  SET PROBABILITY OF DEPOSITION 
            IF(TAUBSED(L).LT.TAUD(NS)) 
     &         PROBDEP=(TAUD(NS)-TAUBSED(L))/TAUD(NS)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
cjmh216            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
            SEDBTMP=SEDB(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
C           IF(SEDBTMP.LT.0.0)THEN
C             SEDF(L,0,NS)=0.
C             SEDBTMP=SEDB1(L,KBT(L),NS)
C             SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
C           ENDIF
C           SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
C     &                        +S2TL*SEDB1(L,KBT(L),NS)
C           SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0)THEN
c1104              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
c1104
              SEDF(L,0,NS)=DELTI*SEDB(L,KBT(L),NS)
c1104
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)
     &                       +(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            ENDIF
cjmh216            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
cjmh216     &                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB1(L,KBT(L),NS)=SEDB(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     &                  +VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
          ENDDO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.EQ.3
C
          IF(ISTOPT(6).EQ.1)THEN
C
            DO K=1,2
              DO L=2,LA
                CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
                GRADSED=(SED(L,K+1,NS)-SED(L,K,NS))/(DZC(K+1)+DZC(K))
                SEDAVG=0.5*(SED(L,K+1,NS)-SED(L,K,NS)+1.E-16)
                WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*
     &                         GRADSED/SEDAVG
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
                TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SED(L,K,NS)
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
                SED(L,K,NS)=TVAR2N(L,K)
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
            SEDF(L,KC-1,NS)=DELTI*DZC(KC)*HP(L)*(SED(L,KC,NS)
     &     -SEDS(L,KC,NS))
          ENDDO  
          DO L=2,LA
            SEDF(L,1,NS)=DELTI*DZC(KC-1)*HP(L)*(SED(L,KC-1,NS)
     &                 -SEDS(L,KC-1,NS))+SEDF(L,KC-1,NS)
          ENDDO  
C
C----------------------------------------------------------------------C
C
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC.GT.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.GT.3)THEN
        DO NS=1,NSED
          DSEDGMM=1./(1.E6*SSG(NS))
C
C----------------------------------------------------------------------C
C
C **  SET SETTLING VELOCITIES
C
          IF(ISEDVW.EQ.0)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=WSEDO(NS)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.1)THEN
            DO K=0,KS
              DO L=2,LA
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),0.0,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.EQ.2)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &             +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
            ENDDO
            DO K=1,KS
              DO L=2,LA
                LN=LNC(L)
                SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     &              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                   +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),SHEAR,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.GE.3)THEN
            K=0
            DO L=2,LA
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &             +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
            ENDDO
            DO K=1,KS
              DO L=2,LA
                LN=LNC(L)
                STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     &             *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                   +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
                WSETA(L,K,NS)=CSEDSET(L,SED(L,K+1,NS),STRESS,ISEDVW)
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISEDVW.GE.3)THEN
            K=0
            L=857
            TAUBC=QQ(L,0)/CTURB2
            UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
            VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
            CURANG=ATAN2(VTMP,UTMP)
            TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &            +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
            TAUB2=MAX(TAUB2,0.)
            STRESSS(K)=SQRT(TAUB2)
            DO K=1,KS
              L=857
              LN=LNC(L)
              STRESSS(K)=AV(L,K)*SQRT( DZIGSD4(K) )
     &          *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     &                +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
            ENDDO
          ENDIF
C
            DO K=0,KS
            DO L=2,LA
              WSETA(L,K,NS)=SPB(L)*WSETA(L,K,NS)
            ENDDO
            ENDDO
C
C      WRITE(11,6111)TIME,(WSETA(857,K,1),K=0,KS)
C      WRITE(41,6111)TIME,(STRESSS(K),K=0,KS)
 6111 FORMAT(F10.2,10E12.4)
C
C----------------------------------------------------------------------C
C
C **  HORIZONTAL LOOPS
C
          K=KC
          DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
          ENDDO
C
          DO K=KS,2,-1
            DO L=2,LA
              WVEL=DELT*HPI(L)*DZIC(K)
              CLEFT=1.+WSETA(L,K-1,NS)*WVEL
              CRIGHT=MAX(SED(L,K,NS),0.)-SEDF(L,K,NS)*WVEL
              SED(L,K,NS)=CRIGHT/CLEFT
              SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
            ENDDO
          ENDDO
C
          DO L=2,LA
            PROBDEP=0.
            WESE=0.
C **  SET MAXIMUM EROSION RATE
c1104            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            WESEMX=DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
c1104            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUBSED(L).GT.TAURB(L,KBT(L)))THEN
C **  MASS EROSION
              WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
              WESE=MIN(WESE,WESEMX)
            ELSE
              TAUE=0.
              IF(TAUBSED(L).GT.TAURS(L,KBT(L)))THEN
C **  SURFACE EROSION
                WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
c1104                WESE=MIN(WESE,WESEMX)
                TAURTMP=TAURS(L,KBT(L))
	          IF(IWRSP(1).GE.2)TAURTMP=TAUR(1)
                IF(IWRSP(1).ge.99)TAURTMP=TAUNS(L,KBT(L))
			  TAUE=(TAUBSED(L)-TAURS(L,KBT(L)))/TAURTMP
                TAUE=MAX(TAUE,0.0)
	          IF(ISLTAUC(NS).EQ.1.AND.TAUBSED(L).LT.TAURTMP) TAUE=0.
                TMPSEDHID=1.0
                IF(ISTRAN(7).GE.1.AND.COSEDHID(1).NE.0.0) 
     &            TMPSEDHID=(FRACCOH(L,KBT(L)))**COSEDHID(1)
                IF(IWRSP(1).LT.99)THEN
                  WESE=TMPSEDHID*WESE*(TAUE**TEXP(NS))
                ELSE
                  WESE=TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
                ENDIF
c1104
                WESE=MIN(WESE,WESEMX)
c1104
              ELSE
C **  NO EROSION 
                WESE=0.0
              ENDIF
            ENDIF
C **  SET PROBABILITY OF DEPOSITION 
            IF(TAUBSED(L).LT.TAUD(NS)) 
     &         PROBDEP=(TAUD(NS)-TAUBSED(L))/TAUD(NS)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
cjmh216            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
            SEDBTMP=SEDB(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
C           IF(SEDBTMP.LT.0.0)THEN
C             SEDF(L,0,NS)=0.
C             SEDBTMP=SEDB1(L,KBT(L),NS)
C             SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
C           ENDIF
C           SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
C     &                      +S2TL*SEDB1(L,KBT(L),NS)
C           SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0)THEN
c1104              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
c1104
              SEDF(L,0,NS)=DELTI*SEDB(L,KBT(L),NS)
c1104
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)
     &                       +(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            ENDIF
cjmh216            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
cjmh216     &                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB1(L,KBT(L),NS)=SEDB(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)+DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)+DSEDGMM*
     &               ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     &                +VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
C
C            IF(L.EQ.857.AND.NS.EQ.1)THEN
C              WRITE(21,6111)TIME,TAUBC,QQWV2(L),TAUB,WESE,WSETMP,
C     &        SED(L,1,1)
C            ENDIF
C
          ENDDO
C
C         WRITE(31,6111)TIME,(SEDF(857,K,1),K=0,KS)
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.GT.3 
C
          IF(ISTOPT(6).EQ.1)THEN
C
            DO K=1,KS
              DO L=2,LA
                CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
                GRADSED=(SED(L,K+1,NS)-SED(L,K,NS))/(DZC(K+1)+DZC(K))
                SEDAVG=0.5*(SED(L,K+1,NS)+SED(L,K,NS)+1.E-16)
                WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*
     &                        GRADSED/SEDAVG
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
C
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
                TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)
     &                    +MAX(WSETA(L,KC-1,NS),0.)
            ENDDO
            DO K=2,KS
              DO L=2,LA
                TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)
     &         +MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)
              ENDDO
            ENDDO
C     TVAR1E=RIGHT HAND SIDE
            DO K=1,KC
              DO L=2,LA
                TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SED(L,K,NS)
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
     &                 TVAR3S(L)
              ENDDO
            ENDDO
            DO K=KS,1,-1
              DO L=2,LA
                TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
              ENDDO
            ENDDO
            DO K=1,KC
              DO L=2,LA
                SED(L,K,NS)=TVAR2N(L,K)
              ENDDO
            ENDDO
C
          ENDIF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC.GT.3
C
          DO L=2,LA
            SEDF(L,KS,NS)=DELTI*DZC(KC)*HP(L)*
     &                   (SED(L,KC,NS)-SEDS(L,KC,NS))
          ENDDO 
C 
          DO K=KS-1,1,-1
            DO L=2,LA
              SEDF(L,K,NS)=DELTI*DZC(K+1)*HP(L)*
     &                   (SED(L,K+1,NS)-SEDS(L,K+1,NS))+SEDF(L,K+1,NS)
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
C **  UNCOMMENT THE FOLLOWING FOR COHESIVE SEDIMENT DIAGNOSTICS
CDIAG           SEDMX=-1.E+12
CDIAG           SEDBMX=-1.E+12
CDIAG           SEDFMX=-1.E+12
CDIAG           SEDMN=1.E+12
CDIAG           SEDBMN=1.E+12
CDIAG           SEDFMN=1.E+12
CDIAG           DO K=1,KC
CDIAG            DO L=2,LA
CDIAG             IF(SED(L,K,NS).GT.SEDMX)THEN
CDIAG               LMX=L
CDIAG               KMX=K
CDIAG               SEDMX=SED(L,K,NS)
CDIAG             ENDIF
CDIAG             IF(SED(L,K,NS).LT.SEDMN)THEN
CDIAG               LMN=L
CDIAG               KMN=K
CDIAG               SEDMN=SED(L,K,NS)
CDIAG             ENDIF
CDIAG            ENDDO
CDIAG           ENDDO
CDIAG           DO L=2,LA
CDIAG            IF(SEDB(L,NS).GT.SEDBMX)THEN
CDIAG              LBMX=L
CDIAG              SEDBMX=SEDB(L,NS)
CDIAG            ENDIF
CDIAG            IF(SEDB(L,NS).LT.SEDBMN)THEN
CDIAG              LBMN=L
CDIAG              SEDBMN=SEDB(L,NS)
CDIAG            ENDIF
CDIAG            IF(SEDF(L,0,NS).GT.SEDFMX)THEN
CDIAG              LFMX=L
CDIAG              SEDFMX=SEDF(L,0,NS)
CDIAG            ENDIF
CDIAG            IF(SEDF(L,0,NS).LT.SEDFMN)THEN
CDIAG              LFMN=L
CDIAG              SEDFMN=SEDF(L,0,NS)
CDIAG            ENDIF
CDIAG           ENDDO
C
CDIAG      L=LMX
CDIAG      K=KMX
CDIAG      WRITE(1,101)N,NS,IL(L),JL(L),K,SED(L,K,NS),SEDS(L,K,NS)
CDIAG      L=LMN
CDIAG      K=KMN
CDIAG      WRITE(1,102)N,NS,IL(L),JL(L),K,SED(L,K,NS),SEDS(L,K,NS)
CDIAG      L=LBMX
CDIAG      WRITE(1,103)N,NS,IL(L),JL(L),SEDB(L,NS),SEDBS(L,NS),SEDF(L,0,NS)
CDIAG      L=LBMN
CDIAG      WRITE(1,104)N,NS,IL(L),JL(L),SEDB(L,NS),SEDBS(L,NS),SEDF(L,0,NS)
CDIAG      L=LFMX
CDIAG      WRITE(1,105)N,NS,IL(L),JL(L),SEDF(L,0,NS)
CDIAG      L=LFMN
CDIAG      WRITE(1,106)N,NS,IL(L),JL(L),SEDF(L,0,NS)
C
CDIAG  101 FORMAT(' N,NS,I,J,K,SEDMX,SEDSMX = ',5I5,4E13.4)       
CDIAG  102 FORMAT(' N,NS,I,J,K,SEDMN,SEDSMN = ',5I5,4E13.4)       
CDIAG  103 FORMAT(' N,NS,I,J,SEDBMX,SEDBSMX = ',4I5,4E13.4)       
CDIAG  104 FORMAT(' N,NS,I,J,SEDBMN,SEDBSMN = ',4I5,4E13.4)       
CDIAG  105 FORMAT(' N,NS,I,J,SEDFMX,SEDFSMX = ',4I5,4E13.4)       
CDIAG  106 FORMAT(' N,NS,I,J,SEDFMN,SEDFSMN = ',4I5,4E13.4)       
C
        OPEN(99,FILE='NEGSEDSND.OUT',POSITION='APPEND')
C
        DO NS=1,NSED
          DO K=1,KC
            DO L=2,LA
              IF(SED(L,K,NS).LT.0.)THEN
                WRITE(99,107)TIME,NS,IL(L),JL(L),K,SED(L,K,NS)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        DO NS=1,NSED
          DO L=2,LA
            IF(SEDB(L,KBT(L),NS).LT.0.)THEN
              WRITE(99,108)TIME,NS,IL(L),JL(L),KBT(L),
     &             SEDB(L,KBT(L),NS),SEDF(L,0,NS)
            ENDIF
          ENDDO
        ENDDO
C
        CLOSE(99)
C
C **  ACCUMULATE NET POSTIVE AND NEGATIVE COHESIVE SEDIMENT FLUXES
C
        DO NS=1,NSED
          DO L=2,LA
            SEDFDTAP(L,NS)=SEDFDTAP(L,NS)+DELT*MAX(SEDF(L,0,NS),0.0)
            SEDFDTAN(L,NS)=SEDFDTAN(L,NS)+DELT*MIN(SEDF(L,0,NS),0.0)
          ENDDO
        ENDDO
C
  107 FORMAT(' TIME,NS,I,J,K,NEGSED = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' TIME,NS,I,J,NEGSEDB,SEDF = ',F12.4,4I5,4E13.4)       
C
C**********************************************************************C
C
      RETURN
      END
