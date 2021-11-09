C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BEDINIT
C
C**********************************************************************C
C
C **  SUBROUTINE BEDINIT INITIALIZES SEDIMENT AND TOXIC VARIABLES
C **  IT SEDIMENT BED FOR HOT AND COLD START CONDITIONS
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 02/19/2002        john hamrick       02/19/2002       john hamrick
C  added additional diagnostic output
C 05/29/2002        john hamrick       05/29/2002       john hamrick
C  moved toxic initializations from ssedtox
C----------------------------------------------------------------------C
C
C**********************************************************************C
C  
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION NSP2(NTXM)
      DIMENSION FRACACT(LCM),FRACPAR(LCM),SEDBALL(LCM,KBM)
	DIMENSION FRACCOH(LCM,KBM),FRACNON(LCM,KBM)
	DIMENSION RADJCOHSEDS(LCM)
C
C**********************************************************************C
C
c **  check initial fractions
c
      DO K=1,KB
	DO L=2,LA
	  SEDBALL(L,K)=0.0
	ENDDO
	ENDDO
C
c
      DO NS=1,NSED
      DO K=1,KB
	DO L=2,LA
	  SEDBALL(L,K)=SEDBALL(L,K)+SEDBINIT(L,K,NS)
	ENDDO
	ENDDO
	ENDDO
C
c
      DO NS=1,NSND
      DO K=1,KB
	DO L=2,LA
	  SEDBALL(L,K)=SEDBALL(L,K)+SNDBINIT(L,K,NS)
	ENDDO
	ENDDO
	ENDDO
C
      OPEN(1,FILE='BEDFRACHK.OUT')
	DO L=2,LA
	  WRITE(1,1492)IL(L),JL(L),(SEDBALL(L,K),K=1,KB)
	ENDDO
	CLOSE(1)
C
 1492 FORMAT(2I6,12F16.6)
C
C**********************************************************************C
C
C **  DETERMINE START UP MODE
C
      IHOTSTRT=0
C
      IF(ISRESTI.NE.0)THEN
        IF(ISCI(6).NE.0.OR.ISCI(7).NE.0)THEN
          IHOTSTRT=1
        ENDIF
      ENDIF
C
C     WRITE(6,*)'INTER BED INITIALIZATION'
C
C**********************************************************************C
C
C **  HOT START INITIALIZATION
C
      IF(IHOTSTRT.NE.0)THEN
C
C----------------------------------------------------------------------C
C
C **  SET POROSITY
C
        DO K=1,KB
          DO L=2,LA
            PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
            PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  SET BULK DENSITY
C
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K)=0.
            SNDBT(L,K)=0.
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
            DO K=1,KB
              DO L=2,LA
                SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        DO K=1,KB
          DO L=2,LA
            IF(HBED(L,K).GT.0.)THEN
              BDENBED(L,K)=1000.*PORBED(L,K)
     &       +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
            ELSE
              BDENBED(L,K)=0.
            ENDIF
          ENDDO
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K)=0.
            SNDBT(L,K)=0.
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                SEDBT(L,K)=SEDBT(L,K)+SEDB1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
            DO K=1,KB
              DO L=2,LA
                SNDBT(L,K)=SNDBT(L,K)+SNDB1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        DO K=1,KB
          DO L=2,LA
            IF(HBED1(L,K).GT.0.)THEN
              BDENBED1(L,K)=1000.*PORBED1(L,K)
     &       +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED1(L,K)
            ELSE
              BDENBED1(L,K)=0.
            ENDIF
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  SET TOP BED LAYER
C
        DO L=2,LA
          KBT(L)=1
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            IF(HBED(L,K).GT.0.)KBT(L)=K
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
        IF(ISTRAN(6).GE.1)THEN
          IF(IWRSP(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=TAUR(1)
                WRSPS(L,K)=WRSPO(1)
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSPB(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=1.E6
                WRSPB(L,K)=0.0
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.1.AND.IWRSP(1).LT.99)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
                WRSPS(L,K)=CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1))
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.99.AND.ISHOUSATONIC.EQ.0)THEN
	      OPEN(1,FILE='SSCOHSEDPMAP.INP')
	      OPEN(2,FILE='SSCOHSEDPMAP.OUT')
	      DO NSKIP=1,6
	        READ(1,100)
	      ENDDO
	      READ(1,*)ISSTYPE
	      IF(ISSTYPE.EQ.0)THEN
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L)
	          RADJCOHSEDS(L)=1.
	        ENDDO
            ELSE
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
	        ENDDO
            ENDIF
            DO L=2,LA
	        LCORE=LSSCOHSED(L)
              TAUDS(L)=TAUDSS(LCORE)
            ENDDO
            DO L=2,LA
	        LCORE=LSSCOHSED(L)
              TAUDS(L)=TAUDSS(LCORE)
            ENDDO
            IF(NCOHSEDL.EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(1,LCORE)
                  TAUNS(L,K)=TAUNSS(1,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
                  TEXPS(L,K)=TEXPSS(1,LCORE)
                ENDDO
              ENDDO
            ELSE
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(K,LCORE)
                  TAUNS(L,K)=TAUNSS(K,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
                  TEXPS(L,K)=TEXPSS(K,LCORE)
                ENDDO
              ENDDO
		  ENDIF
            DO L=2,LA
	        K=KBT(L)
              WRITE(2,222)L,IL(L),JL(L),TAURS(L,K),TAUNS(L,K),
     &                                  WRSPS(L,K),TEXPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)
          ENDIF
          IF(IWRSPB(1).GE.1)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=CSEDTAUB(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
                WRSPB(L,K)=CSEDRESB(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
              ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K)=0.
            BEDDINIT(L,K)=0.
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
                VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
                VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C----------------------------------------------------------------------C
C
C**  INITIALIZE BED BOTTOM ELEVATION
C
        DO L=2,LA
          HBEDA(L)=0.
        ENDDO
C
        DO L=2,LA
          DO K=1,KBT(L)
            HBEDA(L)=HBEDA(L)+HBED(L,K)
          END DO
        ENDDO
C
        DO L=2,LA
          ZELBEDA(L)=BELV(L)-HBEDA(L)
        ENDDO
C
C----------------------------------------------------------------------C
C
C**  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
C
      DO K=1,KB
        DO L=2,LA
          SEDBT(L,K)=0.
          SNDBT(L,K)=0.
        ENDDO
      ENDDO
C
      IF(ISTRAN(6).GE.1)THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
      GOTO 1000
C
      ENDIF
C
C **  END HOT START INITIALIZATION
C
C**********************************************************************C
C
C **  COLD START INITIALIZATION: IBMECH=0
C
      IF(IBMECH.EQ.0)THEN
C
C----------------------------------------------------------------------C
C
C **  SET POROSITY AND VOID RATIO
C
        DO K=1,KB
          DO L=2,LA
            PORBED(L,K)=BEDPORC
            PORBED1(L,K)=BEDPORC
            VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
            VDRBED1(L,K)=PORBED1(L,K)/(1.-PORBED1(L,K))
            HBED(L,K)=0.
            HBED1(L,K)=0.
            KBT(L)=1
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  UNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C
        IF(ISEDINT.LE.1)THEN
C
C         WRITE(6,*)'INTER BED INITIALIZATION OPTION'
C
          DO K=1,KB
            DO L=2,LA
              SEDBT(L,K)=0.
              SNDBT(L,K)=0.
            ENDDO
          ENDDO
C
          IF(ISTRAN(6).GE.1)THEN
            DO NS=1,NSED
              DO K=1,KB
                DO L=2,LA
                  HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
                  SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          IF(ISTRAN(7).GE.1)THEN
            DO NX=1,NSND
              NS=NSED+NX
              DO K=1,KB
                DO L=2,LA
                  HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
                  SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C
          DO K=1,KB
            DO L=2,LA
              HBED(L,K)=(1.+VDRBED(L,K))*HBED(L,K)
              IF(HBED(L,K).GT.0.) KBT(L)=K
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              IF(HBED(L,K).GT.0.)THEN
                BDENBED(L,K)=1000.*PORBED(L,K)
     &         +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
              ELSE
                BDENBED(L,K)=0.
              ENDIF
            ENDDO
          ENDDO
C
          DO K=1,KB
            DO L=2,LA
              HBED1(L,K)=HBED(L,K)
              BDENBED1(L,K)=BDENBED(L,K)
            ENDDO
          ENDDO
C
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  AND INITIAL CONDITIONS ARE IN MASS PER UNIT AREA
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C
        IF(ISEDINT.GE.2)THEN
          IF(ISEDBINT.EQ.0)THEN
C
            DO K=1,KB
              DO L=2,LA
                SEDBT(L,K)=0.
                SNDBT(L,K)=0.
              ENDDO
            ENDDO
C
            IF(ISTRAN(6).GE.1)THEN
              DO NS=1,NSED
                DO K=1,KB
                  DO L=2,LA
                    HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
                    SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
            IF(ISTRAN(7).GE.1)THEN
              DO NX=1,NSND
                NS=NSED+NX
                DO K=1,KB
                  DO L=2,LA
                    HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
                    SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
            DO K=1,KB
              DO L=2,LA
                HBED(L,K)=(1.+VDRBED(L,K))*HBED(L,K)
                IF(HBED(L,K).GT.0.) KBT(L)=K
              ENDDO
            ENDDO
C
            DO K=1,KB
              DO L=2,LA
                IF(HBED(L,K).GT.0.)THEN
                  BDENBED(L,K)=1000.*PORBED(L,K)
     &           +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
                ELSE
                  BDENBED(L,K)=0.
                ENDIF
              ENDDO
            ENDDO
C
            DO K=1,KB
              DO L=2,LA
                HBED1(L,K)=HBED(L,K)
                BDENBED1(L,K)=BDENBED(L,K)
              ENDDO
            ENDDO
C
          ENDIF
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  AND INITIAL CONDITIONS ARE IN MASS FRACTION
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C **  THIS OPTION REQUIRES INITIAL LAYER THICKNESSES
C
        IF(ISEDINT.GE.2)THEN
          IF(ISEDBINT.EQ.1)THEN
C
            IF(IBEDLAYU.EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
                  BEDLINIT(L,K)=0.1*BEDLINIT(L,K)
                ENDDO
              ENDDO
            ENDIF
C
            DO K=1,KB
              DO L=2,LA
                HBED(L,K)=BEDLINIT(L,K)
                HBED1(L,K)=BEDLINIT(L,K)
                IF(HBED(L,K).GT.0.)KBT(L)=K
              ENDDO
            ENDDO
C
            DO K=1,KB
              DO L=2,LA
                BDENBED(L,K)=0.
              ENDDO
            ENDDO
C
            IF(ISTRAN(6).GE.1)THEN
              DO NS=1,NSED
                DO K=1,KB
                  DO L=2,LA
                    BDENBED(L,K)=BDENBED(L,K)
     &             +1000.*SSG(NS)*SEDB(L,K,NS)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
            IF(ISTRAN(7).GE.1)THEN
              DO NX=1,NSND
                NS=NSED+NX
                DO K=1,KB
                  DO L=2,LA
                    BDENBED(L,K)=BDENBED(L,K)
     &             +1000.*SSG(NS)*SNDB(L,K,NX)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
            DO K=1,KB
              DO L=2,LA
                BDENBED(L,K)=1000.*PORBED(L,K)
     &         +(1.-PORBED(L,K))*BDENBED(L,K)
                BDENBED1(L,K)=BDENBED(L,K)
              ENDDO
            ENDDO
C
            DO K=1,KB
              DO L=2,LA
                SEDBT(L,K)=1000.*HBED(L,K)*(BDENBED(L,K)
     &         -1000.*PORBED(L,K))
              ENDDO
            ENDDO
C
            IF(ISTRAN(6).GE.1)THEN
              DO NS=1,NSED
                DO K=1,KB
                  DO L=2,LA
                    SEDB(L,K,NS)=SEDB(L,K,NS)*SEDBT(L,K)
                    SEDB1(L,K,NS)=SEDB(L,K,NS)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
            IF(ISTRAN(7).GE.1)THEN
              DO NX=1,NSND
                DO K=1,KB
                  DO L=2,LA
                    SNDB(L,K,NX)=SNDB(L,K,NX)*SEDBT(L,K)
                    SNDB1(L,K,NX)=SNDB(L,K,NX)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
        IF(ISTRAN(6).GE.1)THEN
          IF(IWRSP(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=TAUR(1)
                WRSPS(L,K)=WRSPO(1)
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSPB(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=1.E6
                WRSPB(L,K)=0.0
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.1.AND.IWRSP(1).LT.99)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
                WRSPS(L,K)=CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1))
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.99.AND.ISHOUSATONIC.EQ.0)THEN
	      OPEN(1,FILE='SSCOHSEDPMAP.INP')
	      OPEN(2,FILE='SSCOHSEDPMAP.OUT')
	      DO NSKIP=1,6
	        READ(1,100)
	      ENDDO
	      READ(1,*)ISSTYPE
	      IF(ISSTYPE.EQ.0)THEN
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L)
	          RADJCOHSEDS(L)=1.
	        ENDDO
            ELSE
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
	        ENDDO
            ENDIF
            DO L=2,LA
	        LCORE=LSSCOHSED(L)
              TAUDS(L)=TAUDSS(LCORE)
            ENDDO
            IF(NCOHSEDL.EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(1,LCORE)
                  TAUNS(L,K)=TAUNSS(1,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
                  TEXPS(L,K)=TEXPSS(1,LCORE)
                ENDDO
              ENDDO
            ELSE
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(K,LCORE)
                  TAUNS(L,K)=TAUNSS(K,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
                  TEXPS(L,K)=TEXPSS(K,LCORE)
                ENDDO
              ENDDO
            ENDIF
            DO L=2,LA
	        K=KBT(L)
              WRITE(2,222)L,IL(L),JL(L),TAURS(L,K),TAUNS(L,K),
     &                                  WRSPS(L,K),TEXPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)
          ENDIF
          IF(IWRSPB(1).GE.1)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=CSEDTAUB(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
                WRSPB(L,K)=CSEDRESB(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
              ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K)=0.
            BEDDINIT(L,K)=0.
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
                VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
                VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C **  END COLD START INITIALIZATION: IBMECH=0
C
C**********************************************************************C
C
C **  COLD START INITIALIZATION: IBMECH.GE.1
C
      IF(IBMECH.GE.1)THEN
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND INITIALIZE BED LAYER THICKNESS AND DEFINE
C **  INITIAL TOP LAYER
C
        IF(IBEDLAYU.EQ.1)TMPCVT=0.001
        IF(IBEDLAYU.EQ.2)TMPCVT=0.01
        IF(IBEDLAYU.EQ.3)TMPCVT=1.0

        IF(IBEDLAYU.GE.1)THEN
          DO K=1,KB
            DO L=2,LA
              BEDLINIT(L,K)=TMPCVT*BEDLINIT(L,K)
            ENDDO
          ENDDO
        ENDIF
C
        DO L=2,LA
          KBT(L)=1
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            HBED(L,K)=BEDLINIT(L,K)
            HBED1(L,K)=BEDLINIT(L,K)
            IF(HBED(L,K).GT.0.) KBT(L)=K
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND INITIALIZE BED BULK DENSITY
C
C**   IBEDBDNU=0 BEDBINIT IS NOT BULK DENSITY
C**   IBEDBDNU=1 BEDBINIT BULK DENSITY IN KG/M**3
C**   IBEDBDNU=3 BEDBINIT BULK DENSITY IN GM/CM**3
C
        IF(IBEDBDNU.GE.1)THEN
C
          IF(IBEDBDNU.EQ.2)THEN
            DO K=1,KB
              DO L=2,LA
                BEDBINIT(L,K)=1000.*BEDBINIT(L,K)
              ENDDO
            ENDDO
          ENDIF
C
          DO K=1,KB
            DO L=2,LA
              BDENBED(L,K)=BEDBINIT(L,K)
              BDENBED1(L,K)=BEDBINIT(L,K)
            ENDDO
          ENDDO
C
          ENDIF
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND DRY DENSITY OF BED 
C **  IBEDDDNU=0,1 ACTUAL DRY DENSITY, =2  POROSITY, =3 VOID RATIO
C
        IF(IBEDDDNU.EQ.1)THEN
          DO K=1,KB
            DO L=2,LA
              BEDDINIT(L,K)=1000.*BEDDINIT(L,K)
            ENDDO
          ENDDO
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  CALCULATE POROSITY AND VOID RATIO
C
        IF(IBEDDDNU.LE.1)THEN
          DO K=1,KB
            DO L=2,LA
              PORBED(L,K)=0.001*(BEDBINIT(L,K)-BEDDINIT(L,K))
              VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
            ENDDO
          ENDDO
        ENDIF
C
        IF(IBEDDDNU.EQ.2)THEN
          DO K=1,KB
            DO L=2,LA
              PORBED(L,K)=BEDDINIT(L,K)
              VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
            ENDDO
          ENDDO
        ENDIF
C
        IF(IBEDDDNU.EQ.3)THEN
          DO K=1,KB
            DO L=2,LA
              VDRBED(L,K)=BEDDINIT(L,K)
              PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
            ENDDO
          ENDDO
        ENDIF
C
        DO K=1,KB
          DO L=2,LA
            VDRBED1(L,K)=VDRBED(L,K)
            PORBED1(L,K)=PORBED(L,K)
          ENDDO
        ENDDO
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE BED SEDIMENT FOR MASS FRACTION INPUT BY CACLUALTING
C **  AND STORING TOTAL MASS OF SED/AREA IN BEDDINIT(L,K)
C
        DO K=1,KB
          DO L=2,LA
            BEDDINIT(L,K)=HBED(L,K)*(BDENBED(L,K)-1000.*PORBED(L,K))
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            IF(ISEDBU(NS).EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
                  SEDB(L,K,NS)=1000.*SEDBINIT(L,K,NS)*BEDDINIT(L,K)
                  SEDB1(L,K,NS)=SEDB(L,K,NS)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
            IF(ISNDBU(NS).EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
                  SNDB(L,K,NS)=1000.*SNDBINIT(L,K,NS)*BEDDINIT(L,K)
                  SNDB1(L,K,NS)=SNDB(L,K,NS)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
        IF(ISTRAN(6).GE.1)THEN
          IF(IWRSP(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=TAUR(1)
                WRSPS(L,K)=WRSPO(1)
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSPB(1).EQ.0)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=1.E6
                WRSPB(L,K)=0.0
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.1.AND.IWRSP(1).LT.99)THEN
            DO K=1,KB
              DO L=2,LA
                TAURS(L,K)=CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
                WRSPS(L,K)=CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSP(1))
              ENDDO
            ENDDO
          ENDIF
          IF(IWRSP(1).GE.99.AND.ISHOUSATONIC.EQ.0)THEN
	      OPEN(1,FILE='SSCOHSEDPMAP.INP')
	      OPEN(2,FILE='SSCOHSEDPMAP.OUT')
	      DO NSKIP=1,6
	        READ(1,100)
	      ENDDO
	      READ(1,*)ISSTYPE
	      IF(ISSTYPE.EQ.0)THEN
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L)
	          RADJCOHSEDS(L)=1.
	        ENDDO
            ELSE
	        DO L=2,LA
	          READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
	        ENDDO
            ENDIF
            DO L=2,LA
	        LCORE=LSSCOHSED(L)
              TAUDS(L)=TAUDSS(LCORE)
            ENDDO
            IF(NCOHSEDL.EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(1,LCORE)
                  TAUNS(L,K)=TAUNSS(1,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
                  TEXPS(L,K)=TEXPSS(1,LCORE)
                ENDDO
              ENDDO
            ELSE
              DO K=1,KB
                DO L=2,LA
	              LCORE=LSSCOHSED(L)
                  TAURS(L,K)=TAURSS(K,LCORE)
                  TAUNS(L,K)=TAUNSS(K,LCORE)
                  WRSPS(L,K)=RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
                  TEXPS(L,K)=TEXPSS(K,LCORE)
                ENDDO
              ENDDO
            ENDIF
            DO L=2,LA
	        K=KBT(L)
              WRITE(2,222)L,IL(L),JL(L),TAURS(L,K),TAUNS(L,K),
     &                                  WRSPS(L,K),TEXPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)
          ENDIF
          IF(IWRSPB(1).GE.1)THEN
            DO K=1,KB
              DO L=2,LA
                TAURB(L,K)=CSEDTAUB(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
                WRSPB(L,K)=CSEDRESB(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                            VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
              ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K)=0.
            BEDDINIT(L,K)=0.
          ENDDO
        ENDDO
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
                VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
                VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
                BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                IF(BEDLINIT(L,K).GT.0.0)THEN
                  VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                ELSE
                  VFRBED(L,K,NS)=0.0
                ENDIF
                IF(BEDDINIT(L,K).GT.0.0)THEN
                  VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
                ELSE
                  VFRBED1(L,K,NS)=0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GE.1)THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                IF(BEDLINIT(L,K).GT.0.0)THEN
                  VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
                ELSE
                  VFRBED(L,K,NS)=0.0
                ENDIF
                IF(BEDDINIT(L,K).GT.0.0)THEN
                  VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
                ELSE
                  VFRBED1(L,K,NS)=0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C **  END COLD START INITIALIZATION: IBMECH.GE.1
C
C**********************************************************************C
C
C**  INITIALIZE BED BOTTOM ELEVATION
C
        DO L=2,LA
          HBEDA(L)=0.
        ENDDO
C
        DO L=2,LA
          DO K=1,KBT(L)
            HBEDA(L)=HBEDA(L)+HBED(L,K)
          END DO
        ENDDO
C
        DO L=2,LA
          ZELBEDA(L)=BELV(L)-HBEDA(L)
        ENDDO
C
C**********************************************************************C
C
C**  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
C
      DO K=1,KB
        DO L=2,LA
          SEDBT(L,K)=0.
          SNDBT(L,K)=0.
        ENDDO
      ENDDO
C
      IF(ISTRAN(6).GE.1)THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN
C **  BED AND REINITILIZE IF NECESSARY
C
C      IF(N.EQ.1.AND.ISTRAN(5).GE.1)THEN
      IF(ISTRAN(5).GE.1)THEN
        IF(ISRESTI.EQ.0.OR.ISCI(5).EQ.0)THEN
c
C
C **  CALCULATE TOTAL SEDIMENT IN THE BED
C
      DO K=1,KB
        DO L=1,LC
          SEDBT(L,K)=0.
          SNDBT(L,K)=0.
          SEDBALL(L,K)=0. 
        ENDDO
      ENDDO
C
      IF(ISTRAN(6).GE.1)THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      DO K=1,KB
        DO L=1,LC
          SEDBALL(L,K)=SEDBT(L,K)+SNDBT(L,K)
        ENDDO
      ENDDO
C  
C
C **  CALCULATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
          DO NT=1,NTOX
            NSP2(NT)=NSED+NSND
            IF(ISTOC(NT).EQ.2) NSP2(NT)=NSP2(NT)+1
            IF(ISTOC(NT).EQ.3) NSP2(NT)=NSP2(NT)+2
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
                   TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*TOXPARB(NS,NT)
                 ENDDO
               ENDDO
             ENDDO
            ENDIF
            IF(ISTRAN(7).GE.1)THEN
              DO NX=1,NSND
                NS=NX+NSED
                DO K=1,KB
                  DO L=2,LA
                    TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*TOXPARB(NS,NT)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            IF(ISTOC(NT).EQ.2)THEN
              NS=1+NSED+NSND
              DO K=1,KB
                DO L=2,LA
                  TOXPFB(L,K,NS,NT)=STDOCB(L,K)*TOXPARBC(1,NT)
                ENDDO
              ENDDO
            ENDIF
            IF(ISTOC(NT).EQ.3)THEN
              NS=2+NSED+NSND
              DO K=1,KB
                DO L=2,LA
                  TOXPFB(L,K,NS,NT)=STPOCB(L,K)*TOXPARBC(2,NT)
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
            DO K=1,KB
              DO L=2,LA
                IF(SEDBALL(L,K).GT.0.0)THEN
                  TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     &                     /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
                ELSE
                  TOXPFTB(L,K,NT)=1.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
C
C **  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
C **  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
C **  TO TOXB UNITS OF OF MG/M**2
C
          DO NT=1,NTOX
            IF(ITXBDUT(NT).EQ.0)THEN
              DO K=1,KB
                DO L=2,LA
                  TOXB(L,K,NT)=HBED(L,K)*TOXB(L,K,NT)
                  TOXB1(L,K,NT)=TOXB(L,K,NT)
                ENDDO
              ENDDO
            ENDIF
            IF(ITXBDUT(NT).EQ.1)THEN
              DO K=1,KB
                DO L=2,LA
                  TOXB(L,K,NT)=0.001*TOXB(L,K,NT)*(SEDBT(L,K)
     &               +SNDBT(L,K))/TOXPFTB(L,K,NT) 
                  TOXB1(L,K,NT)=TOXB(L,K,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
C
C ** DIAGNOSTICS OF INITIALIZATION
C
          IF(ISDTXBUG.EQ.1)THEN
            OPEN(2,FILE='TOXBED.DIA')
            CLOSE(2,STATUS='DELETE')
            OPEN(2,FILE='TOXBED.DIA')
            DO L=2,LA
C             TMP1=-999.
C             TMP2=-999.
C             IF(HBED(L).GT.0.)TMP1=TOXB(L,1)/HBED(L)
C             IF(HBED(L).GT.0.)TMP2=TOXB(L,2)/HBED(L)
C             WRITE(2,2222)IL(L),JL(L),HBED(L),TOXB(L,1),TOXB(L,2),TMP1,TMP2
              TMP1=TOXB(L,1,1)/(HBED(L,1)+1.E-12)
              WRITE(2,2222)IL(L),JL(L),TOXPFTB(L,1,1),TOXB(L,1,1),
     &              TMP1,TOX(L,1,1)
            ENDDO
            CLOSE(2)
	    ENDIF
C
        ENDIF
      ENDIF
C
 2222 FORMAT(2I5,7E13.4)
C
C**********************************************************************C
C
C **  INITIALIZE FRACTION OF PARTICULATE ORGANIC CARBON IN BED
C
      IVAL=0
	DO NT=1,NTOX
        IF(ISTOC(NT).GE.2)IVAL=1
	ENDDO
C
      IF(IVAL.EQ.1.AND.ISTPOCB.EQ.4)THEN
	  CALL SETFPOCB(0)
	ENDIF
C
C**********************************************************************C
C
C **  CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS
C
      DO K=1,KB
      DO L=2,LA
	  IF(K.LE.KBT(L))THEN
	    FVOLSSD=1./(1.+VDRBED(L,K))
          FVOLSED=0.0
	    DO NS=1,NSED
            FVOLSED=FVOLSED+VFRBED(L,K,NS)
	    ENDDO
	    FVOLSND=0.0
          DO NX=1,NSND
            NS=NSED+NX
	      FVOLSND=FVOLSND+VFRBED(L,K,NS)
	    ENDDO
          FVOLSED=FVOLSSD*FVOLSED
          FVOLSND=FVOLSSD*FVOLSND
	    VDRBEDSND(L,K)=SNDVDRD
          VDRBEDSED(L,K)=0.0
          IF(FVOLSED.GT.0.0)
     &    VDRBEDSED(L,K)=((FVOLSED+FVOLSND)*VDRBED(L,K)-FVOLSND*SNDVDRD)
     &                   /FVOLSED
        ELSE
	    VDRBEDSND(L,K)=0.0
	    VDRBEDSED(L,K)=0.0
	  ENDIF
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADD ACTIVE ARMORING LAYER IF NO PRESENT IN INITIAL OR RESTART 
C     CONDITIONS
C
      IF(ISNDAL.EQ.2.AND.IALSTUP.GT.0)THEN
C
      DO L=2,LA
	  KTOPTP=KBT(L)
	  KTOPP1=KBT(L)+1
	  FRACACT(L)=HBEDAL/HBED(L,KTOPTP)
	  FRACPAR(L)=(HBED(L,KTOPTP)-HBEDAL)/HBED(L,KTOPTP)
	  HBED(L,KTOPP1)=FRACACT(L)*HBED(L,KTOPTP)
	  HBED(L,KTOPTP)=FRACPAR(L)*HBED(L,KTOPTP)
        PORBED(L,KTOPP1)=PORBED(L,KTOPTP)
        PORBED1(L,KTOPP1)=PORBED1(L,KTOPTP)
        VDRBED(L,KTOPP1)=VDRBED(L,KTOPTP)
        VDRBED1(L,KTOPP1)=VDRBED1(L,KTOPTP)
        BDENBED(L,KTOPP1)=BDENBED(L,KTOPTP)
        BDENBED1(L,KTOPP1)=BDENBED1(L,KTOPTP)
        SEDBT(L,KTOPP1)=FRACACT(L)*SEDBT(L,KTOPTP)
        SEDBT(L,KTOPTP)=FRACPAR(L)*SEDBT(L,KTOPTP)
        SNDBT(L,KTOPP1)=FRACACT(L)*SNDBT(L,KTOPTP)
        SNDBT(L,KTOPTP)=FRACPAR(L)*SNDBT(L,KTOPTP)
        STDOCB(L,KTOPP1)=STDOCB(L,KTOPTP)
        STPOCB(L,KTOPP1)=STPOCB(L,KTOPTP)
	ENDDO
C
      DO NS=1,NSED
      DO L=2,LA
	  KTOPTP=KBT(L)
	  KTOPP1=KBT(L)+1
	  SEDB(L,KTOPP1,NS)=FRACACT(L)*SEDB(L,KTOPTP,NS)
	  SEDB1(L,KTOPP1,NS)=FRACACT(L)*SEDB1(L,KTOPTP,NS)
	  SEDB(L,KTOPTP,NS)=FRACPAR(L)*SEDB(L,KTOPTP,NS)
	  SEDB1(L,KTOPTP,NS)=FRACPAR(L)*SEDB1(L,KTOPTP,NS)
        STFPOCB(L,KTOPP1,NS)=STFPOCB(L,KTOPTP,NS)
      ENDDO
	ENDDO      
C
      DO NS=1,NSND
	NX=NSED+NS
      DO L=2,LA
	  KTOPTP=KBT(L)
	  KTOPP1=KBT(L)+1
	  SNDB(L,KTOPP1,NS)=FRACACT(L)*SNDB(L,KTOPTP,NS)
	  SNDB1(L,KTOPP1,NS)=FRACACT(L)*SNDB1(L,KTOPTP,NS)
	  SNDB(L,KTOPTP,NS)=FRACPAR(L)*SNDB(L,KTOPTP,NS)
	  SNDB1(L,KTOPTP,NS)=FRACPAR(L)*SNDB1(L,KTOPTP,NS)
        STFPOCB(L,KTOPP1,NX)=STFPOCB(L,KTOPTP,NX)
	ENDDO
	ENDDO      
C
      DO NT=1,NTOX
      DO L=2,LA
	  KTOPTP=KBT(L)
	  KTOPP1=KBT(L)+1
	  TOXB(L,KTOPP1,NT)=FRACACT(L)*TOXB(L,KTOPTP,NT)
	  TOXB1(L,KTOPP1,NT)=FRACACT(L)*TOXB1(L,KTOPTP,NT)
	  TOXB(L,KTOPTP,NT)=FRACPAR(L)*TOXB(L,KTOPTP,NT)
	  TOXB1(L,KTOPTP,NT)=FRACPAR(L)*TOXB1(L,KTOPTP,NT)
        TOXPFTB(L,KTOPP1,NT)=TOXPFTB(L,KTOPTP,NT)
      ENDDO
	ENDDO      
C
      DO NT=1,NTOX
	DO NS=1,NSED+NSND+2
      DO L=2,LA
	  KTOPTP=KBT(L)
	  KTOPP1=KBT(L)+1
	  TOXPFB(L,KTOPP1,NS,NT)=TOXPFB(L,KTOPTP,NS,NT)
      ENDDO
	ENDDO      
	ENDDO 
C
      DO L=2,LA
	  KBT(L)=KBT(L)+1
      ENDDO
     
      ENDIF
C
C
C**********************************************************************C
C
C **  ADJUST POROSITY AND VOID RATIO FOR IBMECH.EQ.99
C
      IF(IBMECH.EQ.99)THEN
C
	  DO K=1,KB
        DO L=2,LA
	    FRACCOH(L,K)=0.0
          FRACNON(L,K)=0.0
        ENDDO
        ENDDO
C
        DO NS=1,NSED
	  DO K=1,KB
        DO L=2,LA
	    IF(K.LE.KBT(L))THEN
	      FRACCOH(L,K)=FRACCOH(L,K)+VFRBED(L,K,NS)
	      FRACNON(L,K)=FRACNON(L,K)+VFRBED1(L,K,NS)
          ENDIF
        ENDDO
        ENDDO
        ENDDO

C
        DO K=1,KB
          DO L=2,LA
            IF(K.LE.KBT(L))THEN
              PORBED(L,K)=BMECH1*(FRACCOH(L,K)**BMECH2)+BMECH3
              PORBED1(L,K)=BMECH1*(FRACNON(L,K)**BMECH2)+BMECH3
              VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
              VDRBED1(L,K)=PORBED1(L,K)/(1.-PORBED1(L,K))
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C     WRITE(6,*)'COMPLETED BED INITIALIZATION'
C
C**  WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION
C
COLD      OPEN(1,FILE='BEDINIT.DIA')
COLD      CLOSE(1,STATUS='DELETE')
COLD      OPEN(1,FILE='BEDINIT.DIA')
COLD      WRITE(1,102)  
COLD      DO L=2,LA
COLDC       WRITE(1,101)IL(L),JL(L),SEDB(L,1,1),SNDB(L,1,1),HBED(L,1),
COLD        WRITE(1,101)IL(L),JL(L),SEDBT(L,1),SNDBT(L,1),HBED(L,1),
COLD     &             PORBED(L,1),VDRBED(L,1),BDENBED(L,1),TAURS(L,1)
COLD      ENDDO
COLD      CLOSE(1)
C
 1000 CONTINUE
C
      OPEN(1,FILE='BEDINIT.SED')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.SED')
      WRITE(1,111)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SEDB(L,K,1),K=1,KB)
        IF(NSED.GT.1) THEN
          DO NX=2,NSED
           WRITE(1,102)(SEDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.SAN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.SAN')
      WRITE(1,112)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SNDB(L,K,1),K=1,KB)
        IF(NSND.GT.1)THEN
          DO NX=2,NSND
           WRITE(1,102)(SNDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.VDR')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.VDR')
      WRITE(1,113)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(VDRBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)

C
      DO K=1,KB
        DO L=2,LA
          IF(K.LE.KBT(L))THEN
              VDRBED2(L,K)=VDRBED(L,K)
              SDENAVG(L,K)=(BDENBED(L,K)-1000.0*PORBED(L,K))
     &                    /(1.0-PORBED(L,K))
          ELSE
              VDRBED2(L,K)=VDRBED(L,KBT(L))
              SDENAVG(L,K)=(BDENBED(L,KBT(L))-1000.0*PORBED(L,KBT(L)))
     &                    /(1.0-PORBED(L,KBT(L)))
          ENDIF
        ENDDO
      ENDDO
C
      OPEN(1,FILE='BEDINIT.POR')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.POR')
      WRITE(1,114)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(PORBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.ZHB')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.ZHB')
      WRITE(1,115)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),(HBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.BDN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.BDN')
      WRITE(1,116)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(BDENBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.ELV')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.ELV')
      WRITE(1,117)  
      DO L=2,LA
        SURF=HP(L)+BELV(L)
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),BELV(L),HP(L),SURF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.TOX')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.TOX')
	DO NT=1,NTOX
      WRITE(1,118)NT  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(TOXB(L,K,NT),K=1,KB)
      ENDDO
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINIT.VRS')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINIT.VRS')
      DO L=2,LA
        WRITE(1,191)IL(L),JL(L),(VDRBED(L,K),K=1,KB)
        WRITE(1,192) (VDRBEDSED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDINITC.VVF')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDINITC.VVF')
      OPEN(2,FILE='BEDINITF.VVF')
      CLOSE(2,STATUS='DELETE')
      OPEN(2,FILE='BEDINITF.VVF')
      DO L=2,LA
	  K=KBT(L)
	  IF(HMP(L).GT.0.05)THEN
          WRITE(1,191)IL(L),JL(L),VFRBED(L,K,1),PORBED(L,K),VDRBED(L,K),
     &        VDRBEDSED(L,K)
	  ELSE
          WRITE(2,191)IL(L),JL(L),VFRBED(L,K,1),PORBED(L,K),VDRBED(L,K),
     &        VDRBEDSED(L,K)
	  ENDIF
      ENDDO
      CLOSE(1)
      CLOSE(2)
C
      OPEN(1,FILE='WATINIT.SED')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WATINIT.SED')
      WRITE(1,128)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SED(L,K,1),K=1,KC)
        IF(NSED.GT.1) THEN
          DO NX=2,NSED
           WRITE(1,102)(SED(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='WATRINIT.SAN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WATRINIT.SAN')
      WRITE(1,129)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SND(L,K,1),K=1,KC)
        IF(NSND.GT.1)THEN
          DO NX=2,NSND
           WRITE(1,102)(SND(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
C
  100 FORMAT(1X)
  222 FORMAT(3I6,5E13.4)
  191 FORMAT(2I5,18F10.3)
  192 FORMAT(10X,18F10.3)
  101 FORMAT(2I5,18E13.5)
  102 FORMAT(10X,18E13.5)
  111 FORMAT('   IL   JL    SEDBT(K=1,KB)')
  112 FORMAT('   IL   JL    SNDBT(K=1,KB)')
  113 FORMAT('   IL   JL    VRDBED(K=1,KB)')
  114 FORMAT('   IL   JL    PORBED(K=1,KB)')
  115 FORMAT('   IL   JL    ZBEDB        HBEDT        HBED(K=1,KB)')
  116 FORMAT('   IL   JL    BDENBED(K=1,KB)')
  117 FORMAT('   IL   JL    ZBEDB        HBEDT        BELV',
     &       '        HWCOL        SELV')
  118 FORMAT('   IL   JL    TOXB(K=1,KB,NT)  NT = ',I5)
  128 FORMAT('   IL   JL    SEDW(K=1,KC)')
  129 FORMAT('   IL   JL    SNDW(K=1,KC)')
C
C**********************************************************************C
C
      RETURN
      END
