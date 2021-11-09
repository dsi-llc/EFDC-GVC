C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SSEDTOX(ISTLX,IS2TLX,CORDTX)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 11/02/2001        John Hamrick       11/02/2001       John Hamrick
C  changed add and remove bed layer algorithm - see sections under
C  ibmech=0,1, and ge 2
C 11/06/2001        John Hamrick       11/06/2001        John Hamrick
C  changed add and remove bed layer algorithm - see sections under
C   ibmech=0,1, and ge 2
C  modified noncohesive resuspension formulation to be consistent
C   with add and remove bed layer modification
C  added noncohesive bedload-suspended load distribution factor
C   calculated by function fsedmode.f
C  added hiding factor correction to bagnold form of bed load transport
C   calculation
C 11/14/2001        John Hamrick       11/14/2001       John Hamrick
C   changed dimensions of QSBDLDX and QSBDLDY
C 11/15/2001        john hamrick       11/15/2001       john hamrick
C  modified calls functions FHYDCN, FDSTRSE,and FSTRSE to added use
C  standard exponential form constitutive relationships 
C 11/16/2001        john hamrick       11/16/2001       john hamrick
C  fixed errors in finite strain consolidation (IBMECH.ge.2) 
C 11/20/2001        john hamrick       11/2/2001       john hamrick
C  added bed load outflow/recirculation boundary condition options
C  and by pass array RBPSBL
C 12/03/2001        john hamrick       12/03/2001       john hamrick
C  increased dimensions of bedload flux SNDFBL
C  corrected calculation of SNDFBL
C  modified function FSBDLD call 
C  added exposure and hiding functions PEXP and PHID
C 12/13/2001        john hamrick       12/13/2001       john hamrick
C  added additional logic to prevent divide by zeros
C 12/15/2001        john hamrick       12/15/2001       john hamrick
C  corrected spelling error SNDDMX (was SNDDMAX)
C  corrected spelling error SEDVRDT (was SEDVRT)
C 01/21/2002        john hamrick       01/21/2002       john hamrick
C  added belv1, fixed errors in general bed load function, and 
C  eliminated errorneous return based on toxic parameter value
C 01/31/2002        john hamrick       01/31/2002       john hamrick
C  modified toxic-organic carbon partitioning options
C 02/19/2002        john hamrick       02/19/2002       john hamrick
C  fix some inconsistencies involving sedb1,sndb1,toxb1. added 
C  adjustment to water column concentrations when morphological mode
C  is activated.  fixed error in pore water advection and diffusion 
C  solution
C 03/05/2002        john hamrick       03/05/2002       john hamrick
C  added by pass of bed load transport for dry cells
C 03/06/2002        john hamrick       03/06/2002       john hamrick
C  added adjustment to toxic pore water advection and diffusion to
C  guarantee mass conservation, including variables derrb,toxbbalo,
C  toxbbaln,toxwbalo,toxwbaln
C 05/22/2002        John Hamrick       05/22/2002       John Hamrick
C  fixed bed load transport of sorbed contaminant. and added bed load 
C  toxic flux TOXFBL(L,NT), bed load toxic flux on out flow boundry
C  TOXBLB(NT) and sed-tox debug flag ISDTXBUG
C 05/29/2002        john hamrick       05/29/2002       john hamrick
C  moved toxic initializations to bedinit.for
C 06/05/2002        John Hamrick       06/06/2002       John Hamrick
C  added roundoff control to sediment-water column exchange of
C  toxics.  moved local array SNDFBL to global common
C 06/06/2002        John Hamrick       06/06/2002       John Hamrick
C  made local arrays TAUB(L),USTAR(L),UCELLCTR(L),VCELLCTR(L) of former
C  scaler variables of same names.  calculated each once at start
C  of sediment transport.
C 06/06/2002        John Hamrick       06/06/2002       John Hamrick
C  rewrote sediment bed flow and recirculation boundary condition
C  implementation in bed load transport section added global arrays
C  QSBDLOT and QSBDLIN
C 06/06/2002        John Hamrick       06/06/2002       John Hamrick
C  added QMORPH(L), the equivalent water column volume source associate
C  with change in bed elevation, for use in mass balance
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE SSEDTOX CALCULATES SETTLING AND WATER COLUMN-BED 
C **  EXCHANGE OF SEDIMENT AND SORBED TOXIC CONTAMINANTS
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
      DIMENSION QCELLAD1SQ(LCM),QCELLAD1ZZ(LCM),SESNFA(LCM),
     &          SESNFP(LCM),DELBED(LCM),TAUBSEDS(LCM),TAUBSNDS(LCM),
     &          ISSBCP(LCM)
C
C**********************************************************************C
C
      ISTL=ISTLX
	IS2TL=IS2TLX
	CORT=CORDTX
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF(ISTL.NE.3)THEN
        DELT=DT
        S3TL=0.0
        S2TL=1.0
        ISUD=0
      ENDIF
      IF(IS2TL.EQ.1)THEN
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
        ELSE
          DELT=DTDYN
        END IF
        S3TL=1.0
        S2TL=0.0
        ISUD=1
      ENDIF
C
      IF(ISEDDT.GT.1) DELT=DTSED
C
      IF(N.EQ.1)THEN
	  WRITE(8,*)'S3TL,S2TL,ISUD'
	  WRITE(8,*)S3TL,S2TL,ISUD
	END IF 
C     
      DELTI=1./DELT
C
      SEDMDGM=SQRT(SEDMDMX*SEDMDMN)
CBEGMOD
      BEDEX=1.
CENDMOD
CJH      BEDEX=0.
      NVAL=MOD(N,2)
CJH      IF(NVAL.EQ.0)THEN
CJH        IF(ISTL.EQ.3) BEDEX=1.
CJH      ENDIF
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
      FOURDPI=4./PI
C
      DO L=2,LA
        CTMPDRY(L)=1.
C       IF(ISCDRY(L).NE.0) CTMPDRY(L)=0.
      ENDDO
C
C**********************************************************************C
C
C **  SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS
C
      IF(ISCORTBC.GE.1) THEN
C
C **  SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS
C
      DO L=1,LC
	  ISSBCP(L)=0
	ENDDO
C
	DO L=2,LA
	  IF(RSSBCE(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCW(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCN(L).GT.1.5)ISSBCP(L)=1
	  IF(RSSBCS(L).GT.1.5)ISSBCP(L)=1
	ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      IF(ISDTXBUG.EQ.1)THEN
      IF(N.EQ.1)THEN
        OPEN(1,FILE='SSEDTOX0.DIA',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='SSEDTOX0.DIA',STATUS='UNKNOWN')
        OPEN(11,FILE='SSEDTOX1.DIA',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        OPEN(11,FILE='SSEDTOX1.DIA',STATUS='UNKNOWN')
        OPEN(21,FILE='SSEDTOX2.DIA',STATUS='UNKNOWN')
        CLOSE(21,STATUS='DELETE')
        OPEN(21,FILE='SSEDTOX2.DIA',STATUS='UNKNOWN')
        OPEN(31,FILE='SSEDTOX3.DIA',STATUS='UNKNOWN')
        CLOSE(31,STATUS='DELETE')
        OPEN(31,FILE='SSEDTOX3.DIA',STATUS='UNKNOWN')
        OPEN(41,FILE='SSEDTOX4.DIA',STATUS='UNKNOWN')
        CLOSE(41,STATUS='DELETE')
        OPEN(41,FILE='SSEDTOX4.DIA',STATUS='UNKNOWN')
      ELSE
        OPEN(1,FILE='SSEDTOX0.DIA',POSITION='APPEND',STATUS='UNKNOWN')
        OPEN(11,FILE='SSEDTOX1.DIA',POSITION='APPEND',STATUS='UNKNOWN')
        OPEN(21,FILE='SSEDTOX2.DIA',POSITION='APPEND',STATUS='UNKNOWN')
        OPEN(31,FILE='SSEDTOX3.DIA',POSITION='APPEND',STATUS='UNKNOWN')
        OPEN(41,FILE='SSEDTOX4.DIA',POSITION='APPEND',STATUS='UNKNOWN')
      ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  IF N=1 CALCULATE INITIAL SEDIMENT BED THICKNESS
C **  MOVED TO SUBROUTINE BEDINIT  IN 8 AUGUST 2001 VERSION
C
C      IF(N.EQ.1)THEN
C
C      DO K=1,KB
C      DO L=2,LA
C       HBED(L,K)=0.
C       SEDBALL(L,K)=0.
C      ENDDO
C      ENDDO
C
C      IF(ISTRAN(6).GE.1)THEN
C      DO NS=1,NSED
C       DO K=1,KB
C       DO L=2,LA
C        HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
C        SEDBALL(L,K)=SEDBALL(L,K)+SEDB(L,K,NS)
C       ENDDO
C       ENDDO
C      ENDDO
C      ENDIF
C
C      IF(ISTRAN(7).GE.1)THEN
C      DO NX=1,NSND
C       NS=NSED+NX
C       DO K=1,KB
C       DO L=2,LA
C        HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
C        SEDBALL(L,K)=SEDBALL(L,K)+SNDB(L,K,NX)
C       ENDDO
C       ENDDO
C      ENDDO
C      ENDIF
C
C      TMPVAL=1./(1.-PORBED(L,K))
C      DO K=1,KB
C      DO L=2,LA
C       HBED(L,K)=TMPVAL*HBED(L,K)
C       HBED(L,K)=MAX(1.E-9,HBED(L,K))
C       VOLBW2(L,K)=HBED(L,K)
C      ENDDO
C      ENDDO
C
C ** DIAGNOSTICS OF INITIALIZATION
C
C       OPEN(2,FILE='DEPBED.DIA')
C       CLOSE(2,STATUS='DELETE')
C       OPEN(2,FILE='DEPBED.DIA')
C       DO L=2,LA
C        WRITE(2,2222)IL(L),JL(L),HP(L),SEDB(L,1,1),SNDB(L,1,1),HBED(L,1)
C       ENDDO
C       CLOSE(2)
C
C      ENDIF
C      
C**********************************************************************C
C
C **   UPDATE SEDIMENT PROCESSES
C
C----------------------------------------------------------------------C
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
CZZDIFF C APPLY CONSOLIDATION TO BULK DENSITY  
CZZDIFF C GM/CC TO BE USED FOR CRITICAL SHEARS  
CZZDIFF C  
CZZDIFF      IF(IBMECH.GT.0)THEN  
CZZDIFF        DO K=1,KB  
CZZDIFF          DO L=2,LA  
CZZDIFF            BDENBED1(L,K)=BDENBED(L,K)  
CZZDIFF            IF(HBED(L,K).GT.0.)THEN  
CZZDIFF C            BB := PORBED*GM/CC+1/1000000*(SEDB+SNDB)*M^3/(HBED*CC)  
CZZDIFF             BDENBED(L,K)=1.*PORBED(L,K)+0.000001*SEDBALL(L,K)/HBED(L,K)  
CZZDIFF            ELSE  
CZZDIFF              BDENBED(L,K)=0.  
CZZDIFF            ENDIF  
CZZDIFF          ENDDO  
CZZDIFF        ENDDO  
CZZDIFF      ENDIF  
CZZDIFF C  
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
      DO NX=1,NSED+NSND
      DO K=1,KB
        DO L=2,LA
          VFRBED(L,K,NX)=0.
          VFRBED1(L,K,NX)=0.
        ENDDO
      ENDDO
	ENDDO

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
      DO L=2,LA
        QWBDTOP(L)=0.
        QSBDTOP(L)=0.
      ENDDO
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
c	    FRACCOH(L,K)=FRACCOH(L,K)+VFRBED(L,KBT(L),NS)
	    FRACCOH(L,K)=FRACCOH(L,K)+VFRBED(L,K,NS)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO NX=1,NSND
	NS=NX+NSED
	DO K=1,KB
      DO L=2,LA
	  IF(K.LE.KBT(L))THEN
c          FRACNON(L,K)=FRACNON(L,K)+VFRBED(L,KBT(L),NS)
          FRACNON(L,K)=FRACNON(L,K)+VFRBED(L,K,NS)
        ENDIF
      ENDDO
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
        IF(IWRSPB(1).GE.1)THEN
          DO K=1,KB
            DO L=2,LA
              TAURB(L,K)=CSEDTAUB(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                          VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
              WRSPB(L,K)=CSEDRESB(BDENBED(L,K),WRSPO(1),VDRRSPO(1),
     &                          VDRBED(L,K),VDRBED(L,K),IWRSPB(1))
            ENDDO
          ENDDO
        ENDIF
C#######################################################################
C     HQI Change, 01/05/04, SO and RM
C     Change to implement HQI Sed-flume analysis based critical shear 
C     stress options
        IF(IWRSP(1).EQ.99)THEN
          DO K=1,KB
            DO L=2,LA
              TAURS(L,K)=CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),
     &                VDRBED(L,K),VDRBEDSED(L,K),IWRSP(1),L)
              TAURB(L,K)=1.0E6
C     HQI change to implement spatially varying coefficient and exponent
C     in resuspension formulation based on Ed G analysis of Sed-flume data
C     RM 12/11/03
c             WRSPS(L,K)=WRSPO(1)
              if(L.LE.265) then  ! Woods Pond
                 WRSPS(L,K) = 7.20
              else               ! North of Woods Pond
                 WRSPS(L,K) = 12.8
              endif
              WRSPB(L,K)=0.0
            ENDDO
          ENDDO
        ENDIF
        IF(IWRSP(1).EQ.999)THEN
          DO K=1,KB
            DO L=2,LA
              TAURS(L,K)=TAUCRCOH(L,K)
              TAURB(L,K)=1.0E6
C     HQI change to implement spatially varying coefficient and exponent
C     in resuspension formulation based on Ed G analysis of Sed-flume data
C     RM 12/11/03
cSO           WRSPS(L,K)=WRSPO(1)
cSO           if(L.LE.265) then  ! Woods Pond
cSO              WRSPS(L,K) = 7.20
cSO           else               ! North of Woods Pond
cSO 05.11.04     WRSPS(L,K) = 12.8
cSO              WRSPS(L,K) = 16.8
cSO           endif

cSO   HQI change to implement new coefficient and exponent in spatially
cSO       varying resuspension formula based on new partition of domain
cSO       using Sedflume data (see CALSED.FOR, line 206)
cSO 05/12/04
c             IF ( L.LE.2042 ) THEN 
c               WRSPS(L,K) = 6.69   ! All of 5C + Woods Pond
c             ELSE
c               WRSPS(L,K) = 7.78   ! North of (All of 5C + Woods Pond)
c             END IF
cSO 05/14/04
CJMH050305              IF ( L.LE.2042 ) THEN 
CJMH050305                WRSPS(L,K) = 9.88   ! All of 5C + Woods Pond
CJMH050305              ELSE
CJMH050305                WRSPS(L,K) = 6.98   ! North of (All of 5C + Woods Pond)
CJMH050305              END IF
CJMH050305
              WRSPS(L,K)=WRSPO(1)
CJMH050305
              WRSPB(L,K)=0.0
            ENDDO
          ENDDO
        ENDIF
C#######################################################################
C
        DO L=2,LA
        DO K=1,KB
	     WRSPS(L,K)=SPB(L)*WRSPS(L,K)
	     WRSPB(L,K)=SPB(L)*WRSPB(L,K)
        ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN
C **  BED AND REINITILIZE IF NECESSARY
C
c      IF(N.EQ.1.AND.ISTRAN(5).GE.1)THEN
c        IF(ISRESTI.EQ.0.OR.ISCI(5).EQ.0)THEN
C
C **  CALCULATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
c          DO NT=1,NTOX
c            NSP2(NT)=NSED+NSND
c            IF(ISTOC(NT).EQ.2) NSP2(NT)=NSP2(NT)+1
c            IF(ISTOC(NT).EQ.3) NSP2(NT)=NSP2(NT)+2
c          END DO
C
c          DO NT=1,NTOX
c            DO NS=1,NSP2(NT)
c              DO K=1,KB
c                DO L=2,LA
c                  TOXPFB(L,K,NS,NT)=0.
c                ENDDO
c              ENDDO
c            ENDDO
c          ENDDO
C
c          DO NT=1,NTOX
c           IF(ISTRAN(6).GE.1)THEN
c             DO NS=1,NSED
c               DO K=1,KB
c                 DO L=2,LA
c                   TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*TOXPARB(NS,NT)
c                 ENDDO
c               ENDDO
c             ENDDO
c            ENDIF
c            IF(ISTRAN(7).GE.1)THEN
c              DO NX=1,NSND
c                NS=NX+NSED
c                DO K=1,KB
c                  DO L=2,LA
c                    TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*TOXPARB(NS,NT)
c                  ENDDO
c                ENDDO
c              ENDDO
c            ENDIF
c            IF(ISTOC(NT).EQ.2)THEN
c              NS=1+NSED+NSND
c              DO K=1,KB
c                DO L=2,LA
c                  TOXPFB(L,K,NS,NT)=STDOCB(L,K)*TOXPARBC(1,NT)
c                ENDDO
c              ENDDO
c            ENDIF
c            IF(ISTOC(NT).EQ.3)THEN
c              NS=2+NSED+NSND
c              DO K=1,KB
c                DO L=2,LA
c                  TOXPFB(L,K,NS,NT)=STPOCB(L,K)*TOXPARBC(2,NT)
c                ENDDO
c              ENDDO
c            ENDIF
c          ENDDO
C
c          DO NT=1,NTOX
c            DO K=1,KB
c              DO L=2,LA
c                TOXPFTB(L,K,NT)=0.
c              ENDDO
c            ENDDO
c            DO NS=1,NSP2(NT)
c              DO K=1,KB
c                DO L=2,LA
c                  TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
c                ENDDO
c              ENDDO
c            ENDDO
c          ENDDO
C
c          DO NT=1,NTOX
c            DO K=1,KB
c              DO L=2,LA
c               IF(SEDBALL(L,K).GT.0.0)THEN
c                  TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
c     &                     /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
c                ELSE
c                  TOXPFTB(L,K,NT)=1.
c                ENDIF
c              ENDDO
c            ENDDO
c          ENDDO
C
C **  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
C **  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
C **  TO TOXB UNITS OF OF MG/M**2
C
c          DO NT=1,NTOX
c            IF(ITXBDUT(NT).EQ.0)THEN
c              DO K=1,KB
c                DO L=2,LA
c                  TOXB(L,K,NT)=HBED(L,K)*TOXB(L,K,NT)
c                  TOXB1(L,K,NT)=TOXB(L,K,NT)
c                ENDDO
c              ENDDO
c            ENDIF
c            IF(ITXBDUT(NT).EQ.1)THEN
c              DO K=1,KB
c                DO L=2,LA
c                  TOXB(L,K,NT)=0.001*TOXB(L,K,NT)*(SEDBT(L,K)
c     &               +SNDBT(L,K))/TOXPFTB(L,K,NT) 
c                  TOXB1(L,K,NT)=TOXB(L,K,NT)
c                ENDDO
c              ENDDO
c            ENDIF
c          ENDDO
C
C ** DIAGNOSTICS OF INITIALIZATION
C
c          IF(ISDTXBUG.EQ.1)THEN
c            OPEN(2,FILE='TOXBED.DIA')
c            CLOSE(2,STATUS='DELETE')
c            OPEN(2,FILE='TOXBED.DIA')
c            DO L=2,LA
C             TMP1=-999.
C             TMP2=-999.
C             IF(HBED(L).GT.0.)TMP1=TOXB(L,1)/HBED(L)
C             IF(HBED(L).GT.0.)TMP2=TOXB(L,2)/HBED(L)
C             WRITE(2,2222)IL(L),JL(L),HBED(L),TOXB(L,1),TOXB(L,2),TMP1,TMP2
c              TMP1=TOXB(L,1,1)/(HBED(L,1)+1.E-12)
c              WRITE(2,2222)IL(L),JL(L),TOXPFTB(L,1,1),TOXB(L,1,1),
c     &              TMP1,TOX(L,1,1)
c            ENDDO
c            CLOSE(2)
c	    ENDIF
C
c        ENDIF
c      ENDIF
C
 2222 FORMAT(2I5,7E13.4)
C      
C**********************************************************************C
C
C **  SAVE OLD VALUES
C
      IF(ISTRAN(5).GE.1)THEN
        DO NT=1,NTOX
          DO K=1,KC
            DO L=2,LA
              TOXS(L,K,NT)=TOX(L,K,NT)
            ENDDO
          ENDDO
        ENDDO
        DO NT=1,NTOX
          DO K=1,KB
            DO L=2,LA
              TOXBS(L,K,NT)=TOXB(L,K,NT)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(6).GE.1)THEN
        DO NS=1,NSED
          DO K=1,KC
            DO L=2,LA
              SEDS(L,K,NS)=SED(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              SEDBS(L,K,NS)=SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTRAN(7).GE.1)THEN
        DO NX=1,NSND
          DO K=1,KC
            DO L=2,LA
              SNDS(L,K,NX)=SND(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
        DO NX=1,NSND
          DO K=1,KB
            DO L=2,LA
              SNDBS(L,K,NX)=SNDB(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  SET MEAN D50 AND D90
C
      IF(ISTRAN(7).GE.1)THEN
C
        DO K=1,KB
          DO L=2,LA
            SEDDIA50(L,K)=0.
            SEDDIA90(L,K)=0.
            SEDGEOSTD(L,K)=0.
            SNDBT(L,K)=0.
          ENDDO
        ENDDO
C
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO L=2,LA
              SEDDIA50(L,K)=SEDDIA50(L,K)+SNDB(L,K,NX)*LOG(SEDDIA(NS))
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            IF(SNDBT(L,K).GT.0.)THEN
              SEDDIA50(L,K)=SEDDIA50(L,K)/SNDBT(L,K)
            ENDIF
          ENDDO
        ENDDO
C
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO L=2,LA
              SEDGEOSTD(L,K)=SEDGEOSTD(L,K)
     &         +SNDB(L,K,NX)*(( LOG(SEDDIA(NS))-SEDDIA50(L,K) )**2)
            ENDDO
          ENDDO
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            IF(SNDBT(L,K).GT.0.)THEN
              SEDGEOSTD(L,K)=SEDGEOSTD(L,K)/SNDBT(L,K)
            ENDIF
          ENDDO
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            IF(SNDBT(L,K).GT.0.)THEN
              SEDDIA50(L,K)=EXP(SEDDIA50(L,K))
              SEDGEOSTD(L,K)=EXP(SEDGEOSTD(L,K))
            ELSE
              SEDDIA50(L,K)=30.*ZBR(L)
              SEDGEOSTD(L,K)=1.0
            ENDIF
          ENDDO
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
            SEDDIA90(L,K)=(SEDGEOSTD(L,K)**1.28)*SEDDIA50(L,K)
          ENDDO
        ENDDO
C
        IF(ISBSDIAM.EQ.0)THEN
        DO K=1,KB
          DO L=2,LA
            SEDDIAGS(L,K)=SEDDIA50(L,K)
          ENDDO
        ENDDO
	  ENDIF
C
        IF(ISBSDIAM.EQ.1)THEN
        DO K=1,KB
          DO L=2,LA
            SEDDIAGS(L,K)=2.*SEDDIA50(L,K)
          ENDDO
        ENDDO
	  ENDIF
C
        IF(ISBSDIAM.EQ.2)THEN
        DO K=1,KB
          DO L=2,LA
            SEDDIAGS(L,K)=SEDDIA90(L,K)
          ENDDO
        ENDDO
	  ENDIF
C
        IF(ISBSDIAM.EQ.3)THEN
        DO K=1,KB
          DO L=2,LA
            SEDDIAGS(L,K)=2.*SEDDIA90(L,K)
          ENDDO
        ENDDO
	  ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  SET CELL CENTER BED STRESS FOR SEDIMENT RESUSPENSION AND DEPOSITION
C
      IF(ISWAVE.GT.0)THEN
        DO L=2,LA
          TAUBC=QQ(L,0)/CTURB2
          UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
          VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
          CURANG=ATAN2(VTMP,UTMP)
          TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     &       +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
          TAUB2=MAX(TAUB2,0.)
          TAUB(L)=SQRT(TAUB2)
          USTAR(L)=SQRT(TAUB(L))
        ENDDO
	ELSE
        DO L=2,LA
          TAUBC=QQ(L,0)/CTURB2
          TAUB(L)=TAUBC
C	    DROUGH=SEDDIA50(L,KBT(L))
C	    TOPTMP=(LOG(2.))**2
C	    BOTTMP=(LOG(0.8*ZBR(L)/DROUGH))**2
C	    TAUB(L)=TOPTMP*TAUB(L)/BOTTMP
          USTAR(L)=SQRT(TAUB(L))
        ENDDO
      ENDIF
C
      DO L=2,LA
        LN=LNC(L)
        UCELLCTR(L)=0.5*(RSSBCW(L)*WCORWST(L)*U(L,1)
     &                  +RSSBCE(L)*WCOREST(L)*U(L+1,1))
        VCELLCTR(L)=0.5*(RSSBCS(L)*WCORSTH(L)*V(L,1)
     &                  +RSSBCN(L)*WCORNTH(L)*V(LN ,1))
        QCELLCTR(L)=SQRT(UCELLCTR(L)*UCELLCTR(L)
     &                  +VCELLCTR(L)*VCELLCTR(L))
        IF(QCELLCTR(L).GT.0.0) THEN
          UCELLCTR(L)=UCELLCTR(L)/QCELLCTR(L)
          VCELLCTR(L)=VCELLCTR(L)/QCELLCTR(L)
	    CBEDTOTAL(L)=TAUB(L)/(QCELLCTR(L)*QCELLCTR(L))
	  ELSE
          UCELLCTR(L)=0.
          VCELLCTR(L)=0.
	    CBEDTOTAL(L)=0.0
        ENDIF
      ENDDO
C
C **  CALCULATE CELL CENTER DEPTH DEVIATION FROM UNIFORM FLOW
C
      IF(ISBSDFUF.GE.1)THEN
      DO L=2,LA
        LN=LNC(L)
        HDFUFXX=0.5*(RSSBCW(L)*WCORWST(L)*HDFUFX(L)
     &                  +RSSBCE(L)*WCOREST(L)*HDFUFX(L+1))
        IF(HDFUFXX.LE.0.0)HDFUFXX=1.0
        HDFUFYY=0.5*(RSSBCS(L)*WCORSTH(L)*HDFUFY(L)
     &                  +RSSBCN(L)*WCORNTH(L)*HDFUFY(LN ))
        IF(HDFUFYY.LE.0.0)HDFUFYY=1.0
        HDFUF(L)=SQRT(HDFUFXX*HDFUFXX+HDFUFYY*HDFUFYY)
CJH060305        HDFUF(L)=MIN(HDFUF(L),1.0)
      ENDDO
      ENDIF
C
      IF(KC.GT.1)THEN
	  TMPEXP=16./7.
	  DO L=2,LA
          LN=LNC(L)
	    UTMP=0.5*(RSSBCW(L)*WCORWST(L)*UHDYE(L)
     &             +RSSBCE(L)*WCOREST(L)*UHDYE(L+1))/(DYP(L)*HP(L))
          VTMP=0.5*(RSSBCS(L)*WCORSTH(L)*VHDXE(L)
     &             +RSSBCN(L)*WCORNTH(L)*VHDXE(LN ))/(DXP(L)*HP(L))
          QCELLCTRA=SQRT(UTMP*UTMP+VTMP*VTMP)
	    IF(QCELLCTR(L).GT.0.0) THEN
	      QCELLAD1SQ(L)=(QCELLCTRA/QCELLCTR(L))**2
	      QCELLAD1ZZ(L)=(QCELLCTRA/QCELLCTR(L))**TMPEXP
	    ELSE
	      QCELLAD1SQ(L)=1.
	      QCELLAD1ZZ(L)=1.
	    ENDIF
	  ENDDO
	ELSE
	  DO L=2,LA
	    QCELLAD1SQ(L)=1.
	    QCELLAD1ZZ(L)=1.
	  ENDDO
	ENDIF
C
      DO L=2,LA  
        TAUBSEDS(L)=TAUBSED(L)
	  TAUBSNDS(L)=TAUBSND(L)
	ENDDO
C
      DO L=2,LA  
        TAUBSED(L)=TAUB(L)
	  TAUBSND(L)=TAUB(L)
        USTARSED(L)=USTAR(L)
	  USTARSND(L)=USTAR(L)
	ENDDO
C
C----------------------------------------------------------------------C
C
C **  PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS
C     ORIGINAL MIXED FORM
C
      IF(ISBEDSTR.EQ.1.OR.ISBEDSTR.EQ.2)THEN  
	
	 TMPEXP=2./7. 
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
C         TVAR3S(L)=0.375E+6*HP(L)*QCELLCTR(L)
C HAMRICK CORRECTED 052504
         TMPVAL=1./(COEFTSBL*VISMUDST)
         TVAR3S(L)=TMPVAL*HP(L)*QCELLCTR(L)
         TVAR3N(L)=SEDDIAGS(L,KBT(L))*HPI(L)
	   HGDH(L)=0.0
         ENDIF
       ENDDO
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
	   IF(TVAR3S(L).GT.0.0)THEN
           TVAR3S(L)=1./TVAR3S(L)
         ENDIF
         ENDIF
       ENDDO
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
	   IF(TAUBSEDS(L).GT.0.0)THEN
           TVAR3E(L)=QCELLCTR(L)/SQRT(TAUBSEDS(L))
         ELSE
	     TVAR3E(L)=0.0
	   ENDIF
         ENDIF
       ENDDO
C  
C      TVAR3S(L)=4*KINIMATIC_VISCOSITY/(DEPTH*VELOCITY_MAGNITUDE)
C      TVAR3E(L)=VELOCITY_MAGNITUDE/COHESIVE_BEDSTRESS
C      TVAR3N(L)=D50/DEPTH
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
         TVAR3W(L)=TVAR3S(L)**0.333333
         TVAR3E(L)=TVAR3E(L)**0.333333
         TVAR3N(L)=TVAR3N(L)**0.333333
         ENDIF
       ENDDO
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
         TVAR3S(L)=TVAR3S(L)**TMPEXP
         ENDIF
       ENDDO
C  
C      TVAR3S(L)=(4*KINIMATIC_VISCOSITY/(DEPTH*VELOCITY_MAGNITUDE))**(1/4)
C      TVAR3W(L)=(4*KINIMATIC_VISCOSITY/(DEPTH*VELOCITY_MAGNITUDE))**(1/3)
C      TVAR3E(L)=(VELOCITY_MAGNITUDE/COHESIVE_BEDSTRESS)**(1/3)
C      TVAR3N(L)=(D50/DEPTH)**1/3
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
	   IF(CBEDTOTAL(L).GT.0.0)THEN
	     K=KBT(L)
C           HGDH(L)=0.025*QCELLAD1SQ(L)*(FRACCOH(L,K)*TVAR3W(L)*TVAR3E(L)
C HAMRICK CORRECTED 052504 and added non-uniform flow factor 
           HGDH(L)=0.014*HDFUF(L)*QCELLAD1SQ(L)
     &                   *(FRACCOH(L,K)*TVAR3W(L)*TVAR3E(L)
     &                    +FRACNON(L,K)*TVAR3N(L))/CBEDTOTAL(L)
	   ENDIF
         ENDIF
       ENDDO
C
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
         HGDH(L)=HGDH(L)**0.75
         ENDIF
       ENDDO
C   
c       DO L=2,LA
c	   IF(LMASKDRY(L))THEN
c         IF(HGDH(L).GT.1.0)THEN
c	     WRITE(8,869)IL(L),JL(L),HGDH(L)
c	   ENDIF
c         ENDIF
c       ENDDO
C   
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
         HGDH(L)=MIN(HGDH(L),1.0)
         ENDIF
       ENDDO
C 
C   convert hgdh to 1/hgdh  ie hdhg
C  
       DO L=2,LA
	   IF(LMASKDRY(L))THEN
	   IF(HGDH(L).GT.0.0)THEN
	     HGDH(L)=1./HGDH(L)
	   ENDIF
         ENDIF
       ENDDO
C
       DO L=2,LA 
	   IF(LMASKDRY(L))THEN
	   IF(TAUB(L).GT.0.0)THEN 
           TAUBSED(L)=HGDH(L)**TMPEXP
	     TAUBSND(L)=HGDH(L)**0.333333
	     TVAR3W(L)=QCELLCTR(L)*QCELLCTR(L)
         ENDIF
         ENDIF
	 ENDDO
C
       DO L=2,LA 
	   IF(LMASKDRY(L))THEN
	   IF(TAUB(L).GT.0.0)THEN
	     TAUBSEDS(L)=TAUBSED(L) 
	     TAUBSNDS(L)=TAUBSND(L) 
c          TAUBSED(L)=0.042*TVAR3S(L)*TAUBSED(L)*TVAR3W(L)*QCELLAD1ZZ(L)
c	     TAUBSND(L)=0.025*TVAR3N(L)*TAUBSND(L)*TVAR3W(L)*QCELLAD1SQ(L)
c hamrick corrected 052504
           TAUBSED(L)=0.026*TVAR3S(L)*TAUBSED(L)*TVAR3W(L)*QCELLAD1ZZ(L)
 	     TAUBSND(L)=0.014*TVAR3N(L)*TAUBSND(L)*TVAR3W(L)*QCELLAD1SQ(L)
         ENDIF
         ENDIF
	 ENDDO
C
C
C       OPEN(1,FILE='GRAINSTR.OUT')
C       CLOSE(1,STATUS='DELETE')
C       OPEN(1,FILE='GRAINSTR.OUT')
C       DO L=2,LA
C	   IF(LMASKDRY(L))THEN
C         WRITE(1,4344)IL(L),JL(L),HGDH(L),TAUB(L),TAUBSED(L),TVAR3S(L),
C     &     TAUBSEDS(L),TVAR3W(L),QCELLAD1ZZ(L),TAUBSND(L),
C     &     TVAR3N(L),TAUBSNDS(L),TVAR3W(L),QCELLAD1SQ(L)
C	   ENDIF
C       ENDDO
C 	CLOSE(1)
C 4344 FORMAT(2I5,12E14.5)
C
C **  IF ISBEDSTR=2, APPLY WEIGHTED AVERAGE TO BOTH SED AND SND
C
      IF(ISBEDSTR.EQ.2.OR.N.EQ.1)THEN 
       DO L=2,LA 
	   IF(LMASKDRY(L))THEN
         K=KBT(L)
         TVAR3E(L)=FRACCOH(L,K)*TAUBSED(L)+FRACNON(L,K)*TAUBSND(L)
         ENDIF
	 ENDDO
       DO L=2,LA 
	   IF(LMASKDRY(L))THEN
         TAUBSED(L)=TVAR3E(L)
         TAUBSND(L)=TVAR3E(L)
         ENDIF
	 ENDDO
	ENDIF     
C
       DO L=2,LA 
	   IF(LMASKDRY(L))THEN
	   USTARSED(L)=SQRT(TAUBSED(L))
	   USTARSND(L)=SQRT(TAUBSND(L))
         ENDIF
	 ENDDO
C
      ENDIF 
C 
C  ENDIF ON GRAIN STRESS PARTITIONING FOR ISBEDSTR.EQ.1.OR.ISBEDSTR.EQ.2
C
C----------------------------------------------------------------------C
C
C **  PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS
C     INDEPENDENTLY SET GRAIN STRESS
C
      IF(ISBEDSTR.EQ.3)THEN  
C
       DO L=2,LA 
         HDZBR=HP(L)/ZBRSED(L)
         TVAR3E(L)=0.16/( (LOG(HDZBR)-1.)**2 )         
	 ENDDO
       DO L=2,LA 
         TAUBSED(L)=TVAR3E(L)*QCELLCTR(L)*QCELLCTR(L)
	 ENDDO
       DO L=2,LA 
         TAUBSND(L)=TAUBSED(L)
	 ENDDO
       DO L=2,LA 
	   USTARSED(L)=SQRT(TAUBSED(L))
	   USTARSND(L)=SQRT(TAUBSND(L))
	 ENDDO
C
      ENDIF 
C
C----------------------------------------------------------------------C
C
C **  PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS
C     USING WEIGHTED ROUGHNESS AND LOG RESISTANCE LAW
C
      IF(ISBEDSTR.EQ.4)THEN 
C
C     RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN 
C 
      DO L=2,LA
	  IF(CBEDTOTAL(L).GT.0.0)THEN
	    TMP=EXP(1.+0.4/SQRT(CBEDTOTAL(L)))
          ZOTOTAL(L)=HP(L)/TMP
        ELSE
	    ZOTOTAL(L)=ZBR(L)
	  ENDIF
        ZOGRAIN(L)=FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
	ENDDO
C 
      DO L=2,LA
	  IF(TAUBSED(L).GT.0.0)THEN
          ZOGRAINCOH=0.041*VISMUDST/SQRT(TAUBSED(L))
          ZOGRAIN(L)=ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
          ZOGRAIN(L)=MAX(ZOGRAIN(L),ZOGRAINCOH)
        ENDIF
	ENDDO
C
C     ITERATE RVAL = SQRT(HG/H)
C 
      DO L=2,LA
        TMPTOP=LOG(HP(L)/ZOTOTAL(L))-1.
        TMPBOT=LOG(HP(L)/ZOGRAIN(L))-1.
	  RVAL=1.
	  DO KK=1,100
	    RVAL=TMPTOP/(LOG(RVAL*RVAL)+TMPBOT)
        ENDDO
        RVAL=MIN(RVAL,1.0)
        TAUBSED(L)=RVAL*RVAL*TAUB(L)
	  TAUBSND(L)=RVAL*RVAL*TAUB(L)
	ENDDO
C
      DO L=2,LA 
	  USTARSED(L)=SQRT(TAUBSED(L))
	  USTARSND(L)=SQRT(TAUBSND(L))
	ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS
C     USING WEIGHTED ROUGHNESS AND POWER RESISTANCE LAW
C
      IF(ISBEDSTR.EQ.5)THEN 
C
C     RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN 
C 
      DO L=2,LA
	  IF(CBEDTOTAL(L).GT.0.0)THEN
	    TMP=CBEDTOTAL(L)/0.04736
          ZOTOTAL(L)=HP(L)*(TMP**3)
        ELSE
	    ZOTOTAL(L)=ZBR(L)
	  ENDIF
        ZOGRAIN(L)=FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
	ENDDO
C 
      DO L=2,LA
	  IF(TAUBSED(L).GT.0.0)THEN
          ZOGRAINCOH=0.041*VISMUDST/SQRT(TAUBSED(L))
          ZOGRAIN(L)=ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
          ZOGRAIN(L)=MAX(ZOGRAIN(L),ZOGRAINCOH)
	  ENDIF
	ENDDO
C
C     CALCULATE GRAIN STRESS DIRECTLY
C 
      DO L=2,LA
        TMP=(ZOGRAIN(L)/ZOTOTAL(L))**0.25
        TMP=MIN(TMP,1.0)
        TAUBSED(L)=TMP*TAUB(L)
	  TAUBSND(L)=TMP*TAUB(L)
	ENDDO
C
      DO L=2,LA 
	  USTARSED(L)=SQRT(TAUBSED(L))
	  USTARSND(L)=SQRT(TAUBSND(L))
	ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
  869 FORMAT(' I,J,HGDH = ',2I5,F10.3)
C
C**********************************************************************C
C
        DO L=2,LA
	    HBEDA(L)=0.0
          DO K=1,KBT(L)
            HBEDA(L)=HBEDA(L)+HBED(L,K)
          END DO
        ENDDO
C
C**********************************************************************C
C
C **  CALCULATE COHESIVE SEDIMENT SETTLING, DEPOSITION AND RESUSPENSION
C
      IF(ISTRAN(6).GE.1) CALL CALSED
C
C**********************************************************************C
C
C **  CALCULATE NONCOHESIVE SEDIMENT BEDLOAD TRANSPORT, SETTLING,
C **  DEPOSITION AND RESUSPENSION
C
      IF(ISTRAN(7).GE.1) CALL CALSND
C
C**********************************************************************C
C
C **  CALCULATE BANK EROSION AND ADJUST SEDIMENT AND WATER VOLUME 
C     FLUXES
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
        IF(ISBKERO.GE.1) THEN
          CALL BANKEROSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            QSBDTOP(LBANK)=QSBDTOP(LBANK)+QSBDTOPBEBKB(NP)
            QWBDTOP(LBANK)=QWBDTOP(LBANK)+QWBDTOPBEBKB(NP)
            QSBDTOP(LCHAN)=QSBDTOP(LCHAN)+QSBDTOPBECHB(NP)
            QWBDTOP(LCHAN)=QWBDTOP(LCHAN)+QWBDTOPBECHB(NP)
          ENDDO
        ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE PARENT TO ACTIVE LAYER SEDIMENT FLUX
C
      DO L=2,LA
        QWATPA(L)=0.0		  
        QSSDPA(L)=0.0
      ENDDO
C
      IF(ISNDAL.EQ.2)THEN
C
c      DO L=2,LA
c	  SESNFA(L)=0.0
c	  SESNFP(L)=0.0
c	ENDDO
c	DO NS=1,NSED
c        DO L=2,LA
c	    SESNFA(L)=SESNFA(L)+SEDF(L,0,NS)
c	  ENDDO
c	ENDDO
c	DO NS=1,NSND
c        DO L=2,LA
c	    SESNFA(L)=SESNFA(L)+SNDF(L,0,NS)
c	  ENDDO
c	ENDDO
C
      IF(IALTYP.EQ.0)THEN
C
C     CONSTANT ACTIVE ARMOR LAYER THICKNESS
C
      DO NS=1,NSED
        DSEDGMM=1./(1.E6*SSG(NS))
        DSEDGMMI=1.E6*SSG(NS)
        DO L=2,LA
	    KTOPTP=KBT(L)
	    KTOPM1=KBT(L)-1
          QSWPOS=(QSBDTOP(L)+QWBDTOP(L))/(1.+VDRBED1(L,KTOPM1))
          QSWNEG=(QSBDTOP(L)+QWBDTOP(L))/(1.+VDRBED1(L,KTOPTP))
          SEDFPA(L,NS)=DSEDGMMI*(VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.)
     &                          +VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.))
          QSSDPA(L)=QSSDPA(L)+DSEDGMM*SEDFPA(L,NS)
          QWATPA(L)=QWATPA(L)+DSEDGMM*
     &                 ( VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.)		  
     &                  +VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.))
	    SEDB(L,KTOPTP,NS)=SEDB(L,KTOPTP,NS)+DELT*SEDFPA(L,NS)
	    SEDB(L,KTOPM1,NS)=SEDB(L,KTOPM1,NS)-DELT*SEDFPA(L,NS)
	  ENDDO
	ENDDO
C		  
      DO NX=1,NSND
	  NS=NX+NSED
        DSEDGMM=1./(1.E6*SSG(NS))
        DSEDGMMI=1.E6*SSG(NS)
        DO L=2,LA
	    KTOPTP=KBT(L)
	    KTOPM1=KBT(L)-1
          QSWPOS=(QSBDTOP(L)+QWBDTOP(L))/(1.+VDRBED1(L,KTOPM1))
          QSWNEG=(QSBDTOP(L)+QWBDTOP(L))/(1.+VDRBED1(L,KTOPTP))
          SNDFPA(L,NX)=DSEDGMMI*(VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.)
     &                          +VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.))
          QSSDPA(L)=QSSDPA(L)+DSEDGMM*SNDFPA(L,NX)
          QWATPA(L)=QWATPA(L)+DSEDGMM*
     &                 ( VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.)		  
     &                  +VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.))
	    SNDB(L,KTOPTP,NX)=SNDB(L,KTOPTP,NX)+DELT*SNDFPA(L,NX)
	    SNDB(L,KTOPM1,NX)=SNDB(L,KTOPM1,NX)-DELT*SNDFPA(L,NX)
	  ENDDO
	ENDDO
C
      ELSE
C     
C     CONSTANT ACTIVE ARMOR LAYER TOTAL SEDIMENT MASS
C
      DO NS=1,NSED
        DSEDGMM=1./(1.E6*SSG(NS))
        DSEDGMMI=1.E6*SSG(NS)
        DO L=2,LA
	    KTOPTP=KBT(L)
	    KTOPM1=KBT(L)-1
          SEDFPA(L,NS)=VFRBED(L,KTOPM1,NS)*MAX(QSBDTOP(L),0.)
     &                +VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
          QSSDPA(L)=QSSDPA(L)+SEDFPA(L,NS)
          QWATPA(L)=QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.)		  
     &                       +VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.)
          SEDFPA(L,NS)=DSEDGMMI*SEDFPA(L,NS)
	    SEDB(L,KTOPTP,NS)=SEDB(L,KTOPTP,NS)+DELT*SEDFPA(L,NS)
	    SEDB(L,KTOPM1,NS)=SEDB(L,KTOPM1,NS)-DELT*SEDFPA(L,NS)
c	    SESNFP(L)=SESNFP(L)+SEDFPA(L,NS)
	  ENDDO
	ENDDO
C		  
      DO NX=1,NSND
	  NS=NX+NSED
        DSEDGMM=1./(1.E6*SSG(NS))
        DSEDGMMI=1.E6*SSG(NS)
        DO L=2,LA
	    KTOPTP=KBT(L)
	    KTOPM1=KBT(L)-1
          SNDFPA(L,NX)=VFRBED(L,KTOPM1,NS)*MAX(QSBDTOP(L),0.)
     &                +VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
          QSSDPA(L)=QSSDPA(L)+SNDFPA(L,NX)
          QWATPA(L)=QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.)		  
     &                       +VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.)
          SNDFPA(L,NX)=DSEDGMMI*SNDFPA(L,NX)
	    SNDB(L,KTOPTP,NX)=SNDB(L,KTOPTP,NX)+DELT*SNDFPA(L,NX)
	    SNDB(L,KTOPM1,NX)=SNDB(L,KTOPM1,NX)-DELT*SNDFPA(L,NX)
c	    SESNFP(L)=SESNFP(L)+SNDFPA(L,NX)
	  ENDDO
	ENDDO
C
c	ERR=SESNFA(11)-SESNFP(11)
c	WRITE(8,8669)N,TIME,ERR,SESNFA(11),SESNFP(11),QSSDPA(11),
c     &             QSBDTOP(11),QWATPA(11),QWBDTOP(11)
c	WRITE(8,8669)KBT(11),TIME,HBED(11,KBT(11)-1),
c     &         (VFRBED(11,KBT(11)-1,NS),NS=1,NSED+NSND)
C
      ENDIF
      ENDIF
C
 8669 FORMAT('PA ERR ',I10,F10.5,8E14.6)
C
C**********************************************************************C
C
C **  UPDATE TOP BED LAYER THICKNESS AND VOID RATIO
C **  FOR DEPOSITION-RESUSPENSION STEP
C
      DO L=2,LA
        K=KBT(L)
        HBED1(L,K)=HBED(L,K)
        VDRBED1(L,K)=VDRBED(L,K)
        QWTRBEDA(L)=QSBDTOP(L)
	  QWTRBEDA1(L)=QWBDTOP(L)
        HBED(L,K)=HBED(L,K)-DELT*(QSBDTOP(L)+QWBDTOP(L))
     &                     +DELT*(QSSDPA(L)+QWATPA(L))
        TMPVAL=HBED1(L,K)/(1.+VDRBED1(L,K))
        TMPVAL=TMPVAL-DELT*(QSBDTOP(L)-QSSDPA(L))
        IF(TMPVAL.GT.0.0) THEN
          VDRBED(L,K)=(HBED(L,K)/TMPVAL)-1.
        ELSE
          VDRBED(L,K)=0.0
        END IF
      ENDDO
C
C **  UPDATE PARENT LAYER BED LAYER THICKNESS AND VOID RATIO
C **  FOR DEPOSITION-RESUSPENSION STEP
C
      IF(ISNDAL.EQ.2)THEN
C
      DO L=2,LA
        K=KBT(L)-1
        HBED1(L,K)=HBED(L,K)
        VDRBED1(L,K)=VDRBED(L,K)
        HBED(L,K)=HBED(L,K)-DELT*(QSSDPA(L)+QWATPA(L))
        TMPVAL=HBED1(L,K)/(1.+VDRBED1(L,K))
        TMPVAL=TMPVAL-DELT*QSSDPA(L)
        IF(TMPVAL.GT.0.0) THEN
          VDRBED(L,K)=(HBED(L,K)/TMPVAL)-1.
        ELSE
          VDRBED(L,K)=0.0
        END IF
      ENDDO
C
      ENDIF
C
      IF(ISNDAL.EQ.0)THEN
C
      DO K=1,KB
        DO L=2,LA
          IF(K.LT.KBT(L))THEN
            HBED1(L,K)=S3TL*HBED(L,K)+S2TL*HBED1(L,K)
            VDRBED1(L,K)=S3TL*VDRBED(L,K)+S2TL*VDRBED1(L,K)
          ENDIF
        ENDDO
      ENDDO
C
      ELSE
C
      DO K=1,KB
        DO L=2,LA
          IF(K.LT.KBT(L)-1)THEN
            HBED1(L,K)=S3TL*HBED(L,K)+S2TL*HBED1(L,K)
            VDRBED1(L,K)=S3TL*VDRBED(L,K)+S2TL*VDRBED1(L,K)
          ENDIF
        ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE TOXIC SETTLING, DEPOSITION AND RESUSPENSION
C
      IF(ISTRAN(5).GT.0)CALL CALTOX
C
C**********************************************************************C
C
C **  UPDATE TOP BED LAYER THICKNESS AND VOID RATIO
C **  FOR DEPOSITION-RESUSPENSION STEP
C **  PRESENTLY ACTIVE BEFORE THE WATER COLUMN-BED TOXICS EXCHANGE
C **  CHECK PLACEMENT THERE AND HERE FOR 
C
C      DO L=2,LA
C        K=KBT(L)
C        HBED(L,K)=HBED1(L,K)-DELT*(QSBDTOP(L)+QWBDTOP(L))
C        TMPVAL=HBED1(L,K)/(1.+VDRBED(L,K))
C        TMPVAL=TMPVAL-DELT*QSBDTOP(L)
C        VDRTMP=(TMPVAL/HBED(L,K))-1.
C        HBED1(L,K)=S3TL*HBED(L,K)+S2TL*HBED1(L,K)
C        VDRBED1(L,K)=S3TL*VDRBED(L,K)+S2TL*VDRBED1(L,K)
C        VDRBED(L,K)=VDRTMP
C      ENDDO
C
C      DO K=1,KB
C      DO L=2,LA
C       IF(K.LT.KBT(L))THEN
C         HBED1(L,K)=S3TL*HBED(L,K)+S2TL*HBED1(L,K)
C         VDRBED1(L,K)=S3TL*VDRBED(L,K)+S2TL*VDRBED1(L,K)
C       ENDIF
C      ENDDO
C      ENDDO
C
C**********************************************************************C
C
C **  UPDATE SEDIMENT BED LAYERING
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)CALL CALBLAY
C
C**********************************************************************C
C
C **  RESET SEDIMENT VOLUME FRACTIONS
C
      DO K=1,KB
        DO L=2,LA
          BEDLINIT(L,K)=0.
        ENDDO
      ENDDO
C
      DO NX=1,NSED+NSND
      DO K=1,KB
        DO L=2,LA
          VFRBED(L,K,NX)=0.
        ENDDO
      ENDDO
	ENDDO

      IF(ISTRAN(6).GE.1)THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
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
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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
        ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO NX=1,NSND
	NS=NX+NSED
	DO K=1,KB
      DO L=2,LA
	  IF(K.LE.KBT(L))THEN
          FRACNON(L,K)=FRACNON(L,K)+VFRBED(L,K,NS)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C **  CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS
C
      DO K=1,KB
      DO L=2,LA
	  IF(K.LE.KBT(L))THEN
	    VDRBEDSND(L,K)=SNDVDRD
          VDRBEDSED(L,K)=0.0
          IF(FRACCOH(L,K).GT.0.0)THEN
            VDRBEDSED(L,K)=
     &        ( (FRACCOH(L,K)+FRACNON(L,K))*VDRBED(L,K)
     &         -FRACNON(L,K)*SNDVDRD )/FRACCOH(L,K)
          ENDIF
        ELSE
	    VDRBEDSND(L,K)=0.0
	    VDRBEDSED(L,K)=0.0
	  ENDIF
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  UPDATE SEDIMENT BED PHYSICAL PROPERTIES
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN
	  IF(IBMECH.EQ.9)THEN
	    CALL CALBED9
	  ELSE
	    CALL CALBED
	  ENDIF
	ENDIF
C
C**********************************************************************C
C
C ++  CHANGE BED MORPHOLOGY
C
      IF(IMORPH.GT.0.OR.ISGWIT.GE.2)THEN
C
        DO L=2,LA
          TVAR3S(L)=HBEDA(L)
          BELV1(L)=BELV(L)
          HTMP(L)=HP(L)
          H1P(L)=HP(L)
          P1(L)=P(L)
        ENDDO
C
        DO L=2,LA
          HBEDA(L)=0.0
          DELBED(L)=TVAR3S(L)
        ENDDO
C
        DO K=1,KB
          DO L=2,LA
	      IF(K.LE.KBT(L))THEN
              HBEDA(L)=HBEDA(L)+HBED(L,K)
            ENDIF
          END DO
        ENDDO
C
        DO L=2,LA
	    IF(HBEDA(L).NE.DELBED(L))THEN
             DELBED(L)=DELBED(L)-HBEDA(L)
          ELSE
	       DELBED(L)=0.0
          ENDIF
        ENDDO
C

        DO L=2,LA
          BELV(L)=ZELBEDA(L)+HBEDA(L)
        ENDDO
C
	  IF(ISGWIT.GE.2)THEN
          DO L=2,LA
C            DELBED(L)=TVAR3S(L)-HBEDA(L)+DELT*QGW(L)*DXYIP(L)
            HP(L)=HP(L)+DELBED(L)+DELT*QGW(L)*DXYIP(L)
            P(L)=P(L)+G*DELT*QGW(L)*DXYIP(L)
	      WSEL1=GI*P1(L)
	      WSEL=GI*P(L)
c	      write(8,8800)L,P1(L),P(L),WSEL1,WSEL,H1P(L),HP(L)
          ENDDO
	  ELSE
          DO L=2,LA
C            DELBED(L)=TVAR3S(L)-HBEDA(L)
            HP(L)=HP(L)+DELBED(L)
          ENDDO
	  ENDIF
C
        DO L=2,LA
	    HPI(L)=1./HP(L)
	    QMORPH(L)=DELTI*DXYP(L)*(HP(L)-H1P(L))
        ENDDO
C
        ITMP=0
        DO L=2,LA
c          IF(JL(L).EQ.233)THEN
c            IF(IL(L).GE.29.AND.IL(L).LE.33)THEN
c	        WRITE(8,2346)IL(L),JL(L),HP(L),HTMP(L),DELBED(L),
c     &                         TVAR3S(L),HBEDA(L)
c	      ENDIF
c	    ENDIF
          IF(HP(L).LT.0.0)THEN
	      ITMP=1
	      WRITE(8,2345)IL(L),JL(L),HBED1(L,KBT(L)),HBED(L,KBT(L)),
     &             BELV1(L),BELV(L),H1P(L),HP(L),ZELBEDA(L),DELBED(L)
	      WRITE(8,2347)L,KBT(L),(HBED(L,K),K=1,KBT(L)),HBEDA(L)
	      WRITE(8,2347)L,KBT(L),(HBED1(L,K),K=1,KBT(L)),HBEDA(L)
	      WRITE(8,2347)L,KBT(L),DELT,QSBDTOP(L),QWBDTOP(L)
	    ENDIF
        ENDDO
C
        IF(ITMP.EQ.1)THEN
          CALL RESTOUT(1)
          IF(NDRYSTP.LT.0) THEN
            OPEN(1,FILE='DRYLOSS.OUT')
            CLOSE(1,STATUS='DELETE')
            OPEN(1,FILE='DRYLOSS.OUT')
            DO L=2,LA
	        IF(VDWASTE(L).GT.0.0)THEN
	          TMPVAL=VDWASTE(L)/DXYP(L)
                WRITE(1,1993)IL(L),JL(L),VDWASTE(L),TMPVAL,QDWASTE(L)
              ENDIF
            ENDDO
            CLOSE(1)
          ENDIF
	    STOP
	  ENDIF
C
      END IF
C
 2345 FORMAT('NEG DEPTH DUE TO MORPH CHANGE', 2I5,12F12.5)
 2347 FORMAT('                             ', 2I5,12F12.5)
 2346 FORMAT('MORP ERR ',2I5,6E15.6)
 1993 FORMAT(2I6,4E14.6)
C
C**********************************************************************C
C
C ++  ADJUST CONCENTRATIONS OF TRANSPORT VARIABLES IN RESPONSE TO
C ++  CHANGE IN BED MORPHOLOGY
C
      IF(IMORPH.GT.0.OR.ISGWIT.GE.2)THEN
C
        IF(ISTRAN(1).GT.0)THEN
          DO K=1,KC
            DO L=2,LA
              SAL(L,K)=HTMP(L)*SAL(L,K)/HP(L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(2).GT.0)THEN
          DO K=1,KC
            DO L=2,LA
              TEM(L,K)=HTMP(L)*TEM(L,K)/HP(L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(3).GT.0)THEN
          DO K=1,KC
            DO L=2,LA
              DYE(L,K)=HTMP(L)*DYE(L,K)/HP(L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(4).GT.0)THEN
          DO K=1,KC
            DO L=2,LA
              SFL(L,K)=HTMP(L)*SFL(L,K)/HP(L)
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(5).GT.0)THEN
	    DO NT=1,NTOX
            DO K=1,KC
              DO L=2,LA
                TOX(L,K,NT)=HTMP(L)*TOX(L,K,NT)/HP(L)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(6).GT.0)THEN
	    DO NS=1,NSED
            DO K=1,KC
              DO L=2,LA
                SED(L,K,NS)=HTMP(L)*SED(L,K,NS)/HP(L)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(ISTRAN(7).GT.0)THEN
	    DO NS=1,NSND
            DO K=1,KC
              DO L=2,LA
                SND(L,K,NS)=HTMP(L)*SND(L,K,NS)/HP(L)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
      END IF
C
C**********************************************************************C
C
C **  POREWATER ADVECTION AND DIFFUSION OF TOXICS
C
      IF(ISTRAN(5).GT.0)CALL CALTOXB
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT REACTIONS
C
      IF(ISTRAN(5).GE.1)CALL TOXCHEM
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT OUTPUT TO FOOD CHAIN MODEL
C
      IF(ISTRAN(5).GE.1.AND.ISFDCH.GE.1)CALL FOODCHAIN(0)
C
C**********************************************************************C
C
 8800 FORMAT(I5,8E14.5)
C
      CLOSE(1)
      CLOSE(11)
      CLOSE(21)
      CLOSE(31)
      CLOSE(41)
C
C**********************************************************************C
C
      RETURN
      END
