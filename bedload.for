C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BEDLOAD(NX,NS)
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

C **  BED LOAD TRANSPORT HORIZONTAL LOOP
C
C **  INITIALIZE BED LOAD TRANSPORTS
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
          IF(ISBDLDBC.EQ.0)RETURN
C
C **  CALCULATE CELL CENTER TRANSPORT RATES USING GENERIC BED LOAD EQUATION
C
          IF(ISBDLD(NS).EQ.0)THEN
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            CSHIELDS=TCSHIELDS(NS)
            IF(ISEDEFF.EQ.2)THEN
              TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
              CSHIELDS=TMPVAL*CSHIELDS
            ENDIF
            BDLDTMPP=SBDLDP(NX)
            BDLDTMP=SQRT(GPDIASED)*DIASED/DSEDGMM 
            IF(BDLDTMPP.GT.0.0)THEN
              FACBEDL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
              SHIELDS=TAUBSND(L)/GPDIASED
              IF(SHIELDS.GT.CSHIELDS)THEN
                IF(SBDLDA(NX).GT.0.0)THEN
                  BDLDTMPA=(SBDLDG1(NX)*SHIELDS
     &                    -SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
                ELSE
                  BDLDTMPA=1.0
                ENDIF
                IF(SBDLDB(NX).GT.0.0)THEN
                  BDLDTMPB=(SBDLDG3(NX)*SQRT(SHIELDS)
     &                   -SBDLDG4(NX)*SQRT(CSHIELDS))**SBDLDB(NX)
                ELSE
                  BDLDTMPB=1.0
                ENDIF
                QSBDLDP(L)=FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*
     &                     BDLDTMPP*BDLDTMPA*BDLDTMPB
C
          IF(ISEDEFF.EQ.1)THEN
              TMPVAL=EXP(-COEHEFF*FRACCOH(L,KBT(L)))
              QSBDLDP(L)=TMPVAL*QSBDLDP(L)
          ENDIF
C
              ENDIF
            ENDIF
          ENDIF
          ENDDO
C 
          ENDIF
C
C **  CALCULATE CELL CENTER TRANSPORT RATES USING VAN RIJN BED LOAD EQUATION
C
          IF(ISBDLD(NS).EQ.1)THEN
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            CSHIELDS=TCSHIELDS(NS)
            IF(ISEDEFF.EQ.2)THEN
              TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
              CSHIELDS=TMPVAL*CSHIELDS
            ENDIF
            IF(ISNDAL.EQ.1)THEN
              TMPVAL=LOG10(19.*DIASED/SEDDIA50(L,KBT(L)))
              TMPVAL=1.66667/(TMPVAL**2)
              CSHIELDSC=CSHIELDS50(L)*TMPVAL
            ELSE
              CSHIELDSC=TCSHIELDS(NS)
            END IF
            IF(ISEDEFF.EQ.2)THEN
              TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
              CSHIELDSC=TMPVAL*CSHIELDS
            ENDIF
            BDLDTMPP=FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),
     &           PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
            IF(ISNDAL.EQ.1)THEN
              BDLDTMPP=((DIASED/SEDDIA50(L,KBT(L)))**0.3)*BDLDTMPP
            ENDIF
            BDLDTMP=SQRT(GPDIASED)*DIASED/DSEDGMM 
            FACBEDL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
            SHIELDS=TAUBSND(L)/GPDIASED
            IF(SHIELDS.GT.CSHIELDSC.OR.ISBDLD(NS).GT.1)THEN
              BDLDTMPA=(SBDLDG1(NX)*SHIELDS
     &                -SBDLDG2(NX)*CSHIELDSC)**SBDLDA(NX)
              QSBDLDP(L)=FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*
     &                     BDLDTMPP*BDLDTMPA
C
C
          IF(ISEDEFF.EQ.1)THEN
              TMPVAL=EXP(-COEHEFF*FRACCOH(L,KBT(L)))
              QSBDLDP(L)=TMPVAL*QSBDLDP(L)
          ENDIF
C
C
            ENDIF
          ENDIF
          ENDDO
C 
          ENDIF
C
C **  CALCULATE CELL CENTER TRANSPORT RATES USING ENGULAND-HANSEN
C
          IF(ISBDLD(NS).EQ.2)THEN
C
          IF(IBLTAUC(NS).EQ.0) CSHIELDS=0.
          IF(IBLTAUC(NS).EQ.1) CSHIELDS=TCSHIELDS(NS)
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
cjh0917
C     ADDED OPTION TO THE CALCULATE CRITICAL SHIELD FOR A SEDIMENT
C     CLASS IN A MIXTURE AND INCREASE CRITICAL STRESS DUE
C     PRESENCE OF COHESIVE MATERIAL
C
            IF(IBLTAUC(NS).EQ.2) 
     &         CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
            IF(IBLTAUC(NS).EQ.3) 
     &         CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
            IF(ISEDEFF.EQ.2)THEN
              TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
              CSHIELDS=TMPVAL*CSHIELDS
            ENDIF
cjh0917
            SHIELDS=TAUBSND(L)/GPDIASED
            IF(SHIELDS.GE.CSHIELDS)THEN
              BDLDTMPP=FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),
     &           PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
cjh0917
C     NOTE HGDH IS ACTUALLY H/HG AND THE CORRECTION MAKES THE
C     EH FORMULATION CONSISTENT WITH THE GRAIN STRESS SEPARATION 
C
              IF(HGDH(L).GT.0.0) BDLDTMPP=BDLDTMPP/(HGDH(L)**0.333)
cjh0917
              BDLDTMP=SQRT(GPDIASED)*DIASED/DSEDGMM 
              FACBEDL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
              BDLDTMPA=(SBDLDG1(NX)*SHIELDS)**SBDLDA(NX)
              QSBDLDP(L)=FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*
     &                     BDLDTMPP*BDLDTMPA
C
C
              IF(ISEDEFF.EQ.1)THEN
              TMPVAL=EXP(-COEHEFF*FRACCOH(L,KBT(L)))
              QSBDLDP(L)=TMPVAL*QSBDLDP(L)
              ENDIF
C
C
            ENDIF
          ENDIF
          ENDDO
C 
          ENDIF
C
C **  CALCULATE CELL CENTER TRANSPORT RATES USING WU, WANG, AND JIA
C
          IF(ISBDLD(NS).EQ.3)THEN
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            CSHIELDS=TCSHIELDS(NS)
            BDLDTMPP=FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),
     &           PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
            CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
            IF(ISEDEFF.EQ.2)THEN
              TMPVAL=1.+(COEHEFF2-1.)
     &                 *( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
              CSHIELDS=TMPVAL*CSHIELDS
            ENDIF
            BDLDTMP=SQRT(GPDIASED)*DIASED/DSEDGMM 
            FACBEDL(L)=FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),
     &                   RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
            SHIELDS=TAUBSND(L)/GPDIASED
            IF(SHIELDS.GT.CSHIELDS)THEN
              BDLDTMPA=(SBDLDG1(NX)*SHIELDS
     &                 -SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
              QSBDLDP(L)=FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*
     &                     BDLDTMPP*BDLDTMPA
C
C
          IF(ISEDEFF.EQ.1)THEN
              TMPVAL=EXP(-COEHEFF*FRACCOH(L,KBT(L)))
              QSBDLDP(L)=TMPVAL*QSBDLDP(L)
          ENDIF
C
C
            ENDIF
          ENDIF
          ENDDO
C 
          ENDIF
C
          DO L=2,LA
            QSBDLDP(L)=SPB(L)*QSBDLDP(L)
          ENDDO
C
C **  CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION
C
          IF(ISBLFUC.EQ.0)THEN
          DO L=2,LA
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            IF(LMASKDRY(L))THEN
              LN=LNC(L)
              IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).GE.0.0)THEN
                QSBDLDX(L+1,NX)=SUB(L+1)*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(LN ,NX)=SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
              ENDIF
              IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).LT.0.0)THEN
                QSBDLDX(L+1,NX)=SUB(L+1)*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(L  ,NX)=SVB(L  )*QSBDLDP(L)*VCELLCTR(L)
              ENDIF
              IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).GE.0.0)THEN
                QSBDLDX(L  ,NX)=SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(LN ,NX)=SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
              ENDIF
              IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
                QSBDLDX(L  ,NX)=SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(L  ,NX)=SVB(L  )*QSBDLDP(L)*VCELLCTR(L)
              ENDIF
            ENDIF
          ENDDO
          ENDIF
C
C **  CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION
C **  WITH CORNER EFFECTS CORRECTION
C
          IF(ISBLFUC.EQ.1)THEN
          DO L=2,LA
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            IF(LMASKDRY(L))THEN
              LN=LNC(L)
              IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).GE.0.0)THEN
                UCELLCTRM=SUB(L+1)*ABS(UCELLCTR(L))
                VCELLCTRM=SVB(LN )*ABS(VCELLCTR(L))
                QSBDLDX(L+1,NX)=SUB(L+1)*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(LN ,NX)=SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
                IF(UCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDY(LN ,NX)=SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
                ENDIF
                IF(VCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDX(L+1,NX)=SUB(L+1)*SIGN(QSBDLDP(L),UCELLCTR(L))
                ENDIF
              ENDIF
              IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).LT.0.0)THEN
                UCELLCTRM=SUB(L+1)*ABS(UCELLCTR(L))
                VCELLCTRM=SVB(L  )*ABS(VCELLCTR(L))
                QSBDLDX(L+1,NX)=SUB(L+1)*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(L  ,NX)=SVB(L  )*QSBDLDP(L)*VCELLCTR(L)
                IF(UCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDY(L  ,NX)=SVB(L  )*SIGN(QSBDLDP(L),VCELLCTR(L))
                ENDIF
                IF(VCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDX(L+1,NX)=SUB(L+1)*SIGN(QSBDLDP(L),UCELLCTR(L))
                ENDIF
              ENDIF
              IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).GE.0.0)THEN
                UCELLCTRM=SUB(L  )*ABS(UCELLCTR(L))
                VCELLCTRM=SVB(LN )*ABS(VCELLCTR(L))
                QSBDLDX(L  ,NX)=SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(LN ,NX)=SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
                IF(UCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDY(LN ,NX)=SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
                ENDIF
                IF(VCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDX(L  ,NX)=SUB(L  )*SIGN(QSBDLDP(L),UCELLCTR(L))
                ENDIF
              ENDIF
              IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
                UCELLCTRM=SUB(L  )*ABS(UCELLCTR(L))
                VCELLCTRM=SVB(L  )*ABS(VCELLCTR(L))
                QSBDLDX(L  ,NX)=SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
                QSBDLDY(L  ,NX)=SVB(L  )*QSBDLDP(L)*VCELLCTR(L)
                IF(UCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDY(L  ,NX)=SVB(L  )*SIGN(QSBDLDP(L),VCELLCTR(L))
                ENDIF
                IF(VCELLCTRM.LT.1.0E-9)THEN
                  QSBDLDX(L  ,NX)=SUB(L  )*SIGN(QSBDLDP(L),UCELLCTR(L))
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          ENDIF
C
C **  CALCULATE CELL FACE TRANSPORT RATES BY AVERAGING 
C **  VECTOR COMPONENTS FROM CELL CENTERS TO FACES
C
          IF(ISBLFUC.EQ.2)THEN
          DO L=2,LA
            QSBDLDOT(L,NX)=0.
            QSBDLDIN(L,NX)=0.
            LS=LSC(L)
            QSBDLDX(L,NX)=0.5*SUB(L)*(QSBDLDP(L)*UCELLCTR(L)
     &                               +QSBDLDP(L-1)*UCELLCTR(L+1))
            QSBDLDY(L,NX)=0.5*SVB(L)*(QSBDLDP(L)*VCELLCTR(L)
     &                               +QSBDLDP(LS )*UCELLCTR(LS ))
          ENDDO
          ENDIF
C
C **  CONVERT TRANSPORT VECTORS TO FACE VECTORS
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            QSBDLDX(L,NX)=SUB(L)*DYU(L)*QSBDLDX(L,NX)
            QSBDLDY(L,NX)=SVB(L)*DXV(L)*QSBDLDY(L,NX)
          ENDIF
          ENDDO
C
C **  ELIMINATE BEDLOAD TRANSPORT UP ADVERSE SLOPES IN DIRECTION OF FLOW
C
          IF(BLBSNT.GT.0.0)THEN
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            IF(QSBDLDX(L,NX).GT.0.0)THEN
              SLOPE=(BELV(L)-BELV(L-1))*DXIU(L)
              IF(SLOPE.GT.BLBSNT) QSBDLDX(L,NX)=0.0
            ENDIF
            IF(QSBDLDX(L,NX).LT.0.0)THEN
              SLOPE=(BELV(L-1)-BELV(L))*DXIU(L)
              IF(SLOPE.GT.BLBSNT) QSBDLDX(L,NX)=0.0
            ENDIF
            IF(QSBDLDY(L,NX).GT.0.0)THEN
              SLOPE=(BELV(L)-BELV(LSC(L)))*DYIV(L)
              IF(SLOPE.GT.BLBSNT) QSBDLDY(L,NX)=0.0
            ENDIF
            IF(QSBDLDY(L,NX).LT.0.0)THEN
              SLOPE=(BELV(LSC(L))-BELV(L))*DYIV(L)
              IF(SLOPE.GT.BLBSNT) QSBDLDY(L,NX)=0.0
            ENDIF
          ENDIF
          ENDDO
          ENDIF
C
C **  INSERT OUTFLOW OR RECIRCULATION BOUNDARY CONDITION FOR BED
C **  BED LOAD TRANSPORT 
c     hamrick fixed to add multiplication by width of out and inflow sections
c     in 0819 code patch
C
          IF(NSBDLDBC.GT.0) THEN 
            QSBLLDXY=0.0
            DO NSB=1,NSBDLDBC
              LUTMP=LSBLBCU(NSB)
              LDTMP=LSBLBCD(NSB)
              IF(LDTMP.GT.0) THEN
C               OUTFLOW ON POSITIVE X FACE RECIRCULATED TO NEGATIVE X FACE
                IF(ISDBLDIR(NSB).EQ.1)THEN
                  IF(UCELLCTR(LUTMP).GT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                DYP(LUTMP)*QSBDLDP(LUTMP)*UCELLCTR(LUTMP)
                    QSBDLDIN(LDTMP,NX)=QSBDLDOT(LUTMP,NX)
                  ENDIF
                ENDIF
C               OUTFLOW ON NEGATIVE X FACE RECIRCULATED TO POSITIVE X FACE
                IF(ISDBLDIR(NSB).EQ.-1)THEN
                  IF(UCELLCTR(LUTMP).LT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                -DYP(LUTMP)*QSBDLDP(LUTMP)*UCELLCTR(LUTMP)
                    QSBDLDIN(LDTMP,NX)=QSBDLDOT(LUTMP,NX)
                  ENDIF
                ENDIF
C               OUTFLOW ON POSITIVE Y FACE RECIRCULATED TO NEGATIVE Y FACE
                IF(ISDBLDIR(NSB).EQ.2)THEN
                  IF(VCELLCTR(LUTMP).GT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                DXP(LUTMP)*QSBDLDP(LUTMP)*VCELLCTR(LUTMP)
                    QSBDLDIN(LDTMP,NX)=QSBDLDOT(LUTMP,NX)
                  ENDIF
                ENDIF
C               OUTFLOW ON NEGATIVE Y FACE RECIRCULATED TO POSITIVE Y FACE
                IF(ISDBLDIR(NSB).EQ.-2)THEN
                  IF(VCELLCTR(LUTMP).LT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                -DXP(LUTMP)*QSBDLDP(LUTMP)*VCELLCTR(LUTMP)
                    QSBDLDIN(LDTMP,NX)=QSBDLDOT(LUTMP,NX)
                  ENDIF
                ENDIF
              ELSE
C               OUTFLOW ON POSITIVE X FACE 
                IF(ISDBLDIR(NSB).EQ.1)THEN
                  IF(UCELLCTR(LUTMP).GT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                DYP(LUTMP)*QSBDLDP(LUTMP)*UCELLCTR(LUTMP)
                  ENDIF
                ENDIF
C               OUTFLOW ON NEGATIVE X FACE 
                IF(ISDBLDIR(NSB).EQ.-1)THEN
                  IF(UCELLCTR(LUTMP).LT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                -DYP(LUTMP)*QSBDLDP(LUTMP)*UCELLCTR(LUTMP)
                  ENDIF
                ENDIF
C               OUTFLOW ON POSITIVE X FACE 
                IF(ISDBLDIR(NSB).EQ.2)THEN
                  IF(VCELLCTR(LUTMP).GT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                DXP(LUTMP)*QSBDLDP(LUTMP)*VCELLCTR(LUTMP)
                  ENDIF
                ENDIF
C               OUTFLOW ON NEGATIVE Y FACE 
                IF(ISDBLDIR(NSB).EQ.-2)THEN
                  IF(VCELLCTR(LUTMP).LT.0.0)THEN
                    QSBDLDOT(LUTMP,NX)=
     &                -DXP(LUTMP)*QSBDLDP(LUTMP)*VCELLCTR(LUTMP)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
C
C **  LIMIT OUTGOING FLUXES IN EACH CELL
C
          DO L=2,LA
          IF(LMASKDRY(L))THEN
            LN=LNC(L)
            IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).GE.0.0)THEN
              SNDFBL(L,NX)=DXYIP(L)*(QSBDLDX(L+1,NX)+QSBDLDY(LN,NX)
     &                               +QSBDLDOT(L,NX))
              ASNDFBL=ABS(SNDFBL(L,NX))
              SNDBTMP=SNDB(L,KBT(L),NX)-DELT*SNDFBL(L,NX)
              IF(SNDBTMP.LT.0.0.AND.ASNDFBL.GT.0.0)THEN
                SNDFBLM=0.5*SNDB(L,KBT(L),NX)/DELT
                FLUXFAC=SNDFBLM/SNDFBL(L,NX)
                QSBDLDX(L+1,NX)=FLUXFAC*QSBDLDX(L+1,NX)
                QSBDLDY(LN ,NX)=FLUXFAC*QSBDLDY(LN ,NX)
                QSBDLDOT(L ,NX)=FLUXFAC*QSBDLDOT(L ,NX)
              ENDIF
            ENDIF
            IF(UCELLCTR(L).GE.0.0.AND.VCELLCTR(L).LT.0.0)THEN
              SNDFBL(L,NX)=DXYIP(L)*(-QSBDLDY(L,NX)+QSBDLDX(L+1,NX)
     &                               +QSBDLDOT(L,NX))
              ASNDFBL=ABS(SNDFBL(L,NX))
              SNDBTMP=SNDB(L,KBT(L),NX)+DELT*SNDFBL(L,NX)
              IF(SNDBTMP.LT.0.0.AND.ASNDFBL.GT.0.0)THEN
                SNDFBLM=0.5*SNDB(L,KBT(L),NX)/DELT
                FLUXFAC=SNDFBLM/SNDFBL(L,NX)
                QSBDLDX(L+1,NX)=FLUXFAC*QSBDLDX(L+1,NX)
                QSBDLDY(L  ,NX)=FLUXFAC*QSBDLDY(L  ,NX)
                QSBDLDOT(L ,NX)=FLUXFAC*QSBDLDOT(L ,NX)
              ENDIF
            ENDIF
            IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).GE.0.0)THEN
              SNDFBL(L,NX)=DXYIP(L)*(-QSBDLDX(L,NX)+QSBDLDY(LN,NX)
     &                               +QSBDLDOT(L,NX))
              ASNDFBL=ABS(SNDFBL(L,NX))
              SNDBTMP=SNDB(L,KBT(L),NX)+DELT*SNDFBL(L,NX)
              IF(SNDBTMP.LT.0.0.AND.ASNDFBL.GT.0.0)THEN
                SNDFBLM=0.5*SNDB(L,KBT(L),NX)/DELT
                FLUXFAC=SNDFBLM/SNDFBL(L,NX)
                QSBDLDX(L ,NX)=FLUXFAC*QSBDLDX(L ,NX)
                QSBDLDY(LN,NX)=FLUXFAC*QSBDLDY(LN,NX)
                QSBDLDOT(L ,NX)=FLUXFAC*QSBDLDOT(L ,NX)
              ENDIF
            ENDIF
            IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
              SNDFBL(L,NX)=DXYIP(L)*(-QSBDLDX(L,NX)-QSBDLDY(L,NX)
     &                               +QSBDLDOT(L,NX))
              ASNDFBL=ABS(SNDFBL(L,NX))
              SNDBTMP=SNDB(L,KBT(L),NX)+DELT*SNDFBL(L,NX)
              IF(SNDBTMP.LT.0.0.AND.ASNDFBL.GT.0.0)THEN
                SNDFBLM=0.5*SNDB(L,KBT(L),NX)/DELT
                FLUXFAC=SNDFBLM/SNDFBL(L,NX)
                QSBDLDX(L,NX)=FLUXFAC*QSBDLDX(L,NX)
                QSBDLDY(L,NX)=FLUXFAC*QSBDLDY(L,NX)
                QSBDLDOT(L ,NX)=FLUXFAC*QSBDLDOT(L ,NX)
              ENDIF
            ENDIF
          ENDIF
          ENDDO
C
          IF(NSBDLDBC.GT.0) THEN 
            DO NSB=1,NSBDLDBC
              LUTMP=LSBLBCU(NSB)
              LDTMP=LSBLBCD(NSB)
              IF(LDTMP.GT.0) THEN
                QSBDLDOT(LUTMP,NX)=QSBDLDIN(LDTMP,NX)
              ENDIF
            ENDDO
          ENDIF
C
C **  CALCULATE EQUIVALENT CONCENTRATIONS
C
          DO L=2,LA
            LN=LNC(L)
            CQBEDLOADX(L,NX)=0.0
            CQBEDLOADY(L,NX)=0.0
            IF(LMASKDRY(L))THEN
              UVARTMP=0.5*(RSSBCW(L)*U(L,1)+RSSBCE(L)*U(L+1,1))
              VVARTMP=0.5*(RSSBCS(L)*V(L,1)+RSSBCN(L)*V(LN ,1))
              VELMAG=SQRT(UVARTMP*UVARTMP+VVARTMP*VVARTMP)
              IF(VELMAG.GT.0.0) THEN
chqi's form
                  CQBEDLOADX(L,NX)=QSBDLDP(L)/(HP(L)*VELMAG)
                  CQBEDLOADY(L,NX)=0.
chqi's form
Chamorg                     IF(UCELLCTR(L).GT.0.0.AND.VCELLCTR(L).GT.0.0)THEN
Chamorg                  CQBEDLOADX(L,NX)=QSBDLDX(L+1,NX)
Chamorg     &                            /(DYU(L+1)*HP(L)*UCELLCTR(L)*VELMAG)
Chamorg                  CQBEDLOADY(L,NX)=QSBDLDY(LN ,NX)
Chamorg     &                            /(DXV(LN )*HP(L)*VCELLCTR(L)*VELMAG)
Chamorg             ENDIF
Chamorg             IF(UCELLCTR(L).GT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
Chamorg                  CQBEDLOADX(L,NX)=QSBDLDX(L+1,NX)
Chamorg     &                            /(DYU(L+1)*HP(L)*UCELLCTR(L)*VELMAG)
Chamorg                  CQBEDLOADY(L,NX)=QSBDLDY(L  ,NX)
Chamorg     &                            /(DXV(L  )*HP(L)*VCELLCTR(L)*VELMAG)
Chamorg             ENDIF
Chamorg             IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).GT.0.0)THEN
Chamorg                  CQBEDLOADX(L,NX)=QSBDLDX(L  ,NX)
Chamorg     &                            /(DYU(L  )*HP(L)*UCELLCTR(L)*VELMAG)
Chamorg                  CQBEDLOADY(L,NX)=QSBDLDY(LN ,NX)
Chamorg     &                            /(DXV(LN )*HP(L)*VCELLCTR(L)*VELMAG)
Chamorg             ENDIF
Chamorg             IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
Chamorg                  CQBEDLOADX(L,NX)=QSBDLDX(L  ,NX)
Chamorg     &                            /(DYU(L  )*HP(L)*UCELLCTR(L)*VELMAG)
Chamorg                  CQBEDLOADY(L,NX)=QSBDLDY(L  ,NX)
Chamorg     &                            /(DXV(L  )*HP(L)*VCELLCTR(L)*VELMAG)
Chamorg             ENDIF
chammod                     IF(UCELLCTR(L).GT.0.0.AND.VCELLCTR(L).GT.0.0)THEN
chammod                  CQBEDLOADX(L,NX)=(QSBDLDX(L+1,NX)/DYU(L+1))**2
chammod     &                            +(QSBDLDY(LN ,NX)/DXV(LN ))**2
chammod                  CQBEDLOADX(L,NX)=SQRT(CQBEDLOADX(L,NX))/(HP(L)*VELMAG)
chammod                  CQBEDLOADY(L,NX)=0.0
chammod             ENDIF
chammod             IF(UCELLCTR(L).GT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
chammod                  CQBEDLOADX(L,NX)=(QSBDLDX(L+1,NX)/DYU(L+1))**2
chammod     &                            +(QSBDLDY(L  ,NX)/DXV(L  ))**2
chammod                  CQBEDLOADX(L,NX)=SQRT(CQBEDLOADX(L,NX))/(HP(L)*VELMAG)
chammod                  CQBEDLOADY(L,NX)=0.0
chammod             ENDIF
chammod             IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).GT.0.0)THEN
chammod                  CQBEDLOADX(L,NX)=(QSBDLDX(L  ,NX)/DYU(L  ))**2
chammod     &                            +(QSBDLDY(LN ,NX)/DXV(LN ))**2
chammod                  CQBEDLOADX(L,NX)=SQRT(CQBEDLOADX(L,NX))/(HP(L)*VELMAG)
chammod                  CQBEDLOADY(L,NX)=0.0
chammod             ENDIF
chammod             IF(UCELLCTR(L).LT.0.0.AND.VCELLCTR(L).LT.0.0)THEN
chammod                  CQBEDLOADX(L,NX)=(QSBDLDX(L  ,NX)/DYU(L  ))**2
chammod     &                            +(QSBDLDY(L  ,NX)/DXV(L  ))**2
chammod                  CQBEDLOADX(L,NX)=SQRT(CQBEDLOADX(L,NX))/(HP(L)*VELMAG)
chammod                  CQBEDLOADY(L,NX)=0.0
chammod             ENDIF
              ENDIF
            ENDIF
          ENDDO
C
C **  INSERT OUTFLOW OR RECIRCULATION BOUNDARY CONDITION FOR BED 
C **  LOAD TRANSPORT
C
          QSBLLDXY=0.0
          IF(NSBDLDBC.GT.0) THEN 
            DO NSB=1,NSBDLDBC
              LUTMP=LSBLBCU(NSB)
              LDTMP=LSBLBCD(NSB)
              IF(LDTMP.EQ.0) THEN
                QSBLLDXY=QSBLLDXY+QSBDLDOT(LUTMP,NX)
              ENDIF
            ENDDO
          ENDIF

C **  CALCULATE MASS PER UNIT AREA CHANGE IN BED CONCENTRATION DUE TO
C **  TO NET BED LOAD
C
          SNDFBLTOT=0.0
          DO L=2,LA
            IF(LMASKDRY(L))THEN
              LN=LNC(L)
              SNDFBL(L,NX)=DXYIP(L)*(QSBDLDX(L+1,NX)
     &                  -QSBDLDX(L,NX)+QSBDLDY(LN,NX)-QSBDLDY(L,NX)
     &                  +QSBDLDOT(L,NX)-QSBDLDIN(L,NX))
              SNDFBLTOT=SNDFBLTOT+DXYP(L)*SNDFBL(L,NX)
            ENDIF
          ENDDO
C
          DO L=2,LA
              SNDFBL(L,NX)=SPB(L)*SNDFBL(L,NX)
          ENDDO
C
c      IF(LA.LT.16)THEN
c      DO L=2,LA
c        TMPVAL=DXYP(L)*SNDFBL(L,NX)
c        WRITE(8,8999)N,L,NX,QSBDLDX(L,NX),QSBDLDX(L+1,NX),
c     &               QSBDLDOT(L,NX),QSBDLDIN(L,NX),TMPVAL
c      ENDDO
c      ENDIF
C
c
c            WRITE(8,8862)N,NX,SNDFBLTOT,QSBLLDXY
c
CXX          CORSNDBL=QSBLLDXY-SNDFBLTOT
CXX          SNDFBLTOT=0.0
CXX          DO L=2,LA
CXX            IF(LMASKDRY(L))THEN
CXX              SNDFBL(L,NX)=SNDFBL(L,NX)+CORSNDBL*DXYIP(L)
CXX              SNDFBLTOT=SNDFBLTOT+DXYP(L)*SNDFBL(L,NX)
CXX            ENDIF
CXX          ENDDO
c
c            WRITE(8,8862)N,NX,SNDFBLTOT,QSBLLDXY
C
c      IF(LA.LT.16)THEN
c      DO L=2,LA
c        TMPVAL=DXYP(L)*SNDFBL(L,NX)
c        WRITE(8,8999)N,L,NX,QSBDLDX(L,NX),QSBDLDX(L+1,NX),
c     &               QSBDLDOT(L,NX),QSBDLDIN(L,NX),TMPVAL
c      ENDDO
c       ENDIF
C
 8999 FORMAT(' BL ',3I5,5E14.5)
 8862 FORMAT(' SNDFBLTOT,QSBLLDXY',3I5,5E14.5)
C
C**********************************************************************C
C
      RETURN
      END
