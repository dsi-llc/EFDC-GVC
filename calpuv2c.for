C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV2C
C
C **  PREVIOUS NAME WAS CALPUV2TC
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
C 02/28/2002        john hamrick       02/28/2002       john hamrick
C  modified drying and wetting scheme. the old formulation remains
C  see (isdry.gt.0.and.isdry.lt.98). the new formulation is activated 
C  by (isdry.eq.99). also added option to waste water from essentially
C  dry cells having water depths greater than hdry.  ie the high and
C  wet cells blocked by dry cells. this is actived by a negative value
C  of ndrystp parameter is the efdc.inp file
C 03/05/2002        john hamrick       03/05/2002       john hamrick
C  added save of old values of horizontal flow face switches sub1 & svb1
C  and transport bypass mask, imaskdry for dry cells. add variable
C  idrydwn to mark wasting from blocked cells
C 06/06/2002        john hamrick       06/06/2002       john hamrick
C  added QDWASTE(L) to save source equivalent of volume loss rate
C  for reducing depth of high/dry cells.  also added concentration 
C  adjustment
C**********************************************************************C
C
C ** SUBROUTINE CALPUV2TC CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOWS WITH PROVISIONS FOR WETTING
C ** AND DRYING OF CELLS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
c
      DIMENSION QCHANUT(NCHANM),QCHANVT(NCHANM)
      DIMENSION QSUMTMP(LCM)
C      DIMENSION XSMLC(LCM-2),XLSMLC(LCM-2),XUSMLC(LCM-2)
C      DIMENSION WSMLC(LCM*LCM)
C      DIMENSION IWSMLC(2*LCM)
C      DIMENSION ASMLC(1,LCM)
C      DIMENSION BSMLC(1)
      DIMENSION IQDRYDWN(LCM)
      DIMENSION IACTIVE(NCHANM)
C
C**********************************************************************C
C
C      WRITE(6,6000)N
C 6000 FORMAT(' CALLED CALPUV9, N = ',I10)
C
      IF(N.EQ.1)THEN
        OPEN(1,FILE='MODCHAN.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(N.EQ.1.AND.ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP1.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(N.EQ.2.AND.ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP2.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      IF(ISDSOLV.EQ.1)THEN
        OPEN(1,FILE='FUV.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FP.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
        DELTD2=0.5*DT
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      ENDIF
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      ISTL=2  
C
      RLAMN=QCHERR
      RLAMO=1.-RLAMN
C
C**********************************************************************C
C
C **  SET SWITCHES FOR DRYING AND WETTING  
C
C----------------------------------------------------------------------C
C     
      ITERHP=0
      NCORDRY=0
      ICORDRY=0
      NEWDRY=0  
      DO L=1,LC
        IQDRYDWN(L)=0
        ISCDRY(L)=0
      ENDDO
      DO L=1,LC
        SUB1(L)=SUB(L)
        SVB1(L)=SVB(L)
      ENDDO
C
C**********************************************************************C
C
C **  INITIALIZE SUBGRID SCALE CHANNEL INTERACTIONS
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          QCHANUT(NMD)=QCHANU(NMD)
          QCHANVT(NMD)=QCHANV(NMD)
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
C
C
      IF(BSC.GT.1.E-6) CALL CALEBI
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL PRESSURE GRADIENTS 
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C **  SNLPX=SNLPX*GID2*DYU & SNLPY=SNLPY*GID2*DXV
C
C----------------------------------------------------------------------C
C
C      IF(BSC.GT.1.E-6)THEN
C        DO L=2,LA
C          FPGXE(L)=-SBX(L)*HU(L)*((BI2(L)+BI2(L-1))*(HP(L)-HP(L-1))
C     &           +2.*HU(L)*(BI1(L)-BI1(L-1))
C     &           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1)))
C        ENDDO
C        DO L=2,LA
C          LS=LSC(L)      
C          FPGYE(L)=-SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
C     &           +2.*HV(L)*(BI1(L)-BI1(LS))
C     &           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS)))
C       ENDDO
C      ENDIF
C
C
      IF(IS2TLPG.EQ.1)THEN
	  ROLDPG=-0.5
	  RNEWPG=1.5
	ELSE
	  ROLDPG=0.0
	  RNEWPG=1.0
	ENDIF
C

      IF(BSC.GT.1.E-6)THEN
        DO L=2,LA
          FPGXE(L)=ROLDPG*FPGXE(L)+RNEWPG*(
     &           -SBX(L)*HU(L)*((BI2(L)+BI2(L-1))*(HP(L)-HP(L-1))
     &           +2.*HU(L)*(BI1(L)-BI1(L-1))
     &           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1)))  )
        ENDDO
        DO L=2,LA
          LS=LSC(L)      
          FPGYE(L)=ROLDPG*FPGYE(L)+RNEWPG*(
     &           -SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
     &           +2.*HV(L)*(BI1(L)-BI1(LS))
     &           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS)))  )
       ENDDO
      ENDIF
C
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      HUTMP(L)=HU(L)
      HVTMP(L)=HV(L)
      H2P(L)=HP(L)
      ENDDO
C
      IF(ISHOUSATONIC.EQ.1)THEN
      IF(HWET.EQ.HDRY)THEN
      DO L=2,LA
        LS=LSC(L)
        IDEPFLGU=0
        IDEPFLGV=0
        IF(IMASKDRY(L).EQ.1)IDEPFLGU=IDEPFLGU+1
        IF(IMASKDRY(L-1).EQ.1)IDEPFLGU=IDEPFLGU+1
        IF(IMASKDRY(L).EQ.1)IDEPFLGV=IDEPFLGV+1
        IF(IMASKDRY(LS).EQ.1)IDEPFLGV=IDEPFLGV+1
        IF(IDEPFLGU.GE.1)THEN
          TMPVAL=MIN(HP(L),HP(L-1))
          HUTMP(L)=MIN(HUTMP(L),TMPVAL)
        ENDIF
        IF(IDEPFLGV.GE.1)THEN
          TMPVAL=MIN(HP(L),HP(LS))
          HVTMP(L)=MIN(HVTMP(L),TMPVAL)
        ENDIF
      ENDDO
      ENDIF
      ENDIF
C
      DO L=2,LA
	  TVAR3S(L)=P(LSC(L))
	ENDDO
C
      DO L=2,LA
        FUHDYE(L)=UHDYE(L)
     &         -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P(L)-P(L-1))
     &         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L))
     &         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
        FVHDXE(L)=VHDXE(L)
     &         -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P(L)-TVAR3S(L))
     &         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L))
     &         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      ENDDO
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='FUV.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='FUV1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='FUV2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HUTMP(L),P1(L),
     &         P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
C
C**********************************************************************C
C
C **  SET IMPLICIT BOTTOM AND VEGETATION DRAG AS APPROPRIATE
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      RCX(L)=1./(1.+FBODYFXI(L))
      RCY(L)=1./(1.+FBODYFYI(L))
      ENDDO
C    
      RCX(1)=0.
      RCY(1)=0.
      RCX(LC)=0.
      RCY(LC)=0.
C
C * SINGLE LAYER NO VEGETATION
C
      IF(KC.EQ.1.AND.ISVEG.EQ.0)THEN
        DO L=2,LA
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * SINGLE LAYER WITH VEGETATION
C
      IF(KC.EQ.1.AND.ISVEG.GE.1)THEN
        DO L=2,LA
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1))
     &    +DELT*FXVEGE(L) )
C     &    +RITB*DELT*FXVEGE(L) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
     &    +DELT*FYVEGE(L) )
C     &    +RITB*DELT*FYVEGE(L) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C
C * MULTIPLE LAYERS NO VEGETATION
C
      IF(KC.GT.1.AND.ISVEG.EQ.0)THEN
        DO L=2,LA
	  TMPX=1.0
	  TMPY=1.0
	  IF(UHE(L).NE.0.0) TMPX=U(L,1)*HU(L)/UHE(L)
	  IF(VHE(L).NE.0.0) TMPY=V(L,1)*HV(L)/VHE(L)
        RCX(L)=1./( 1.+FBODYFXI(L)
     &  +TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &  +TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * MULTIPLE LAYERS WITH VEGETATION
C
      IF(KC.GT.1.AND.ISVEG.GE.1)THEN
        DO L=2,LA
	  TMPX=1.0
	  TMPY=1.0
	  IF(UHE(L).NE.0.0) TMPX=U(L,1)*HU(L)/UHE(L)
	  IF(VHE(L).NE.0.0) TMPY=V(L,1)*HV(L)/VHE(L)
        RCX(L)=1./( 1.+FBODYFXI(L)
     &    +TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1))
     &    +DELT*FXVEGE(L) )
        RCY(L)=1./( 1.+FBODYFYI(L)
     &    +TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
     &    +DELT*FYVEGE(L) )
        FUHDYE(L)=FUHDYE(L)*RCX(L)
        FVHDXE(L)=FVHDXE(L)*RCY(L)
        ENDDO
      ENDIF
C
C * MULTIPLE LAYERS WITH VEGETATION  ORIGINAL
C
CORG      IF(KC.GT.1.AND.ISVEG.GE.1)THEN
CORG        DO L=2,LA
C        RCX(L)=1./( 1.+RITB*DELT*FXVEGE(L) )
C        RCY(L)=1./( 1.+RITB*DELT*FYVEGE(L) )
CORG        RCX(L)=1./( 1.+DELT*FXVEGE(L) )
CORG        RCY(L)=1./( 1.+DELT*FYVEGE(L) )
CORG        FUHDYE(L)=FUHDYE(L)*RCX(L)
CORG        FVHDXE(L)=FVHDXE(L)*RCY(L)
CORG        ENDDO
CORG      ENDIF
C
C**********************************************************************C
C
C **  RESET BOUNDARY CONDITIONS SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
        SUB(L)=SUBO(L)
        SVB(L)=SVBO(L)
        SBX(L)=SBXO(L)
        SBY(L)=SBYO(L)
        SUB(L+1)=SUBO(L+1)
        SBX(L+1)=SBXO(L+1)
      ENDDO
C
      DO L=2,LA
        LN=LNC(L)
        SVB(LN)=SVBO(LN)
        SBY(LN)=SBYO(LN)
      ENDDO
C

      DO L=1,LC
      FP(L)=0.
      FP1(L)=0.
      ENDDO
C
C**********************************************************************C
C
C **  BLOCK DRY CELLS WITH NO POTENTIAL FOR WETTING
C
C----------------------------------------------------------------------C
C
C      DO L=2,LA
C       LN=LNC(L)
C       LS=LSC(L)
C       IF(H1P(L).LT.HWET)THEN
C         IF(SUBO(L).GT.0.5)THEN
C           IF(P(L-1).LT.P(L).OR.HP(L-1).LT.HDRY) SUB(L)=0.0
C           IF(P(L-1).LT.P(L)) SUB(L)=0.0
C         ENDIF
C         IF(SUBO(L+1).GT.0.5)THEN
C           IF(P(L+1).LT.P(L).OR.HP(L+1).LT.HDRY) SUB(L+1)=0.0
C           IF(P(L+1).LT.P(L)) SUB(L+1)=0.0
C         ENDIF
C         IF(SVBO(L).GT.0.5)THEN
C           IF(P(LS).LT.P(L).OR.HP(LS).LT.HDRY) SVB(L)=0.0
C           IF(P(LS).LT.P(L)) SVB(L)=0.0
C         ENDIF
C         IF(SVBO(LN).GT.0.5)THEN
C           IF(P(LN).LT.P(L).OR.HP(LN).LT.HDRY) SVB(LN)=0.0
C           IF(P(LN).LT.P(L)) SVB(LN)=0.0
C         ENDIF
C       ENDIF
C      ENDDO
C
c      DO L=2,LA
c       LN=LNC(L)
c       LS=LSC(L)
c	 ITMPVAL=0
c       IF(IMASKDRY(L).GE.1)THEN
c         IF(SUBO(L).LT.0.5.OR.IMASKDRY(L-1).GE.1)
c     &     ITMPVAL=ITMPVAL+1
c         IF(SUBO(L+1).LT.0.5.OR.IMASKDRY(L+1).GE.1)
c     &     ITMPVAL=ITMPVAL+1
c         IF(SVBO(L).LT.0.5.OR.IMASKDRY(LS).GE.1)
c     &     ITMPVAL=ITMPVAL+1
c         IF(SVBO(LN).LT.0.5.OR.IMASKDRY(LN).GE.1)
c     &     ITMPVAL=ITMPVAL+1
c	   IF(ITMPVAL.EQ.4.AND.QSUME(L).LE.0.0)THEN
c	     SUB(L)=0.
c           SVB(L)=0.
c           SBX(L)=0.
c           SBY(L)=0.
c           SUB(L+1)=0.
c           SVB(LN)=0.
c           SBX(L+1)=0.
c           SBY(LN)=0.
c	   ENDIF
c       ENDIF
c      ENDDO
C
C**********************************************************************C
C
C **  SET OPEN BOUNDARY SURFACE ELEVATIONS 
C
      IVAL=NPBW+NPBE+NPBS+NPBN
	IF(IVAL.GT.0) CALL SETOBC2T(DELT,DELTD2,DELTI)
C
C**********************************************************************C
C
C **  ADJUST VOLUME SOURCE AND SINKS
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.EQ.0)THEN
C
      DO L=2,LA
      IF(QSUME(L).LE.0.)THEN
        IF(H1P(L).LE.HDRY)THEN
          QSUMTMP(L)=0.
         ELSE
          QSUMTMP(L)=-(H1P(L)-HDRY)*DXYP(L)*DELTI
          QSUMTMP(L)=MAX(QSUMTMP(L),QSUME(L))
        ENDIF
       ELSE
        QSUMTMP(L)=QSUME(L)
      ENDIF
      ENDDO
C 
      DO L=2,LA
       DIFQVOL=QSUME(L)-QSUMTMP(L)
       DO K=1,KC
       QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
       ENDDO
       QSUME(L)=QSUMTMP(L)
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C 
C **  ADJUST SOURCES AND SINKS ESTIMATING SURFACE AND GROUNDWATER
C **  AVAILABLE FOR EVAPOTRANSPIRATON AND INFILTRATION
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.GE.1)THEN
C
      DO L=2,LA
      RIFTR(L)=0.
      EVAPSW(L)=0.
      EVAPGW(L)=0.
      IF(H1P(L).GT.HDRY)THEN
C       APPLY MAXIMUM ET
        IF(EVAPCVT.LT.0.)THEN 
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/     
     &              (1.+0.00412*TEM(L,KC))))            
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)       
        ENDIF                                            
        EVAPSW(L)=EVAPT(L)*DXYP(L)
        RIFTR(L)=0.
C       CALCULATE DEPTH OF ACTIVE GROUNDWATER ELEV BELOW SURFACE
        DTAGW=BELV(L)-AGWELV(L)
        IF(DTAGW.GT.0.0)THEN
C         INFLITRATION CAN OCCUR, CALCULATE LIMITING RATE TO BRING
C         GW ELEV TO SOIL SURFACE
          RIFTRL=RNPOR*DTAGW*DELTI
C         SET RIFTRL TO MIN OF LIMITING RATE OR ACTUAL RATE
          RIFTRL=MIN(RIFTRM,RIFTRL)
C         ESTIMATE RATE BASED ON AVAILABLE SURFACE WATER 
          RAVAIL=(H1P(L)-HDRY)*DELTI-EVAPT(L)
C         SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
          RIFTRL=MIN(RAVAIL,RIFTRL)
C         CONVERT TO VOLUME FLOW UNITS
          RIFTR(L)=RIFTRL*DXYP(L)         
        ENDIF
C       ADJUST VOLUME OUTFLOWS OF WET CELLS
        IF(QSUME(L).LT.0.0)THEN
          QSUMIET=RIFTR(L)+EVAPSW(L)
          QEAVAIL=DXYP(L)*(H1P(L)-HDRY)*DELTI-QSUMIET
          QEAVAIL=MAX(QEAVAIL,0.0)
          QEAVAIL=-QEAVAIL
          QSUMTMP(L)=MAX(QSUME(L),QEAVAIL)
         ELSE
          QSUMTMP(L)=QSUME(L)
        ENDIF         
       ELSE
        RIFTR(L)=0.
        EVAPSW(L)=0.
        QSUMTMP(L)=MAX(QSUME(L),0.0)       
      ENDIF
      ENDDO
C
      DO L=2,LA
      DIFQVOL=QSUME(L)-QSUMTMP(L)
      DO K=1,KC
      QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
      ENDDO
      QSUME(L)=QSUMTMP(L)
      ENDDO
C       
      ENDIF
C
C**********************************************************************C
C
C **  ADVANCE EXTERNAL VARIABLES
C
C----------------------------------------------------------------------C
C
        DO L=2,LA
        UHDY2E(L)=UHDY1E(L)
        VHDX2E(L)=VHDX1E(L)
        UHDY1E(L)=UHDYE(L)
        VHDX1E(L)=VHDXE(L)
        U1V(L)=UV(L)
        V1U(L)=VU(L)
        P1(L)=P(L)
        H1U(L)=HU(L)
        H1V(L)=HV(L)
        H1UI(L)=HUI(L)
        H1VI(L)=HVI(L)
        H2P(L)=H1P(L)
        H1P(L)=HP(L)
        AGWELV2(L)=AGWELV1(L)
        AGWELV1(L)=AGWELV(L)
        ENDDO
C
C**********************************************************************C
C
C **  SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR 
C **  NON BOUNDARY POINTS
C **  HRU=HMU*DYU/DXU & HRV=HMV*DXV/DYV 
C **  DXYIP=1/(DXP*DYP)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
        TVAR3N(L)=VHDXE(LNC(L))
	ENDDO
      DO L=2,LA
      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*P(L)
     &      -0.5*G*(UHDYE(L+1)-UHDYE(L)
     &             +TVAR3N(L)-VHDXE(L)) )
      ENDDO
C
C**********************************************************************C
C
C **  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
C **  HOST-GUEST CHANNAL INTERACTION FOR NON BOUNDARY POINTS 
C
C **  REENTER AT 1000 FOR WETTING-DRYING CORRECTION AND CHANNEL 
C **  INTERACTION
C
C----------------------------------------------------------------------C
C
 1000 CONTINUE
C
C
      DO L=2,LA
        LN=LNC(L)
        TVAR3N(L)=SVB(LN )*FVHDXE(LN)
	ENDDO
      DO L=2,LA
      FP(L)=FP1(L)-0.5*G*SPB(L)*
     &      ( SUB(L+1)*FUHDYE(L+1)-SUB(L)*FUHDYE(L)
     &       +TVAR3N(L)-SVB(L)*FVHDXE(L)
     &       -2.0*QSUME(L) )
CC      P(L)=0.
      ENDDO
C
      IF(ISGWIE.GE.1)THEN
        DO L=2,LA
        FP(L)=FP(L)-G*SPB(L)*(RIFTR(L)+EVAPSW(L))
        ENDDO
      ENDIF
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='FP.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='FP1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='FP2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),FP1(L),FUHDYE(L),FUHDYE(L+1),
     &          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      C1=-0.5*DELTD2*G*SPB(L)
      CS(L)=C1*SVB(L  )*HRVO(L  )*RCY(L  )*HVTMP(L  )
C    &     +(1.-SPB(L))*CS(L)
      CW(L)=C1*SUB(L  )*HRUO(L  )*RCX(L  )*HUTMP(L  )
C    &     +(1.-SPB(L))*CW(L)
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
C    &     +(1.-SPB(L))*CE(L)
      ENDIF
      ENDDO
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      LN=LNC(L)
      C1=-0.5*DELTD2*G*SPB(L)
      CN(L)=C1*SVB(LN )*HRVO(LN )*RCY(LN )*HVTMP(LN )
C    &     +(1.-SPB(L))*CN(L)
      ENDIF
      ENDDO
C
      DO L=2,LA
      IF(SPB(L).GT.0.)THEN
      CC(L)=SPB(L)*(DELTI*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
C    &     +(1.-SPB(L))*CC(L)
      ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  INSERT IMPLICT SUB-GRID SCALE CHANNEL INTERACTIONS
C
      IF(MDCHH.GE.1)CALL SUBCHAN(QCHANUT,QCHANVT,IACTIVE,RLAMN,RLAMO,
     &    DELT,IACTALL)
C
C**********************************************************************C
C
C **  SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
C
C----------------------------------------------------------------------C
C
      CCMNM=1.E+18
      DO L=2,LA
        CCMNM=MIN(CCMNM,CC(L))
        FPTMP(L)=FP(L)
      ENDDO
C
      CCMNMI=1./CCMNM
C     WRITE(6,666) CCMNM,CCMNMI
C 666 FORMAT(' CCMIN, CCMINI = ',E12.4,1X,E12.4)
C
      DO LL=1,NPBW
      IF(ISPBW(LL).EQ.0)THEN
        L=LPBW(LL)
        CW(L+1)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBE
      IF(ISPBE(LL).EQ.0)THEN
        L=LPBE(LL)
        CE(L-1)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBS
      IF(ISPBS(LL).EQ.0)THEN
        L=LPBS(LL)
        LN=LNC(L)
        CS(LN)=0.
      ENDIF
      ENDDO
C
      DO LL=1,NPBN
      IF(ISPBN(LL).EQ.0)THEN
        L=LPBN(LL)
        LS=LSC(L)
        CN(LS)=0.
      ENDIF
      ENDDO
C
      CC(1)=1.
      CC(LC)=1.
C
C----------------------------------------------------------------------C
C
C **  SCALE BY MINIMUM DIAGONAL
C
      IF(IRVEC.EQ.9)THEN
C
      DO L=2,LA
      CCS(L)=CS(L)*CCMNMI
      CCW(L)=CW(L)*CCMNMI
      CCE(L)=CE(L)*CCMNMI
      CCN(L)=CN(L)*CCMNMI
      CCC(L)=CC(L)*CCMNMI
      FPTMP(L)=FPTMP(L)*CCMNMI
      CCCI(L)=1./CCC(L)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          CCCCHH(NMD)=CCCCHH(NMD)*CCMNMI
        ENDDO
      ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  SCALE TO NORMAL FORM - NOTE NORMAL FORM CURRENTLY CANNOT BE USED
C **  WITH SUB-GRID SCALE CHANNEL MODEL
C
C      IF(IRVEC.EQ.99)THEN
C
C      DO L=2,LA
C        CCS(L)=CS(L)/SQRT( CC(L)*CC(LSC(L)) )
C        CCW(L)=CW(L)/SQRT( CC(L)*CC(L-1   ) )
C        CCE(L)=CE(L)/SQRT( CC(L)*CC(L+1   ) )
C        CCN(L)=CN(L)/SQRT( CC(L)*CC(LNC(L)) )
C        CCC(L)=1.
C        FPTMP(L)=FPTMP(L)/SQRT( CC(L) )
C        P(L)=P(L)*SQRT( CC(L) )
C        CCCI(L)=1.
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
C **  CALL EQUATION SOLVER
C
      IF(MDCHH.EQ.0) CALL CONGRAD (ISTL)
      IF(MDCHH.GE.1) CALL CONGRADC (ISTL)
C
C----------------------------------------------------------------------C
C
C BEGIN CONSTRAINED MINIMIZATION SOLUTION
C
C      DO L=2,LA
C       XSMLC(L-1)=P(L)
CJH       XLSMLC(L-1)=G*(BELV(L)+HDRY)
C       XLSMLC(L-1)=G*BELV(L)
C       XUSMLC(L-1)=1.E4
C      ENDDO
C
C      IF(IRVEC.EQ.99)THEN
C        DO L=2,LA
C          XLSMLC(L-1)=XLSMLC(L-1)*SQRT( CC(L) )
C          XUSMLC(L-1)=XUSMLC(L-1)*SQRT( CC(L) )
C        ENDDO
C      ENDIF
C
C      DO L=2,LA
C       XSMLC(L-1)=MAX(XSMLC(L-1),XLSMLC(L-1))
C      ENDDO
C
C      NTMP=LC-2
C      MTMP=0
C      MEQ=0
CJH     ASMLC=ASMLC(1,1)
C      LDA=1
CJH     BSMLC=BSMLC(1)
C      IPRINT=0
C      NFMAX=ITERM
CJH     IWSMLC=IWSM(LIW)
C      LIWSM=4+MTMP+2*NTMP
C      LWSM=3+MTMP+NTMP*(NTMP+16)
CJH     WSMLC=IWSM(LIW)
C      CALL SMLC01(NTMP,MTMP,MEQ,ASMLC,LDA,BSMLC,XLSMLC,
C     &     XUSMLC,XSMLC,RSQM,IPRINT,NFMAX,IWSMLC,LIWSM,WSMLC,LWSM)
C      DO L=2,LA
C       P(L)=XSMLC(L-1)
C      ENDDO
C
C END CONSTRAINED MINIMIZATION SOLUTION
C
C----------------------------------------------------------------------C
C
C **  REVERSE SCALE IF NORMAL FORM SOLUTION IS USED
C
C      IF(IRVEC.EQ.99)THEN
C
C      DO L=2,LA
C      P(L)=P(L)/SQRT( CC(L) )
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
C ** DIAGNOSTICS
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
C        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
        WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
        ENDDO
        IF(MDCHH.GE.1)THEN
          DO NMD=1,MDCHH
            WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                        CCCCHV(NMD),QCHANUT(NMD),QCHANVT(NMD)
          ENDDO
        ENDIF
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='EQCOEF1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          SURFTMP=GI*P(L)
C          WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
          WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
          ENDDO
          IF(MDCHH.GE.1)THEN
            DO NMD=1,MDCHH
             WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                        CCCCHV(NMD),QCHANUT(NMD),QCHANVT(NMD)
            ENDDO
          ENDIF
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='EQCOEF2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          SURFTMP=GI*P(L)
C          WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
C     &               FPTMP(L),SURFTMP
          WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     &               FP(L),SURFTMP
          ENDDO
          IF(MDCHH.GE.1)THEN
            DO NMD=1,MDCHH
             WRITE(1,1001)NMD,MDCHTYP(NMD),CCCCHH(NMD),CCCCHU(NMD),
     &                        CCCCHV(NMD),QCHANUT(NMD),QCHANVT(NMD)
            ENDDO
          ENDIF
          CLOSE(1)
        ENDIF
      ENDIF
C
      IF(ISDSOLV.GE.1)THEN
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
        IF(N.EQ.1)THEN
          OPEN(1,FILE='EQTERM1.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
          ENDDO
          CLOSE(1)
        ENDIF
        IF(N.EQ.2)THEN
          OPEN(1,FILE='EQTERM2.OUT',POSITION='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF
C
 1001 FORMAT(2I5,10(1X,E12.4))
 1002 FORMAT(3I4,10(1X,E9.2))
C
C**********************************************************************C
C
C **  CALCULATE UHEX AND VHEX AND TOTAL DEPTHS AT TIME LEVEL (N+1)
C **  HRU=SUB*DYU/DXU & HRV=SVB*DXV/DYV 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
	  TVAR3S(L)=P(LSC(L))
	ENDDO
      DO L=2,LA
        UHDYE(L)=SUB(L)*( FUHDYE(L)
     &            -DELTD2*HRUO(L)*RCX(L)*HUTMP(L)*(P(L)-P(L-1)) )
        VHDXE(L)=SVB(L)*( FVHDXE(L)
     &            -DELTD2*HRVO(L)*RCY(L)*HVTMP(L)*(P(L)-TVAR3S(L)) )
      ENDDO
      DO L=2,LA
        UHE(L)=UHDYE(L)*DYIU(L)
        VHE(L)=VHDXE(L)*DXIV(L)
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE NEW SUB-GRID SCALE CHANNEL EXCHANGE FLOWS
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          IF (IACTIVE(NMD).GT.0)THEN 
            LHOST=LMDCHH(NMD)
            LCHNU=LMDCHU(NMD)
            LCHNV=LMDCHV(NMD)
            IF(MDCHTYP(NMD).EQ.1)THEN
C             QCHANU(NMD)=0.
C             IF(IACTIVE(NMD).EQ.1)THEN
              QCHANU(NMD)=CCCCHU(NMD)*QCHANUT(NMD)
     &             -RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNU))
     &             -RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNU))
C             ENDIF
              QCHANUN(NMD)=QCHANUT(NMD)
              QCHANV(NMD)=0.
              QCHANVN(NMD)=QCHANVT(NMD)
            ENDIF
            IF(MDCHTYP(NMD).EQ.2)THEN
C             QCHANV(NMD)=0.
C             IF(IACTIVE(NMD).EQ.1)THEN
              QCHANV(NMD)=CCCCHU(NMD)*QCHANVT(NMD)
     &             -RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNV))
     &             -RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNV))
C             ENDIF
              QCHANVN(NMD)=QCHANVT(NMD)
              QCHANU(NMD)=0.
              QCHANUN(NMD)=QCHANUT(NMD)
            ENDIF
          ELSE
            QCHANV(NMD)=0.
            QCHANVN(NMD)=0.  
            QCHANU(NMD)=0.  
            QCHANUN(NMD)=0.  
          ENDIF
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL 
C **  TRANSPORTS AT (N+1)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
        LN=LNC(L)
	  TVAR3N(L)=VHDXE(LN)+VHDX1E(LN)
	ENDDO
      DO L=2,LA
        TVAR3C(L)=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     &       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &       +TVAR3N(L)-VHDXE(L)-VHDX1E(L)))
      ENDDO
      IF(ISGWIE.GE.1)THEN
	  DO L=2,LA
          TVAR3C(L)=TVAR3C(L)-DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        ENDDO
      ENDIF
	DO L=2,LA
	  HP(L)=SPB(L)*TVAR3C(L)+(1.-SPB(L))*(GI*P(L)-BELV(L))
      ENDDO
C
C**********************************************************************C
C
C **  ADD CHANNEL INTERACTION EXCHANGES
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
          IF(IACTIVE(NMD).GT.0)THEN  
            LHOST=LMDCHH(NMD)
            LCHNU=LMDCHU(NMD)
            LCHNV=LMDCHV(NMD)
            IF(MDCHTYP(NMD).EQ.1)THEN
C             IF(IACTIVE(NMD).EQ.1)THEN
              TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUT(NMD))
              HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
              HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
C             ENDIF            
            ENDIF            
            IF(MDCHTYP(NMD).EQ.2)THEN
C             IF(IACTIVE(NMD).EQ.1)THEN
              TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVT(NMD))
              HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
              HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
C             ENDIF            
            ENDIF            
          ENDIF            
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      ENDDO
C
C
      DO L=2,LA
	  TVAR3C(L)=DXP(LSC(L) )*DYP(LSC(L) )
	  TVAR3S(L)=HP(LSC(L) )
	ENDDO
C
      IF(ISHOUSATONIC.EQ.0)THEN
      DO L=2,LA
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+TVAR3C(L)*TVAR3S(L))
     &           /(DXV(L)*DYV(L))
      ENDDO
	ENDIF
C
      IF(ISHOUSATONIC.EQ.1)THEN
      DO L=2,LA
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
cjah 
       IF(ABS(HP(L)).LT.HDRY/1000.) THEN
         HP(L)=HDRY/1000.
C        write(6,*) ' *** WARNING CELL RESET TO HYDR/100:',L,'***'
       ENDIF
cjah   IF(HP(L-1).LE.0.or.HP(L).LE.0) then
cjah      write(6,*) ' ERROR L, HP(L),L-1=',HP(L),HP(L-1)
cjajah    stop
cjah   endif
cjah end
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       IF(ABS(HU(L)).LT.HDRY/1000.) THEN
       HU(L)=HDRY/1000.
C        write(6,*) ' *** WARNING CELL RESET TO HYDR/100:',L,'***'
       ENDIF
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+TVAR3C(L)*TVAR3S(L))
     &           /(DXV(L)*DYV(L))
       IF(ABS(HV(L)).LT.HDRY/1000.) THEN
       HV(L)=HDRY/1000.
C        write(6,*) ' *** WARNING CELL RESET TO HYDR/100:',L,'***'
       ENDIF
      ENDDO
	ENDIF
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      ENDDO
C
C
C**********************************************************************C
C
C **  CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
C     OLD PRE 28 FEB 02 VERSION
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.GT.0.AND.ISDRY.LT.98)THEN
      OPEN(1,FILE='DRYWET.LOG',POSITION='APPEND',STATUS='UNKNOWN')
C
      ICORDRY=0
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
C
      IF(HP(L).LE.HDRY)THEN
        IF(ISCDRY(L).EQ.0)THEN
          ISCDRY(L)=1
          ICORDRY=1
C         WRITE(1,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(6,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(8,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
        ENDIF
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
      ENDIF
C
      ENDDO
C
      IF(ICORDRY.EQ.1)THEN
        NCORDRY=NCORDRY+1
        GOTO 1000
      ENDIF
C
      CLOSE(1)
      ENDIF
C
C      WRITE(6,6960)NCORDRY
CTMP      WRITE(8,6960)NCORDRY
C
 6960 FORMAT(' NCORDRY =', I5)
 6961 FORMAT(' UNSTABLE, NCORDRY =', I5)
C
 9999 CONTINUE
C
C**********************************************************************C
C
C **  CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
C     NEW VERSION 28 FEB 02  ISDRY=99
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.EQ.99)THEN
      OPEN(1,FILE='DRYWET.LOG',POSITION='APPEND')
C
      HDRY2=2.*HDRY
      ICORDRY=0
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
C
      IF(HP(L).LE.HDRY)THEN
	  SUBW=SUB(L)
	  SUBE=SUB(L+1)
	  SVBS=SVB(L)
	  SVBN=SVB(LN)
	  DHPDT=DELTI*(HP(L)-H1P(L))
	  IF(DHPDT.GT.0.0)THEN
CZZ
          SUB(L)=0.0  
          SUB(L+1)=0.0  
          SVB(L)=0.0  
          SVB(LN)=0.0  
          SBX(L)=0.0  
          SBX(L+1)=0.0  
          SBY(L)=0.0  
          SBY(LN)=0.0  
CZZ
	    IF(SUBO(L).GT.0.5)THEN
	      IF(UHDYE(L).GT.0.0.AND.HP(L-1).GT.HDRY2)THEN
	        SUB(L)=1.
              SBX(L)=1.
            ENDIF
          ENDIF
	    IF(SUBO(L+1).GT.0.5)THEN
	      IF(UHDYE(L+1).LT.0.0.AND.HP(L+1).GT.HDRY2)THEN
	        SUB(L+1)=1.
              SBX(L+1)=1.
            ENDIF
          ENDIF
	    IF(SVBO(L).GT.0.5)THEN
	      IF(VHDXE(L).GT.0.0.AND.HP(LS).GT.HDRY2)THEN
	        SVB(L)=1.
              SBY(L)=1.
            ENDIF
          ENDIF
	    IF(SVBO(LN).GT.0.5)THEN
	      IF(VHDXE(LN).LT.0.0.AND.HP(LN).GT.HDRY2)THEN
	        SVB(LN)=1.
              SBY(LN)=1.
            ENDIF
          ENDIF
	    RDRY=SUB(L)+SUB(L+1)+SVB(L)+SVB(LN)
	    IF(RDRY.LT.0.5)THEN
            ISCDRY(L)=1
          ELSE
            ISCDRY(L)=0
          ENDIF
	    TMPVAL=ABS(SUB(L)-SUBW)
	    IF(TMPVAL.GT.0.5)ICORDRY=1
	    TMPVAL=ABS(SUB(L+1)-SUBE)
	    IF(TMPVAL.GT.0.5)ICORDRY=1
	    TMPVAL=ABS(SVB(L)-SVBS)
	    IF(TMPVAL.GT.0.5)ICORDRY=1
	    TMPVAL=ABS(SVB(LN)-SVBN)
	    IF(TMPVAL.GT.0.5)ICORDRY=1
	  ELSE
	    SUB(L)=0.0
	    SUB(L+1)=0.0
	    SVB(L)=0.0
	    SVB(LN)=0.0
	    SBX(L)=0.0
	    SBX(L+1)=0.0
	    SBY(L)=0.0
	    SBY(LN)=0.0
          IF(ISCDRY(L).EQ.0)THEN
		  ISCDRY(L)=1
            ICORDRY=1
          ENDIF
	  ENDIF
	ENDIF
C
      ENDDO
C
      IF(ICORDRY.EQ.1)THEN
        NCORDRY=NCORDRY+1
        GOTO 1000
      ENDIF
C
      CLOSE(1)
      ENDIF
C
C     WRITE(6,6960)NCORDRY
      WRITE(8,6960)NCORDRY
C
C**********************************************************************C
C
C **  COUNT THE NUMBER TO TIME STEPS A CELL IS DRY, AND IF IT HAS BEEN
C **  DRY FOR MORE THAN ABS(NDRYSTP), AND ITS BOTTOM ELEVATION IS HIGHER
C **  THAN THE SURROUNDING DRY CELLS, THEN REDUCE ITS DEPTH BELOW THE 
C **  DRYING DEPTH IF NECESSARY.  SAVE VOLUME REDUCTION RATE AS QDWASTE
C **  DEFINED AS POSITIVE OUT. THEN ADJUST CONCENTRATIONS
C
      IF(ISDRY.GT.0) THEN
        IF(NDRYSTP.LT.0) THEN
	    NTMP=ABS(NDRYSTP)
	    DO L=2,LA
	      LN=LNC(L)
	      LS=LSC(L)
            QDWASTE(L)=0.
		  IQDRYDWN(L)=0
	      RDRY=SUB(L)+SUB(L+1)+SVB(L)+SVB(LN)
	      IF(RDRY.GT.0.5)NATDRY(L)=0
            IF(RDRY.LT.0.5)NATDRY(L)=NATDRY(L)+1
	      IF(NATDRY(L).GT.NTMP)THEN	
	        IF(HP(L).GE.HDRY)THEN
	          BELVAVG=0.0
	          RVAL=0.0
	          IF(HP(L+1).LT.HDRY.AND.SUBO(L+1).GT.0.5)THEN
	            BELVAVG=BELVAVG+BELV(L+1)
                  RVAL=RVAL+1.
	          ENDIF
	          IF(HP(L-1).LT.HDRY.AND.SUBO(L).GT.0.5)THEN
	            BELVAVG=BELVAVG+BELV(L-1)
                  RVAL=RVAL+1.
	          ENDIF
	          IF(HP(LN).LT.HDRY.AND.SVBO(LN).GT.0.5)THEN
	            BELVAVG=BELVAVG+BELV(LN)
                  RVAL=RVAL+1.
	          ENDIF
	          IF(HP(LS).LT.HDRY.AND.SVBO(L).GT.0.5)THEN
	            BELVAVG=BELVAVG+BELV(LS)
                  RVAL=RVAL+1.
	          ENDIF
	          IF(BELV(L).GT.BELVAVG)THEN
	            HOLDTMP=HP(L)
			    IQDRYDWN(L)=1
	            HP(L)=0.90*HDRY
	            NATDRY(L)=0
	            QDWASTE(L)=DELTI*DXYP(L)*(HOLDTMP-HP(L))
	            VDWASTE(L)=VDWASTE(L)+DXYP(L)*(HOLDTMP-HP(L))
	            TMPVAL=HOLDTMP/HP(L)
	            DO K=1,KC
	              SAL(L,K)=TMPVAL*SAL(L,K)
	              TEM(L,K)=TMPVAL*SAL(L,K)
	              DYE(L,K)=TMPVAL*SAL(L,K)
	              DO NT=1,NTOX
	                TOX(L,K,NT)=TMPVAL*TOX(L,K,NT)
	              ENDDO
	              DO NS=1,NSED
	                SED(L,K,NS)=TMPVAL*SED(L,K,NS)
	              ENDDO
	              DO NX=1,NSND
	                SND(L,K,NX)=TMPVAL*SND(L,K,NX)
	              ENDDO
	            ENDDO
	          ENDIF
              END IF
	      ENDIF
	      IF(QDWASTE(L).GT.0.0)THEN
	        TMPVAL=QDWASTE(L)/DXYP(L)
	        WRITE(8,8888)IL(L),JL(L),TIME,RDRY,HOLDTMP,HP(L),
     &             QDWASTE(L),TMPVAL
	        WRITE(6,8888)IL(L),JL(L),TIME,RDRY,HOLDTMP,HP(L),
     &             QDWASTE(L),TMPVAL
	      ENDIF
          ENDDO
	  END IF
	END IF	      
C
 8888 FORMAT(' QDW ',2I6,6E14.6)
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      ENDDO
C
      DO L=2,LA
	  TVAR3C(L)=DXP(LSC(L) )*DYP(LSC(L) )
	  TVAR3S(L)=HP(LSC(L) )
	ENDDO

      DO L=2,LA
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+TVAR3C(L)*TVAR3S(L))
     &           /(DXV(L)*DYV(L))
      ENDDO
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      ENDDO
C
C**********************************************************************C
C
C **  SET TRANSPORT MASK FOR DRY CELLS
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.GT.0)THEN
C
      DO L=2,LA
        IMASKDRY(L)=0
        LMASKDRY(L)=.TRUE.
      END DO
C
      IF(IDRYTBP.EQ.1)THEN
        DO L=2,LA
	    LN=LNC(L)
          IUW=0
          IUE=0
          IVS=0
          IVN=0
	      IF(SUB1(L).LT.0.5.AND.SUB(L).LT.0.5)IUE=1
	      IF(SUB1(L+1).LT.0.5.AND.SUB(L+1).LT.0.5)IUW=1
	      IF(SVB1(L).LT.0.5.AND.SVB(L).LT.0.5)IVS=1
	      IF(SVB1(LN).LT.0.5.AND.SVB(LN).LT.0.5)IVN=1
	      IFACE=IUW+IUE+IVS+IVN
	      IF(IFACE.EQ.4)THEN
	        IMASKDRY(L)=1
		     LMASKDRY(L)=.FALSE.
	        IF(H1P(L).EQ.HP(L))IMASKDRY(L)=2
            END IF
		  IF(IQDRYDWN(L).EQ.1)THEN
		    IMASKDRY(L)=0
		    LMASKDRY(L)=.TRUE.
            ENDIF
        END DO
      END IF
C
      END IF
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.GT.0)THEN
C
C     COUNT NUMBER OF WET CELLS
C
      NWETCELLSNEW=0
	DO L=2,LA
	IF(LMASKDRY(L))NWETCELLSNEW=NWETCELLSNEW+1
	ENDDO
C
      IF(NWETCELLSNEW.NE.NWETCELLSOLD)THEN
	  OPEN(1,FILE='WETDRYCHG.OUT',POSITION='APPEND')
	  IF(ISDYNSTP.EQ.0)THEN
          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIME=TIMESEC/86400.
        ENDIF
	  WRITE(1,1234)TIME,NWETCELLSNEW,NWETCELLSOLD
	ENDIF
C
      NWETCELLSOLD=NWETCELLSNEW
C
      ENDIF
C
 1234 FORMAT(F12.4,2I8)
C
C**********************************************************************C
C
C **  OUTPUT DIAGNOSTICS FOR 2 GRID INTERATCTION
C
        IF(MDCHH.GT.0)THEN
        IF(MDCHHD.GT.0)THEN
        IVAL=MOD(N,MDCHHD2)
        IF(IVAL.EQ.0)THEN
        IF(IACTALL.GT.0)THEN
        OPEN(1,FILE='MODCHAN.OUT',POSITION='APPEND')
C
        DO NMD=1,MDCHH
          WRITE(1,8000)
          LHOST=LMDCHH(NMD)
          IHOST=IL(LHOST)        
          JHOST=JL(LHOST)        
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
          ICHNU=IL(LCHNU)
          JCHNU=JL(LCHNU)
          SRFCHAN=HP(LCHNU)+BELV(LCHNU)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNU)+BELV(LCHNU)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(1,8001)N,NMD,MDCHTYP(NMD),ICHNU,JCHNU,ISCDRY(LCHNU),
     &                 SRFCHAN,HP(LCHNU),SRFCHAN1,H1P(LCHNU)
          WRITE(1,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),SRFHOST1,H1P(LHOST)
          WRITE(1,8003)QCHANU(NMD),QCHANUT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
          ICHNV=IL(LCHNV)
          JCHNV=JL(LCHNV)
          SRFCHAN=HP(LCHNV)+BELV(LCHNV)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNV)+BELV(LCHNV)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(1,8001)N,NMD,MDCHTYP(NMD),ICHNV,JCHNV,ISCDRY(LCHNV),
     &                 SRFCHAN,HP(LCHNV),SRFCHAN1,H1P(LCHNV)
          WRITE(1,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),SRFHOST1,H1P(LHOST)
          WRITE(1,8003)QCHANV(NMD),QCHANVT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
          WRITE(1,8004)
        ENDDO
C
        CLOSE(1)
        ENDIF
        ENDIF
        ENDIF
        ENDIF
C
C**********************************************************************C
C
C **  PERFORM UPDATE ON GROUNDWATER ELEVATION
C
C----------------------------------------------------------------------C
C
      IF(ISGWIE.GE.1)THEN 
C
        DO L=2,LA
        QSUM(L,KC)=QSUM(L,KC)-EVAPSW(L)
        QSUM(L,1 )=QSUM(L,1 )-RIFTR(L)
        ENDDO
C
C       INFILTRATION STEP
C
        RNPORI=1./RNPOR
        IF(ISTL.EQ.3)THEN
          DO L=2,LA
          AGWELV(L)=AGWELV2(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          ENDDO
         ELSE
          DO L=2,LA
          AGWELV(L)=AGWELV1(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          ENDDO
        ENDIF
        DO L=2,LA
        AGWELV(L)=MIN(AGWELV(L),BELV(L))
        ENDDO
C
C       ET STEP
C
        DO L=2,LA
        IF(EVAPCVT.LT.0.)THEN                    
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/ 
     &              (1.+0.00412*TEM(L,KC))))      
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)       
        ENDIF                                    
        ETGWTMP=EVAPT(L)-EVAPSW(L)*DXYIP(L)
        ETGWTMP=MAX(ETGWTMP,0.0)
        ETGWAVL=RNPOR*DELTI*(AGWELV(L)-BELAGW(L))
        ETGWAVL=MAX(ETGWAVL,0.0)
        ETGWTMP=MIN(ETGWTMP,ETGWAVL)
        EVAPGW(L)=ETGWTMP*DXYP(L)
        ENDDO
        DO L=2,LA
        AGWELV(L)=AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)
        ENDDO
        DO L=2,LA
        AGWELV(L)=MAX(AGWELV(L),BELAGW(L))
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
       IF(N.EQ.NTS)THEN
C
        IF(MDCHH.GT.0)THEN
        DO NMD=1,MDCHH
          WRITE(8,8000)
          LHOST=LMDCHH(NMD)
          IHOST=IL(LHOST)        
          JHOST=JL(LHOST)        
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
          ICHNU=IL(LCHNU)
          JCHNU=JL(LCHNU)
          SRFCHAN=HP(LCHNU)+BELV(LCHNU)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNU)+BELV(LCHNU)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(8,8001)N,NMD,MDCHTYP(NMD),ICHNU,JCHNU,ISCDRY(LCHNU),
     &                 SRFCHAN,HP(LCHNU),P1(LCHNU),H1P(LCHNU)
          WRITE(8,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),P1(LHOST),H1P(LHOST)
          WRITE(8,8003)QCHANU(NMD),QCHANUT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
          ICHNV=IL(LCHNV)
          JCHNV=JL(LCHNV)
          SRFCHAN=HP(LCHNV)+BELV(LCHNV)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNV)+BELV(LCHNV)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(8,8001)N,NMD,MDCHTYP(NMD),ICHNV,JCHNV,ISCDRY(LCHNV),
     &                 SRFCHAN,HP(LCHNV),SRFCHAN1,H1P(LCHNV)
          WRITE(8,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),SRFHOST1,H1P(LHOST)
          WRITE(8,8003)QCHANV(NMD),QCHANVT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
          WRITE(8,8004)
        ENDDO
        ENDIF
C
       ENDIF
C
C**********************************************************************C
C
C **  CHECK FOR NEGATIVE DEPTHS
C
      IF(ISNEGH.GE.1.AND.ISHOUSATONIC.EQ.0)
     &  CALL NEGDEP(QCHANUT,QCHANVT,2)
      IF(ISNEGH.GE.1.AND.ISHOUSATONIC.EQ.1)
     &  CALL NEGDEPHOUS(QCHANUT,QCHANVT,2)
C
C**********************************************************************C
C
 6910 FORMAT('  DRYING AT N,I,J =',I10,2I6,'  HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6911 FORMAT('  DRY W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6912 FORMAT('  DRY E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6913 FORMAT('  DRY S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6914 FORMAT('  DRY N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6920 FORMAT('  WETTING AT N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6921 FORMAT('  WET S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6922 FORMAT('  WET W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6923 FORMAT('  WET E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6924 FORMAT('  WET N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6930 FORMAT('  WET BY VOL  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6940 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6941 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HUE,HP,H1P ='
     &         ,3(2X,E12.4))
 6942 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HUW,HP,H1P ='
     &         ,3(2X,E12.4))
 6943 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HVS,HP,H1P ='
     &         ,3(2X,E12.4))
 6944 FORMAT('  RESOLVE,  N,I,J =',I10,2I6,' HVN,HP,H1P ='
     &         ,3(2X,E12.4))
 6945 FORMAT('  RESOLVE NEG,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 6950 FORMAT('  RESOLVE, NEG DEP N,I,J =',I10,2I6,' HP,H1P,H2P ='
     &         ,3(2X,E12.4))
 8001 FORMAT(I7,5I5,4E13.4)
 8002 FORMAT(17X,3I5,4E13.4)
 8003 FORMAT(32X,4E13.4)
 8000 FORMAT('    N    NMD  MTYP   I    J  IDRY      P           H',
     &       '           P1           H1')
 8004 FORMAT('                                     QCHANU',
     &       '       QCHANUT      CCCCHU       CCCCHV ')   
C
C**********************************************************************C
C
C **  CALCULATE THE EXTERNAL DIVERGENCE 
C
C----------------------------------------------------------------------C
C
      IF(ISDIVEX.EQ.1)THEN
C
      DIVEXMX=0.
      DIVEXMN=1000000.
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      IF(SPB(L).NE.0)THEN
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI
     &     +0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     &     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L)
     &     +RIFTR(L)+EVAPSW(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      ENDIF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      ENDIF
      ENDIF
      ENDDO
C
      IMAX=IL(LMAX)
      JMAX=JL(LMAX)
      IMIN=IL(LMIN)
      JMIN=JL(LMIN)
C
      WRITE(6,6628)DIVEXMX,IMAX,JMAX
      WRITE(6,6629)DIVEXMN,IMIN,JMIN
C
C----------------------------------------------------------------------C
C
      ENDIF
C
C----------------------------------------------------------------------C
C
  566 FORMAT('  I=',I5,3X,'J=',I5,3X,'HP=',F12.4)
 6628 FORMAT('  DIVEXMX=',E13.5,5X,2I10)
 6629 FORMAT('  DIVEXMN=',E13.5,5X,2I10)
C
C**********************************************************************C
C
C **  UPDATE ZERO DIMENSION VOLUME BALANCE
C 
C----------------------------------------------------------------------C
C
      ISTL=2
      IF(ISDRY.GE.1.AND.ISTL.EQ.3)THEN
        VOLADD=0.
        DO L=2,LA
        IF(SPB(L).NE.0)THEN
          VOLADD=VOLADD+QSUME(L)-RIFTR(L)-EVAPSW(L)
        ENDIF
        ENDDO
        VOLADD=VOLADD*DT
        VOLZERD=VOLZERD+VOLADD
        VETZERD=VETZERD+VOLADD+DT*EVAPSW(L)
      ENDIF
 
 5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E13.5))
C                         
C**********************************************************************C
C
      RETURN
      END
