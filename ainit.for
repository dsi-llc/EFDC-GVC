C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE AINIT
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 03/05/2002        john hamrick       03/05/2002       john hamrick
C  added transport bypass mask, imaskdry for dry cells
C 03/19/2002        john hamrick       03/19/2002       john hamrick
C  added transport bypass mask, lmaskdry for dry cells
C  modified definition of chanlen in initialization rather than 
C  in subs caltbxy and calpuv2c and calpuv9c
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  INITIALIZE ARRAYS  
C
C----------------------------------------------------------------------C
C
      ZBR(1)=ZBRADJ
      ZBRE(1)=ZBRADJ
      HMP(1)=HMIN
      HMU(1)=HMIN
      HMV(1)=HMIN
      HMC(1)=HMIN
      H1C(1)=HMIN
      HWQ(1)=HMIN
      H2WQ(1)=HMIN
      DXP(1)=DX
      DYP(1)=DY
      DXU(1)=DX
      DYU(1)=DY
      DXV(1)=DX
      DYV(1)=DY
      DXYP(1)=DX*DY
      MVEGL(1)=1
      BELV(1)=BELV(2)
      ZBR(LC)=ZBRADJ
      ZBRE(LC)=ZBRADJ
      HMP(LC)=HMIN
      HMU(LC)=HMIN
      HMV(LC)=HMIN
      HMC(LC)=HMIN
      H1C(LC)=HMIN
      HWQ(LC)=HMIN
      H2WQ(LC)=HMIN
      DXP(LC)=DX
      DYP(LC)=DY
      DXU(LC)=DX
      DYU(LC)=DY
      DXV(LC)=DX
      DYV(LC)=DY
      DXYP(LC)=DX*DY
      MVEGL(LC)=1
      BELV(LC)=BELV(LA)
C
      IF(ISGWIE.EQ.0) DAGWZ=0.
      DO L=2,LA
      I=IL(L)
      J=JL(L)
C     HMP(L)=.25*(H(I,J)+H(I+1,J)+H(I,J+1)+H(I+1,J+1))
C     HMU(L)= .5*(H(I,J)+H(I,J+1))                    
C     HMU(L)=0.5*(H(I-1,J)+H(I,J))
C     HMV(L)= .5*(H(I,J)+H(I+1,J))                    
C     HMV(L)=0.5*(H(I,J-1)+H(I,J))
      BELAGW(L)=BELV(L)-DAGWZ
      ZBRE(L)=ZBR(L)
C     DLON(L)=-77.0+(0.25*FLOAT(I)-1.125)/60.
      DLON(L)=CDLON1+(CDLON2*FLOAT(I)+CDLON3)/60.
C     DLAT(L)=36.0+(0.2*FLOAT(J)+48.9)/60.
      DLAT(L)=CDLAT1+(CDLAT2*FLOAT(J)+CDLAT3)/60.
      CUE(L)=1.
      CVE(L)=0.
      CUN(L)=0.
      CVN(L)=1.
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA      
      LS=LSC(L)
      DXU(L)=0.5*(DXP(L)+DXP(L-1))
      DYU(L)=0.5*(DYP(L)+DYP(L-1))
      DXV(L)=0.5*(DXP(L)+DXP(LS))
      DYV(L)=0.5*(DYP(L)+DYP(LS))
      ENDDO
C
      DO L=2,LA      
      LS=LSC(L)
C      HMU(L)=0.5*(HMP(L)+HMP(L-1))
C      HMV(L)=0.5*(HMP(L)+HMP(LS))
       HMU(L)=0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(L-1)*DYP(L-1)*HMP(L-1))
     &           /(DXU(L)*DYU(L))
       HMV(L)=0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(LS )*DYP(LS)*HMP(LS ))
     &           /(DXV(L)*DYV(L))
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO L=1,LC
      AMCP(L)=0.
      AMSP(L)=0.
      BE(L)=0.
      BI1(L)=0.
      BI2(L)=0.
      CS(L)=0.
      CW(L)=0.
      CC(L)=1.
      CE(L)=0.
      CN(L)=0.
      CCS(L)=0.
      CCW(L)=0.
      CCC(L)=1.
      CCE(L)=0.
      CCN(L)=0.
      FP(L)=0.
      FCAXE(L)=0.
      FCAYE(L)=0.
      FCAX1E(L)=0.
      FCAY1E(L)=0.
      FPGXE(L)=0.
      FPGYE(L)=0.
      FXE(L)=0.
      FYE(L)=0.
      FUHDYE(L)=0.
      FVHDXE(L)=0.
      P(L)=G*(HMP(L)+BELV(L))
      P1(L)=G*(HMP(L)+BELV(L))
      PRED(L)=0.
      PBLK(L)=0.
      PAM(L)=0.
      PPH(L)=0.
      QSUME(L)=0.
      QSUMELPF(L)=0.
      HP(L)=HMP(L)+PDGINIT
      HU(L)=HMU(L)+PDGINIT
      HV(L)=HMV(L)+PDGINIT
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      HWQ(L)=HMP(L)+PDGINIT
      H1P(L)=HMP(L)+PDGINIT
      H2P(L)=HMP(L)+PDGINIT
      H1U(L)=HMU(L)+PDGINIT
      H1V(L)=HMV(L)+PDGINIT
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
      H2WQ(L)=HMP(L)+PDGINIT
      HLPF(L)=0.
      HRU(L)=0.
      HRV(L)=0.
      SCB(L)=1.
      SPB(L)=1.
      SUB(L)=1.
      SVB(L)=1.
      SWB(L)=1.
      STCUV(L)=1.
      STCAP(L)=1.
      STBX(L)=1.
      STBY(L)=1.
      SAAX(L)=1.
      SAAY(L)=1.
      SNLPX(L)=1.
      SNLPY(L)=1.
      SCAX(L)=1.
      SCAY(L)=1.
      SBX(L)=1.
      SBY(L)=1.
      SDX(L)=1.
      SDY(L)=1.
      TBX(L)=0.
      TBY(L)=0.
      TSX(L)=0.
      TSY(L)=0.
      TBX1(L)=0.
      TBY1(L)=0.
      TSX1(L)=0.
      TSY1(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      UV(L)=0.
      U1V(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      VU(L)=0.
      V1U(L)=0.
      SFLSBOT(L)=0.
      TVAR3S(L)=0.
      TVAR3W(L)=0.
      TVAR3E(L)=0.
      TVAR3N(L)=0.
      RIFTR(L)=0.
      EVAPSW(L)=0.
      EVAPGW(L)=0.
      QQL(L,0)=0.
      QQL1(L,0)=0.
      QQL2(L,0)=0.
      QQL(L,KC)=0.
      QQL1(L,KC)=0.
      QQL2(L,KC)=0.
      DML(L,0)=0.
      DML(L,KC)=0.
	QGW(L)=0.0
      QDWASTE(L)=0.0
      VDWASTE(L)=0.0
c
      QSBDTOPBEBKB(L)=0.
      QSBDTOPBECHB(L)=0.
      QSBDTOPBECHW(L)=0.
      QWBDTOPBEBKB(L)=0.
      QWBDTOPBECHB(L)=0.
      QWBDTOPBECHW(L)=0.
c
      IMASKDRY(L)=0
      LMASKDRY(L)=.TRUE.
	GVCSCLP(L)=1.0
	GVCSCLU(L)=1.0
	GVCSCLV(L)=1.0
	GVCSCLPI(L)=1.0
	GVCSCLUI(L)=1.0
	GVCSCLVI(L)=1.0
	KGVCP(L)=1
	KGVCU(L)=1
	KGVCV(L)=1
	KGVCW(L)=1
      ENDDO
C
      DO NS=1,NSED
      DO L=1,LC
        SEDFBEBKB(L,NS)=0.
        SEDFBECHB(L,NS)=0.
        SEDFBECHW(L,NS)=0.
      ENDDO
      ENDDO
C
      DO NS=1,NSND
      DO L=1,LC
        SNDFBEBKB(L,NS)=0.
        SNDFBECHB(L,NS)=0.
        SNDFBECHW(L,NS)=0.
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO L=1,LC
        TOXFBEBKB(L,NT)=0.
        TOXFBECHB(L,NT)=0.
        TOXFBECHW(L,NT)=0.
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=1,KB
      DO L=1,LC
       TOXB(L,K,NT)=TOXBINIT(L,K,NT)
       TOXB1(L,K,NT)=TOXBINIT(L,K,NT)
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSED
       DO K=1,KB
       DO L=1,LC
       SEDB(L,K,NS)=SEDBINIT(L,K,NS)
       SEDB1(L,K,NS)=SEDBINIT(L,K,NS)
       ENDDO
       ENDDO
      ENDDO
      DO NS=1,NSND
      NX=NS+NSED
       DO K=1,KB
       DO L=1,LC
       SNDB(L,K,NS)=SNDBINIT(L,K,NS)
       SNDB1(L,K,NS)=SNDBINIT(L,K,NS)
       ENDDO
       ENDDO
      ENDDO

C
      DO NS=1,NSED
        DO L=1,LC
          SEDFDTAP(L,NS)=0.0
          SEDFDTAN(L,NS)=0.0
        ENDDO
      ENDDO
C
      DO NS=1,NSND
        DO L=1,LC
          SNDFDTAP(L,NS)=0.0
          SNDFDTAN(L,NS)=0.0
        ENDDO
      ENDDO
C
      DO L=1,LC
      CSR(L)=0.
      CWR(L)=0.
      CER(L)=0.
      CNR(L)=0.
      CSB(L)=0.
      CWB(L)=0.
      CEB(L)=0.
      CNB(L)=0.
      CCSR(L)=0.
      CCWR(L)=0.
      CCER(L)=0.
      CCNR(L)=0.
      CCSB(L)=0.
      CCWB(L)=0.
      CCEB(L)=0.
      CCNB(L)=0.
      FPR(L)=0.
      FPB(L)=0.
      ISCDRY(L)=0
      NATDRY(L)=0
      ENDDO
C
      IF(IS1DCHAN.EQ.1)THEN
      DO L=1,LC
      FADYP(L)=1.
      FADYP1(L)=1.
      FADYP2(L)=1.
      WPDYP(L)=1.
      WPDYP1(L)=1.
      FADXP(L)=1.
      FADXP1(L)=1.
      FADXP2(L)=1.
      WPDXP(L)=1.
      WPDXP1(L)=1.
      FADYU(L)=1.
      FADYU1(L)=1.
      WPDYU(L)=1.
      WPDYU1(L)=1.
      FADXV(L)=1.
      FADXV1(L)=1.
      WPDXV(L)=1.
      WPDXV1(L)=1.
      DADH(L)=1.
      DADH1(L)=1.
      SRFXP(L)=0.
      SRFYP(L)=0.
      SRFXP1(L)=0.
      SRFYP1(L)=0.
      SRFXV(L)=0.
      SRFYU(L)=0.
      SRFXV1(L)=0.
      SRFYU1(L)=0.
      ENDDO
      ENDIF
C
      DO L=1,NLRPD
      NLRPDL(L)=1
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=1,LC
C
      UUU(L,K)=0.
      VVV(L,K)=0.
C     AMCW(L,K)=0.
C     AMSW(L,K)=0.
      AV(L,K)=AVO
      AVVI(L,K)=1./AVO
      AVUI(L,K)=1./AVO
      AB(L,K)=ABO
      ABLPF(L,K)=0.
      ABEFF(L,K)=0.
C     AMCAB(L,K)=0.
C     AMSAB(L,K)=0.
C     AMC2AB(L,K)=0.
C     AMS2AB(L,K)=0.
      QQL(L,K)=QQLMIN
      QQL1(L,K)=QQLMIN
      QQL2(L,K)=QQLMIN
      DML(L,K)=DMLMIN
      WIRT(L,K)=0.
      WTLPF(L,K)=0.
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=1,LC
C
      AH(L,K)=AHO
      AHU(L,K)=AHO
      AHULPF(L,K)=AHO
      AHV(L,K)=AHO
      AHVLPF(L,K)=AHO
      AHC(L,K)=AHO
      AQ(L,K)=AVO
      AMCU(L,K)=0.
      AMSU(L,K)=0.
      AMCV(L,K)=0.
      AMSV(L,K)=0.
      B(L,K)=0.
      SALLPF(L,K)=0.
      TEMLPF(L,K)=0.
      SFLLPF(L,K)=0.
      DYELPF(L,K)=0.
      B1(L,K)=0.
      CAC(L,K)=0.
      CU1(L,K)=0.
      CU2(L,K)=0.
      DXU1(L,K)=0.
      DYU1(L,K)=0.
      DXV1(L,K)=0.
      DYV1(L,K)=0.
      DU(L,K)=0.
      DV(L,K)=0.
      FX(L,K)=0.
      FY(L,K)=0.
      FBBX(L,K)=0.
      FBBY(L,K)=0.
      FCAX(L,K)=0.
      FCAY(L,K)=0.
      FCAX1(L,K)=0.
      FCAY1(L,K)=0.
      FMDUX(L,K)=0.
      FMDUY(L,K)=0.
      FMDVY(L,K)=0.
      FMDVX(L,K)=0.
      FUHU(L,K)=0.
      FUHV(L,K)=0.
      FVHU(L,K)=0.
      FVHV(L,K)=0.
      FWQQ(L,K)=0.
      U(L,K)=0.     
      U1(L,K)=0.     
      U2(L,K)=0.     
      UHDY(L,K)=0.     
      UHDY1(L,K)=0.     
      UHDY2(L,K)=0.
      UHDYWQ(L,K)=0.
      VHDXWQ(L,K)=0.
      UWQ(L,K)=0.
      VWQ(L,K)=0.     
      V(L,K)=0.     
      V1(L,K)=0.     
      V2(L,K)=0.     
      VHDX(L,K)=0.
      VHDX1(L,K)=0.
      VHDX2(L,K)=0.
      UVPT(L,K)=0.
      VVPT(L,K)=0.
      VPZ(L,K)=0.
      UHLPF(L,K)=0.
      UIRT(L,K)=0.
      ULPF(L,K)=0.
      UTLPF(L,K)=0.
      VHLPF(L,K)=0.
      VIRT(L,K)=0.
      VLPF(L,K)=0.
      VTLPF(L,K)=0.
      UUU(L,K)=0.
      VVV(L,K)=0.
C     UDBDXI(L,K)=0.   REPLACED BY TEMPORARY VARIABLE UUU
C     VDBDYI(L,K)=0.   REPLACED BY TEMPORARY VARIABLE VVV
      SAL(L,K)=0.
      SAL1(L,K)=0.
      SFL(L,K)=0.
      SFL2(L,K)=0.
      CWQ(L,K)=0.
      CWQ2(L,K)=0.
      TEM(L,K)=TEMO
      TEM1(L,K)=TEMO
      DYE(L,K)=0.
      DYE1(L,K)=0.
      QSUM(L,K)=0.
      QSUMLPF(L,K)=0.
	QWSEDA(L,K)=0.
      TVAR1S(L,K)=0.
      TVAR1W(L,K)=0.
      TVAR1E(L,K)=0.
      TVAR1N(L,K)=0.
      TVAR2S(L,K)=0.
      TVAR2W(L,K)=0.
      TVAR2E(L,K)=0.
      TVAR2N(L,K)=0.
      CTURBB1(L,K)=CTURB
      CTURBB2(L,K)=CTURB2B
      SUB3D(L,K)=SUB(L)
      SUB3DO(L,K)=SUB(L)
      SVB3D(L,K)=SVB(L)
      SVB3DO(L,K)=SVB(L)
      SBX3D(L,K)=SBX(L)
      SBY3D(L,K)=SBY(L)
      SDX3D(L,K)=SDX(L)
      SDY3D(L,K)=SDY(L)
	LGVCP(L,K)=.TRUE.
	LGVCU(L,K)=.TRUE.
	LGVCV(L,K)=.TRUE.
C
      ENDDO
      ENDDO
C
      NTMPC=MAX(NSED,1)
      DO NS=1,NTMPC
      DO K=1,KC
      DO L=1,LC
       SED(L,K,NS)=SEDO(NS)
       SED1(L,K,NS)=SEDO(NS)
       SEDLPF(L,K,NS)=0.
      ENDDO
      ENDDO
      ENDDO
C
       NS=1
	 DO K=1,KC
	 DO L=2,LA
c	 FLOCDIA(L,K)=0.001
       FLOCDIA(L,K)=0.0004  !  initialized in cm equivalent to 4 um
       WSETFLOC(L,K)=80.E-6
       SEDFLOCDIA(L,K)=172.*SED(L,K,NS)/
     &                        (FLOCDIA(L,K)*WSETFLOC(L,K))
	 ENDDO
	 ENDDO
C
      NTMPN=MAX(NSND,1)
      DO NX=1,NTMPN
      NS=NX+NTMPC
      DO K=1,KC
      DO L=1,LC
       SND(L,K,NX)=SEDO(NS)
       SND1(L,K,NX)=SEDO(NS)
       SNDLPF(L,K,NX)=0.
      ENDDO
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=1,KC
      DO L=1,LC
       TOX(L,K,NT)=TOXINTW(NT)
       TOX1(L,K,NT)=TOXINTW(NT)
       TOXLPF(L,K,NT)=0.
      ENDDO
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=0,KC
      DO L=1,LC
C
      W(L,K)=0.
      W1(L,K)=0.
      W2(L,K)=0.
      WWQ(L,K)=0.
      WVPT(L,K)=0.
      WLPF(L,K)=0.
      FWU(L,K)=0.
      FWV(L,K)=0.
      QQ(L,K)=QQMIN
      QQ1(L,K)=QQMIN
      QQ2(L,K)=QQMIN
      VPX(L,K)=0.
      VPY(L,K)=0.
      WWW(L,K)=0.
      BBT(L,K)=0.
      BBT1(L,K)=0.
      BBT2(L,K)=0.
	SWB3D(L,K)=SWB(L)
      SWB3DO(L,K)=SWB(L)
C     WDBDZI(L,K)=0.  REPLACED BY TEMPORARY VARIABLE WWW
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KB
      DO L=1,LC
       VOLBW2(L,K)=0.
       VOLBW3(L,K)=0.
       PARTMIXZ(L,K)=0.
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
C
        DO NMD=1,MDCHH
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         SET HOST DRYING DEPTH
          IF(PMDCH(NMD).LT.0.0) PMDCH(NMD)=HWET
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            IF(CHANLEN(NMD).LT.0.0)THEN
              CHANLEN(NMD)=0.25*DYP(LHOST)
            ELSE
	        CHANLEN(NMD)=CHANLEN(NMD)-0.5*DYP(LCHNU)
	      ENDIF
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            IF(CHANLEN(NMD).LT.0.0)THEN
		    CHANLEN(NMD)=0.25*DXP(LHOST)
	      ELSE
	        CHANLEN(NMD)=CHANLEN(NMD)-0.5*DXP(LCHNV)
          ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** INITIALIZE ORGANIC CARBON VARIABLES OF SEDIMENT-TOXICS
C
      IVAL=0
	DO NT=1,NTOX
	  IF(ISTOC(NT).GT.0)IVAL=1
	ENDDO
C
	IF(IVAL.EQ.0)THEN
	  DO K=1,KB
	    DO L=1,LC
            STDOCB(L,K)=0.0
            STPOCB(L,K)=0.0
	    ENDDO
	  ENDDO
	  DO NS=1,NSED+NSND
	    DO K=1,KB
	      DO L=1,LC
              STFPOCB(L,K,NS)=1.0
            ENDDO
	    ENDDO
	  ENDDO
	  DO K=1,KC
	    DO L=1,LC
            STDOCW(L,K)=0.0
            STPOCW(L,K)=0.0
	    ENDDO
	  ENDDO
	  DO NS=1,NSED+NSND
	    DO K=1,KC
	      DO L=1,LC
              STFPOCW(L,K,NS)=1.0
            ENDDO
	    ENDDO
	  ENDDO
	ENDIF
C
      IF(IVAL.EQ.1)THEN
C     
      IF(ISTDOCB.EQ.0)THEN      
        DO K=1,KB
          DO L=1,LC
            STDOCB(L,K)=STDOCBC
          ENDDO
        ENDDO
	ENDIF
C
      IF(ISTPOCB.EQ.0)THEN      
        DO K=1,KB
          DO L=1,LC
            STPOCB(L,K)=STPOCBC
          ENDDO
        ENDDO
	ENDIF
C
      IF(ISTPOCB.EQ.2)THEN      
	  DO NS=1,NSED+NSND
          DO K=1,KB
            DO L=1,LC
              STFPOCB(L,K,NS)=FPOCBST(NS,1)
            ENDDO
          ENDDO
        ENDDO
	ENDIF
C
      IF(ISTDOCW.EQ.0)THEN      
        DO K=1,KC
          DO L=1,LC
            STDOCW(L,K)=STDOCWC
          ENDDO
        ENDDO
	ENDIF
C
      IF(ISTPOCW.EQ.0)THEN      
        DO K=1,KC
          DO L=1,LC
            STPOCW(L,K)=STPOCWC
          ENDDO
        ENDDO
	ENDIF
C
      IF(ISTPOCW.EQ.2)THEN      
	  DO NS=1,NSED+NSND
          DO K=1,KC
            DO L=1,LC
              STFPOCW(L,K,NS)=FPOCWST(NS,1)
            ENDDO
          ENDDO
        ENDDO
	ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C INITIALIZATION OF FOODCHAIN VARIABLES FOR HOUSATONIC
C
cjah FORTRAN 90 only
C
      ATOXPFWX=0.
      ATOXPFWY=0.
      ATOXPFBW=0.
      ATOXPFBWP=0.
      ATOXPFBWN=0.
      ATOXPFBLX=0.
      ATOXPFBLY=0.
      ATOXFDFWX=0.
      ATOXFDFWY=0.
      ATOXFDFBW=0.
      ATOXFDFBWP=0.
      ATOXFDFBWN=0.
      ATOXCDFWX=0.
      ATOXCDFWY=0.
      ATOXCDFBW=0.
      ATOXCDFBWP=0.
      ATOXCDFBWN=0.
      ATOXFBEBKB=0.
      ATOXFBECHB=0.
      ATOXFBECHW=0.
      ATOXVOL=0.
      ATOXSOUR=0.
      ATOXSINK=0.
      ATOXCDFBWAD=0.
      ATOXCDFBWADP=0.
      ATOXCDFBWADN=0.
      ATOXCDFBWAD=0.
      ATOXCDFBWADP=0.
      ATOXCDFBWADN=0.
      ASEDFWX=0.
      ASEDFWY=0.
      ASEDFBLX=0.
      ASEDFBLY=0.
C      
cjah      DO NT=1,NTOX
cjah        DO NS=1,NSED+NSND+2
cjah          DO L=1,LC
cjah            ATOXPFWX(L,NS,NT)=0.0
cjah            ATOXPFWY(L,NS,NT)=0.0
cjah            ATOXPFBW(L,NS,NT)=0.0
cjah            ATOXPFBWP(L,NS,NT)=0.0
cjah            ATOXPFBWN(L,NS,NT)=0.0
cjah          ENDDO
cjah        ENDDO
cjah      ENDDO
cjah
C
C
cjah      DO NT=1,NTOX
cjah        DO NS=1,NSND
cjah          DO L=1,LC
cjah            ATOXPFBLX(L,NS,NT)=0.0
cjah            ATOXPFBLY(L,NS,NT)=0.0
cjah          ENDDO
cjah        ENDDO
cjah      ENDDO
C
cjah      DO NT=1,NTOX
cjah        DO L=1,LC
cjah          ATOXFDFWX(L,NT)=0.0
cjah          ATOXFDFWY(L,NT)=0.0
cjah          ATOXFDFBW(L,NT)=0.0
cjah          ATOXFDFBWP(L,NT)=0.0
cjah          ATOXFDFBWN(L,NT)=0.0
cjah          ATOXCDFWX(L,NT)=0.0
cjah            ATOXCDFWY(L,NT)=0.0
cjah          ATOXCDFBW(L,NT)=0.0
cjah          ATOXCDFBWP(L,NT)=0.0
cjah          ATOXCDFBWN(L,NT)=0.0
cjah          ATOXFBEBKB(L,NT)=0.0
cjah            ATOXFBECHB(L,NT)=0.0
cjah          ATOXFBECHW(L,NT)=0.0
cjah          ATOXVOL(L,NT)=0.0
cjah            ATOXSOUR(L,NT)=0.0
cjah            ATOXSINK(L,NT)=0.0
cjah          ATOXCDFBWAD(L,NT)=0.0
cjah          ATOXCDFBWADP(L,NT)=0.0
cjah          ATOXCDFBWADN(L,NT)=0.0
cjah          ATOXCDFBWAD(L,NT)=0.0
cjah          ATOXCDFBWADP(L,NT)=0.0
cjah          ATOXCDFBWADN(L,NT)=0.0
cjah        ENDDO
cjah      ENDDO
cjahC
cjah        DO NS=1,NSED+NSND
cjah          DO L=1,LC
cjah           ASEDFWX(L,NS)=0.0
cjah             ASEDFWY(L,NS)=0.0
cjah          ENDDO
cjah        ENDDO
C
cjah        DO NS=1,NSND
cjah          DO L=1,LC
cjah           ASEDFBLX(L,NS)=0.0
cjah             ASEDFBLY(L,NS)=0.0
cjah          ENDDO
cjah        ENDDO
C
C
C**********************************************************************C
C
      RETURN 
      END
