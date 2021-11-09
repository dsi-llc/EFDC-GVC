C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WAVEBL
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C----------------------------------------------------------------------C
cmlt +++++ modification for SWAN 40.11 input to be verified.
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER*9 FNWAVE
      CHARACTER*1 CFNWAVE(0:9)
C
C**********************************************************************C
C
C **  INITIALIZE AND INPUT WAVE INFORMATION
C
      IF(JSWAVE.EQ.1) GOTO 100
C
      JSWRPH=1
C
      DO L=1,LC
      HMPW(L)=0.
      HMCW(L)=0.
      HMUW(L)=0.
      HMVW(L)=0.
      WVWHA(L)=0.
C     WVACOS(L)=0.
C     WVASIN(L)=0.
      WVKHP(L)=0.
      WVKHC(L)=0.
      WVKHU(L)=0. 
      WVKHV(L)=0. 
      WVTMP1(L)=0.
      WVTMP2(L)=0.
      WVTMP3(L)=0.
      WVTMP4(L)=0.
      UWVMAG(L)=0.
      VWVMAG(L)=0.
      WVENEP(L)=0.
      UWVSQ(L)=0.
      QQWC(L)=1.E-12
      QQWCR(L)=1.E-12
      QQWV1(L)=1.E-12
      QQWV2(L)=1.E-12
      QQWV3(L)=1.E-12
      WACCWE(L)=0.
      ENDDO
C
      DO K=1,KC
      DO L=1,LC
      WVHUU(L,K)=0.
      WVHVV(L,K)=0.
      WVHUV(L,K)=0.
      WVPP(L,K)=0.
      WVPU(L,K)=0. 
      WVPV(L,K)=0.
      WVDISP(L,K)=0.
      UWVRE(L,K)=0.
      UWVIM(L,K)=0.
      VWVRE(L,K)=0.
      VWVIM(L,K)=0.
      FXWAVE(L,K)=0.
      FYWAVE(L,K)=0.
      ENDDO
      ENDDO
C
      OPEN(1,FILE='WAVEBL.INP',STATUS='UNKNOWN')
C
      DO NSKIP=1,11
      READ(1,1,IOSTAT=ISO)
      IF(ISO.GT.0) GOTO 1081
      ENDDO
C
      READ(1,*,IOSTAT=ISO)NWVDAT,CVTWHA,ISWCBL,NWUPDT,NTSWV,ISDZBR
      IF(ISO.GT.0) GOTO 1082
C
      CLOSE(1)
C
      FNWAVE='WV001.INP'
C
      OPEN(1,FILE=FNWAVE)
C
      DO L=2,LA
	 READ(1,*,IOSTAT=ISO) wwvh, wwdir, wwprdp
cmlt      READ(1,*,IOSTAT=ISO)IWVH,IWDIR,IWPRD1,IWPRD2
      IF(ISO.GT.0) GOTO 1083
      HMPW(L)=HMP(L)
	WVWHA(L)= wwvh
cmlt      WVWHA(L)=0.001*FLOAT(IWVH)
	WACCWE(L)= wwdir/57.29578
cmlt      WACCWE(L)=FLOAT(IWDIR)/57.29578
	IF(wwprdp.GT.0.) then
		WVFRQL(L)= wwprdp
	else
		WVFRQL(L)= 1.
		wvwha(l) = 0.
		waccwe(l) = 0.
	endif
cmlt      IF(IWPRD2.GT.0) WVFRQL(L)=0.01*FLOAT(IWPRD2)
cmlt      IF(IWPRD2.LE.0) WVFRQL(L)=1.
      WVFRQL(L)=2.*PI/WVFRQL(L)
      ENDDO
C
      CLOSE(1)
C
      OPEN(1,FILE='WAVEBL.DIA',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
C

      NTMP=0
      WRITE(6,666)NTMP,FNWAVE
      WRITE(8,666)NTMP,FNWAVE
C
      JSWAVE=1
      ITWCBL1=1
      ITWCBL2=0
      ITWCBL3=0
      NWCUNT=1
C
      GOTO 200
C
C**********************************************************************C
C
C ** UPDATE WAVE FIELD AT NWUPDT INTERVALS
C 
  100 CONTINUE
C
      JSWRPH=0
C
      IF(JSWAVE.EQ.1)THEN
        IF(NWCUNT.EQ.NWUPDT)THEN
          NWCUNT=1
          ITWCBL1=ITWCBL1+1
          IF(ITWCBL1.GT.9)THEN
           ITWCBL1=0
           ITWCBL2=ITWCBL2+1
          ENDIF
          IF(ITWCBL2.GT.9)THEN
           ITWCBL2=0
           ITWCBL3=ITWCBL3+1
          ENDIF
          CFNWAVE(0)='0'
          CFNWAVE(1)='1'
          CFNWAVE(2)='2'
          CFNWAVE(3)='3'
          CFNWAVE(4)='4'
          CFNWAVE(5)='5'
          CFNWAVE(6)='6'
          CFNWAVE(7)='7'
          CFNWAVE(8)='8'
          CFNWAVE(9)='9'
          FNWAVE='WV' // CFNWAVE(ITWCBL3) // CFNWAVE(ITWCBL2)
     &                // CFNWAVE(ITWCBL1)// '.INP'
        ELSE
          NWCUNT=NWCUNT+1
          RETURN
        ENDIF
      ENDIF
C
      OPEN(1,FILE=FNWAVE)
C
      DO L=2,LA
	 READ(1,*,IOSTAT=ISO) wwvh, wwdir, wwprdp
cmlt      READ(1,*,IOSTAT=ISO)IWVH,IWDIR,IWPRD1,IWPRDP
      IF(ISO.GT.0) GOTO 1083
      HMPW(L)=HMP(L)
	WVWHA(L)= 0.5 * (wwvh + WVWHA(L))
cmlt      WVWHA(L)=0.5*(0.001*FLOAT(IWVH)+WVWHA(L))
	WACCWE(L)= wwdir/57.29578
cmlt      WACCWE(L)=FLOAT(IWDIR)/57.29578
	IF(wwprdp.GT.0.) then
		WVFRQL(L)= wwprdp
	else
		WVFRQL(L)= 1.
		wvwha(l) = 0.
		waccwe(l) = 0.
	endif
cmlt      IF(IWPRDP.GT.0) WVFRQL(L)=0.01*FLOAT(IWPRDP)
cmlt      IF(IWPRDP.LE.0) WVFRQL(L)=1.
      WVFRQL(L)=2.*PI/WVFRQL(L)
      ENDDO
C
      CLOSE(1)
C
      WRITE(6,666)N,FNWAVE
      WRITE(8,666)N,FNWAVE
C
  666 FORMAT(' UPDATED WAVE FIELD N,FNWAVE = ',I12,A12)
C
C**********************************************************************C
C
C **  GENERATE WAVE TABLE
C
  200 CONTINUE
C
      HMXTMP=0.
      WVFRQM=0.
      DO L=2,LA
       HMXTMP=MAX(HMXTMP,HMPW(L))
       IF(WVWHA(L).GT.0.) WVFRQM=MAX(WVFRQM,WVFRQL(L))
      ENDDO
C
      FKHMAX=1.5*GI*WVFRQM*WVFRQM*HMXTMP
      RKHTMP=0.001
   10 CONTINUE
      FKHTMP=RKHTMP*TANH(RKHTMP)
      IF(FKHTMP.LT.FKHMAX)THEN
        RKHTMP=2.*RKHTMP
        GOTO 10
       ELSE
        DKH=RKHTMP/1000.
      ENDIF
C
      RKHTAB(1)=0.
      FUNKH(1)=0.
      DO NKH=2,1001
       RKHTAB(NKH)=RKHTAB(NKH-1)+DKH
       FUNKH(NKH)=RKHTAB(NKH)*TANH(RKHTAB(NKH))
      ENDDO
C
      DO L=2,LA
       HFFDG=GI*WVFRQL(L)*WVFRQL(L)*HMPW(L)
       WVKHP(L)=1.
       IF(WVWHA(L).GT.0.) WVKHP(L)=VALKH(HFFDG)
      ENDDO
C
      IF(JSWRPH.EQ.1)THEN
        OPEN(1,FILE='WVTAB.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='WVTAB.OUT',STATUS='UNKNOWN')
        DO NKH=1,1001
         WRITE(1,111)RKHTAB(NKH),FUNKH(NKH)
        ENDDO
        CLOSE(1)
      ENDIF
C
      GOTO 400
C
 1081 WRITE(6,1091)
      STOP
 1082 WRITE(6,1092)
      STOP
 1083 WRITE(6,1093) NWVDAT
      STOP
 1084 WRITE(6,1094) IWVH
      STOP
C
    1 FORMAT(120X)
 1091 FORMAT('  READ ERROR ON FILE WAVE.INP , HEADER')
 1092 FORMAT('  READ ERROR ON FILE WAVE.INP , 1ST DATA')
 1093 FORMAT('  READ ERROR ON FILE WAVE.INP , 2ND DATA, NWV = ',I5)
 1094 FORMAT('  READ ERROR ON FILE WAVE.INP , 3RD DATA, NWV = ',I5)
  111 FORMAT(2E14.4)
C
C**********************************************************************C
C
  400 CONTINUE
C
      DO L=2,LA
       IF(HMP(L).LT.0.55) WVWHA(L)=0. 
       IF(MVEGL(L).NE.MVEGOW) WVWHA(L)=0. 
       IWVRDC=0 
       IF(MVEGL(L).EQ.MVEGOW)THEN
         IF(MVEGL(L-1).NE.MVEGOW) IWVRDC=1 
         IF(MVEGL(L+1).NE.MVEGOW) IWVRDC=1 
         IF(MVEGL(LSC(L)).NE.MVEGOW) IWVRDC=1 
         IF(MVEGL(LNC(L)).NE.MVEGOW) IWVRDC=1 
       ENDIF
       IF(IWVRDC.GT.0) WVWHA(L)=0.5*WVWHA(L)
      ENDDO
C 
C**********************************************************************C
C
C **  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING 
C **  THE WAVE TURBULENT INTENSITY, QQWV
C **  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE
C
       OPEN(1,FILE='WAVEBL.DIA')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WAVEBL.DIA')
C
      DO L=2,LA
       AEXTMP=0.5*WVWHA(L)/SINH(WVKHP(L))
       UWORBIT=AEXTMP*WVFRQL(L)
       UWVSQ(L)=UWORBIT*UWORBIT
       VISMUD=1.E-6
       IF(ISMUD.GE.1) VISMUD=CSEDVIS(SED(L,1,1))
       DELWAVE=SQRT(VISMUD/WVFRQL(L))
       REYWAVE=UWORBIT*AEXTMP/VISMUD
       QQWV1(L)=0.
       QQWV2(L)=0.
C
C ** LAMINAR WAVE BOUNDARY LAYER
C
cjh       IF(REYWAVE.LT.5.E5)THEN
cjh         CDTMP=0.
cjh         IF(REYWAVE.GT.0.) CDTMP=1./SQRT(REYWAVE)
cjh         QQWV2(L)=CDTMP*UWORBIT*UWORBIT
cjh       ENDIF
C
C ** TURBULENT WAVE BOUNDARY LAYER
C
cjh       IF(REYWAVE.GE.5.E5)THEN
C
C ** TURBULENT SMOOTH WAVE BOUNDARY LAYER
C
         IF(ZBR(L).LE.0.)THEN
           CDTMP=0.012/(REYWAVE**0.123)
           QQWV1(L)=CDTMP*UWORBIT*UWORBIT
C
C ** TURBULENT ROUGH WAVE BOUNDARY LAYER
C
          ELSE
CJMH BEGIN JMH MODIFICATION
	     IF(WVWHA(L).GT.0.0)THEN
             TMPVAL=30.*ZBR(L)/AEXTMP
             TMPVAL=5.5*(TMPVAL**0.2)-6.3
             CDTMP=0.5*EXP(TMPVAL)
	     ELSE
	       CDTMP=0.0
	     ENDIF
CJMH END JMH MODIFICATION
           QQWV1(L)=CDTMP*UWORBIT*UWORBIT
           ZBRE(L)=ZBR(L)
           IF(QQ(L,0).GT.0.)THEN
             TMPVAL=UWORBIT*SQRT( AEXTMP/(30.*ZBR(L)) )
             USTARC=SQRT(QQ(L,0)/CTURB2)
             TMPVAL=TMPVAL/USTARC
             ZBRE(L)=ZBR(L)*(1.+0.19*TMPVAL) 
           ENDIF
         ENDIF
cjh       ENDIF
C
C       TVAL1=SINH(WVKHP(L))
C       TVAL2=CSEDVIS(SED(L,1,1))
        WRITE(1,600)L,IL(L),JL(L),WVWHA(L),WVFRQL(L),AEXTMP,UWORBIT,
     &             VISMUD,REYWAVE,CDTMP,QQWV1(L),QQWV2(L),ZBR(L),ZBRE(L)
      ENDDO
C
       CLOSE(1)
C
  600 FORMAT(3I5,11E12.4)
C
C      DO L=2,LA
C      AEXTMP=0.5*WVWHA(L)/SINH(WVKHP(L))
C      UWVSQ(L)=AEXTMP*AEXTMP*WVFRQ*WVFRQ
C      CDRGTMP=(30.*ZBR(L)/AEXTMP)**0.2
C      CDRGTMP=5.57*CDRGTMP-6.13
C      CDRGTMP=EXP(CDRGTMP)
C      CDRGTMP=MIN(CDRGTMP,0.22)
C      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
C      QQWV1(L)=CTURB2*TMPVAL
C      TAUTMP=TMPVAL/TAUR(NSED+1)
C      CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
C      ZBRE(L)=CORZBR*ZBR(L)
C      CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2
C      CDRGTMP=5.57*CDRGTMP-6.13
C      CDRGTMP=EXP(CDRGTMP)
C      CDRGTMP=MIN(CDRGTMP,0.22)
C      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
C      QQWV2(L)=CTURB2*TMPVAL
C      IF(ISTRAN(7).EQ.0) QQWV2(L)=0. 
C      ENDDO
C
C      IF(ISRESTI.NE.0)THEN
C       OPEN(1,FILE='WVQWCP.INP',STATUS='UNKNOWN')
C       DO L=2,LA
C       READ(1,*)IDUM,JDUM,QQWV1(L),QQWV2(L),QQWV2(L),QQWC(L),QQWCR(L)
C       ENDDO
C      ENDIF     
C
C**********************************************************************C
C
C      IF(ISWAVE.GE.1) CALL WRSPLTH
C
C**********************************************************************C
C
      RETURN
      END
