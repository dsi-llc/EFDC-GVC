C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WAVESXY
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
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  INPUT WAVE INFORMATION
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
      OPEN(1,FILE='WAVE.INP',STATUS='UNKNOWN')
C
      DO NSKIP=1,40
      READ(1,1,IOSTAT=ISO)
      IF(ISO.GT.0) GOTO 1081
      ENDDO
C
      READ(1,*,IOSTAT=ISO)NWVDAT,WVPRD,CVTWHA,ISWCBL,ISWRSR,ISWRSI,
     &          NWUPDT,NTSWV,WVDISV,WVDISH,WVLSH,WVLSX,ISWVSD,ISDZBR 
      IF(ISO.GT.0) GOTO 1082
      WVFRQ=2.*PI/WVPRD
      NWCUNT=NWUPDT-1
      JSWAVE=1
      RSWRSR=FLOAT(ISWRSR)
      RSWRSI=FLOAT(ISWRSI)
C
      DO NWV=1,NWVDAT
      READ(1,*,IOSTAT=ISO)IWV,JWV,HMPWV,HMCWV,ENETMP,SXXTMP,SYYTMP,
     &                     SXYTMP,DISPTMP
      IF(ISO.GT.0) GOTO 1083
      L=LIJ(IWV,JWV)
      HMPW(L)=HMPWV
      HMCW(L)=HMCWV
      WVENEP(L)=ENETMP
      WVHUU(L,KC)=SXXTMP
      WVHVV(L,KC)=SYYTMP
      WVHUV(L,KC)=SXYTMP
      WVDISP(L,KC)=DISPTMP
      ENDDO
C
      DO NWV=1,NWVDAT
      READ(1,*,IOSTAT=ISO)IWV,JWV,HMUWV,HMVWV,UWVRET,UWVIMT,VWVRET,
     &                     VWVIMT
      IF(ISO.GT.0) GOTO 1084
      L=LIJ(IWV,JWV)
      HMUW(L)=HMUWV
      HMVW(L)=HMVWV
      UWVRE(L,KC)=UWVRET
      UWVIM(L,KC)=UWVIMT
      VWVRE(L,KC)=VWVRET
      VWVIM(L,KC)=VWVIMT
      ENDDO
C
      CLOSE(1)
C
C**********************************************************************C
C
C **  DETERMINE ORTIBAL AMPLITUDES AND APPOXIMATE ANGLE
C
      DO L=2,LA
      LN=LNC(L)
      LE=L+1
      UWVMAG(L)=0.5*SQRT(UWVRE(L ,KC)*UWVRE(L ,KC)
     &                  +UWVIM(L ,KC)*UWVIM(L ,KC))
     &         +0.5*SQRT(UWVRE(LE,KC)*UWVRE(LE,KC)
     &                  +UWVIM(LE,KC)*UWVIM(LE,KC))+1.E-12
      VWVMAG(L)=0.5*SQRT(VWVRE(L ,KC)*VWVRE(L ,KC)
     &                  +VWVIM(L ,KC)*VWVIM(L ,KC))
     &         +0.5*SQRT(VWVRE(LN,KC)*VWVRE(LN,KC)
     &                  +VWVIM(LN,KC)*VWVIM(LN,KC))
      WACCWE(L)=ATAN2(VWVMAG(L),UWVMAG(L))
      WVWHA(L)=SQRT(2.*WVENEP(L)/G)
      ENDDO
C
C**********************************************************************C
C
C **  INITIALIZE VERTICAL DISTRIBUTION OF WAVE DISSIPATION AS SOURCE
C **  TO VERTICAL TKE CLOSURE
C
      IF(KC.EQ.2)THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=WVDISV
      ENDIF
C
      IF(KC.EQ.3)THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=0.5*WVDISV
      WVDTKEM(2)=0.5*WVDISV
      WVDTKEP(2)=WVDISV
      ENDIF
C
      IF(KC.GE.4)THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=0.5*WVDISV
      WVDTKEM(KS)=0.5*WVDISV
      WVDTKEP(KS)=WVDISV
      DO K=2,KS-1
       WVDTKEM(K)=0.5*WVDISV
       WVDTKEP(K)=0.5*WVDISV
      ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  GENERATE WAVE TABLE
C
      HMXTMP=0.
      DO L=2,LA
      HMXTMP=MAX(HMXTMP,HMPW(L))
      ENDDO
C
      FKHMAX=1.5*GI*WVFRQ*WVFRQ*HMXTMP
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
      OPEN(1,FILE='WVTAB.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WVTAB.OUT',STATUS='UNKNOWN')
C
      DO NKH=1,1001
      WRITE(1,111)RKHTAB(NKH),FUNKH(NKH)
      ENDDO
C
      CLOSE(1)
C
      GOTO 100
C
 1081 WRITE(6,1091)
      STOP
 1082 WRITE(6,1092)
      STOP
 1083 WRITE(6,1093) NWV
      STOP
 1084 WRITE(6,1094) NWV
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
C **  INITIALIZE OR UPDATE WAVE FIELD
C
  100 CONTINUE
C
      NWCUNT=NWCUNT+1
      IF(NWCUNT.LT.NWUPDT) RETURN
C
      NWCUNT=0
      GDFRQ=G/WVFRQ
C
C**********************************************************************C
C
C **  DISTRIBUTE WVHUU, WVHVV, AND WVDISP OVER DEPTH 
C
      IF(ISWAVE.GE.2)THEN
C
      DO L=1,LC
      HFFDG=GI*WVFRQ*WVFRQ*HMPW(L)
      WVKHP(L)=VALKH(HFFDG)
      ENDDO
C
C     WVTMP1(L)=0.5*G*WVWHA*WVWHA
C     WVTMP2(L)=SINH(2KH)
C     WVTMP3(L)=COSH(KH)
C     WVTMP4(L)=GROUP VEL/PHASE VEL
C
      DO L=2,LA
C     WVTMP1(L)=0.5*G*WVWHA(L)*WVWHA(L)
      WVTMP2(L)=SINH(2.*WVKHP(L))
      WVTMP3(L)=COSH(WVKHP(L))
      WVTMP4(L)=0.5+( WVKHP(L)/WVTMP2(L) )
      ENDDO
C
      DO K=1,KC
      ZTOP=Z(K)
      ZBOT=Z(K-1)
      DO L=2,LA
      RKHM1=WVKHP(L)
      RKHM2=2.*WVKHP(L)
      SINHTOP=SINH(RKHM2*ZTOP)
      SINHBOT=SINH(RKHM2*ZBOT)
      COSHTOP=COSH(RKHM1*ZTOP)
      COSHBOT=COSH(RKHM1*ZBOT)
      TMPVAL=(RKHM2*(ZTOP-ZBOT)+SINH(RKHM2*ZTOP)-SINH(RKHM2*ZBOT))
     &      /(RKHM2+WVTMP2(L))
      TMPP1=-0.5*(ZTOP-ZBOT)+(ZTOP*COSHTOP-ZBOT*COSHBOT)/WVTMP3(L)
      TMPP2=(WVTMP4(L)-1.)*(SINHTOP-SINHBOT-2.*(ZTOP-ZBOT))
     &      /(WVTMP2(L)-2.)
      WVHUU(L,K)=TMPVAL*WVHUU(L,KC)
      WVHVV(L,K)=TMPVAL*WVHVV(L,KC)
      WVDISP(L,K)=TMPVAL*WVDISP(L,KC)
      WVPP(L,K)=WVENEP(L)*(TMPP1+TMPP2)
      ENDDO 
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING 
C **  THE WAVE TURBULENT INTENSITY, QQWV
C **  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE
C
      IF(ISWAVE.GE.1)THEN
C
      DO L=2,LA
      AEXTMP=WVWHA(L)/SINH(WVKHP(L))
      UWVSQ(L)=AEXTMP*AEXTMP*WVFRQ*WVFRQ
      CDRGTMP=(30.*ZBR(L)/AEXTMP)**0.2
      CDRGTMP=5.57*CDRGTMP-6.13
      CDRGTMP=EXP(CDRGTMP)
      CDRGTMP=MIN(CDRGTMP,0.22)
      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
      QQWV1(L)=CTURB2*TMPVAL
      TAUTMP=TMPVAL/TAUR(NSED+1)
      CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
      ZBRE(L)=CORZBR*ZBR(L)
      CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2
      CDRGTMP=5.57*CDRGTMP-6.13
      CDRGTMP=EXP(CDRGTMP)
      CDRGTMP=MIN(CDRGTMP,0.22)
      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
      QQWV2(L)=CTURB2*TMPVAL
      ENDDO
C
      IF(ISRESTI.NE.0)THEN
       OPEN(1,FILE='WVQWCP.INP',STATUS='UNKNOWN')
       DO L=2,LA
       READ(1,*)IDUM,JDUM,QQWV1(L),QQWV2(L),QQWV2(L),QQWC(L),QQWCR(L)
       ENDDO
      ENDIF     
C
      ENDIF
C
C**********************************************************************C
C
C **  COMPUTE CELL CORNER QUANTITY WVHUV 
C
      IF(ISWAVE.GE.1)THEN
C
      DO L=2,LA
      HFFDG=GI*WVFRQ*WVFRQ*HMCW(L)
      WVKHC(L)=VALKH(HFFDG)
      ENDDO
C
C     WVTMP1(L)=0.5*G*WVWHA*WVWHA
C     WVTMP2(L)=SINH(2KH)
C     WVTMP3(L)=COSH(KN)
C     WVTMP4(L)=GROUP VEL/PHASE VEL
C
      DO L=2,LA
      WVTMP2(L)=SINH(2.*WVKHC(L))
      WVTMP4(L)=0.5+( WVKHC(L)/WVTMP2(L) )
      ENDDO
C
      DO K=1,KC
      ZTOP=Z(K)
      ZBOT=Z(K-1)
      DO L=2,LA
      RKHM2=2.*WVKHC(L)
      SINHTOP=SINH(RKHM2*ZTOP)
      SINHBOT=SINH(RKHM2*ZBOT)
      TMPVAL=(RKHM2*(ZTOP-ZBOT)+SINH(RKHM2*ZTOP)-SINH(RKHM2*ZBOT))
     &      /(RKHM2+WVTMP2(L))
      WVHUV(L,K)=TMPVAL*WVHUV(L,KC)
      ENDDO 
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  COMPUTE CELL FACE QUANTITIES WVPU,WVPV
C
      IF(ISWAVE.GE.1)THEN
C
C     WVTMP1(L)=SINH(WVKHU(L))      
C     WVTMP2(L)=0.5*GI*WVWHA*WVWHA*WVFRQ*WVFRQ/(WVTMP1*WVTMP1)
C     WVTMP3(L)=SINH(WVKHV(L))      
C     WVTMP4(L)=0.5*GI*WVWHA*WVWHA*WVFRQ*WVFRQ/(WVTMP3*WVTMP3)
C
      TMPVAL=0.5*WVFRQ*WVFRQ
C
      DO L=2,LA
      HFFDG=GI*WVFRQ*WVFRQ*HMUW(L)
      WVKHU(L)=VALKH(HFFDG)
      HFFDG=GI*WVFRQ*WVFRQ*HMVW(L)
      WVKHV(L)=VALKH(HFFDG)
      ENDDO
C
      DO L=2,LA
      LS=LSC(L)
      WVTMP1(L)=SINH(WVKHU(L))
      WVWHAUT=(WVWHA(L)+SUB(L)*WVWHA(L-1))/(1.+SUB(L))
      WVTMP2(L)=TMPVAL*WVWHAUT*WVWHAUT
     &         /(WVTMP1(L)*WVTMP1(L))
      WVWHAVT=(WVWHA(L)+SVB(L)*WVWHA(LS ))/(1.+SVB(L))
      WVTMP3(L)=SINH(WVKHV(L))
      WVTMP4(L)=TMPVAL*WVWHAVT*WVWHAVT
     &         /(WVTMP3(L)*WVTMP3(L))
      ENDDO
C      
      DO K=1,KC
      ZTOP=Z(K)
      ZBOT=Z(K-1)
      DO L=2,LA
      SNHTOPU=SINH(WVKHU(L)*ZTOP)
      SNHBOTU=SINH(WVKHU(L)*ZBOT)
      SNHTOPV=SINH(WVKHV(L)*ZTOP)
      SNHBOTV=SINH(WVKHV(L)*ZBOT)
      TMPPU=(1.-ZTOP)*SNHTOPU*(ZTOP*WVTMP1(L)-SNHTOPU)
     &     -(1.-ZBOT)*SNHBOTU*(ZBOT*WVTMP1(L)-SNHBOTU)
      TMPPV=(1.-ZTOP)*SNHTOPV*(ZTOP*WVTMP3(L)-SNHTOPV)
     &     -(1.-ZBOT)*SNHBOTV*(ZBOT*WVTMP3(L)-SNHBOTV)
      WVPU(L,K)=WVTMP2(L)*TMPPU
      WVPV(L,K)=WVTMP4(L)*TMPPV
      ENDDO
      ENDDO    
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE THE NET X AND  Y WAVE REYNOLDS STRESS FORCINGS
C
      IF(ISWAVE.GE.1)THEN
C
      DO K=1,KC
      DZITMP=1./DZC(K)
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      FXWAVE(L,K)=DZITMP*SUB(L)*SPB(L)
     &           *( RSWRSI*(DYU(L)*(WVPP(L,K)-WVPP(L-1,K))
     &                     +DYU(L)*WVPU(L,K)*(HMPW(L)-HMPW(L-1)))
     &           +RSWRSR*(DYP(L)*WVHUU(L,K)-DYP(L-1)*WVHUU(L-1,K)
     &           +0.5*(DXV(LN )+DXV(LNW))*WVHUV(LN,K)
     &           -0.5*(DXV(L  )+DXV(L-1))*WVHUV(L ,K)) )
C    &           -0.5*(DXV(L  )+DXV(L-1))*WVHUV(L ,K)
C    &           -0.5*(DYU(L  )-DYU(L-1))*WVHVV(L-1,K)
C    &           -0.5*(DYU(L+1)-DYU(L  ))*WVHVV(L  ,K)
C    &           +0.5*(DXV(LNW)-DXV(L-1))*WVHUVP(L-1,K)
C    &           +0.5*(DXV(LN )-DXV(L  ))*WVHUVP(L  ,K) )
      FYWAVE(L,K)=DZITMP*SVB(L)*SPB(L)
     &           *( RSWRSI*(DXV(L)*(WVPP(L,K)-WVPP(LS ,K))
     &                     +DXV(L)*WVPV(L,K)*(HMPW(L)-HMPW(LS )))
     &           +RSWRSR*(DXP(L)*WVHVV(L,K)-DXP(LS )*WVHVV(LS ,K)
     &           +0.5*(DYU(L+1)+DYU(LSE))*WVHUV(L+1,K)
     &           -0.5*(DYU(L  )+DYU(LS ))*WVHUV(L  ,K)) )
C    &           -0.5*(DXV(L  )+DXV(LS ))*WVHUV(L  ,K)
C    &           +0.5*(DYU(LSE)-DYU(LS ))*WVHUVP(LS ,K)
C    &           +0.5*(DYU(L+1)-DYU(L  ))*WVHUVP(L  ,K)
C    &           -0.5*(DXV(L  )-DXV(LS ))*WVHUU(LS ,K)
C    &           -0.5*(DXV(LN )-DXV(L  ))*WVHUU(L  ,K) )  
      ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISPGNS.GE.2)THEN
C
      DO K=1,KC  
      DO NPNS=1,NPNSBP
      L=LIJ(ISPNS(NPNS),JSPNS(NPNS))
      FXWAVE(L,K)=0.
      FYWAVE(L,K)=0.
      L=LIJ(INPNS(NPNS),JNPNS(NPNS))
      FXWAVE(L,K)=0.
      FYWAVE(L,K)=0.
      ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE THE NONDIVERGENT OR VECTOR POTENTIAL COMPONENT OF THE
C **  WAVE STOKES DRIFT.
C
C     IF(ISWAVE.GE.1)THEN
C
C     DO K=1,KS
C     DO L=2,LA
C     LT=L
C     SINHB=SINH(WVKHP(LT))
C     COSHR=COSH(WVKHP(LT)*Z(K))/SINHB
C     SINHR=SINH(WVKHP(LT)*Z(K))/SINHB
C     VPMUL=0.25*WVFRQ*WVWHA(LT)*WVWHA(LT)*(COSHR*SINHR-Z(K)*COSHR)
C     VPXTMP=VPMUL*WVASIN(LT)
C     VPYTMP=VPMUL*WVACOS(LT)
C     LT=LSC(L)
C     SINHB=SINH(WVKHP(LT))
C     COSHR=COSH(WVKHP(LT)*Z(K))/SINHB
C     SINHR=SINH(WVKHP(LT)*Z(K))/SINHB
C     VPXTMP=0.25*WVFRQ*WVWHA(LT)*WVWHA(LT)*WVASIN(LT)*( COSHR*SINHR
C    &      -Z(K)*COSHR )+VPXTMP
C     VPX(L,K)=SVB(L)*VPXTMP
C     LT=L-1
C     SINHB=SINH(WVKHP(LT))
C     COSHR=COSH(WVKHP(LT)*Z(K))/SINHB
C     SINHR=SINH(WVKHP(LT)*Z(K))/SINHB
C     VPYTMP=0.25*WVFRQ*WVWHA(LT)*WVWHA(LT)*WVACOS(LT)*( COSHR*SINHR
C    &      -Z(K)*COSHR )+VPYTMP
C     VPY(L,K)=-SUB(L)*VPYTMP
C     ENDDO
C     ENDDO
C
C     DO K=1,KC
C     DO L=2,LA
C     UVPT(L,K)=-DZIC(K)*(VPY(L,K)-VPY(L,K-1))
C     VVPT(L,K)=DZIC(K)*(VPX(L,K)-VPX(L,K-1))
C     ENDDO
C     ENDDO      
C
C     DO K=1,KS
C     DO L=2,LA
C     LN=LNC(L)
C     WVPT(L,K)=( DYU(L+1)*VPY(L+1,K)-DYU(L)*VPY(L,K)
C    &           -DXV(LN )*VPX(LN ,K)+DXV(L)*VPX(L,K) )/DXYP(L)
C     ENDDO
C     ENDDO      
C
C     ENDIF
C
C**********************************************************************C
C
      IF(ISWAVE.GE.1) CALL WRSPLTH
C
C**********************************************************************C
C
      RETURN
      END
