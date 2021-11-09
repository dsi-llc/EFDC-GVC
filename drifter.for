C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE DRIFTER
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
C **  SUBROUTINE DRIFTER CALCULATES THREE DIMENSIONAL TRAJECTORIES 
C **  OF NEUTRALLY BUOYANT PARTICLES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      CHARACTER*80 TITLE
      TITLE='DRIFTER DATA'
      IF(JSPD.NE.1) GOTO 5
C
      OPEN(96,FILE='DRIFTER.OUT',STATUS='UNKNOWN')
      CLOSE(96,STATUS='DELETE')
      OPEN(96,FILE='DRIFTER.OUT',STATUS='UNKNOWN')
C
      JSPD=0
      WRITE(96,99) TITLE
      WRITE(96,98) NPD,KC
   98 FORMAT(2I10)
   99 FORMAT(A80)
  200 FORMAT(I10)
  201 FORMAT(2I4,1X,12F10.6)
C
C----------------------------------------------------------------------C
C
      WRITE(96,200)N
      DO NP=1,NPD
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      L=LIJ(I,J)
      ZS=(RK(NP)-0.5)/FLOAT(KC)
      DLOND=CDLON1+(CDLON2*RI(NP)+CDLON3)/60.
      DLATD=CDLAT1+(CDLAT2*RJ(NP)+CDLAT3)/60.
      WRITE(96,201)I,J,DLOND,DLATD,ZS,HP(L)
      ENDDO
C
      CLOSE(96)
C
C----------------------------------------------------------------------C
C
    5 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  LOOP OVER THE NUMBER OF PARTICLES
C
      DO NP=1,NPD
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      ENDIF
C
      TOLD=0.
C
  100 CONTINUE
C
      DTMIN=DELT
      TX=DELT
      TY=DELT
      TZ=DELT
      TNEW=TOLD+DTMIN
C
C----------------------------------------------------------------------C
C
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      K=NINT(RK(NP))
      L=LIJ(I,J)      
      LN=LNC(L)
C
      U1L=DTI*DXIU(L)*U1(L,K)
      U1LE=DTI*DXIU(L+1)*U1(L+1,K)
      V1L=DTI*DYIV(L)*V1(L,K)
      V1LN=DTI*DYIV(LN)*V1(LN,K)
      W1K=DTI*DZI*W1(L,K)/H1P(L)
      W1KM=DTI*DZI*W1(L,K-1)/H1P(L)
C
      UL=DTI*DXIU(L)*U(L,K)
      ULE=DTI*DXIU(L+1)*U(L+1,K)
      VL=DTI*DYIV(L)*V(L,K)
      VLN=DTI*DYIV(LN)*V(LN,K)
      WK=DTI*DZI*W(L,K)*HPI(L)
      WKM=DTI*DZI*W(L,K-1)*HPI(L)
C
      RIM=FLOAT(I)-0.5
      RJM=FLOAT(J)-0.5
      RKM=FLOAT(K)-0.5
      RIP=FLOAT(I)+0.5
      RJP=FLOAT(J)+0.5
      RKP=FLOAT(K)+0.5
C
      UOLD=U1L*(RI(NP)-RIP)*(TOLD-DT)-U1LE*(RI(NP)-RIM)*(TOLD-DT)
     &     -UL*(RI(NP)-RIP)*TOLD+ULE*(RI(NP)-RIM)*TOLD
C
      VOLD=V1L*(RJ(NP)-RJP)*(TOLD-DT)-V1LN*(RJ(NP)-RJM)*(TOLD-DT)
     &     -VL*(RJ(NP)-RJP)*TOLD+VLN*(RJ(NP)-RJM)*TOLD
C
      WOLD=W1KM*(RK(NP)-RKP)*(TOLD-DT)-W1K*(RK(NP)-RKM)*(TOLD-DT)
     &     -WKM*(RK(NP)-RKP)*TOLD+WK*(RK(NP)-RKM)*TOLD
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((U1LE*RIM-U1L*RIP)*(TNEW-DT)
     &   +(UL*RIP-ULE*RIM)*TNEW+UOLD)+RI(NP)
C
      BOT=(U1LE-U1L)*(TNEW-DT)+(UL-ULE)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RINEW=TOP/BOT
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((V1LN*RJM-V1L*RJP)*(TNEW-DT)
     &   +(VL*RJP-VLN*RJM)*TNEW+VOLD)+RJ(NP)
C
      BOT=(V1LN-V1L)*(TNEW-DT)+(VL-VLN)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RJNEW=TOP/BOT
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((W1K*RKM-W1KM*RKP)*(TNEW-DT)
     &   +(WKM*RKP-WK*RKM)*TNEW+WOLD)+RK(NP)
C
      BOT=(W1K-W1KM)*(TNEW-DT)+(WKM-WK)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RKNEW=TOP/BOT
C
C----------------------------------------------------------------------C
C
C **  DETERMINE TIME INCREMENT IN X DIRECTION
C
      IF(RINEW.GT.RIP)THEN
C
      AAQ=ULE-U1LE
      BQ=(ULE*TOLD-U1LE*(TOLD-DT)+UOLD)
      CQ=-2.*(RIP-RI(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TX=TQ
      ELSE
       BQ=(ULE*TOLD-U1LE*(TOLD-DT)+UOLD)/AAQ
       CQ=-2.*(RIP-RI(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 10
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT) 
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TX=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TX=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   10 CONTINUE
C
      IF(RINEW.LT.RIM)THEN
C
      AAQ=UL-U1L
      BQ=(UL*TOLD-U1L*(TOLD-DT)+UOLD)
      CQ=-2.*(RIM-RI(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TX=TQ
      ELSE
       BQ=(UL*TOLD-U1L*(TOLD-DT)+UOLD)/AAQ
       CQ=-2.*(RIM-RI(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 20
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TX=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TX=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   20 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  DETERMINE TIME INCREMENT IN Y DIRECTION
C
      IF(RJNEW.GT.RJP)THEN
C
      AAQ=VLN-V1LN
      BQ=(VLN*TOLD-V1LN*(TOLD-DT)+VOLD)
      CQ=-2.*(RJP-RJ(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TY=TQ
      ELSE
       BQ=(VLN*TOLD-V1LN*(TOLD-DT)+VOLD)/AAQ
       CQ=-2.*(RJP-RJ(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 30
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN 
        TY=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TY=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   30 CONTINUE
C
      IF(RJNEW.LT.RJM)THEN
C
      AAQ=VL-V1L
      BQ=(VL*TOLD-V1L*(TOLD-DT)+VOLD)
      CQ=-2.*(RJM-RJ(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TY=TQ
      ELSE
       BQ=(VL*TOLD-V1L*(TOLD-DT)+VOLD)/AAQ
       CQ=-2.*(RJM-RJ(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 40
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TY=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TY=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   40 CONTINUE
C
C----------------------------------------------------------------------C
C
C **  DETERMINE TIME INCREMENT IN Z DIRECTION
C
      IF(RKNEW.GT.RKP)THEN
      AAQ=WK-W1K
      BQ=(WK*TOLD-W1K*(TOLD-DT)+WOLD)
      CQ=-2.*(RKP-RK(NP))
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TZ=TQ
      ELSE
       BQ=(WK*TOLD-W1K*(TOLD-DT)+WOLD)/AAQ
       CQ=-2.*(RKP-RK(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 50
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TZ=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TZ=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   50 CONTINUE
C
      IF(RKNEW.LT.RKM)THEN
C
      AAQ=WKM-W1KM
      BQ=(WKM*TOLD-W1KM*(TOLD-DT)+WOLD)
      CQ=-2.*(RKM-RK(NP))
C
      IF(AAQ.EQ.0.)THEN
       TQ=-CQ/BQ
       IF(TQ.GT.0.) TZ=TQ
      ELSE
       BQ=(WKM*TOLD-W1KM*(TOLD-DT)+WOLD)/AAQ
       CQ=-2.*(RKM-RK(NP))/AAQ
       RSQR=BQ*BQ-4.*CQ
       IF(RSQR.LT.0.0) GOTO 60
       SQROOT=SQRT(RSQR)
       TQP=0.5*(-BQ+SQROOT)
       TQM=0.5*(-BQ-SQROOT)
       IF(TQP.GT.0.0.AND.TQM.GT.0.0)THEN
        TZ=MIN(TQP,TQM)
       ELSE
        TQ=MAX(TQP,TQM)
        IF(TQ.GT.0.) TZ=TQ
       ENDIF
      ENDIF
C
      ENDIF
C
   60 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(TX.LT.DTMIN) DTMIN=TX
      IF(TY.LT.DTMIN) DTMIN=TY
      IF(TZ.LT.DTMIN) DTMIN=TZ
C
      TNEW=TOLD+DTMIN
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((U1LE*RIM-U1L*RIP)*(TNEW-DT)
     &   +(UL*RIP-ULE*RIM)*TNEW+UOLD)+RI(NP)
C
      BOT=(U1LE-U1L)*(TNEW-DT)+(UL-ULE)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RI(NP)=TOP/BOT
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((V1LN*RJM-V1L*RJP)*(TNEW-DT)
     &   +(VL*RJP-VLN*RJM)*TNEW+VOLD)+RJ(NP)
C
      BOT=(V1LN-V1L)*(TNEW-DT)+(VL-VLN)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RJ(NP)=TOP/BOT
C
C----------------------------------------------------------------------C
C
      TOP=0.5*DTMIN*((W1K*RKM-W1KM*RKP)*(TNEW-DT)
     &   +(WKM*RKP-WK*RKM)*TNEW+WOLD)+RK(NP)
C
      BOT=(W1K-W1KM)*(TNEW-DT)+(WKM-WK)*TNEW
      BOT=1.0+0.5*DTMIN*BOT
C
      RK(NP)=TOP/BOT
C
C----------------------------------------------------------------------C
C
      DELT=DELT-DTMIN
      TOLD=TNEW
C
      IF(DELT.GT.0.) GOTO 100
C
C----------------------------------------------------------------------C
C
      ENDDO
C
C**********************************************************************C
C
      IF(NCPD.EQ.NWPD)THEN
C
      OPEN(96,FILE='DRIFTER.OUT',POSITION='APPEND',STATUS='UNKNOWN')
C
      NCPD=1
      WRITE(96,200)N
      DO NP=1,NPD
      I=NINT(RI(NP))
      J=NINT(RJ(NP))
      L=LIJ(I,J)
      ZS=(RK(NP)-0.5)/FLOAT(KC)
      DLOND=CDLON1+(CDLON2*RI(NP)+CDLON3)/60.
      DLATD=CDLAT1+(CDLAT2*RJ(NP)+CDLAT3)/60.
      WRITE(96,201)I,J,DLOND,DLATD,ZS,HP(L)
      ENDDO
C
      ELSE
      NCPD=NCPD+1
C
      CLOSE(96)
C
      ENDIF
C
C**********************************************************************C
C
      RETURN 
      END
