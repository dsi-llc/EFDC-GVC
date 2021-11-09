C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP9 (ISTL)
C
C **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
C **  THIS EXPERIMENTAL VERSION IMPLEMENTS COSMIC ADVECTION
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
      DIMENSION UCTR(LCM,KCM),VCTR(LCM,KCM),WCTR(LCM,KCM),
     &          UCOR(LCM,KCM),VCOR(LCM,KCM),WCOR(LCM,KCM)
      DIMENSION UABU(LCM,0:KCM),VABU(LCM,0:KCM),WABU(LCM,0:KCM),
     &          UABV(LCM,0:KCM),VABV(LCM,0:KCM),WABV(LCM,0:KCM)
C
C**********************************************************************C
C
       DELT=DT2
       DELTD2=DT
      IF(ISTL.EQ.2)THEN
       DELT=DT
       DELTD2=0.5*DT
      ENDIF
C
      DELTI=1./DELT
C
      IF(N.EQ.1)THEN
       OPEN(1,FILE='CFLMOM.OUT')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='CFLMOM.OUT')
      ELSE
       OPEN(1,FILE='CFLMOM.OUT',POSITION='APPEND')
      ENDIF
C
C**********************************************************************C
C
C **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
C
C----------------------------------------------------------------------C
C
      FCAXE(1)=0.
      FCAYE(1)=0.
      FCAX1E(1)=0.
      FCAY1E(1)=0.
      FXE(1)=0.
      FYE(1)=0.
C
      FCAXE(LC)=0.
      FCAYE(LC)=0.
      FCAX1E(LC)=0.
      FCAY1E(LC)=0.
      FXE(LC)=0.
      FYE(LC)=0.
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        FCAXE(L)=0.
        FCAYE(L)=0.
        FCAX1E(L)=0.
        FCAY1E(L)=0.
        FXE(L)=0.
        FYE(L)=0.
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3) GOTO 300
C
C**********************************************************************C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
C
  200 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=0.5*(UHDY1(L,K)+UHDY(L,K))
        VHDX2(L,K)=0.5*(VHDX1(L,K)+VHDX(L,K))
        W2(L,K)=0.5*(W1(L,K)+W(L,K))
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        UCTR(L,K)=0.
        VCTR(L,K)=0.
        WCTR(L,K)=0.
        UCOR(L,K)=0.
        VCOR(L,K)=0.
        WCOR(L,K)=0.
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FVHV(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=0,KC
       DO L=1,LC
        UABU(L,K)=0.
        VABU(L,K)=0.
        WABU(L,K)=0.
        UABV(L,K)=0.
        VABV(L,K)=0.
        WABV(L,K)=0.
        FWU(L,K)=0.
        FWV(L,K)=0.
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      TMP=0.25*DELT
C
      DO K=1,KC
       DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        UCTR(L,K)=TMP*(U(L,K)+U1(L,K))/DXU(L)
     &           +TMP*(U(L+1,K)+U1(L+1,K))/DXU(L+1)
        VCTR(L,K)=TMP*(V(L,K)+V1(L,K))/DYV(L)
     &           +TMP*(V(LN,K)+V1(LN,K))/DYV(LN)
        WCTR(L,K)=TMP*(W(L,K)+W(L,K-1))/(HP(L)*DZC(K))
     &           +TMP*(W1(L,K)+W1(L,K-1))/(H1P(L)*DZC(K))
        UCOR(L,K)=TMP*(U(L,K)+U1(L,K))/DXU(L)
     &           +TMP*(U(LS,K)+U1(LS,K))/DXU(LS)
        VCOR(L,K)=TMP*(V(L,K)+V1(L,K))/DYV(L)
     &           +TMP*(V(L-1,K)+V1(L-1,K))/DYV(L-1)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        LSW=LSWC(L)
        WCOR(L,K)=0.25*(WCTR(L,K)+WCTR(LS,K)+WCTR(LSW,K)+WCTR(L-1,K))
       ENDDO
      ENDDO
C
      DO K=1,KC-1
       DO L=1,LC
        LS=LSC(L)
        UABU(L,K)=TMP*(U(L,K)+U1(L,K))/DXU(L)
     &           +TMP*(U(L,K+1)+U1(L,K+1))/DXU(L)
        VABV(L,K)=TMP*(V(L,K)+V1(L,K))/DYV(L)
     &           +TMP*(V(L,K+1)+V1(L,K+1))/DYV(L)
        WABU(L,K)=TMP*W(L,K)/(HP(L)*DZC(K))
     &           +TMP*W1(L,K)/(H1P(L)*DZC(K))
     &           +TMP*W(L-1,K)/(HP(L-1)*DZC(K))
     &           +TMP*W1(L-1,K)/(H1P(L-1)*DZC(K))
        WABV(L,K)=TMP*W(L,K)/(HP(L)*DZC(K))
     &           +TMP*W1(L,K)/(H1P(L)*DZC(K))
     &           +TMP*W(LS,K)/(HP(LS)*DZC(K))
     &           +TMP*W1(LS,K)/(H1P(LS)*DZC(K))
       ENDDO
      ENDDO
C
      DO K=1,KC-1
       DO L=1,LC
        LN=LNC(L)
        LS=LSC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        UABV(L,K)=0.25*(UABU(L,K)+UABU(LS,K)+UABU(L+1,K)+UABU(LSE,K))
        VABU(L,K)=0.25*(VABV(L,K)+VABV(L-1,K)+VABV(LN,K)+VABV(LNW,K))
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      IF(ABS(UCTR(L,K)).GT.1.) WRITE(1,801)N,IL(L),JL(L),K,UCTR(L,K)
      IF(ABS(VCTR(L,K)).GT.1.) WRITE(1,802)N,IL(L),JL(L),K,VCTR(L,K)
      IF(ABS(WCTR(L,K)).GT.1.) WRITE(1,803)N,IL(L),JL(L),K,WCTR(L,K)
      IF(ABS(UCOR(L,K)).GT.1.) WRITE(1,804)N,IL(L),JL(L),K,UCOR(L,K)
      IF(ABS(VCOR(L,K)).GT.1.) WRITE(1,805)N,IL(L),JL(L),K,VCOR(L,K)
      IF(ABS(WCOR(L,K)).GT.1.) WRITE(1,806)N,IL(L),JL(L),K,WCOR(L,K)
      IF(ABS(UABU(L,K)).GT.1.) WRITE(1,807)N,IL(L),JL(L),K,UABU(L,K)
      IF(ABS(VABU(L,K)).GT.1.) WRITE(1,808)N,IL(L),JL(L),K,VABU(L,K)
      IF(ABS(WABU(L,K)).GT.1.) WRITE(1,809)N,IL(L),JL(L),K,WABU(L,K)
      IF(ABS(UABV(L,K)).GT.1.) WRITE(1,810)N,IL(L),JL(L),K,UABV(L,K)
      IF(ABS(VABV(L,K)).GT.1.) WRITE(1,811)N,IL(L),JL(L),K,VABV(L,K)
      IF(ABS(WABV(L,K)).GT.1.) WRITE(1,812)N,IL(L),JL(L),K,WABV(L,K)
      ENDDO
      ENDDO
C
  801 FORMAT('CRN: N,I,J,K,UCTR = ',I10,3I5,F10.2)
  802 FORMAT('CRN: N,I,J,K,VCTR = ',I10,3I5,F10.2)
  803 FORMAT('CRN: N,I,J,K,WCTR = ',I10,3I5,F10.2)
  804 FORMAT('CRN: N,I,J,K,UCOR = ',I10,3I5,F10.2)
  805 FORMAT('CRN: N,I,J,K,VCOR = ',I10,3I5,F10.2)
  806 FORMAT('CRN: N,I,J,K,WCOR = ',I10,3I5,F10.2)
  807 FORMAT('CRN: N,I,J,K,UABU = ',I10,3I5,F10.2)
  808 FORMAT('CRN: N,I,J,K,VABU = ',I10,3I5,F10.2)
  809 FORMAT('CRN: N,I,J,K,WABU = ',I10,3I5,F10.2)
  810 FORMAT('CRN: N,I,J,K,UABV = ',I10,3I5,F10.2)
  811 FORMAT('CRN: N,I,J,K,VABV = ',I10,3I5,F10.2)
  812 FORMAT('CRN: N,I,J,K,WABV = ',I10,3I5,F10.2)
C
C----------------------------------------------------------------------C
C
C ** FUHU
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        FUHU(L,K)=0.5*(UHDY2(L,K)+UHDY2(L+1,K))*UFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).LT.0.)THEN
          KW=K+1
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
          DZDYU=U1(LV,K+1)-U1(LVS,K+1)-U1(LV,K)+U1(LVS,K)
          UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        ENDIF
        FUHU(L,K)=0.5*(UHDY2(L,K)+UHDY2(L+1,K))*UFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        KW=K
        IF(WCTR(L,K).LT.0.) KW=K+1 
        UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDYU=U1(LV,K)-U1(LVS,K)-U1(LV,K-1)+U1(LVS,K-1)
         ELSE
          DZDYU=U1(LV,K+1)-U1(LVS,K+1)-U1(LV,K)+U1(LVS,K)
        ENDIF
        UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
        FUHU(L,K)=0.5*(UHDY2(L,K)+UHDY2(L+1,K))*UFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
          DZDYU=U1(LV,K)-U1(LVS,K)-U1(LV,K-1)+U1(LVS,K-1)
          UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        ENDIF
        FUHU(L,K)=0.5*(UHDY2(L,K)+UHDY2(L+1,K))*UFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FVHU
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        FVHU(L,K)=0.5*(VHDX2(L,K)+VHDX2(L-1,K))*UFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCOR(L,K).LT.0.)THEN
          KW=K+1
          UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
          DZDXU=U1(LU,K+1)-U1(LW,K+1)-U1(LU,K)+U1(LW,K)
          UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        ENDIF
        FVHU(L,K)=0.5*(VHDX2(L,K)+VHDX2(L-1,K))*UFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        KW=K
        IF(WCOR(L,K).LT.0.) KW=K+1 
        UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDXU=U1(LU,K)-U1(LW,K)-U1(LU,K-1)+U1(LW,K-1)
         ELSE
          DZDXU=U1(LU,K+1)-U1(LW,K+1)-U1(LU,K)+U1(LW,K)
        ENDIF
        UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
        FVHU(L,K)=0.5*(VHDX2(L,K)+VHDX2(L-1,K))*UFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
          DZDXU=U1(LU,K)-U1(LW,K)-U1(LU,K-1)+U1(LW,K-1)
          UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        ENDIF
        FVHU(L,K)=0.5*(VHDX2(L,K)+VHDX2(L-1,K))*UFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FWU
C
      IF(KC.GE.2)THEN
C
      DO K=1,KS
       DO L=2,LA
        KW=K
        IF(WABU(L,K).LT.0.) KW=K+1
        LU=L
        IF(UABU(L,K).LT.0.) LU=L+1
        LW=LU-1
        LV=L
        IF(VABU(L,K).LT.0.) LV=LNC(L)
        LVS=LSC(LV)
        UFACE=U1(L,KW)
     &       -0.5*UABU(L,K)*(U1(LU,KW)-U1(LW,KW))
     &       -0.5*VABU(L,K)*(U1(LV,KW)-U1(LVS,KW))
        IF(VABU(L,K).GT.0.)THEN
          DXDYU=U1(LU,KW)-U1(LW,KW)-U1(LSC(LU),KW)+U1(LSC(LW),KW)
         ELSE
          DXDYU=U1(LNC(LU),KW)-U1(LNC(LW),KW)-U1(LU,KW)+U1(LW,KW)
        ENDIF
        UFACE=UFACE+UABU(L,K)*VABU(L,K)*DXDYU/3.
        FWU(L,K)=0.5*(DXYP(L)*W2(L,K)+DXYP(L-1)*W2(L-1,K))*UFACE
       ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FUHV
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        FUHV(L,K)=0.5*(UHDY2(L,K)+UHDY2(LSC(L),K))*VFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCOR(L,K).LT.0.)THEN
          KW=K+1
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
          DZDYV=V1(LV,K+1)-V1(LVS,K+1)-V1(LV,K)+V1(LVS,K)
          VFACE=VFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
         ELSE
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        ENDIF
        FUHV(L,K)=0.5*(UHDY2(L,K)+UHDY2(LSC(L),K))*VFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        KW=K
        IF(WCOR(L,K).LT.0.) KW=K+1 
        VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDYV=V1(LV,K)-V1(LVS,K)-V1(LV,K-1)+V1(LVS,K-1)
         ELSE
          DZDYV=V1(LV,K+1)-V1(LVS,K+1)-V1(LV,K)+V1(LVS,K)
        ENDIF
        VFACE=VFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
        FUHV(L,K)=0.5*(UHDY2(L,K)+UHDY2(LSC(L),K))*VFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
          DZDYV=V1(LV,K)-V1(LVS,K)-V1(LV,K-1)+V1(LVS,K-1)
          VFACE=UFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
         ELSE
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        ENDIF
        FUHV(L,K)=0.5*(UHDY2(L,K)+UHDY2(LSC(L),K))*VFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FVHV
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
        FVHV(L,K)=0.5*(VHDX2(L,K)+VHDX2(LNC(L),K))*VFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).LT.0.)THEN
          KW=K+1
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
          DZDXV=V1(LU,K+1)-V1(LW,K+1)-V1(LU,K)+V1(LW,K)
          VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
         ELSE
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
        ENDIF
        FVHV(L,K)=0.5*(VHDX2(L,K)+VHDX2(LNC(L),K))*VFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        KW=K
        IF(WCTR(L,K).LT.0.) KW=K+1 
        VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDXV=V1(LU,K)-V1(LW,K)-V1(LU,K-1)+V1(LW,K-1)
         ELSE
          DZDXV=V1(LU,K+1)-V1(LW,K+1)-V1(LU,K)+V1(LW,K)
        ENDIF
        VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
        FVHV(L,K)=0.5*(VHDX2(L,K)+VHDX2(LNC(L),K))*VFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
          DZDXV=V1(LU,K)-V1(LW,K)-V1(LU,K-1)+V1(LW,K-1)
          VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
         ELSE
          VFACE=V1(LV,K)
     &       -0.5*WCTR(L,K)*(V1(LU,K)-V1(LW,K))
        ENDIF
        FVHV(L,K)=0.5*(VHDX2(L,K)+VHDX2(LNC(L),K))*VFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FWV
C
      IF(KC.GE.2)THEN
C
      DO K=1,KS
       DO L=2,LA
        LS=LSC(L)
        KW=K
        IF(WABV(L,K).LT.0.) KW=K+1
        LU=L
        IF(UABV(L,K).LT.0.) LU=L+1
        LW=LU-1
        LV=L
        IF(VABV(L,K).LT.0.) LV=LNC(L)
        LVS=LSC(LV)
        VFACE=V1(L,KW)
     &       -0.5*UABV(L,K)*(V1(LU,KW)-V1(LW,KW))
     &       -0.5*VABV(L,K)*(V1(LV,KW)-V1(LVS,KW))
        IF(VABV(L,K).GT.0.)THEN
          DXDYV=V1(LU,KW)-V1(LW,KW)-V1(LSC(LU),KW)+V1(LSC(LW),KW)
         ELSE
          DXDYV=V1(LNC(LU),KW)-V1(LNC(LW),KW)-V1(LU,KW)+V1(LW,KW)
        ENDIF
        VFACE=VFACE+UABV(L,K)*VABV(L,K)*DXDYV/3.
        FWV(L,K)=0.5*(DXYP(L)*W2(L,K)+DXYP(LS)*W2(LS,K))*VFACE
       ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      GOTO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        UCTR(L,K)=0.
        VCTR(L,K)=0.
        WCTR(L,K)=0.
        UCOR(L,K)=0.
        VCOR(L,K)=0.
        WCOR(L,K)=0.
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FVHV(L,K)=0.
       ENDDO
      ENDDO
C
      DO K=0,KC
       DO L=1,LC
        UABU(L,K)=0.
        VABU(L,K)=0.
        WABU(L,K)=0.
        UABV(L,K)=0.
        VABV(L,K)=0.
        WABV(L,K)=0.
        FWU(L,K)=0.
        FWV(L,K)=0.
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      TMP=0.5*DELT
C
      DO K=1,KC
       DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        UCTR(L,K)=TMP*U(L,K)/DXU(L)+TMP*U(L+1,K)/DXU(L+1)
        VCTR(L,K)=TMP*V(L,K)/DYV(L)+TMP*V(LN,K)/DYV(LN)
        WCTR(L,K)=TMP*(W(L,K)+W(L,K-1))/(HP(L)*DZC(K))
        UCOR(L,K)=TMP*U(L,K)/DXU(L)+TMP*U(LS,K)/DXU(LS)
        VCOR(L,K)=TMP*V(L,K)/DYV(L)+TMP*V(L-1,K)/DYV(L-1)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        LS=LSC(L)
        LSW=LSWC(L)
        WCOR(L,K)=0.25*(WCTR(L,K)+WCTR(LS,K)+WCTR(LSW,K)+WCTR(L-1,K))
       ENDDO
      ENDDO
C
      DO K=1,KC-1
       DO L=1,LC
        LS=LSC(L)
        UABU(L,K)=TMP*(U(L,K)+U(L,K+1))/DXU(L)
        VABV(L,K)=TMP*(V(L,K)+V(L,K+1))/DYV(L)
        WABU(L,K)=TMP*W(L,K)/(HP(L)*DZC(K))
     &           +TMP*W(L-1,K)/(HP(L-1)*DZC(K))
        WABV(L,K)=TMP*W(L,K)/(HP(L)*DZC(K))
     &           +TMP*W(LS,K)/(HP(LS)*DZC(K))
       ENDDO
      ENDDO
C
      DO K=1,KC-1
       DO L=1,LC
        LN=LNC(L)
        LS=LSC(L)
        LNW=LNWC(L)
        LSE=LSEC(L)
        UABV(L,K)=0.25*(UABU(L,K)+UABU(LS,K)+UABU(L+1,K)+UABU(LSE,K))
        VABU(L,K)=0.25*(VABV(L,K)+VABV(L-1,K)+VABV(LN,K)+VABV(LNW,K))
       ENDDO
      ENDDO
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      IF(ABS(UCTR(L,K)).GT.1.) WRITE(1,801)N,IL(L),JL(L),K,UCTR(L,K)
      IF(ABS(VCTR(L,K)).GT.1.) WRITE(1,802)N,IL(L),JL(L),K,VCTR(L,K)
      IF(ABS(WCTR(L,K)).GT.1.) WRITE(1,803)N,IL(L),JL(L),K,WCTR(L,K)
      IF(ABS(UCOR(L,K)).GT.1.) WRITE(1,804)N,IL(L),JL(L),K,UCOR(L,K)
      IF(ABS(VCOR(L,K)).GT.1.) WRITE(1,805)N,IL(L),JL(L),K,VCOR(L,K)
      IF(ABS(WCOR(L,K)).GT.1.) WRITE(1,806)N,IL(L),JL(L),K,WCOR(L,K)
      IF(ABS(UABU(L,K)).GT.1.) WRITE(1,807)N,IL(L),JL(L),K,UABU(L,K)
      IF(ABS(VABU(L,K)).GT.1.) WRITE(1,808)N,IL(L),JL(L),K,VABU(L,K)
      IF(ABS(WABU(L,K)).GT.1.) WRITE(1,809)N,IL(L),JL(L),K,WABU(L,K)
      IF(ABS(UABV(L,K)).GT.1.) WRITE(1,810)N,IL(L),JL(L),K,UABV(L,K)
      IF(ABS(VABV(L,K)).GT.1.) WRITE(1,811)N,IL(L),JL(L),K,VABV(L,K)
      IF(ABS(WABV(L,K)).GT.1.) WRITE(1,812)N,IL(L),JL(L),K,WABV(L,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C ** FUHU
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        FUHU(L,K)=0.5*(UHDY(L,K)+UHDY(L+1,K))*UFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).LT.0.)THEN
          KW=K+1
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
          DZDYU=U1(LV,K+1)-U1(LVS,K+1)-U1(LV,K)+U1(LVS,K)
          UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        ENDIF
        FUHU(L,K)=0.5*(UHDY(L,K)+UHDY(L+1,K))*UFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        KW=K
        IF(WCTR(L,K).LT.0.) KW=K+1 
        UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDYU=U1(LV,K)-U1(LVS,K)-U1(LV,K-1)+U1(LVS,K-1)
         ELSE
          DZDYU=U1(LV,K+1)-U1(LVS,K+1)-U1(LV,K)+U1(LVS,K)
        ENDIF
        UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
        FUHU(L,K)=0.5*(UHDY(L,K)+UHDY(L+1,K))*UFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LU=L
        IF(UCTR(L,K).LT.0.) LU=L+1
        LV=LU
        IF(VCTR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
     &       -0.5*WCTR(L,K)*(U1(LU,KW)-U1(LU,KW-1))
          DZDYU=U1(LV,K)-U1(LVS,K)-U1(LV,K-1)+U1(LVS,K-1)
          UFACE=UFACE+VCTR(L,K)*WCTR(L,K)*DZDYU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*VCTR(L,K)*(U1(LV,K)-U1(LVS,K))
        ENDIF
        FUHU(L,K)=0.5*(UHDY(L,K)+UHDY(L+1,K))*UFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FVHU
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        FVHU(L,K)=0.5*(VHDX(L,K)+VHDX(L-1,K))*UFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCOR(L,K).LT.0.)THEN
          KW=K+1
          UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
          DZDXU=U1(LU,K+1)-U1(LW,K+1)-U1(LU,K)+U1(LW,K)
          UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        ENDIF
        FVHU(L,K)=0.5*(VHDX(L,K)+VHDX(L-1,K))*UFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        KW=K
        IF(WCOR(L,K).LT.0.) KW=K+1 
        UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDXU=U1(LU,K)-U1(LW,K)-U1(LU,K-1)+U1(LW,K-1)
         ELSE
          DZDXU=U1(LU,K+1)-U1(LW,K+1)-U1(LU,K)+U1(LW,K)
        ENDIF
        UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
        FVHU(L,K)=0.5*(VHDX(L,K)+VHDX(L-1,K))*UFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LV=LSC(L)
        IF(VCOR(L,K).LT.0.) LV=L
        LU=LV
        IF(UCOR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          UFACE=U1(LV,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
     &       -0.5*WCOR(L,K)*(U1(LV,KW)-U1(LV,KW-1))
          DZDXU=U1(LU,K)-U1(LW,K)-U1(LU,K-1)+U1(LW,K-1)
          UFACE=UFACE+UCOR(L,K)*WCOR(L,K)*DZDXU/3.
         ELSE
          UFACE=U1(LU,K)
     &       -0.5*UCOR(L,K)*(U1(LU,K)-U1(LW,K))
        ENDIF
        FVHU(L,K)=0.5*(VHDX(L,K)+VHDX(L-1,K))*UFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FWU
C
      IF(KC.GE.2)THEN
C
      DO K=1,KS
       DO L=2,LA
        KW=K
        IF(WABU(L,K).LT.0.) KW=K+1
        LU=L
        IF(UABV(L,K).LT.0.) LU=L+1
        LW=LU-1
        LV=L
        IF(VABU(L,K).LT.0.) LV=LNC(L)
        LVS=LSC(LV)
        UFACE=U1(L,KW)
     &       -0.5*UABU(L,K)*(U1(LU,KW)-U1(LW,KW))
     &       -0.5*VABU(L,K)*(U1(LV,KW)-U1(LVS,KW))
        IF(VABU(L,K).GT.0.)THEN
          DXDYU=U1(LU,KW)-U1(LW,KW)-U1(LSC(LU),KW)+U1(LSC(LW),KW)
         ELSE
          DXDYU=U1(LNC(LU),KW)-U1(LNC(LW),KW)-U1(LU,KW)+U1(LW,KW)
        ENDIF
        UFACE=UFACE+UABU(L,K)*VABU(L,K)*DXDYU/3.
        FWU(L,K)=0.5*(DXYP(L)*W(L,K)+DXYP(L-1)*W(L-1,K))*UFACE
       ENDDO
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FUHV
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        FUHV(L,K)=0.5*(UHDY(L,K)+UHDY(LSC(L),K))*VFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCOR(L,K).LT.0.)THEN
          KW=K+1
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
          DZDYV=V1(LV,K+1)-V1(LVS,K+1)-V1(LV,K)+V1(LVS,K)
          VFACE=VFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
         ELSE
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        ENDIF
        FUHV(L,K)=0.5*(UHDY(L,K)+UHDY(LSC(L),K))*VFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        KW=K
        IF(WCOR(L,K).LT.0.) KW=K+1 
        VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDYV=V1(LV,K)-V1(LVS,K)-V1(LV,K-1)+V1(LVS,K-1)
         ELSE
          DZDYV=V1(LV,K+1)-V1(LVS,K+1)-V1(LV,K)+V1(LVS,K)
        ENDIF
        VFACE=VFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
        FUHV(L,K)=0.5*(UHDY(L,K)+UHDY(LSC(L),K))*VFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LU=L-1
        IF(UCOR(L,K).LT.0.) LU=L
        LV=LU
        IF(VCOR(L,K).LT.0.) LV=LNC(LU)
        LVS=LSC(LV)
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
     &       -0.5*WCOR(L,K)*(V1(LU,KW)-V1(LU,KW-1))
          DZDYV=V1(LV,K)-V1(LVS,K)-V1(LV,K-1)+V1(LVS,K-1)
          VFACE=UFACE+VCOR(L,K)*WCOR(L,K)*DZDYV/3.
         ELSE
          VFACE=V1(LU,K)
     &       -0.5*VCOR(L,K)*(V1(LV,K)-V1(LVS,K))
        ENDIF
        FUHV(L,K)=0.5*(UHDY(L,K)+UHDY(LSC(L),K))*VFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FVHV
C
      IF(KC.EQ.1)THEN
       K=1
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
        FVHV(L,K)=0.5*(VHDX(L,K)+VHDX(LNC(L),K))*VFACE
       ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
       K=1 
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).LT.0.)THEN
          KW=K+1
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
          DZDXV=V1(LU,K+1)-V1(LW,K+1)-V1(LU,K)+U1(LW,K)
          VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
         ELSE
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
        ENDIF
        FVHV(L,K)=0.5*(VHDX(L,K)+VHDX(LNC(L),K))*VFACE
       ENDDO
C
      IF(KC.GT.2)THEN
      DO K=2,KC-1
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        KW=K
        IF(WCTR(L,K).LT.0.) KW=K+1 
        VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
        IF(WCTR(L,K).GT.0.)THEN
          DZDXV=V1(LU,K)-V1(LW,K)-V1(LU,K-1)+V1(LW,K-1)
         ELSE
          DZDXV=V1(LU,K+1)-V1(LW,K+1)-V1(LU,K)+V1(LW,K)
        ENDIF
        VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
        FVHV(L,K)=0.5*(VHDX(L,K)+VHDX(LNC(L),K))*VFACE
       ENDDO
      ENDDO
      ENDIF
C
       K=KC
       DO L=2,LA
        LV=L
        IF(VCTR(L,K).LT.0.) LV=LNC(L)
        LU=LV
        IF(UCTR(L,K).LT.0.) LU=LV+1
        LW=LU-1
        IF(WCTR(L,K).GT.0.)THEN
          KW=K
          VFACE=V1(LV,K)
     &       -0.5*UCTR(L,K)*(V1(LU,K)-V1(LW,K))
     &       -0.5*WCTR(L,K)*(V1(LV,KW)-V1(LV,KW-1))
          DZDXV=V1(LU,K)-V1(LW,K)-V1(LU,K-1)+V1(LW,K-1)
          VFACE=VFACE+UCTR(L,K)*WCTR(L,K)*DZDXV/3.
         ELSE
          VFACE=V1(LV,K)
     &       -0.5*WCTR(L,K)*(V1(LU,K)-V1(LW,K))
        ENDIF
        FVHV(L,K)=0.5*(VHDX(L,K)+VHDX(LNC(L),K))*VFACE
       ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** FWV
C
      IF(KC.GE.2)THEN
C
      DO K=1,KS
       DO L=2,LA
        LS=LSC(L)
        KW=K
        IF(WABV(L,K).LT.0.) KW=K+1
        LU=L
        IF(UABV(L,K).LT.0.) LU=L+1
        LW=LU-1
        LV=L
        IF(VABV(L,K).LT.0.) LV=LNC(L)
        LVS=LSC(LV)
        VFACE=V1(L,KW)
     &       -0.5*UABV(L,K)*(V1(LU,KW)-V1(LW,KW))
     &       -0.5*VABV(L,K)*(V1(LV,KW)-V1(LVS,KW))
        IF(VABV(L,K).GT.0.)THEN
          DXDYV=V1(LU,KW)-V1(LW,KW)-V1(LSC(LU),KW)+V1(LSC(LW),KW)
         ELSE
          DXDYV=V1(LNC(LU),KW)-V1(LNC(LW),KW)-V1(LU,KW)+V1(LW,KW)
        ENDIF
        VFACE=VFACE+UABV(L,K)*VABV(L,K)*DXDYV/3.
        FWV(L,K)=0.5*(DXYP(L)*W(L,K)+DXYP(LS)*W(LS,K))*VFACE
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
  500 CONTINUE
C
      DO K=1,KC
       DO L=1,LA
       FUHU(L,K)=STCUV(L)*FUHU(L,K)
       FVHV(L,K)=STCUV(L)*FVHV(L,K)
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      IF(ISDCCA.EQ.0)THEN
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       ENDDO
      ENDDO
C
      ELSE
C
      CFMAX=CF
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)          
       CAC(L,K)=( FCORC(L)*DXYP(L)
     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
       CFMAX=MAX(CFMAX,CFEFF)
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        FCAX1(L,K)=0.
        FCAY1(L,K)=0.     
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.3)THEN
C      IF(IRVEC.EQ.10.OR.IRVEC.EQ.11)THEN
      IF(IRVEC.EQ.14)THEN
C
      DO L=1,LC
       TVAR3E(L)=0.
       TVAR3N(L)=0.
       TVAR3W(L)=0.
       TVAR3S(L)=0.
       TVAR1S(L,1)=0.
       TVAR2S(L,1)=0.
      ENDDO
C
      DO K=1,KC
      DZTMP=DZC(K) 
      DO L=2,LA
       TVAR3E(L)=TVAR3E(L)+DZTMP*U(L,K)
       TVAR3N(L)=TVAR3N(L)+DZTMP*V(L,K)
       TVAR3W(L)=TVAR3W(L)+DZTMP*U1(L,K)
       TVAR3S(L)=TVAR3S(L)+DZTMP*V1(L,K)
      ENDDO
      ENDDO
C
C       CAC(L,K)=( FCORC(L)*DXYP(L)
C     &        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
C     &        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
C
C       TVAR1S=CAC @ N
C       TVAR2S=CAC @ N-1
C
C       TVAR3E=U @ N
C       TVAR3N=V @ N
C       TVAR3W=U @ N-1
C       TVAR3S=V @ N-1
C
      DO L=2,LA
       LN=LNC(L)          
       TVAR1S(L,1)=( FCORC(L)*DXYP(L)
     &      +0.5*SNLT*(TVAR3N(LN)+TVAR3N(L))*(DYU(L+1)-DYU(L))
     &      -0.5*SNLT*(TVAR3E(L+1)+TVAR3E(L))*(DXV(LN)-DXV(L)) )*HP(L)
       TVAR2S(L,1)=( FCORC(L)*DXYP(L)
     &      +0.5*SNLT*(TVAR3S(LN)+TVAR3S(L))*(DYU(L+1)-DYU(L))
     &      -0.5*SNLT*(TVAR3W(L+1)+TVAR3W(L))*(DXV(LN)-DXV(L)) )*H1P(L)
      ENDDO
C
C       FCAX(L,K)=ROLD*FCAX(L,K)
C     &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
C     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
C       FCAY(L,K)=ROLD*FCAY(L,K)
C     &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
C     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
C
C
C       TVAR1S=CAC @ N
C       TVAR2S=CAC @ N-1
C
C       TVAR3E=U @ N
C       TVAR3N=V @ N
C       TVAR3W=U @ N-1
C       TVAR3S=V @ N-1
C
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=0.0625*SCAX(L)*
     &            ( TVAR2S(L,1)*(TVAR3S(LN)+TVAR3S(L))
     &             +TVAR2S(L-1,1)*(TVAR3S(LNW)+TVAR3S(L-1)) )
       FCAY1E(L)=0.0625*SCAY(L)*
     &            ( TVAR2S(L,1)*(TVAR3W(L+1)+TVAR3W(L))
     &             +TVAR2S(LS,1)*(TVAR3W(LSE)+TVAR3W(LS)) )
      ENDDO
C
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=FCAX1E(L)-0.125*SCAX(L)*
     &            ( TVAR1S(L,1)*(TVAR3N(LN)+TVAR3N(L))
     &             +TVAR1S(L-1,1)*(TVAR3N(LNW)+TVAR3N(L-1)) )
       FCAY1E(L)=FCAY1E(L)-0.125*SCAY(L)*
     &            ( TVAR1S(L,1)*(TVAR3E(L+1)+TVAR3E(L))
     &             +TVAR1S(LS,1)*(TVAR3E(LSE)+TVAR3E(LS)) )
      ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=ROLD*FCAX(L,K)
     &          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     &                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
       FCAY(L,K)=ROLD*FCAY(L,K)
     &          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
     &                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K))
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K))
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF(ISHDMF.GE.1)THEN
C
      DO K=1,KC
       DO L=2,LA
       LS=LSC(L)
       LN=LNC(L)
       FX(L,K)=FX(L,K)
     &        -SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
       FY(L,K)=FY(L,K)
     &        -SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.2)THEN
C
      IF(N.LT.NTSWV)THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
       ELSE
        WVFACT=1.0
      ENDIF        
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FXE(L)=FXE(L)+WVFACT*SAAX(L)*FXWAVE(L,K)*DZC(K)
         FYE(L)=FYE(L)+WVFACT*SAAY(L)*FYWAVE(L,K)*DZC(K)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF(ISWAVE.GE.2)THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
         FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
C **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*
     &            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(L-1,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     &            (HMP(L)-HMP(L-1)-Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*
     &            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
     &                     +(B(L,K)-B(LS,K))*DZC(K) )
     &           +(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     &            (HMP(L)-HMP(LS)-Z(K)*(HP(L)-HP(LS))) )
      ENDDO
      ENDDO
C
C     IF(N.EQ.1)THEN
C       OPEN(1,FILE='BUOY.DIA',STATUS='UNKNOWN')
C       DO L=2,LA
C        DO K=1,KS
C        TMP3D(K)=SUBO(L)*FBBX(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C        DO K=1,KS
C        TMP3D(K)=SVBO(L)*FBBY(L,K)
C        ENDDO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C       ENDDO
C       CLOSE(1)
C     ENDIF
C
 1111 FORMAT(2I5,2X,8E12.4)       
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
        DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.2)THEN
C 
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       ENDDO
      ENDDO
C
      ENDIF
C
      CLOSE(1)
C
C**********************************************************************C
C
      RETURN
      END
