C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTRANI (ISTL,M,CON,CON1)
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
C **  SUBROUTINE CALTRANI CALCULATES THE ADVECTIVE AND DIFFUSIVE 
C **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
C **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES 
C **  THE NUMBER OF TIME LEVELS IN THE STEP. THE SOLUTION IS IMPLICIT
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION CON(LCM,KCM), CON1(LCM,KCM)
C     DIMENSION  CH(LCM,KCM), CMAX(LCM,KCM), CMIN(LCM,KCM),
C    &         CONT(LCM,KCM), CON2(LCM,KCM) 
C
C**********************************************************************C
C
      BSMALL=1.0E-6
C
      DELT=DT2
      DELTA=DT2
      IF(ISCDCA(M).EQ.2) DELTA=DT
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF(ISTL.NE.3)THEN
       DELT=DT
       DELTA=DT
       DELTD2=0.5*DT
       S3TL=0.0
       S2TL=1.0
       ISUD=0
      ENDIF
C
      IF(ISLTMT.GE.1) ISUD=1
C
C**********************************************************************C
C
      DO K=1,KC
      DO L=1,LC
      CONT(L,K)=0.
      CMAX(L,K)=0.
      CMIN(L,K)=0.
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT (N+1)
C
C----------------------------------------------------------------------C
C
      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
      FUHU(L,K)=MAX(UHDY2(L,K),0.)
      FUHV(L,K)=MIN(UHDY2(L,K),0.)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)
      FVHV(L,K)=MIN(VHDX2(L,K),0.)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)
      FWV(L,K)=MIN(W2(L,K),0.)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUX COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      CACTMP=FUHU(L+1,K)+FVHU(LN,K)-FUHV(L,K)-FVHV(L,K)
     &        +(FWU(L,K)-FWV(L,K-1))*DZIC(K)
      CACNEW=1./( CACTMP+2.*DXYP(L)*HP(L)/DELT )
      HTTMP=S3TL*H2P(L)+S2TL*H1P(L)
      CONT(L,K)=(CACTMP-2.*DXYP(L)*HTTMP/DELT)*CACNEW*CON1(L,K)      
      CAS(L,K)=-FVHU(L,K)*CACNEW
      CAW(L,K)=-FUHU(L,K)*CACNEW
      CAE(L,K)=FUHV(L+1,K)*CACNEW 
      CAN(L,K)=FVHV(LN,K)*CACNEW
      CAM(L,K)=-FWU(L,K-1)*DZIC(K)*CACNEW
      CAP(L,K)=FWV(L,K)*DZIC(K)*CACNEW
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  MODIFY ADVECTION COEFFICIENTS TO
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
      CAS(L,K)=0.
      CAW(L,K)=0.
      CAE(L,K)=0.
      CAN(L,K)=0.
      CAM(L,K)=0.
      CAP(L,K)=0.
      CON(L,K)=0.
C
      IF(VHDX2(LN,K).LT.0.)THEN
       HTTMP=S3TL*H2P(L)+S2TL*H1P(L)
       CACTMP=-FVHV(LN,K)-DXYP(L)*(HP(L)-HTTMP)/DELT
       CACNEW=1./( CACTMP+2.*DXYP(L)*HP(L)/DELT )
       CONT(L,K)=(CACTMP-2.*DXYP(L)*HTTMP/DELT)*CACNEW*CON1(L,K)      
       CAN(L,K)=FVHV(LN,K)*CACNEW
C      CLOS(LL,K,M)=CON(L,K)
       NLOS(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOS(LL,K,M)
       IF(NMNLO.GE.NTSCRS(LL))THEN
        CONT(L,K)=-CBT
       ELSE
        CONT(L,K)=-CLOS(LL,K,M)
     &         -(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)
      CAS(L,K)=0.
      CAW(L,K)=0.
      CAE(L,K)=0. 
      CAN(L,K)=0.
      CAM(L,K)=0.
      CAP(L,K)=0.      
C
      IF(UHDY2(L+1,K).LT.0.)THEN
       HTTMP=S3TL*H2P(L)+S2TL*H1P(L)
       CACTMP=-FUHV(L+1,K)-DXYP(L)*(HP(L)-HTTMP)/DELT
       CACNEW=1./( CACTMP+2.*DXYP(L)*HP(L)/DELT )
       CONT(L,K)=(CACTMP-2.*DXYP(L)*HTTMP/DELT)*CACNEW*CON1(L,K)      
       CAE(L,K)=FUHV(L+1,K)*CACNEW 
C      CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOW(LL,K,M)
       IF(NMNLO.GE.NTSCRW(LL))THEN
        CONT(L,K)=CBT
       ELSE
        CONT(L,K)=-CLOW(LL,K,M)
     &         -(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)
      CAS(L,K)=0.
      CAW(L,K)=0.
      CAE(L,K)=0.
      CAN(L,K)=0.
      CAM(L,K)=0.
      CAP(L,K)=0.      
C
      IF(UHDY2(L,K).GT.0.)THEN
       HTTMP=S3TL*H2P(L)+S2TL*H1P(L)
       CACTMP=FUHU(L,K)-DXYP(L)*(HP(L)-HTTMP)/DELT
       CACNEW=1./( CACTMP+2.*DXYP(L)*HP(L)/DELT )
       CONT(L,K)=(CACTMP-2.*DXYP(L)*HTTMP/DELT)*CACNEW*CON1(L,K)      
       CAW(L,K)=-FUHU(L,K)*CACNEW
C      CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOE(LL,K,M)
       IF(NMNLO.GE.NTSCRE(LL))THEN
        CONT(L,K)=-CBT
       ELSE
        CONT(L,K)=-CLOE(LL,K,M)
     &         -(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
      CAS(L,K)=0.
      CAW(L,K)=0.
      CAE(L,K)=0.
      CAN(L,K)=0.
      CAM(L,K)=0.
      CAP(L,K)=0.      
C
      IF(VHDX2(L,K).GT.0.)THEN
       HTTMP=S3TL*H2P(L)+S2TL*H1P(L)
       CACTMP=FVHU(L,K)-DXYP(L)*(HP(L)-HTTMP)/DELT
       CACNEW=1./( CACTMP+2.*DXYP(L)*HP(L)/DELT )
       CONT(L,K)=(CACTMP-2.*DXYP(L)*HTTMP/DELT)*CACNEW*CON1(L,K)      
       CAS(L,K)=-FVHU(L,K)*CACNEW
C      CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLON(LL,K,M)
       IF(NMNLO.GE.NTSCRN(LL))THEN
        CONT(L,K)=-CBT
       ELSE
        CONT(L,K)=-CLON(LL,K,M)
     &         -(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       ENDIF
      ENDIF
C
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  SET EXPLICIT TERMS
C
C----------------------------------------------------------------------C
C
      IF(KC.EQ.1)THEN
      DO L=1,LA
      LS=LSC(L)
      LN=LNC(L)
      CONT(L,K)=CONT(L,1)
     &         +CAS(L,1)*CON1(LS,1)+CAW(L,1)*CON1(L-1,1)
     &         +CAE(L,1)*CON1(L+1,1)+CAN(L,1)*CON1(LN,1)
      ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
      DO L=1,LA
      LS=LSC(L)
      LN=LNC(L)
      CONT(L,1)=CONT(L,1)
     &         +CAS(L,1)*CON1(LS,1)+CAW(L,1)*CON1(L-1,1)
     &         +CAE(L,1)*CON1(L+1,1)+CAN(L,1)*CON1(LN,1)
     &         +CAP(L,1)*CON1(L,2)
      CONT(L,KC)=CONT(L,KC)
     &         +CAS(L,KC)*CON1(LS,KC)+CAW(L,KC)*CON1(L-1,KC)
     &         +CAE(L,KC)*CON1(L+1,KC)+CAN(L,KC)*CON1(LN,KC)
     &         +CAM(L,KC)*CON1(L,KS)
      ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN       
      DO K=2,KS
      DO L=1,LA
      LS=LSC(L)
      LN=LNC(L)
      CONT(L,K)=CONT(L,K)
     &         +CAS(L,K)*CON1(LS,K)+CAW(L,K)*CON1(L-1,K)
     &         +CAE(L,K)*CON1(L+1,K)+CAN(L,K)*CON1(LN,K)
     &         +CAM(L,K)*CON1(L,K-1)+CAP(L,K)*CON1(L,K+1)
      ENDDO
      ENDDO
      ENDIF
C
      IF(ISUD.EQ.1)THEN
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=CON(L,K)
      ENDDO
      ENDDO
      ENDIF
C
C   
C**********************************************************************C
C
C **  SOLVE EQUATONS BY ZEBRA RELAXATION
C
C----------------------------------------------------------------------C
C
      ITER=1
C
  300 CONTINUE
      RSQ=0.
C
C----------------------------------------------------------------------C
C
C **  RED CELL LOOP
C
      IF(KC.EQ.1)THEN
      DO LL=1,NRC
      L=LRC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,1)=CONT(L,1)+CON(L,1)
     &      +CAS(L,1)*CON(LS,1)+CAW(L,1)*CON(L-1,1)
     &      +CAE(L,1)*CON(L+1,1)+CAN(L,1)*CON(LN,1)
      ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
      DO LL=1,NRC
      L=LRC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,1)=CONT(L,1)+CON(L,1)
     &      +CAS(L,1)*CON(LS,1)+CAW(L,1)*CON(L-1,1)
     &      +CAE(L,1)*CON(L+1,1)+CAN(L,1)*CON(LN,1)
     &      +CAP(L,1)*CON(L,2)
      RSDZ(L,KC)=CONT(L,KC)+CON(L,KC)
     &      +CAS(L,KC)*CON(LS,KC)+CAW(L,KC)*CON(L-1,KC)
     &      +CAE(L,KC)*CON(L+1,KC)+CAN(L,KC)*CON(LN,KC)
     &      +CAM(L,KC)*CON(L,KS)
      ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
      DO K=2,KS
      DO LL=1,NRC
      L=LRC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,K)=CONT(L,K)+CON(L,K)
      RSDZ(L,K)=CONT(L,K)+CON(L,K)
     &      +CAS(L,K)*CON(LS,K)+CAW(L,K)*CON(L-1,K)
     &      +CAE(L,K)*CON(L+1,K)+CAN(L,K)*CON(LN,K)
     &      +CAM(L,K)*CON(L,K-1)+CAP(L,K)*CON(L,K+1)
      ENDDO
      ENDDO
      ENDIF 
C
      IF(KC.EQ.1)THEN
      DO LL=1,NRC
      L=LRC(LL)
      RSQ=RSQ+RSDZ(L,1)*RSDZ(L,1)
      CON(L,1)=CON(L,1)-RPIA*RSDZ(L,1)
      ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
      DO K=1,KC
      DO LL=1,NRC
      L=LRC(LL)
      RSQ=RSQ+RSDZ(L,K)*RSDZ(L,K)
      RSDZ(L,K)=-RPIA*RSDZ(L,K)
      ENDDO
      ENDDO
C
      DO LL=1,NRC
      L=LRC(LL)
      CU1(L,1)=CAP(L,1)
      ENDDO
C
      DO K=2,KS
      DO LL=1,NRC
      L=LRC(LL)
      EB=1./(1.-CAM(L,K)*CU1(L,K-1))
      CU1(L,K)=CAP(L,K)*EB
      RSDZ(L,K)=(RSDZ(L,K)-CAM(L,K)*RSDZ(L,K-1))*EB
      ENDDO
      ENDDO
C
      K=KC
      DO LL=1,NRC
      L=LRC(LL)
      CMBTMP=1
      EB=1./(1.-CAM(L,K)*CU1(L,K-1))
      RSDZ(L,K)=(RSDZ(L,K)-CAM(L,K)*RSDZ(L,K-1))*EB
      ENDDO
C
      DO K=KC-1,1,-1
      DO LL=1,NRC
      L=LRC(LL)
      RSDZ(L,K)=RSDZ(L,K)-CU1(L,K)*RSDZ(L,K+1)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NRC
      L=LRC(LL)
      CON(L,K)=CON(L,K)+RSDZ(L,K)
      ENDDO
      ENDDO
C
      ENDIF 
C
C----------------------------------------------------------------------C
C
C **  BLACK CELL LOOP
C
      IF(KC.EQ.1)THEN
      DO LL=1,NBC
      L=LBC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,1)=CONT(L,1)+CON(L,1)
     &      +CAS(L,1)*CON(LS,1)+CAW(L,1)*CON(L-1,1)
     &      +CAE(L,1)*CON(L+1,1)+CAN(L,1)*CON(LN,1)
      ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
      DO LL=1,NBC
      L=LBC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,1)=CONT(L,1)+CON(L,1)
     &      +CAS(L,1)*CON(LS,1)+CAW(L,1)*CON(L-1,1)
     &      +CAE(L,1)*CON(L+1,1)+CAN(L,1)*CON(LN,1)
     &      +CAP(L,1)*CON(L,2)
      RSDZ(L,KC)=CONT(L,KC)+CON(L,KC)
     &      +CAS(L,KC)*CON(LS,KC)+CAW(L,KC)*CON(L-1,KC)
     &      +CAE(L,KC)*CON(L+1,KC)+CAN(L,KC)*CON(LN,KC)
     &      +CAM(L,KC)*CON(L,KS)
      ENDDO
      ENDIF
C
      IF(KC.GE.3)THEN
      DO K=2,KS
      DO LL=1,NBC
      L=LBC(LL)
      LN=LNC(L)
      LS=LSC(L)
      RSDZ(L,K)=CONT(L,K)+CON(L,K)
     &      +CAS(L,K)*CON(LS,K)+CAW(L,K)*CON(L-1,K)
     &      +CAE(L,K)*CON(L+1,K)+CAN(L,K)*CON(LN,K)
     &      +CAM(L,K)*CON(L,K-1)+CAP(L,K)*CON(L,K+1)
      ENDDO
      ENDDO
      ENDIF 
C
      IF(KC.EQ.1)THEN
      DO LL=1,NBC
      L=LBC(LL)
      RSQ=RSQ+RSDZ(L,1)*RSDZ(L,1)
      CON(L,1)=CON(L,1)-RPIA*RSDZ(L,1)
      ENDDO
      ENDIF
C
      IF(KC.GE.2)THEN
C
      DO K=1,KC
      DO LL=1,NBC
      L=LBC(LL)
      RSQ=RSQ+RSDZ(L,K)*RSDZ(L,K)
      RSDZ(L,K)=-RPIA*RSDZ(L,K)
      ENDDO
      ENDDO
C
      DO LL=1,NBC
      L=LBC(LL)
      CU1(L,1)=CAP(L,1)
      ENDDO
C
      DO K=2,KS
      DO LL=1,NBC
      L=LBC(LL)
      EB=1./(1.-CAM(L,K)*CU1(L,K-1))
      CU1(L,K)=CAP(L,K)*EB
      RSDZ(L,K)=(RSDZ(L,K)-CAM(L,K)*RSDZ(L,K-1))*EB
      ENDDO
      ENDDO
C
      K=KC
      DO LL=1,NBC
      L=LBC(LL)
      EB=1./(1.-CAM(L,K)*CU1(L,K-1))
      RSDZ(L,K)=(RSDZ(L,K)-CAM(L,K)*RSDZ(L,K-1))*EB
      ENDDO
C
      DO K=KC-1,1,-1
      DO LL=1,NBC
      L=LBC(LL)
      RSDZ(L,K)=RSDZ(L,K)-CU1(L,K)*RSDZ(L,K+1)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NBC
      L=LBC(LL)
      CON(L,K)=CON(L,K)+RSDZ(L,K)
      ENDDO
      ENDDO
C
      ENDIF 
C
C----------------------------------------------------------------------C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      IF(RSQ .LE. RSQMIA) GOTO 900
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF(ITER .GE. ITRMIA) GOTO 900
C  
      ITER=ITER+1
      GOTO 300
C
  900 CONTINUE
C
      WRITE(6,601)ITER,RSQ 
C
  601 FORMAT('  CALTRANI VALUES, ITER=',I5,3X,'RSQ=',E12.4)
C
C**********************************************************************C
C
C **  RECORD LAST OUTFLOWING CONCENTRATION
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
C
      IF(VHDX2(LN,K).LT.0.)THEN
       CLOS(LL,K,M)=CON(L,K)
       NLOS(LL,K,M)=N
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)
C
      IF(UHDY2(L+1,K).LT.0.)THEN
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
C
      IF(UHDY2(L,K).GT.0.)THEN
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ENDIF
C
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
C
      IF(VHDX2(L,K).GT.0.)THEN
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ENDIF
C
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
      IF(ISADAC(M).EQ.0) RETURN
      IF(ISCDCA(M).EQ.1) RETURN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)          
      UHU=ABS(UHDY2(L,K))*(CON(L,K)-CON(L-1,K))/
     &       (CON(L,K)+CON(L-1,K)+BSMALL)
      VHV=ABS(VHDX2(L,K))*(CON(L,K)-CON(LS,K))/
     &       (CON(L,K)+CON(LS,K)+BSMALL)
      FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)
     &         +MIN(UHU,0.)*CON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)
     &         +MIN(VHV,0.)*CON(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      WW=ABS(W2(L,K))*(CON(L,K+1)-CON(L,K))/
     &      (CON(L,K+1)+CON(L,K)+BSMALL)
      FWU(L,K)=MAX(WW,0.)*CON(L,K)
     &        +MIN(WW,0.)*CON(L,K+1)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DO K=1,KC
      FVHU(LN,K)=0.0
      ENDDO
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DO K=1,KC
      FUHU(L+1,K)=0.0
      ENDDO
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DO K=1,KC
      FUHU(L,K)=0.0
      ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      DO K=1,KC
      FVHU(L,K)=0.0
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
C
      IF(ISFCT(M).EQ.0) GOTO 600
C
C----------------------------------------------------------------------C
C
C **  DETERMINE MAX AND MIN CONCENTRATIONS
C
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=MAX(CON(L,K),CON2(L,K))
      ENDDO
      ENDDO
C
      DO L=2,LA
      CMAX(L,1)=MAX(CONT(L,1),CONT(L,2))
      CMAX(L,KC)=MAX(CONT(L,KS),CONT(L,KC))
      ENDDO
C
      DO K=2,KS
      DO L=2,LA
      CMAX(L,K)=MAX(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB(L)*CONT(L-1,K) 
      CEMAX=SUB(L+1)*CONT(L+1,K)
      CSMAX=SVB(L)*CONT(LS,K)
      CNMAX=SVB(LN)*CONT(LN,K)
      CMAX(L,K)=MAX(CMAX(L,K),CNMAX,CEMAX,CSMAX,CWMAX)
      ENDDO
      ENDDO
C
      DO L=2,LA
      DO K=1,KC
      CONT(L,K)=MIN(CON(L,K),CON2(L,K))
      ENDDO
      ENDDO
C
      DO L=2,LA
      CMIN(L,1)=MIN(CONT(L,1),CONT(L,2))
      CMIN(L,KC)=MIN(CONT(L,KS),CONT(L,KC))
      ENDDO
C
      DO K=2,KS
      DO L=2,LA
      CMIN(L,K)=MIN(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMIN=SUB(L)*CONT(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB(L+1)*CONT(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB(L)*CONT(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB(LN)*CONT(LN,K)+1.E+6*(1.-SVB(LN))
      CMIN(L,K)=MIN(CMIN(L,K),CNMIN,CEMIN,CSMIN,CWMIN)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  SEPARATE POSITIVE AND NEGATIVE FLUXES PUTTING NEGATIVE FLUXES 
C **  INTO FUHV, FVHV, AND FWV
C
      DO K=1,KC
      DO L=2,LA
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS AND LOAD
C **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)
     &                               +FVHU(L,K)-FVHV(LN,K))
     &                      +DZIC(K)*(FWU(L,K-1)-FWV(L,K)) )*HPI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)
     &                               +FVHU(LN,K)-FVHV(L,K))
     &                      +DZIC(K)*(FWU(L,K)-FWV(L,K-1)) )*HPI(L)
      ENDDO
      ENDDO
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DO K=1,KC
      DU(LN,K)=0.
      DV(LN,K)=0.
      ENDDO
      ENDDO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DO K=1,KC
      DU(L+1,K)=0.
      DV(L+1,K)=0.
      ENDDO
      ENDDO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DO K=1,KC
      DU(L-1,K)=0.
      DV(L-1,K)=0.
      ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      DO K=1,KC
      DU(LS,K)=0.
      DV(LS,K)=0.
      ENDDO
      ENDDO
C
      DO NS=1,NQSIJ
       L=LQS(NS)
       NQSTMP=NQSERQ(NS)
       DO K=1,KC
        QQQTMP=ABS(QSS(K,NS)+QSERT(K,NQSTMP))
        IF(QQQTMP.GE.1.E-12)THEN
          DU(L,K)=0.
          DV(L,K)=0.
C         DU(L+1,K)=0.
C         DV(LNC(L),K)=0.
        ENDIF       
       ENDDO
      ENDDO
C
      DO NCTL=1,NQCTL
       IU=IQCTLU(NCTL)
       JU=JQCTLU(NCTL)
       LU=LIJ(IU,JU)
       ID=IQCTLD(NCTL)
       JD=JQCTLD(NCTL)
       IF(ID.EQ.0.AND.JD.EQ.0)THEN
         LD=LC
        ELSE
         LD=LIJ(ID,JD)
       ENDIF
       DO K=1,KC
        QQQTMP=ABS(QCTLT(K,NCTL))
        IF(QQQTMP.GE.1.E-12)THEN
          DU(LU,K)=0.
          DV(LU,K)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
          DU(LD,K)=0.
          DV(LD,K)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
        ENDIF
       ENDDO
      ENDDO
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
        NQSTMP=NQWRSERQ(NWR)
        QQQTMP=ABS(QWR(NWR)+QWRSERT(NQSTMP))
        IF(QQQTMP.GE.1.E-12)THEN
          DU(LU,KU)=0.
          DV(LU,KU)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
          DU(LD,KD)=0.
          DV(LD,KD)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
        ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DU AND DV
C
      DO K=1,KC
      DO L=2,LA
      IF(DU(L,K).GT.0.) DU(L,K)=(CMAX(L,K)-CON(L,K))/(DU(L,K)+BSMALL)
      IF(DV(L,K).GT.0.) DV(L,K)=(CON(L,K)-CMIN(L,K))/(DV(L,K)+BSMALL)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      DU(L,K)=MIN(DU(L,K),1.)
      DV(L,K)=MIN(DV(L,K),1.)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  LIMIT FLUXES
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)
     &         +MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)
     &         +MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      ENDDO
      ENDDO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MIN(DV(L,K),DU(L,K+1))*FWU(L,K)
     &        +MIN(DU(L,K),DV(L,K+1))*FWV(L,K)
      ENDDO
      ENDDO  
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTION CALCULATION
C
  600 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      CH(L,K)=CON(L,K)*HP(L)
     &       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     &              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     &             +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  DIAGNOSE FCT SCHEME
C
      IF(ISFCT(M).EQ.2)THEN
      WRITE(6,6010)N
C
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF(CCMAX.GT.0.)THEN
       WRITE(6,6011)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      ENDIF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF(CCMIN.GT.0.)THEN
       WRITE(6,6012)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C
 6010 FORMAT('  FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT('  CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT('  CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
      RETURN
      END
