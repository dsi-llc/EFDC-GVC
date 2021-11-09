C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SETOBC2T(DELT,DELTD2,DELTI)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 02/15/2002        John Hamrick       02/11/2002       John Hamrick
C  added alternate sor equation solver relax2t
C----------------------------------------------------------------------C  
C
C
C**********************************************************************C
C
C ** SUBROUTINE SETOBC2T SETS OPEN BOUNDARY CONDITIONS FOR CALPUV2T AND
C ** CALPUV2C
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  SET OPEN BOUNDARY SURFACE ELEVATIONS 
C
C----------------------------------------------------------------------C
C
      IF(ISDYNSTP.EQ.0)THEN
        TN=DT*FLOAT(N)+TCON*TBEGIN
      ELSE
        TN=TIMESEC
      ENDIF
C
      DO M=1,MTIDE
      TM=MOD(TN,TCP(M))
      TM=PI2*TM/TCP(M)
      CCCOS(M)=COS(TM)
      SSSIN(M)=SIN(TM)
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBW
      L=LPBW(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      FP1(L)=PSERT(NPSERW(LL))+0.5*PSERZDF(NPSERW(LL))
     &      +PSERST(NPSERW(LL))+0.5*PSERZDS(NPSERW(LL))
	IF(NPFORT.GE.1.AND.NPSERW1(LL).GT.0)THEN
        TMPVAL=PSERT(NPSERW1(LL))+0.5*PSERZDF(NPSERW1(LL))
     &      +PSERST(NPSERW1(LL))+0.5*PSERZDS(NPSERW1(LL))
        FP1(L)=FP1(L)+TPCOORDW(LL)*(TMPVAL-FP1(L))
      ENDIF
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
      ENDDO
      CET=0.5*DELTD2*G*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HUTMP(L+1))*DXIU(L+1)
C       TMP=DELTD2*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/TMP
       CE(L)=-CET
       FP1(L)=CET*(2.*FP1(L)
     & -SQRT(G*HUTMP(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/TMP
C     & -SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/TMP
      ELSE
       FP1(L+1)=CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.      
      FP1(L)=PSERT(NPSERE(LL))+0.5*PSERZDF(NPSERE(LL))
     &      +PSERST(NPSERE(LL))+0.5*PSERZDS(NPSERE(LL))
	IF(NPFORT.GE.1.AND.NPSERE1(LL).GT.0)THEN
        TMPVAL=PSERT(NPSERE1(LL))+0.5*PSERZDF(NPSERE1(LL))
     &      +PSERST(NPSERE1(LL))+0.5*PSERZDS(NPSERE1(LL))
        FP1(L)=FP1(L)+TPCOORDE(LL)*(TMPVAL-FP1(L))
      ENDIF
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
      ENDDO
      CWT=0.5*DELTD2*G*HRUO(L  )*RCX(L)*HUTMP(L  )
      IF(ISPBE(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HUTMP(L))*DXIU(L)
C       TMP=DELTD2*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/TMP
       CW(L)=-CWT
       FP1(L)=CWT*(2.*FP1(L)
     & +SQRT(G*HUTMP(L))*FUHDYE(L)*DYIU(L)*HUI(L))/TMP
C     & +SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/TMP
      ELSE
       FP1(L-1)=CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LN=LNC(L)
      FP1(L)=PSERT(NPSERS(LL))+0.5*PSERZDF(NPSERS(LL))
     &      +PSERST(NPSERS(LL))+0.5*PSERZDS(NPSERS(LL))
	IF(NPFORT.GE.1.AND.NPSERS1(LL).GT.0)THEN
        TMPVAL=PSERT(NPSERS1(LL))+0.5*PSERZDF(NPSERS1(LL))
     &      +PSERST(NPSERS1(LL))+0.5*PSERZDS(NPSERS1(LL))
        FP1(L)=FP1(L)+TPCOORDS(LL)*(TMPVAL-FP1(L))
      ENDIF
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
      ENDDO
      CNT=0.5*DELTD2*G*HRVO(LN )*RCY(LN)*HVTMP(LN )
      IF(ISPBS(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HVTMP(LN))*DYIV(LN)
C       TMP=DELTD2*SQRT(G*HMV(LN))*DYIV(LN)
       CC(L)=CNT*(1.+TMP)/TMP
       CN(L)=-CNT
       FP1(L)=CNT*(2.*FP1(L)
     & -SQRT(G*HVTMP(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/TMP
C     & -SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/TMP
      ELSE
       FP1(LN)=CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LS=LSC(L)
      FP1(L)=PSERT(NPSERN(LL))+0.5*PSERZDF(NPSERN(LL))
     &      +PSERST(NPSERN(LL))+0.5*PSERZDS(NPSERN(LL))
	IF(NPFORT.GE.1.AND.NPSERN1(LL).GT.0)THEN
        TMPVAL=PSERT(NPSERN1(LL))+0.5*PSERZDF(NPSERN1(LL))
     &      +PSERST(NPSERN1(LL))+0.5*PSERZDS(NPSERN1(LL))
        FP1(L)=FP1(L)+TPCOORDN(LL)*(TMPVAL-FP1(L))
      ENDIF
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
      ENDDO
      CST=0.5*DELTD2*G*HRVO(L  )*RCY(L)*HVTMP(L  )
      IF(ISPBN(LL).GE.1)THEN
       TMP=DELTD2*SQRT(G*HVTMP(L))*DYIV(L)
C       TMP=DELTD2*SQRT(G*HMV(L))*DYIV(L)
       CC(L)=CST*(1.+TMP)/TMP
       CS(L)=-CST
       FP1(L)=CST*(2.*FP1(L)
     & +SQRT(G*HVTMP(L))*FVHDXE(L)*DXIV(L)*HVI(L))/TMP
C     & +SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/TMP
      ELSE
       FP1(LS)=CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      ENDIF
      ENDDO
C
C**********************************************************************C
C
      RETURN
      END
