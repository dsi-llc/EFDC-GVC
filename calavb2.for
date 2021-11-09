C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALAVB2 (ISTL)
C
C **  SUBROUTINE CALAV CALCULATES VERTICAL VISCOSITY AND DIFFUSIVITY
C **  USING GLAPERIN ET AL'S MODIFICATION OF THE MELLOR-YAMADA MODEL
C **  (NOTE AV, AB, AND AQ ARE ACTUALLY DIVIDED BY H)
C **  IF ISGA=1 VALUES ARE GEOMETRIC AVERAGES WITH THE PREVIOUS VALUES
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
C     DIMENSION QQI(LCM)
C
      DATA ATURB1,ATURB2,TURBC1/0.92,0.74,0.08/
C
C**********************************************************************C 
C
C   SHTOP    =      0.4939
C   SHBOT    =     34.6764
C   SMTOP1   =      0.3933
C   SMTOP2   =      7.8464
C   SMBOT1   =     34.6764
C   SMBOT2   =      6.1272
C   RLIMIT   =      0.0233
C   SHMIN    =      0.0934
C   SMMIN    =      0.1099
C   SHMAX    =      5.2073
C   SMMAX    =      4.9639
C
      AVMAX=AVO
      ABMAX=ABO
      AVMIN=10.
      ABMIN=10.
C     RIQMIN=-1./44.
C030705      RIQMIN=-0.023
      RAVBTMP=1.
      IF(ISAVBMN.GE.1) RAVBTMP=0.
C
C----------------------------------------------------------------------C
C
C
      DO K=1,KS
       DO L=2,LA
         CTURBB1(L,K)=CTURB
         CTURBB2(L,K)=CTURB2B
       ENDDO
      ENDDO
C
      DO K=1,KS
       DO L=2,LA
         DELBTMP=(B(L,K+1)-B(L,K))*DZIG(K)
         RITMP=-GP*HP(L)*DELBTMP/QQ(L,K)
         IF(RITMP.GT.0.)THEN
           RIQ=DML(L,K)*DML(L,K)*RITMP
           BFUN=EXP(-3.11*RIQ)
           CTURBB1(L,K)=CTURB/(BFUN+1.E-16)
           IF(BBT(L,K).GT.0.)THEN
             TMPVAL=DELBTMP*DELBTMP/(RITMP*BBT(L,K))
             CTURBB2(L,K)=CTURB2B/(1.+0.61*(1.-BFUN)*TMPVAL)
           ENDIF
         ENDIF
       ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C      IF(ISTL.EQ.3)THEN
      IF(ISFAVB.EQ.0)THEN
C
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &    *(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      AB(L,K)=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AV(L,K)=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AV(L,K))
       ABMAX=MAX(ABMAX,AB(L,K))
       AVMIN=MIN(AVMIN,AV(L,K))
       ABMIN=MIN(ABMIN,AB(L,K))
      AV(L,K)=AV(L,K)*HPI(L)
      AB(L,K)=SCB(L)*AB(L,K)*HPI(L)
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
      IF(ISFAVB.EQ.1)THEN
C
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &    *(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=0.5*(AV(L,K)+AVTMP*HPI(L))
      AB(L,K)=SCB(L)*0.5*(AB(L,K)+ABTMP*HPI(L))
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
      IF(ISFAVB.EQ.2)THEN
C
      DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=LF+LDM-1
      DO K=1,KS
       DO L=LF,LL
       QQI(L)=1./QQ(L,K)
       ENDDO
      DO L=LF,LL
      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &    *(B(L,K+1)-B(L,K))*QQI(L)
      RIQ=MAX(RIQ,RIQMIN)
      TMPVAL1=1.-( 6.*ATURB1/CTURBB1(L,K) )
      SBTOP=ATURB2*TMPVAL1
      SBBOT=3.*ATURB2*( 6.*ATURB1+CTURBB2(L,K) )
      SVTOP=ATURB1*( TMPVAL1-3.*TURBC1)
      TMPVAL2=TMPVAL1*( CTURBB2(L,K)-3.*ATURB2 )
      TMPVAL3=-3.*TURBC1*( CTURBB2(L,K)+6.*ATURB1)
      SVTOP2=3.*ATURB2*(TMPVAL2+TMPVAL3)/SVTOP
      SVBOT=9.*ATURB1*ATURB2
      SFAV=SVTOP*(1.+SVTOP2*RIQ)/((1.+SVBOT*RIQ)*(1.+SBBOT*RIQ))
      SFAB=SBTOP/(1.+SBBOT*RIQ)
      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
       AVMAX=MAX(AVMAX,AVTMP)
       ABMAX=MAX(ABMAX,ABTMP)
       AVMIN=MIN(AVMIN,AVTMP)
       ABMIN=MIN(ABMIN,ABTMP)
      AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))
      AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C      ENDIF
C
C      IF(ISTL.EQ.2)THEN
C
C      DO ND=1,NDM
C      LF=2+(ND-1)*LDM
C      LL=LF+LDM-1
C      DO K=1,KS
C       DO L=LF,LL
C       QQI(L)=1./(QQMIN+QQ(L,K))
C       ENDDO
C      DO L=LF,LL
C      RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
C     &    *(B(L,K+1)-B(L,K))*QQI(L)
C      RIQ=MAX(RIQ,RIQMIN)
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
C      SFAB=0.5/(1.+36.*RIQ)
C      ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
C      AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
C       AVMAX=MAX(AVMAX,AVTMP)
C       ABMAX=MAX(ABMAX,ABTMP)
C       AVMIN=MIN(AVMIN,AVTMP)
C       ABMIN=MIN(ABMIN,ABTMP)
C      AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))
C      AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))
C      ENDDO
C      ENDDO
C      ENDDO
C
C      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISAVBMN.GE.1)THEN
        DO K=1,KS
        DO L=2,LA
         AVTMP=AVMN*HPI(L)
         ABTMP=ABMN*HPI(L)
         AV(L,K)=MAX(AV(L,K),AVTMP)
         AB(L,K)=MAX(AB(L,K),ABTMP)
        ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)      
      AVUI(L,K)=2./(AV(L,K)+AV(L-1,K))
      AVVI(L,K)=2./(AV(L,K)+AV(LS,K))
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
      DO K=2,KS
      DO L=2,LA
C      AQ(L,K)=0.255*(AV(L,K-1)+AV(L,K))
      AQ(L,K)=0.205*(AV(L,K-1)+AV(L,K))
      ENDDO
      ENDDO
C
      DO L=2,LA
C      AQ(L,1)=0.255*AV(L,1)
C      AQ(L,KC)=0.255*AV(L,KS)
      AQ(L,1)=0.205*AV(L,1)
      AQ(L,KC)=0.205*AV(L,KS)
      ENDDO
C
      ELSE
C
      DO K=2,KS
      DO L=2,LA
C      AQTMP=0.255*(AV(L,K-1)+AV(L,K))
      AQTMP=0.205*(AV(L,K-1)+AV(L,K))
C      AQ(L,K)=SQRT(AQ(L,K)*AQTMP)
      AQ(L,K)=AQTMP
      ENDDO
      ENDDO
C
      DO L=2,LA
C      AQTMP=0.255*AV(L,1)
      AQTMP=0.205*AV(L,1)
C      AQ(L,1)=SQRT(AQ(L,1)*AQTMP)
      AQ(L,1)=AQTMP
C      AQTMP=0.255*AV(L,KS)
      AQTMP=0.205*AV(L,KS)
C      AQ(L,KC)=SQRT(AQ(L,KC)*AQTMP)
      AQ(L,KC)=AQTMP
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
      RETURN
      END
