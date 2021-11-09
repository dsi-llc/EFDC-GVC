c
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALAVBGVC(ISTL)
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
C 03/19/2002        John Hamrick       03/19/2002       John Hamrick
C  added drycell bypass and consistent initialization of dry values
C----------------------------------------------------------------------C
C
C**********************************************************************C 
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C     DIMENSION QQI(LCM)
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
c
c  GALPERIN
C
      IF(ISTOPT(0).EQ.1)THEN
C
      SFAV0= 0.392010
	SFAV1= 7.760050
	SFAV2=34.676440
	SFAV3= 6.127200
	SFAB0= 0.493928
	SFAB1=34.676440
	RIQMIN=-0.999/SFAB1
C
      ENDIF
C
C  KANTHA AND CLAYSON (1994)
C
      IF(ISTOPT(0).EQ.2)THEN

      SFAV0= 0.392010
	SFAV1= 8.679790
	SFAV2=30.192000
	SFAV3= 6.127200
	SFAB0= 0.493928
	SFAB1=30.192000
	RIQMIN=-0.999/SFAB1
C
      ENDIF
C
C  KANTHA (2003)
C
      IF(ISTOPT(0).EQ.3)THEN
C
      SFAV0= 0.392010
	SFAV1=14.509100
	SFAV2=24.388300
	SFAV3= 3.236400
	SFAB0= 0.490025
	SFAB1=24.388300
	RIQMIN=-0.999/SFAB1
C
      ENDIF
C
      QQIMAX=1./QQMIN
      AVMAX=AVO
      ABMAX=ABO
      AVMIN=10.
      ABMIN=10.
C     RIQMIN=-1./44.
c      RIQMIN=-0.023
c      RIQMAX=0.28
      RAVBTMP=1.
      IF(ISAVBMN.GE.1) RAVBTMP=0.
C
      DO K=1,KC
	DO L=1,LC
	IF(IMASKDRY(L).EQ.1)THEN
      AV(L,K)=AVO*GVCSCLPI(L)*HPI(L)
      AB(L,K)=ABO*GVCSCLPI(L)*HPI(L)
      ENDIF
	ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C      IF(ISTL.EQ.3)THEN
      IF(ISFAVB.EQ.0.OR.N.EQ.1)THEN
C
      DO K=1,KS
        DO L=2,LA
	    IF(LMASKDRY(L))THEN
	      IF(LGVCP(L,K))THEN
C             QQI(L)=1./(QQMIN+QQ(L,K))
              QQI(L)=1./QQ(L,K)
              QQI(L)=MIN(QQI(L),QQIMAX)
            ENDIF
	    ENDIF
        ENDDO
        DO L=2,LA
	    IF(LMASKDRY(L))THEN
	    IF(LGVCP(L,K))THEN
            RIQ=-GP*GVCSCLPI(L)*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &          *(B(L,K+1)-B(L,K))*QQI(L)
            RIQ=MAX(RIQ,RIQMIN)
            IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
C      SFAB=0.5/(1.+36.*RIQ)
c      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
c      SFAB=0.4939/(1.+34.6764*RIQ)
c sfav and sfab above replace by kanatha and clayson
c
c            SFAV=0.3920*(1.+8.6736*RIQ)
c     &          /((1.+30.192*RIQ)*(1.+6.1272*RIQ))
c            SFAB=0.4939/(1.+30.192*RIQ)
c
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)
     &                *(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
C
            AB(L,K)=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
            AV(L,K)=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
            AVMAX=MAX(AVMAX,AV(L,K))
            ABMAX=MAX(ABMAX,AB(L,K))
            AVMIN=MIN(AVMIN,AV(L,K))
            ABMIN=MIN(ABMIN,AB(L,K))
            AV(L,K)=AV(L,K)*GVCSCLPI(L)*HPI(L)
            AB(L,K)=SCB(L)*AB(L,K)*GVCSCLPI(L)*HPI(L)
          ELSE
            AV(L,K)=0.0
            AB(L,K)=0.0
          ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISFAVB.EQ.1.AND.N.GT.1)THEN
C
      DO K=1,KS
        DO L=2,LA
	    IF(LMASKDRY(L))THEN
	      IF(LGVCP(L,K))THEN
C             QQI(L)=1./(QQMIN+QQ(L,K))
              QQI(L)=1./QQ(L,K)
              QQI(L)=MIN(QQI(L),QQIMAX)
	      ENDIF
	    ENDIF
        ENDDO
        DO L=2,LA
          IF(LMASKDRY(L))THEN
	    IF(LGVCP(L,K))THEN
            RIQ=-GP*GVCSCLPI(L)*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &         *(B(L,K+1)-B(L,K))*QQI(L)
            RIQ=MAX(RIQ,RIQMIN)
            IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
C      SFAB=0.5/(1.+36.*RIQ)
c      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
c      SFAB=0.4939/(1.+34.6764*RIQ)
c sfav and sfab above replace by kanatha and clayson
c
c            SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)
c     &          *(1.+6.1272*RIQ))
c            SFAB=0.4939/(1.+30.192*RIQ)
c
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)
     &                *(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
C
            ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
            AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
            AVMAX=MAX(AVMAX,AVTMP)
            ABMAX=MAX(ABMAX,ABTMP)
            AVMIN=MIN(AVMIN,AVTMP)
            ABMIN=MIN(ABMIN,ABTMP)
            AV(L,K)=0.5*(AV(L,K)+AVTMP*GVCSCLPI(L)*HPI(L))
            AB(L,K)=SCB(L)*0.5*(AB(L,K)+GVCSCLPI(L)*ABTMP*HPI(L))
          ELSE
            AV(L,K)=0.0
            AB(L,K)=0.0
          ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISFAVB.EQ.2.AND.N.GT.1)THEN
C
      DO K=1,KS
        DO L=2,LA
	    IF(LMASKDRY(L))THEN
	      IF(LGVCP(L,K))THEN
C             QQI(L)=1./(QQMIN+QQ(L,K))
              QQI(L)=1./QQ(L,K)
              QQI(L)=MIN(QQI(L),QQIMAX)
	      ENDIF
	    ENDIF
        ENDDO
        DO L=2,LA
	    IF(LMASKDRY(L))THEN
	    IF(LGVCP(L,K))THEN
            RIQ=-GP*GVCSCLPI(L)*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)
     &         *(B(L,K+1)-B(L,K))*QQI(L)
            RIQ=MAX(RIQ,RIQMIN)
            IF(ISLLIM.GE.1) RIQ=MIN(RIQ,RIQMAX)
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))
C      SFAB=0.5/(1.+36.*RIQ)
c      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))
c      SFAB=0.4939/(1.+34.6764*RIQ)
c sfav and sfab above replace by kanatha and clayson
c
c            SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)
c     &          *(1.+6.1272*RIQ))
c            SFAB=0.4939/(1.+30.192*RIQ)
c
            SFAV=SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)
     &                *(1.+SFAV3*RIQ))
            SFAB=SFAB0/(1.+SFAB1*RIQ)
C
            ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*ABO
            AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*SQRT(QQ(L,K))+RAVBTMP*AVO
            AVMAX=MAX(AVMAX,AVTMP)
            ABMAX=MAX(ABMAX,ABTMP)
            AVMIN=MIN(AVMIN,AVTMP)
            ABMIN=MIN(ABMIN,ABTMP)
	      STUFF=AV(L,K)*AVTMP*GVCSCLPI(L)*HPI(L)
	      IF(STUFF.LT.0.0)
     &      WRITE(6,6999)IL(L),JL(L),K,AV(L,K),AVTMP,GVCSCLPI(L),
     &                        HPI(L),HP(L)
            AV(L,K)=SQRT(AV(L,K)*AVTMP*GVCSCLPI(L)*HPI(L))
            AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*GVCSCLPI(L)*HPI(L))
          ELSE
            AV(L,K)=0.
            AB(L,K)=0.
          ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      ENDIF
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
CGVCDIAG
C      IF(N.LE.2)THEN
C      DO L=2,LA
C	  WRITE(8,881)N,IL(L),JL(L),(AV(L,K),K=1,KS)
C      ENDDO
C      ENDIF
C
C  881 FORMAT('AV GVC ',3I5,10E14.6)
CGVCDIAG
c----------------------------------------------------------------------C
C
      IF(ISAVBMN.GE.1)THEN
        DO K=1,KS
          DO L=2,LA
	      IF(LGVCP(L,K))THEN
              AVTMP=AVMN*HPI(L)
              ABTMP=ABMN*HPI(L)
              AV(L,K)=MAX(AV(L,K),AVTMP)
              AB(L,K)=MAX(AB(L,K),ABTMP)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C
C **  INTERPOLATE TO U AND V POINTS
C
      DO K=1,KS
        DO L=2,LA
	    AVUI(L,K)=0.0
	    AVVI(L,K)=0.0
        ENDDO
      ENDDO
C
      DO K=1,KS
        DO L=2,LA
	    IF(LGVCU(L,K)) AVUI(L,K)=2./(AV(L,K)+AV(L-1,K))
        ENDDO
      ENDDO
C
      DO K=1,KS
        DO L=2,LA
          LS=LSC(L)      
          IF(LGVCV(L,K)) AVVI(L,K)=2./(AV(L,K)+AV(LS,K))
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C      IF(ISTL.EQ.3)THEN
C
C
      IF(ISTOPT(0).LE.2.OR.ISSQL.EQ.0)THEN
C
      DO K=2,KS
        DO L=2,LA
	    IF(LGVCP(L,K)) THEN
C           AQ(L,K)=0.255*(AV(L,K-1)+AV(L,K))
            AQ(L,K)=0.205*(AV(L,K-1)+AV(L,K))
	    ELSE
            AQ(L,K)=0.0
          ENDIF
	  ENDDO
      ENDDO
C
      DO L=2,LA
C      AQ(L,1)=0.255*AV(L,1)
C      AQ(L,KC)=0.255*AV(L,KS)
      IF(KGVCP(L).LT.KC) AQ(L,KGVCP(L))=0.205*AV(L,KGVCP(L))
      AQ(L,KC)=0.205*AV(L,KS)
      ENDDO
C
      ENDIF
C
      IF(ISTOPT(0).EQ.3.OR.ISSQL.EQ.1)THEN
C
C     COEF = 0.5*0.628
C
      IF(ISFAVB.EQ.1.AND.N.GT.1)THEN
C
      DO K=2,KS
      DO L=2,LA
	IF(LGVCP(L,K)) THEN
        AQ(L,K)=0.314*( DML(L,K-1)*SQRT(QQ(L,K-1))
     &                     +DML(L,K  )*SQRT(QQ(L,K  )) )
     &       +RAVBTMP*AVO*HPI(L)
	ELSE
        AQ(L,K)=0.
	ENDIF
      ENDDO
      ENDDO
C
      DO L=2,LA
      IF(KGVCP(L).LT.KC) 
     & AQ(L,KGVCP(L))=0.314*DML(L,KGVCP(L))**SQRT(QQ(L,KGVCP(L)))
     &       +RAVBTMP*AVO*HPI(L)
      AQ(L,KC)=0.314*DML(L,KC)*SQRT(QQ(L,KC))
     &       +RAVBTMP*AVO*HPI(L)
      ENDDO
C
      ELSE
C
      DO K=2,KS
      DO L=2,LA
	IF(LGVCP(L,K)) THEN
        AQ(L,K)=0.5*( 0.314*( DML(L,K-1)*SQRT(QQ(L,K-1))
     &                     +DML(L,K  )*SQRT(QQ(L,K  )) )
     &       +RAVBTMP*AVO*HPI(L) ) +0.5*AQ(L,K)
	ELSE
        AQ(L,K)=0.
	ENDIF
      ENDDO
      ENDDO
C
      DO L=2,LA
      IF(KGVCP(L).LT.KC) 
     & AQ(L,KGVCP(L))=0.5*( 0.314*DML(L,KGVCP(L))**SQRT(QQ(L,KGVCP(L)))
     &       +RAVBTMP*AVO*HPI(L) ) +0.5*AQ(L,KGVCP(L))
      AQ(L,KC)=0.5*( 0.314*DML(L,KC)*SQRT(QQ(L,KC))
     &       +RAVBTMP*AVO*HPI(L) ) +0.5*AQ(L,KC)
      ENDDO
C
      ENDIF
C
      ENDIF
C
C
C      ELSE
C
C      DO K=2,KS
C      DO L=2,LA
C      AQTMP=0.255*(AV(L,K-1)+AV(L,K))
C      AQTMP=0.205*(AV(L,K-1)+AV(L,K))
C      AQ(L,K)=SQRT(AQ(L,K)*AQTMP)
C      AQ(L,K)=AQTMP
C      ENDDO
C      ENDDO
C
C      DO L=2,LA
C      AQTMP=0.255*AV(L,1)
C      AQTMP=0.205*AV(L,1)
C      AQ(L,1)=SQRT(AQ(L,1)*AQTMP)
C      AQ(L,1)=AQTMP
C      AQTMP=0.255*AV(L,KS)
C      AQTMP=0.205*AV(L,KS)
C      AQ(L,KC)=SQRT(AQ(L,KC)*AQTMP)
C      AQ(L,KC)=AQTMP
C      ENDDO
C
C      ENDIF
C
  800 FORMAT(2I6,6F10.3)
 6999 FORMAT('NEG VIS ',3I5,5E14.6)

C
C----------------------------------------------------------------------C
C
      RETURN
      END
