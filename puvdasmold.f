C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE PUVDASM(ISTL,ICALL)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C **  PERFORMCE WATER SURFACE ELEVATION AND HORIZONTAL VELOCITY DATA
C **  ASSIMILATION
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
      DIMENSION UROTTMP(KCM),VROTTMP(KCM)
C
C**********************************************************************C
C

      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT2
        DELTD2=DT
        IF(ISTL.EQ.2)THEN
          DELT=DT
          DELTD2=0.5*DT
        ENDIF
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      ENDIF
C
      NDAYA=MOD(N,NTSPTC)
      NDAYA=1+(N-NDAYA)/NTSPTC
	RVAL=-998.0
C
C**********************************************************************C
C
C **  INTERPOLATE INPUT TIME SERIES OF OBSERVED HORIZONTAL VELOCITIES
C
C----------------------------------------------------------------------C
C
      IF(ICALL.EQ.1)THEN
      IF(NUVSER.GT.0)THEN
C
      DO K=1,KC
        USERT(K,0)=0.
        VSERT(K,0)=0.
      ENDDO
C
      DO NS=1,NUVSER
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)/TCUVSER(NS)+TBEGIN*(TCON/TCUVSER(NS))
      ELSE
        TIME=TIMESEC/TCUVSER(NS)
      ENDIF
C
      M1=MUVTLAST(NS)
  100 CONTINUE
      M2=M1+1
      IF(TIME.GT.TUVSER(M2,NS))THEN
       M1=M2
       GOTO 100
      ELSE
       MUVTLAST(NS)=M1
      ENDIF      
C
      TDIFF=TUVSER(M2,NS)-TUVSER(M1,NS)
      WTM1=(TUVSER(M2,NS)-TIME)/TDIFF
      WTM2=(TIME-TUVSER(M1,NS))/TDIFF
      DO K=1,KC
        USERT(K,NS)=WTM1*USER(M1,K,NS)+WTM2*USER(M2,K,NS)
	  IF(USER(M1,K,NS).LT.RVAL)USERT(K,NS)=USER(M1,K,NS)
 	  IF(USER(M2,K,NS).LT.RVAL)USERT(K,NS)=USER(M2,K,NS)
        VSERT(K,NS)=WTM1*VSER(M1,K,NS)+WTM2*VSER(M2,K,NS)
	  IF(VSER(M1,K,NS).LT.RVAL)VSERT(K,NS)=VSER(M1,K,NS)
	  IF(VSER(M2,K,NS).LT.RVAL)VSERT(K,NS)=VSER(M2,K,NS)
	ENDDO
C
C      WRITE(7,*)TIME,NS,USER(M1,KC,NS),USER(M2,KC,NS),
C     &                  VSER(M1,KC,NS),VSER(M2,KC,NS)
C
      ENDDO
C
      ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION DATA ASSIMILATION
C
C     CALCULATIONS A VOLUME SOURCE OR SINK NECESSARY TO FORCE THE 
C     COMPUTE WATER SURFACE ELEVATION TOWARD THE OBSERVED WATER SURFACE
C     ELEVATION. THE VOLUME SOURCE OR SINK IS INSERTED INTO THE
C     INTERNAL AND EXTERNAL CONTINUITY EQUATIONS ON THE NEXT TIME STEP
C     THE COMPUTED SOURCE SINK FLOW IS WRITTEN TO VOLUME SOURCE
C     TIME SERIES OUTPUT AND A DAILY SUMMARY OUTPUT FILE.
C
C----------------------------------------------------------------------C
C
      IF(ICALL.EQ.1)THEN
      IF(ISWSEDA.GT.0)THEN
C
        DO NL=1,NLWSEDA
C
          I=ICWSEDA(NL)
          J=JCWSEDA(NL)
          L=LIJ(I,J)
          TFACTOR=1.0
          IF(IGRIDV.EQ.1) TFACTOR=GVCSCLP(L)
          NS=NWSESERA(NL)
          WSEMODEL=HP(L)+BELV(L)
          WSEOBSER=GI*PSERT(NS)
          IF(WSEOBSER.GT.RVAL)THEN
            QSSWSE=DELTI*TSWSEDA(NL)*DXYP(L)*(WSEOBSER-WSEMODEL)
            DO K=1,KC
              QWSEDA(L,K)=DZC(K)*TFACTOR*QSSWSE
            ENDDO
            QWSEASM(NDAYA,NL)=QWSEASM(NDAYA,NL)+QSSWSE
          ELSE
            DO K=1,KC
              QWSEDA(L,K)=0.0
            ENDDO
          ENDIF
C
        ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  HORIZONTAL VELOCITY DATA ASSIMILATION
C
C     CALCULATIONS VECTOR BODY FORCES (HAVING KINEMATIC STRESS UNITS)
C     NECESSARY TO FORCE THE THE COMPUTED HORIZONTAL VELOCITY COMPONENTS
C     TOWARD THE OBSERVED HORIZONTAL VELOCITY COMPONENTS.  THE
C     BODY FORCES ARE INSERTED INTO THE EXPLICIT MOMENTUM EQUATION
C     TERMS (CALEXP, CALEXP2T, CALEXPGVC) ON THE NEXT TIME STEP
C     THE COMPUTED BODY FORCES (STRESSES) ARE WRITTEN TO THE EXTERNAL 
C     MODE VELOCITY TIME SERIES OUTPUT FILES
C
C----------------------------------------------------------------------C
C
C ** ASSIMILATE ALL LAYERS
C
      IF(ICALL.EQ.1)THEN
      IF(ISUVDA.GE.1)THEN
C
      DO NL=1,NLUVDA
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
        NS=NUVSERA(NL)
	  DO K=1,KC
	    IF(USERT(K,NS).GT.RVAL) UROTTMP(K)=USERT(K,NS)
	    IF(VSERT(K,NS).GT.RVAL) VROTTMP(K)=VSERT(K,NS)
        ENDDO
	  DO K=1,KC
	    IF(USERT(K,NS).GT.RVAL) 
     &      USERT(K,NS)=WINDSXX(L)*UROTTMP(K)+WINDSXY(L)*VROTTMP(K)
	    IF(VSERT(K,NS).GT.RVAL) 
     &      VSERT(K,NS)=WINDSYX(L)*UROTTMP(K)+WINDSYY(L)*VROTTMP(K)
        ENDDO
      ENDDO
C
      DO NL=1,NLUVDA
C
        TMPIMP=1.-FSUVDA(NL)
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
	  LN=LNC(L)
	  LE=L+1
        NS=NUVSERA(NL)
        IF(ISTL.EQ.2)THEN
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUVDA(NL)*HU(L)*(USERT(K,NS)
     &                                          -TMPIMP*U(L,K))
              FBODYFXI(L)=FSUVDA(NL)*TSUVDA(NL)
 	        FBODYFX(LE,K)=DELTI*TSUVDA(NL)*HU(LE)*(USERT(K,NS)
     &                                          -TMPIMP*U(LE,K))
              FBODYFXI(LE)=FSUVDA(NL)*TSUVDA(NL)
            ELSE
	        FBODYFX(L,K)=0.0
              FBODYFXI(L)=0.0
	        FBODYFX(LE,K)=0.0
              FBODYFXI(LE)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSUVDA(NL)*HV(L)*(VSERT(K,NS)
     &                                          -TMPIMP*V(L,K))
              FBODYFYI(L)=FSUVDA(NL)*TSUVDA(NL)
 	        FBODYFY(LN,K)=DELTI*TSUVDA(NL)*HV(LN)*(VSERT(K,NS)
     &                                          -TMPIMP*V(LN,K))
              FBODYFYI(LN)=FSUVDA(NL)*TSUVDA(NL)
            ELSE
	        FBODYFY(L,K)=0.0
              FBODYFYI(L)=0.0
	        FBODYFY(LN,K)=0.0
              FBODYFYI(LN)=0.0
            ENDIF
          ENDDO
        ELSE
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUVDA(NL)*HU(L)*(USERT(K,NS)
     &                                          -TMPIMP*U1(L,K))
              FBODYFXI(L)=FSUVDA(NL)*TSUVDA(NL)
 	        FBODYFX(LE,K)=DELTI*TSUVDA(NL)*HU(LE)*(USERT(K,NS)
     &                                          -TMPIMP*U1(LE,K))
              FBODYFXI(LE)=FSUVDA(NL)*TSUVDA(NL)
            ELSE
	        FBODYFX(L,K)=0.0
              FBODYFXI(L)=0.0
	        FBODYFX(LE,K)=0.0
              FBODYFXI(LE)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSUVDA(NL)*HV(L)*(VSERT(K,NS)
     &                                          -TMPIMP*V1(L,K))
              FBODYFYI(L)=FSUVDA(NL)*TSUVDA(NL)
 	        FBODYFY(LN,K)=DELTI*TSUVDA(NL)*HV(LN)*(VSERT(K,NS)
     &                                          -TMPIMP*V1(LN,K))
              FBODYFYI(LN)=FSUVDA(NL)*TSUVDA(NL)
            ELSE
	        FBODYFY(L,K)=0.0
              FBODYFYI(L)=0.0
	        FBODYFY(LN,K)=0.0
              FBODYFYI(LN)=0.0
            ENDIF
          ENDDO
        ENDIF
C
        ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** CONSTRAIN ASSIMILATION TO DEPTH AVERAGE
C
      IF(ICALL.EQ.1)THEN
      IF(ISUVDA.EQ.1)THEN
C
      IF(IGRIDV.EQ.0)THEN
        DO NL=1,NLUVDA
          I=ICUVDA(NL)
          J=JCUVDA(NL)
          L=LIJ(I,J)
          LN=LNC(L)
          LE=L+1
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*FBODYFY(L,K)
          ENDDO
          DO K=1,KC
            FBODYFX(L,K)=AVGX
	      FBODYFY(L,K)=AVGY
          ENDDO
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*FBODYFX(LE,K)
	      AVGY=AVGY+DZC(K)*FBODYFY(LN,K)
          ENDDO
          DO K=1,KC
            FBODYFX(LE,K)=AVGX
	      FBODYFY(LN,K)=AVGY
          ENDDO
         ENDDO
      ENDIF
C
      IF(IGRIDV.EQ.1)THEN
        DO NL=1,NLUVDA
          I=ICUVDA(NL)
          J=JCUVDA(NL)
          L=LIJ(I,J)
          LN=LNC(L)
          LE=L+1
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*GVCSCLU(L)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*GVCSCLV(L)*FBODYFY(L,K)
          ENDDO
          DO K=KGVCU(L),KC
            FBODYFX(L,K)=AVGX
          ENDDO
          DO K=KGVCV(L),KC
	      FBODYFY(L,K)=AVGY
          ENDDO
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*GVCSCLU(LE)*FBODYFX(LE,K)
	      AVGY=AVGY+DZC(K)*GVCSCLV(LN)*FBODYFY(LN,K)
          ENDDO
          DO K=KGVCU(LE),KC
            FBODYFX(LE,K)=AVGX
          ENDDO
          DO K=KGVCV(LN),KC
	      FBODYFY(LN,K)=AVGY
          ENDDO
        ENDDO
      ENDIF
C
      ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  HORIZONTAL VELOCITY DATA ASSIMILATION
C
C      DIAGNOSTIC FOR OUTPUT
C
C----------------------------------------------------------------------C
C
C ** ASSIMILATE ALL LAYERS
C
      IF(ICALL.EQ.2)THEN
      IF(ISUVDA.GE.1)THEN
C
      DO NL=1,NLUVDA
C
        TMPIMP=1.-FSUVDA(NL)
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
	  LN=LNC(L)
	  LE=L+1
        NS=NUVSERA(NL)
        IF(ISTL.EQ.2)THEN
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUVDA(NL)*HU(L)*(USERT(K,NS)
     &                       -FSUVDA(NL)*U(L,K)-TMPIMP*U1(L,K))
            ELSE
	        FBODYFX(L,K)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSUVDA(NL)*HV(L)*(VSERT(K,NS)
     &                       -FSUVDA(NL)*V(L,K)-TMPIMP*V1(L,K))
            ELSE
	        FBODYFY(L,K)=0.0
            ENDIF
          ENDDO
        ELSE
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUVDA(NL)*HU(L)*(USERT(K,NS)
     &                       -FSUVDA(NL)*U(L,K)-TMPIMP*U2(L,K))
            ELSE
	        FBODYFX(L,K)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSUVDA(NL)*HV(L)*(VSERT(K,NS)
     &                       -FSUVDA(NL)*V(L,K)-TMPIMP*V2(L,K))
            ELSE
	        FBODYFY(L,K)=0.0
            ENDIF
          ENDDO
        ENDIF
C
        ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** CONSTRAIN ASSIMILATION TO DEPTH AVERAGE
C
      IF(ICALL.EQ.2)THEN
      IF(ISUVDA.EQ.1)THEN
C
      IF(IGRIDV.EQ.0)THEN
        DO NL=1,NLUVDA
          I=ICUVDA(NL)
          J=JCUVDA(NL)
          L=LIJ(I,J)
	    LN=LNC(L)
	    LE=L+1
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*FBODYFY(L,K)
          ENDDO
          DO K=1,KC
            FBODYFX(L,K)=AVGX
	      FBODYFY(L,K)=AVGY
          ENDDO
        ENDDO
      ENDIF
C
      IF(IGRIDV.EQ.1)THEN
        DO NL=1,NLUVDA
          I=ICUVDA(NL)
          J=JCUVDA(NL)
          L=LIJ(I,J)
	    LN=LNC(L)
	    LE=L+1
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*GVCSCLU(L)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*GVCSCLV(L)*FBODYFY(L,K)
          ENDDO
          DO K=KGVCU(L),KC
            FBODYFX(L,K)=AVGX
          ENDDO
          DO K=KGVCV(L),KC
	      FBODYFY(L,K)=AVGY
          ENDDO
        ENDDO
      ENDIF
C
      ENDIF
	ENDIF
C
C**********************************************************************C
C
      RETURN
      END
