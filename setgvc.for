c
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SETGVC
C
C **  THIS SUBROUTINE IS PART OF EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 08/01/2004        John Hamrick       01/11/2002       John Hamrick
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE SETGVC READS INFORMATION FOR THE GENERALIZED VERTICAL 
C **  COORDINATE OPTION AND SETS REAL AND LOCICAL MASK
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION IJKLAY(ICM,JCM),IHEADER(ICM)
      DIMENSION BELVCOR(LCM)
	CHARACTER*4 CV3D(ICM,JCM,KCM)
	CHARACTER*4 CUP3D(ICM,JCM,KCM),CW3D(ICM,JCM,KCM)
C
C**********************************************************************C
C
      DO I=1,IC
	  IHEADER(I)=I
	ENDDO
C
      DO J=1,JC
	  DO I=1,IC
	    IJKLAY(I,J)=0
	  ENDDO
	ENDDO
C
      DO K=1,KC
	  DO J=1,JC
	    DO I=1,IC
	      CUP3D(I,J,K)='|   '
	      CW3D(I,J,K)='|   '
	      CV3D(I,J,K)='+ - '
          ENDDO
	  ENDDO
	ENDDO
C
      GVCSCLP(1)=1.0
	GVCSCLU(1)=1.0
	GVCSCLV(1)=1.0
      GVCSCLP(LC)=0.0
	GVCSCLU(LC)=1.0
	GVCSCLV(LC)=1.0
      GVCSCLPI(1)=1.0
	GVCSCLUI(1)=1.0
	GVCSCLVI(1)=1.0
      GVCSCLPI(LC)=1.0
	GVCSCLUI(LC)=1.0
	GVCSCLVI(LC)=1.0
C
      DO L=1,LC
	  KGVCP(L)=1
	  KGVCU(L)=1
	  KGVCV(L)=1
	ENDDO
C
C**********************************************************************C
C
C **  READ VERTICAL GRID CELL MAPPING FROM CELLGVC.INP
C
      OPEN(1,FILE='CELLGVC.INP',STATUS='UNKNOWN')
C
      DO IS=1,4
      READ(1,100)
      ENDDO
C
      IF(IC.GT.120)THEN
        IACROSS=120
        DO IT=1,IC,IACROSS
          IFIRST=IT
          ILAST=IT+IACROSS-1
          IF(ILAST.GT.IC) ILAST=IC
          DO J=JC,1,-1
            READ(1,101,IOSTAT=ISO)JDUMY,(IJCTV(I,J),I=IFIRST,ILAST)
            IF(ISO.GT.0) THEN
	        WRITE(6,102)J
	        STOP
	      ENDIF
          ENDDO
        ENDDO
      ELSE
        IFIRST=1
        ILAST=IC
        DO J=JC,1,-1
          READ(1,101,IOSTAT=ISO)JDUMY,(IJCTV(I,J),I=IFIRST,ILAST)
          IF(ISO.GT.0) THEN
	      WRITE(6,102)J
	      STOP
	    ENDIF
        ENDDO
      ENDIF
C
      CLOSE(1)
C
C**********************************************************************C
C
C **  SET LCTV(L) USING IJCTV(I,J)
C
C     LCTV AND IJCTV EQUAL 1 FOR FOR RESCALED HEIGHT CELLS AND 
C     2 FOR SIGMA CELLS
C
      DO L=1,LC
	  LCTV(L)=0
	ENDDO
C
      DO J=1,JC
	  DO I=1,IC
	    IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LE.5)THEN
	      L=LIJ(I,J)
	      LCTV(L)=IJCTV(I,J)
	      IF(LCTV(L).EQ.0)THEN
	        WRITE(6,103)I,J
	        STOP
            ENDIF
          ENDIF
	  ENDDO
	ENDDO
C
C**********************************************************************C
C
C **  READ NUMBER OF ACTIVE VERTICAL LAYERS FROM GVCLAYER.INP
C **  AND WRITE SUMMARY TO GVCLAYER.OUT FOR ISETGVC=0
C
      IF(ISETGVC.EQ.0) THEN
C
      OPEN(1,FILE='GVCLAYER.INP',STATUS='UNKNOWN')
C
      DO IS=1,4
      READ(1,100)
      ENDDO
C
      DO LL=2,LA
	  READ(1,*)I,J,KLTMP
	  IJKLAY(I,J)=KLTMP
	  L=LIJ(I,J)
	  IF(LCTV(L).EQ.2) THEN
	    KGVCP(L)=KC-KSIG+1
	    GVCSCLP(L)=FLOAT(KC)/FLOAT(KSIG)
	    GVCSCLPI(L)=1./GVCSCLP(L)
        ENDIF
	  IF(LCTV(L).EQ.1) THEN
	    KGVCP(L)=KC-KLTMP+1 
	    GVCSCLP(L)=FLOAT(KC)/FLOAT(KLTMP)
	    GVCSCLPI(L)=1./GVCSCLP(L)
        ENDIF
	ENDDO
C
      CLOSE(1)
C
      OPEN(1,FILE='GVCLAYER.OUT',STATUS='UNKNOWN')
C
      DO L=2,LA
	  IF(LCTV(L).EQ.2) THEN
          WRITE(1,110)IL(L),JL(L),KSIG,KGVCP(L),GVCSCLP(L),
     &                BELV(L),BELV(L)
        ENDIF
	  IF(LCTV(L).EQ.1) THEN
	    KLTMP=KC-KGVCP(L)+1
          WRITE(1,110)IL(L),JL(L),KLTMP,KGVCP(L),GVCSCLP(L),
     &                BELV(L),BELV(L)
        ENDIF
	ENDDO
C
      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  AUTOMATICALLY SET BOTTOM LAYER K INDEX USING SELVREF, SELVREF
C **  AND BELV (IN DXDY.INP), ADJUST BELV, AND WRITE RESULTS 
C **  TO GVCLAYER.OUT FOR ISETGVC=1
C
      IF(ISETGVC.EQ.1) THEN
C
      DO L=2,LA
	  IF(LCTV(L).EQ.2) THEN
	    KGVCP(L)=KC-KSIG+1
	    GVCSCLP(L)=FLOAT(KC)/FLOAT(KSIG)
	    GVCSCLPI(L)=1./GVCSCLP(L)
		IJKLAY(IL(L),JL(L))=KSIG
        ENDIF
	ENDDO
C
	TMPSCL=SELVREF-BELVREF
C
      DO L=2,LA
	  IF(LCTV(L).EQ.1) THEN
	    TMPVAL=(SELVREF-BELV(L))/TMPSCL
          TMPVAL=TMPVAL*FLOAT(KC)
		KLTMP=NINT(TMPVAL)
	    GVCSCLP(L)=FLOAT(KC)/FLOAT(KLTMP)
	    GVCSCLPI(L)=1./GVCSCLP(L)
	    BELVCOR(L)=SELVREF-TMPSCL/GVCSCLP(L)
	    KGVCP(L)=KC-KLTMP+1 
 		IJKLAY(IL(L),JL(L))=KLTMP
       ENDIF
	ENDDO
C
      OPEN(1,FILE='GVCLAYER.OUT',STATUS='UNKNOWN')
C
      DO L=2,LA
	  IF(LCTV(L).EQ.2) THEN
          WRITE(1,110)IL(L),JL(L),KSIG,KGVCP(L),GVCSCLP(L),
     &                BELV(L),BELV(L)
        ENDIF
	  IF(LCTV(L).EQ.1) THEN
	    KLTMP=KC-KGVCP(L)+1
          WRITE(1,110)IL(L),JL(L),KLTMP,KGVCP(L),GVCSCLP(L),
     &                BELV(L),BELVCOR(L)
        ENDIF
	ENDDO
C
      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
      IF(KC.LT.10)THEN
C
      OPEN(1,FILE='CELLLAY.OUT')
C
      DO J=JC,1,-1
	WRITE(1,101)J,(IJKLAY(I,J),I=1,IC)
	ENDDO
C
      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  INITIALIZE LOGICAL MASK
C
      DO K=1,KC
      DO L=1,LC
	  LGVCP(L,K)=.FALSE.
	  LGVCU(L,K)=.FALSE.
	  LGVCV(L,K)=.FALSE.
	  LGVCW(L,K)=.FALSE.
	ENDDO
	ENDDO
C
C **  SET CELL CENTER LOGICAL MASK
C
      DO K=1,KC
      DO L=2,LA
	  IF(K.GE.KGVCP(L))LGVCP(L,K)=.TRUE.
	ENDDO
	ENDDO
C
C **  SET U FACE LOGICAL AND REAL MASKS
C
      DO L=2,LA
	  IF(SUB(L).LT.0.5)KGVCU(L)=KGVCP(L)
	  IF(SVB(L).LT.0.5)KGVCV(L)=KGVCP(L)
	ENDDO
C
      DO L=2,LA
	  GVCSCLU(L)=1.0
	  IF(SUB(L).GT.0.5)THEN
	    KGVCU(L)=MAX(KGVCP(L),KGVCP(L-1))   ! *** Maximum layer number of active bottom layer
	    KLTMP=KC-KGVCU(L)+1                 ! *** Number of active Layers
	    GVCSCLU(L)=FLOAT(KC)/FLOAT(KLTMP)   ! *** Ratio of Total/Minimum Active Layers
          DO K=KGVCU(L),KC 
            LGVCU(L,K)=.TRUE.
	    ENDDO
		  DO K=1,KC 
            IF(LGVCU(L,K))THEN
	        SUB3D(L,K)=1.0
	        SUB3DO(L,K)=1.0
	        SBX3D(L,K)=1.0
	        SDX3D(L,K)=1.0
	      ELSE
	        SUB3D(L,K)=0.0
	        SUB3DO(L,K)=0.0
	        SBX3D(L,K)=0.0
	        SDX3D(L,K)=0.0
	      ENDIF
	    ENDDO
	  ENDIF
	ENDDO
C
C **  SET V FACE LOGICAL AND REAL MASKS
C
      DO L=2,LA
	  GVCSCLV(L)=1.0
	  LS=LSC(L)
	  IF(SVB(L).GT.0.5)THEN
	    KGVCV(L)=MAX(KGVCP(L),KGVCP(LS))
	    KLTMP=KC-KGVCV(L)+1 
	    GVCSCLV(L)=FLOAT(KC)/FLOAT(KLTMP)
		DO K=KGVCV(L),KC 
            LGVCV(L,K)=.TRUE.
	    ENDDO
		DO K=1,KC 
            IF(LGVCV(L,K))THEN
	        SVB3D(L,K)=1.0
	        SVB3DO(L,K)=1.0
	        SBY3D(L,K)=1.0
	        SDY3D(L,K)=1.0
	      ELSE
	        SVB3D(L,K)=0.0
	        SVB3DO(L,K)=0.0
	        SBY3D(L,K)=0.0
	        SDY3D(L,K)=0.0
	      ENDIF
	    ENDDO
	  ENDIF
	ENDDO
C
      DO L=2,LA
	  GVCSCLUI(L)=1./GVCSCLU(L)
	  GVCSCLVI(L)=1./GVCSCLV(L)
	ENDDO
C
C **  SET W FACE LOGICAL AND REAL MASKS
C
      DO L=2,LA
        SWB3D(L,0)=0.0
	  SWB3DO(L,0)=0.0
        SWB3D(L,KC)=0.0
	  SWB3DO(L,KC)=0.0
	ENDDO
C
      DO K=1,KS
      DO L=2,LA
	  IF(LGVCP(L,K))THEN
	    IF(LGVCP(L,K+1))THEN
            LGVCW(L,K)=.TRUE.
	    ENDIF
	  ENDIF
	ENDDO
	ENDDO
C
      DO K=1,KS
      DO L=2,LA
	  IF(LGVCW(L,K))THEN
          SWB3D(L,K)=SWB(L)
          SWB3DO(L,K)=SWB(L)
	  ELSE
          SWB3D(L,K)=0.0
          SWB3DO(L,K)=0.0
	  ENDIF
	ENDDO
	ENDDO
C
C**********************************************************************C
C
C **  RESET CELL CENTER DEPTHS
C
      DO L=2,LA
      LS=LSC(L)      
       IF(SUB(L).GT.0.5)THEN
	   HU(L)=0.5*(DXP(L)*DYP(L)*GVCSCLP(L)*HP(L)
     &           +DXP(L-1)*DYP(L-1)*GVCSCLP(L-1)*HP(L-1))
     &           /(GVCSCLU(L)*DXU(L)*DYU(L))
	   HUI(L)=1./HU(L)
	   H1U(L)=0.5*(DXP(L)*DYP(L)*GVCSCLP(L)*H1P(L)
     &           +DXP(L-1)*DYP(L-1)*GVCSCLP(L-1)*H1P(L-1))
     &           /(GVCSCLU(L)*DXU(L)*DYU(L))
	   H1UI(L)=1./H1U(L)
	 ELSE
	   HU(L)=HP(L)
	   HUI(L)=1./HU(L)
	   H1U(L)=H1P(L)
	   H1UI(L)=1./H1U(L)
	 ENDIF
        IF(SVB(L).GT.0.5)THEN
          HV(L)=0.5*(DXP(L)*DYP(L)*GVCSCLP(L)*HP(L)
     &           +DXP(LS )*DYP(LS)*GVCSCLP(LS)*HP(LS ))
     &           /(GVCSCLV(L)*DXV(L)*DYV(L))
	    HVI(L)=1./HV(L)
          H1V(L)=0.5*(DXP(L)*DYP(L)*GVCSCLP(L)*H1P(L)
     &           +DXP(LS )*DYP(LS)*GVCSCLP(LS)*H1P(LS ))
     &           /(GVCSCLV(L)*DXV(L)*DYV(L))
	    H1VI(L)=1./H1V(L)
        ELSE
	   HV(L)=HP(L)
	   HVI(L)=1./HV(L)
	   H1V(L)=H1P(L)
	   H1VI(L)=1./H1V(L)
	 ENDIF
      ENDDO
C
C
C**********************************************************************C
C
C **  SET INTERNAL MODE SOLUTION COEFFICIENTS
C
c           CDZR(1)=DZC(1)-1.
c
c           DO K=2,KS
c           CDZR(K)=DZC(K)+CDZR(K-1)
c           ENDDO
c
c           DO K=1,KS
c           CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
c           ENDDO
C
      DO K=1,KC
	 DO L=1,LC
	   CDZRGVCU(L,K)=0.0
	   CDZRGVCV(L,K)=0.0
	 ENDDO
	ENDDO
C
      DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBPLUS=KBTMPU+1
	  CDZRGVCU(L,KBTMPU)=GVCSCLU(L)*DZC(KBTMPU)-1.
        DO K=KBPLUS,KS
	    CDZRGVCU(L,K)=GVCSCLU(L)*DZC(K)+CDZRGVCU(L,K-1)
        ENDDO
        DO K=KBTMPU,KS
	    CDZRGVCU(L,K)=CDZRGVCU(L,K)*DZG(K)*CDZL(KBTMPU)
        ENDDO
      ENDDO
  800 FORMAT('CDZR,I,J,K',4I5,2E14.6)
C
      DO L=2,LA
	  KBTMPV=KGVCV(L)
	  KBPLUS=KBTMPV+1
	  CDZRGVCV(L,KBTMPV)=GVCSCLV(L)*DZC(KBTMPV)-1.
        DO K=KBPLUS,KS
	    CDZRGVCV(L,K)=GVCSCLV(L)*DZC(K)+CDZRGVCV(L,K-1)
        ENDDO
        DO K=KBTMPV,KS
	    CDZRGVCV(L,K)=CDZRGVCV(L,K)*DZG(K)*CDZL(KBTMPV)
        ENDDO
      ENDDO
C
C        CDZDGVCU(L,K)=GVCSCLU(L)*(CDZD(K) FIXED FOR GVC)
C	   CDZDGVCU(L,K)=GVCSCLV(L)*(CDZD(K) FIXED FOR GVC)
C
      DO K=1,KC
	 DO L=1,LC
	   CDZDGVCU(L,K)=0.0
	   CDZDGVCV(L,K)=0.0
	 ENDDO
	ENDDO
C
      DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBPLUS=KBTMPU+1
	  CDZDGVCU(L,KBTMPU)=GVCSCLU(L)*DZC(KBTMPU)
        DO K=KBPLUS,KS
	    CDZDGVCU(L,K)=GVCSCLU(L)*DZC(K)+CDZDGVCU(L,K-1)
        ENDDO
      ENDDO
C
      DO L=2,LA
	  KBTMPV=KGVCV(L)
	  KBPLUS=KBTMPV+1
	  CDZDGVCV(L,KBTMPV)=GVCSCLV(L)*DZC(KBTMPV)
        DO K=KBPLUS,KS
	    CDZDGVCV(L,K)=GVCSCLV(L)*DZC(K)+CDZDGVCV(L,K-1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  RESET SCALAR VARIABLES
C
       IF(ISTRAN(1).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) SAL(L,K)=0.0
          IF(K.LT.KGVCP(L)) SAL1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(2).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) TEM(L,K)=0.0
          IF(K.LT.KGVCP(L)) TEM1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(3).GE.1)THEN
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) DYE(L,K)=0.0
          IF(K.LT.KGVCP(L)) DYE1(L,K)=0.0
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(5).GE.1)THEN
         DO NT=1,NTOX      
          DO K=1,KC
          DO L=1,LC
           IF(K.LT.KGVCP(L)) TOX(L,K,NT)=0.0
          ENDDO
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(6).GE.1)THEN
         DO NS=1,NSED
          DO K=1,KC
          DO L=1,LC
           IF(K.LT.KGVCP(L)) SED(L,K,NS)=0.0
          ENDDO
         ENDDO
         ENDDO
       ENDIF
       IF(ISTRAN(7).GE.1)THEN
        DO NS=1,NSND
         DO K=1,KC
         DO L=1,LC
          IF(K.LT.KGVCP(L)) SND(L,K,NS)=0.0
         ENDDO
        ENDDO
         ENDDO
       ENDIF
C
C**********************************************************************C
C
C **  SET TURBULENCE MODEL PROXIMITY FUNCTION
C
      DO K=0,KC
	DO L=2,LA
	  FPROXGVC(L,K)=0.0
	ENDDO
	ENDDO
C
      IF(IFPROX.EQ.1)THEN
	  DO K=1,KS
	  DO L=2,LA
		KBOT=KGVCP(L)-1
	    IF(K.GE.KGVCP(L))THEN
	      FPROXGVC(L,K)=1./(VKC*(Z(K)-Z(KBOT))*(1.-Z(K)))**2
	    ENDIF
        ENDDO
        ENDDO
      ENDIF
C
      IF(IFPROX.EQ.2)THEN
	  DO K=1,KS
	  DO L=2,LA
		KBOT=KGVCP(L)-1
	    IF(K.GE.KGVCP(L))THEN
	      FPROXGVC(L,K)=(1./(VKC*(Z(K)-Z(KBOT)))**2)
     &          +CTE5*(1./(VKC*(1.-Z(K)))**2)/(CTE4+0.00001)
	    ENDIF
        ENDDO
        ENDDO
      ENDIF
C
      OPEN(1,FILE='GVCPROX.OUT',STATUS='UNKNOWN')
	DO L=2,LA
	  WRITE(1,1947)IL(L),JL(L),KGVCP(L),(FPROXGVC(L,K),K=0,KC)
      ENDDO
	CLOSE(1)
C
 1947 FORMAT(3I6,28F12.5)
C
C**********************************************************************C
C
C **  ADJUST STEADY INFLOWS TO SUM TO TOTAL
C

      DO LL=1,NQSIJ
        L=LQS(LL)
        DO K=1,KC
	    IF(K.GE.KGVCP(L))THEN
            QSS(K,LL)=GVCSCLP(L)*QSS(K,LL)
          ELSE
            QSS(K,LL)=0.0
          ENDIF
        ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  OUTPUT MASKING RESULTS
C

      OPEN(1,FILE='GVCMASK.OUT',STATUS='UNKNOWN')
C
      DO K=1,KC
	  DO J=1,JC
	    DO I=1,IC
            IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LE.5)THEN
	        L=LIJ(I,J)
	        IF(LGVCP(L,K)) THEN
	          IF(LGVCU(L,K)) THEN
		        CUP3D(I,J,K)='U O '
                ELSE
		        CUP3D(I,J,K)='| O '
	          ENDIF
	        ENDIF
	        IF(LGVCV(L,K)) CV3D(I,J,K)='+ V '
	      ENDIF
          ENDDO
	  ENDDO
	ENDDO
C
      DO K=2,KC
	  DO J=1,JC
	    DO I=1,IC
            IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LE.5)THEN
	        L=LIJ(I,J)
	        IF(LGVCU(L,K)) THEN
	          IF(LGVCW(L,K-1)) THEN
		        CW3D(I,J,K)='U W '
                ELSE
		        CW3D(I,J,K)='U   '
	          ENDIF
	        ELSE
	          IF(LGVCW(L,K-1)) THEN
		        CW3D(I,J,K)='| W '
                ELSE
		        CW3D(I,J,K)='|   '
	          ENDIF
	        ENDIF
	      ENDIF
          ENDDO
	  ENDDO
	ENDDO
C
        K=1
	  DO J=1,JC
	    DO I=1,IC
            IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LE.5)THEN
	        L=LIJ(I,J)
	        IF(LGVCU(L,K)) THEN
		      CW3D(I,J,K)='U   '
	        ELSE
		      CW3D(I,J,K)='|   '
	        ENDIF
	      ENDIF
          ENDDO
	  ENDDO
C
      DO K=KC,1,-1
	  WRITE(1,120)
	  WRITE(1,121)K
	  WRITE(1,124)(IHEADER(I),I=1,IC)
	  WRITE(1,123)(CV3D(I,1,K),I=1,IC)
	  DO J=JC,1,-1
	    WRITE(1,122)J,(CUP3D(I,J,K),I=1,IC)
	    WRITE(1,123)(CV3D(I,J,K),I=1,IC)
        ENDDO
	  WRITE(1,120)
	  WRITE(1,124)(IHEADER(I),I=1,IC)
	  WRITE(1,123)(CV3D(I,1,K),I=1,IC)
	  DO J=JC,1,-1
	    WRITE(1,122)J,(CW3D(I,J,K),I=1,IC)
	    WRITE(1,123)(CV3D(I,J,K),I=1,IC)
        ENDDO
	ENDDO	   
C
      CLOSE(1)
C
C**********************************************************************C
C
  100 FORMAT(1X)
  101 FORMAT(I3,2X,120I1)
  102 FORMAT(' READ ERROR FOR FILE CELLGVC.INP ON J = ',I5)
  103 FORMAT(' INCONSISTENCY BETWEEN IJCT AND IJCTV AT I,J = ',2I5)
  110 FORMAT(4I6,3F10.4)
  111 FORMAT(3I6,4F6.1)
  120 FORMAT(/)
  121 FORMAT('GVCMASK FOR LAYER  ',I4,/)
  122 FORMAT(I3,2X,200A4)
  123 FORMAT(5X,200A4)
  124 FORMAT(4X,200I4)
C
C**********************************************************************C
C
      RETURN
      END
