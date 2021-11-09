C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTSXY
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
C ** SUBROUTINE CALTSXY UPDATES TIME VARIABLE SURFACE WIND STRESS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION  PATMTT(NASERM),TATMTT(NASERM),TWETTT(NASERM),
     &           RAINTT(NASERM),EVAPTT(NASERM),SOLSWRTT(NASERM),
     &           CLOUDTT(NASERM),SVPAT(NASERM),VPAT(NASERM),
     &           RHAT(NASERM),CLEVAPT(NASERM),CCNHTTT(NASERM),
     &           WINDE(NWSERM),WINDN(NWSERM)
C
C**********************************************************************C
C
C INITIALIZE WIND SHELTERED SURFACE GAS TRANSFER
C
      IF(N.EQ.-1)THEN
      OPEN(1,FILE='WINDSHELT.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WINDSHELT.OUT')
C
      DO L=2,LA
       I=IL(L)
       J=JL(L)
C ** IF WINDSTKA > 0 BOTH X AND Y COMPONENTS ARE APPLIED
       IF(WINDSTKA(L).GT.0.0)THEN
         WINDSXX(L)=CVN(L)
         WINDSXY(L)=-CVE(L)
         WINDSYX(L)=-CUN(L)
         WINDSYY(L)=CUE(L)
        ELSE
C ** IF WINDSTKA < 0 SLECTIVELY APPLY X AND Y COMPONENTS
C ** FIRST CASE IS FULLY OPEN WATER
         WINDSXX(L)=CVN(L)
         WINDSXY(L)=-CVE(L)
         WINDSYX(L)=-CUN(L)
         WINDSYY(L)=CUE(L)
         LS=LSC(L)
         LN=LNC(L)
C ** SECOND CASE IS 1D CHANNEL IN COMP X DIRECTION
         IF(SVB(L).LT.0.5.AND.IJCT(I,J-1).NE.5)THEN
         IF(SVB(LN).LT.0.5.AND.IJCT(I,J+1).NE.5)THEN
           WINDSXX(L)=CVN(L)
           WINDSXY(L)=-CVE(L)
           WINDSYX(L)=-1000.
           WINDSYY(L)=0.
         ENDIF
         ENDIF
C ** THIRD CASE IS 1D CHANNEL IN COMP Y DIRECTION
         IF(SUB(L).LT.0.5.AND.IJCT(I-1,J).NE.5)THEN
         IF(SUB(L+1).LT.0.5.AND.IJCT(I+1,J).NE.5)THEN
           WINDSXX(L)=0.
           WINDSXY(L)=-1000.
           WINDSYX(L)=-CUN(L)
           WINDSYY(L)=CUE(L)
         ENDIF
         ENDIF
       ENDIF
       WRITE(1,1111)IL(L),JL(L),WINDSTKA(L),WINDSXX(L),WINDSXY(L),
     &           WINDSYX(L),WINDSYY(L)
      ENDDO
C
      CLOSE(1)
      ENDIF
C
 1111 FORMAT(2I5,10F13.6)
C
C**********************************************************************C
C
C  NEW MINCURVATURE SPLINE FLOAT(I) AND FLOAT(J) ARE X AND Y COORDINATES
C  CALCULATE KERNAL FOR 8 WIND STATIONS
C
      IF(NWSER.LE.-1.AND.N.EQ.1)THEN
C
      OPEN(1,FILE='WNDKERNAL.OUT')
C
      DO M=1,8
	  XS=XWFKER(M)
	  YS=YWFKER(M)
        DO L=2,LA
	    X=XWFMUL*(FLOAT(IL(L))+XWFADD)
          Y=YWFMUL*(FLOAT(JL(L))+YWFADD)
          DELX=X-XS
          DELY=Y-YS
          DEL=SQRT(DELX*DELX+DELY*DELY)
          WFKER(L,M)=DEL**1.5
	    WRITE(1,1492)M,IL(L),JL(L),XS,YS,X,Y,WFKER(L,M)
	  ENDDO
	ENDDO
C
C      OPEN(1,FILE='WNDKERNAL.OUT')
C	DO L=2,LA
C	  WRITE(1,1492)IL(L),JL(L),(WFKER(L,M),M=1,8)
C	ENDDO
C
      WRITE(6,*)'INITIALIZED WIND FORCING KERNALS'
      WRITE(8,*)'INITIALIZED WIND FORCING KERNALS'
C
      CLOSE(1)
C
      ENDIF
C
 1492 FORMAT(3I5,8F14.4)
C
C**********************************************************************C
C
      IF(NWSER.GT.0)THEN
C
      DO NA=1,NWSER
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)/TCWSER(NA)+TBEGIN*(TCON/TCWSER(NA))
        ELSE
          TIME=TIMESEC/TCWSER(NA)
        ENDIF
        M1=MWTLAST(NA)
        MSAVE=M1    
  200   CONTINUE
        M2=M1+1
        IF(TIME.GT.TWSER(M2,NA))THEN
          M1=M2
          GOTO 200
        ELSE
          MWTLAST(NA)=M1
        ENDIF      
        TDIFF=TWSER(M2,NA)-TWSER(M1,NA)
        WTM1=(TWSER(M2,NA)-TIME)/TDIFF
        WTM2=(TIME-TWSER(M1,NA))/TDIFF
        IF(ISWDINT(NA).LE.1)THEN
          DEGM1=90.-WINDD(M1,NA)
          DEGM2=90.-WINDD(M2,NA)
          WINDS1=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
          WINDS2=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
          WINDE1=WINDS(M1,NA)*COS(DEGM1/57.29578)
          WINDN1=WINDS(M1,NA)*SIN(DEGM1/57.29578) 
          WINDE2=WINDS(M2,NA)*COS(DEGM2/57.29578)
          WINDN2=WINDS(M2,NA)*SIN(DEGM2/57.29578) 
          WINDE(NA)=WTM1*WINDE1+WTM2*WINDE2
          WINDN(NA)=WTM1*WINDN1+WTM2*WINDN2
        ELSE
          WINDS1=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
          WINDS2=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
          WINDE1=WINDS(M1,NA)
          WINDN1=WINDD(M1,NA) 
          WINDE2=WINDS(M2,NA)
          WINDN2=WINDD(M2,NA)
          WINDE(NA)=WTM1*WINDE1+WTM2*WINDE2
          WINDN(NA)=WTM1*WINDN1+WTM2*WINDN2
	  ENDIF
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------c
C
C OLD CONTINUITY CONSTRAINTED CUBIC
C
      IF(NWSER.LE.-98)THEN
C
      IF(NWSER.EQ.-2)THEN
        OPEN(1,FILE='WNDFLDCK.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='WNDFLDCK.OUT')
	ENDIF
C
C ** BI-QUADRATIC WIND FIELD FUNCTION WITH CONTINUITY CONSTRAINT
C
        NA=1
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)/TCWSER(NA)+TBEGIN*(TCON/TCWSER(NA))
        ELSE
          TIME=TIMESEC/TCWSER(NA)
        ENDIF
        M1=MWTLAST(NA)
        MSAVE=M1    
  300   CONTINUE
        M2=M1+1
        IF(TIME.GT.TWSER(M2,NA))THEN
          M1=M2
          GOTO 300
        ELSE
          MWTLAST(NA)=M1
        ENDIF      
        TDIFF=TWSER(M2,NA)-TWSER(M1,NA)
        WTM1=(TWSER(M2,NA)-TIME)/TDIFF
        WTM2=(TIME-TWSER(M1,NA))/TDIFF
C
        DO L=2,LA
	    X=XWFMUL*(DLON(L)+XWFADD)
          Y=YWFMUL*(DLAT(L)+YWFADD)
          WINDE1=WFCOEF(M1,1)+WFCOEF(M1,2)*X+WFCOEF(M1,3)*Y
     &          +WFCOEF(M1,4)*X*X+WFCOEF(M1,5)*X*Y+WFCOEF(M1,6)*Y*Y
          WINDN1=WFCOEF(M1,7)+WFCOEF(M1,8)*X-WFCOEF(M1,2)*Y
     &      +WFCOEF(M1,9)*X*X-2.*WFCOEF(M1,4)*X*Y-0.5*WFCOEF(M1,5)*Y*Y
          WINDE2=WFCOEF(M2,1)+WFCOEF(M2,2)*X+WFCOEF(M2,3)*Y
     &          +WFCOEF(M2,4)*X*X+WFCOEF(M2,5)*X*Y+WFCOEF(M2,6)*Y*Y
          WINDN2=WFCOEF(M2,7)+WFCOEF(M2,8)*X-WFCOEF(M2,2)*Y
     &      +WFCOEF(M2,9)*X*X-2.*WFCOEF(M2,4)*X*Y-0.5*WFCOEF(M2,5)*Y*Y
          WNDVELE(L)=WTM1*WINDE1+WTM2*WINDE2
          WNDVELN(L)=WTM1*WINDN1+WTM2*WINDN2
          WINDST(L)=SQRT( WNDVELE(L)*WNDVELE(L)
     &                 +WNDVELN(L)*WNDVELN(L) )
        ENDDO
C
        DO L=2,LA
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSEAST=1.2E-3*CD10*WINDST(L)*WNDVELE(L)
          TSNORT=1.2E-3*CD10*WINDST(L)*WNDVELN(L)
          TSX(L)=WINDSXX(L)*TSEAST+WINDSXY(L)*TSNORT
          TSY(L)=WINDSYX(L)*TSEAST+WINDSYY(L)*TSNORT
        ENDDO
C
C      IF(NWSER.EQ.-2)THEN
C	  WRITE(1,1949)M1,M2,WTM1,WTM2,TIME
C	  DO L=2,LA
C	    WRITE(1,1949)IL(L),JL(L),WNDVELE(L),WNDVELN(L),WINDST(L),
C     &       TSX(L),TSY(L)      
C        ENDDO
C	  CLOSE(1)
C	ENDIF   

      ENDIF
C
C----------------------------------------------------------------------c
C
C  NEW MINCURVATURE SPLINE FLOAT(I) AND FLOAT(J) ARE X AND Y COORDINATES
C  WIND IS PRE-ROTATED TO CARTESIAN GRID
C
      IF(NWSER.LE.-1)THEN
C
      IF(NWSER.EQ.-2)THEN
        OPEN(1,FILE='WNDFLDCK.OUT')
	  CLOSE(1,STATUS='DELETE')
	  OPEN(1,FILE='WNDFLDCK.OUT')
	ENDIF
C
C ** MIN CURVATURE SPLINE HARD WIRED FOR 8 WIND FIELDS
C
        NA=1
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)/TCWSER(NA)+TBEGIN*(TCON/TCWSER(NA))
        ELSE
          TIME=TIMESEC/TCWSER(NA)
        ENDIF
        M1=MWTLAST(NA)
        MSAVE=M1    
  400   CONTINUE
        M2=M1+1
        IF(TIME.GT.TWSER(M2,NA))THEN
          M1=M2
          GOTO 400
        ELSE
          MWTLAST(NA)=M1
        ENDIF      
        TDIFF=TWSER(M2,NA)-TWSER(M1,NA)
        WTM1=(TWSER(M2,NA)-TIME)/TDIFF
        WTM2=(TIME-TWSER(M1,NA))/TDIFF
C
        DO L=2,LA
C	    X=XWFMUL*(DLON(L)+XWFADD)
C          Y=YWFMUL*(DLAT(L)+YWFADD)
	    X=XWFMUL*(FLOAT(IL(L))+XWFADD)
          Y=YWFMUL*(FLOAT(JL(L))+YWFADD)
          WINDE1=WFCOEF(M1,1)+WFCOEF(M1,2)*X+WFCOEF(M1,3)*Y
	    DO K=4,11
	      KK=K-3
            WINDE1=WINDE1+WFCOEF(M1,K)*WFKER(L,KK)
          ENDDO
          WINDE2=WFCOEF(M2,1)+WFCOEF(M2,2)*X+WFCOEF(M2,3)*Y
	    DO K=4,11
	      KK=K-3
            WINDE2=WINDE2+WFCOEF(M2,K)*WFKER(L,KK)
          ENDDO
          WINDN1=WFCOEF(M1,12)+WFCOEF(M1,13)*X-WFCOEF(M1,14)*Y
	    DO K=15,22
	      KK=K-14
            WINDN1=WINDN1+WFCOEF(M1,K)*WFKER(L,KK)
          ENDDO
          WINDN2=WFCOEF(M2,12)+WFCOEF(M2,13)*X-WFCOEF(M2,14)*Y
	    DO K=15,22
	      KK=K-14
            WINDN2=WINDN2+WFCOEF(M2,K)*WFKER(L,KK)
          ENDDO
          WNDVELE(L)=WTM1*WINDE1+WTM2*WINDE2
          WNDVELN(L)=WTM1*WINDN1+WTM2*WINDN2
          WINDST(L)=SQRT( WNDVELE(L)*WNDVELE(L)
     &                 +WNDVELN(L)*WNDVELN(L) )
        ENDDO
C
        DO L=2,LA
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSEAST=1.2E-3*CD10*WINDST(L)*WNDVELE(L)
          TSNORT=1.2E-3*CD10*WINDST(L)*WNDVELN(L)
C          TSX(L)=WINDSXX(L)*TSEAST+WINDSXY(L)*TSNORT
C          TSY(L)=WINDSYX(L)*TSEAST+WINDSYY(L)*TSNORT
          TSX(L)=TSEAST
          TSY(L)=TSNORT
        ENDDO
C
      IF(NWSER.EQ.-2)THEN
	  WRITE(1,1949)M1,M2,WTM1,WTM2,TIME
	  DO L=2,LA
	    WRITE(1,1949)IL(L),JL(L),WNDVELE(L),WNDVELN(L),WINDST(L),
     &       TSX(L),TSY(L)      
        ENDDO
	  CLOSE(1)
	ENDIF   

      ENDIF
C
 1949 FORMAT(2I5,3F12.5,2E13.5)
C
C**********************************************************************C
C
      IF(NWSER.GT.0)THEN
C
      NWMP=1
	IF(NWNDMAP.GT.1)THEN
	  DO NTMP=1,NWNDMAP
	    IF(TIME.GE.TWNDMAPBEG(NTMP).AND.TIME.LT.TWNDMAPEND(NTMP))THEN
	      NWMP=NTMP
          ENDIF
        ENDDO
	ENDIF
C
      DO L=2,LA
        WNDVELE(L)=0.
        WNDVELN(L)=0.
      ENDDO
C
      DO NA=1,NWSER
        DO L=2,LA
          WNDVELE(L)=WNDVELE(L)+WNDWHT(L,NA,NWMP)*WINDE(NA)
          WNDVELN(L)=WNDVELN(L)+WNDWHT(L,NA,NWMP)*WINDN(NA)
        ENDDO
      ENDDO
C
      DO L=2,LA
C ** CASE 0 MAGNITUDE SHELTERING AND NO DIRECTIONAL SHELTERING
        IF(WINDSTKA(L).GT.0.0)THEN
          WNDFAC=ABS(WINDSTKA(L))
          WNDVELE(L)=WNDFAC*WNDVELE(L)
          WNDVELN(L)=WNDFAC*WNDVELN(L)
          WINDST(L)=SQRT( WNDVELE(L)*WNDVELE(L)
     &                 +WNDVELN(L)*WNDVELN(L) )
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSEAST=1.2E-3*CD10*WINDST(L)*WNDVELE(L)
          TSNORT=1.2E-3*CD10*WINDST(L)*WNDVELN(L)
C          TSX(L)= CVN(L)*TSEAST-CVE(L)*TSNORT
C          TSY(L)=-CUN(L)*TSEAST+CUE(L)*TSNORT
          TSX(L)=WINDSXX(L)*TSEAST+WINDSXY(L)*TSNORT
          TSY(L)=WINDSYX(L)*TSEAST+WINDSYY(L)*TSNORT
        ENDIF
C ** CASE 1 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, OPEN WATER
        IF(WINDSTKA(L).LT.0.0)THEN
        IF(WINDSYX(L).GT.-99.0.AND.WINDSXY(L).GT.-99.0)THEN
          WNDFAC=ABS(WINDSTKA(L))
          WNDVELE(L)=WNDFAC*WNDVELE(L)
          WNDVELN(L)=WNDFAC*WNDVELN(L)
          WINDST(L)=SQRT( WNDVELE(L)*WNDVELE(L)
     &                 +WNDVELN(L)*WNDVELN(L) )
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSEAST=1.2E-3*CD10*WINDST(L)*WNDVELE(L)
          TSNORT=1.2E-3*CD10*WINDST(L)*WNDVELN(L)
C          TSX(L)= CVN(L)*TSEAST-CVE(L)*TSNORT
C          TSY(L)=-CUN(L)*TSEAST+CUE(L)*TSNORT
          TSX(L)=WINDSXX(L)*TSEAST+WINDSXY(L)*TSNORT
          TSY(L)=WINDSYX(L)*TSEAST+WINDSYY(L)*TSNORT
        ENDIF
        ENDIF
C ** CASE 2 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, X CHANNEL
        IF(WINDSTKA(L).LT.0.0)THEN
        IF(WINDSYX(L).LT.-99.0)THEN
          WINDXX=WINDSXX(L)*WNDVELE(L)+WINDSXY(L)*WNDVELN(L)
          WNDFAC=ABS(WINDSTKA(L))
          WINDXX=WNDFAC*WNDVELE(L)
          WINDST(L)=ABS(WINDXX)
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSX(L)=1.2E-3*CD10*WINDST(L)*WINDST(L)*WINDXX
          TSY(L)=0.
        ENDIF
        ENDIF
C ** CASE 3 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, Y CHANNEL
        IF(WINDSTKA(L).LT.0.0)THEN
        IF(WINDSXY(L).LT.-99.0)THEN
          WINDYY=WINDSYX(L)*WNDVELE(L)+WINDSYY(L)*WNDVELN(L)
          WNDFAC=ABS(WINDSTKA(L))
          WINDYY=WNDFAC*WINDYY
          WINDST(L)=ABS(WINDYY)
	    WTMP=MAX(WINDST(L),0.25)
          IF(WTMP.LT.5.)THEN
	      WTMP=1./WTMP
            CD10=3.83111e-005*(WTMP**3)-0.000308715*(WTMP**2) 
     &          +0.00116012*WTMP+0.000899602
	    ENDIF
          IF(WTMP.GE.5.0.and.WTMP.LE.7.)THEN
            CD10=-5.37642e-006*(WTMP**3)+0.000112556*(WTMP**2)
     &           -0.000721203*WTMP+0.00259657
	    ENDIF
          IF(WTMP.GE.7.)THEN
            CD10=-3.99677e-007*(WTMP**2)+7.32937e-005*WTMP+0.000726716
	    ENDIF
          TSX(L)=0
          TSY(L)=1.2E-3*CD10*WINDST(L)*WINDST(L)*WINDYY
        ENDIF
        ENDIF
      ENDDO
C
      ENDIF
C
c      IF(ISVEG.GE.1.OR.ISDRY.GE.1)THEN
c        DO L=2,LA
c          CFTSX=1.
c          CFTSY=1.
c          HHUU=2.*HUWET(L)
c          HHVV=2.*HVWET(L)
c          IF(HU(L).LT.0.3) CFTSX=(HU(L)-HHUU)/(0.3-HHUU)
c          IF(HU(L).LE.HHUU) CFTSX=0.
c          IF(HV(L).LT.0.3) CFTSY=(HV(L)-HHVV)/(0.3-HHVV)
c          IF(HV(L).LE.HHVV) CFTSY=0.
c          TSX(L)=CFTSX*TSX(L)
c          TSY(L)=CFTSY*TSY(L)
c          IF(MVEGL(L).GT.1) TSX(L)=0. 
c          IF(MVEGL(L).GT.1) TSY(L)=0. 
c        ENDDO
c      ENDIF
C
C**********************************************************************C
C
      IF(NASER.GT.0)THEN
C
      DO NA=1,NASER
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)/TCASER(NA)+TBEGIN*(TCON/TCASER(NA))
        ELSE
          TIME=TIMESEC/TCASER(NA)
        ENDIF
        M1=MATLAST(NA)    
  100   CONTINUE
        M2=M1+1
        IF(TIME.GT.TASER(M2,NA))THEN
          M1=M2
          GOTO 100
        ELSE
          MATLAST(NA)=M1
        ENDIF      
        TDIFF=TASER(M2,NA)-TASER(M1,NA)
        WTM1=(TASER(M2,NA)-TIME)/TDIFF
        WTM2=(TIME-TASER(M1,NA))/TDIFF
        PATMTT(NA)=WTM1*PATM(M1,NA)+WTM2*PATM(M2,NA)
        TATMTT(NA)=WTM1*TDRY(M1,NA)+WTM2*TDRY(M2,NA)
        TWETTT(NA)=WTM1*TWET(M1,NA)+WTM2*TWET(M2,NA)
        RAINTT(NA)=WTM1*RAIN(M1,NA)+WTM2*RAIN(M2,NA)
        EVAPTT(NA)=WTM1*EVAP(M1,NA)+WTM2*EVAP(M2,NA)
        SOLSWRTT(NA)=WTM1*SOLSWR(M1,NA)+WTM2*SOLSWR(M2,NA)
        CLOUDTT(NA)=WTM1*CLOUD(M1,NA)+WTM2*CLOUD(M2,NA)
        SVPAT(NA)=
     &  (10.**((0.7859+0.03477*TATMTT(NA))/(1.+0.00412*TATMTT(NA))))
C    &  *(1+1.E-6*PATMTT(NA)*(4.5+0.0006*TATMTT(NA)*TATMT(NA)))
        IF(IRELH(NA).EQ.0)THEN
          RHAT(NA)=1.
     &         -0.00066*(PATMTT(NA)/SVPAT(NA))*(TATMTT(NA)-TWETTT(NA))
        ELSE
          RHAT(NA)=TWETTT(NA)
        ENDIF
        VPAT(NA)=RHAT(NA)*SVPAT(NA)
        CLEVAPT(NA)=1.E-3
        CCNHTTT(NA)=1.E-3
        IF(REVC.GT.1.E-6) CLEVAPT(NA)=REVC*1.E-3
        IF(RCHC.GT.1.E-6) CCNHTTT(NA)=RCHC*1.E-3
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      IF(NASER.GT.0)THEN
C
      DO L=2,LA
        PATMT(L)=0.
        TATMT(L)=0.
        RAINT(L)=0.
        EVAPT(L)=0.
        SOLSWRT(L)=0.
        CLOUDT(L)=0.
        SVPA(L)=0.
        RHA(L)=0.
        VPA(L)=0.
        CLEVAP(L)=0.
        CCNHTT(L)=0.
      ENDDO
C
      DO NA=1,NASER
        DO L=2,LA
          PATMT(L)=PATMT(L)+ATMWHT(L,NA)*PATMTT(NA)
          TATMT(L)=TATMT(L)+ATMWHT(L,NA)*TATMTT(NA)
          RAINT(L)=RAINT(L)+ATMWHT(L,NA)*RAINTT(NA)
          EVAPT(L)=EVAPT(L)+ATMWHT(L,NA)*EVAPTT(NA)
          SOLSWRT(L)=SOLSWRT(L)+ATMWHT(L,NA)*SOLSWRTT(NA)
          CLOUDT(L)=CLOUDT(L)+ATMWHT(L,NA)*CLOUDTT(NA)
          SVPA(L)=SVPA(L)+ATMWHT(L,NA)*SVPAT(NA)
          RHA(L)=RHA(L)+ATMWHT(L,NA)*RHAT(NA)
          VPA(L)=VPA(L)+ATMWHT(L,NA)*VPAT(NA)
          CLEVAP(L)=CLEVAP(L)+ATMWHT(L,NA)*CLEVAPT(NA)
          CCNHTT(L)=CCNHTT(L)+ATMWHT(L,NA)*CCNHTTT(NA)
        ENDDO
      ENDDO
C
      IF(REVC.LT.0.)THEN
        DO L=2,LA
          CLEVAP(L)=1.E-3*(0.8+0.065*WINDST(L))
	    CLEVAPTMP=0.001*ABS(REVC)
          CLEVAP(L)=MAX(CLEVAP(L),CLEVAPTMP)
        ENDDO
	ELSE
        DO L=2,LA
	    CLEVAP(L)=0.001*ABS(REVC)
        ENDDO
      ENDIF
C
      IF(RCHC.LT.0.)THEN
        DO L=2,LA
          CCNHTT(L)=1.E-3*(0.8+0.065*WINDST(L))
	    CCNHTTTMP=0.001*ABS(RCHC)
          CCNHTT(L)=MAX(CCNHTT(L),CCNHTTTMP)
        ENDDO
	ELSE
        DO L=2,LA
	    CCNHTT(L)=0.001*ABS(RCHC)
        ENDDO
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
      IF(NISER.GT.0)THEN
C
      DO NI=1,NISER
        IF(ISDYNSTP.EQ.0)THEN
          TIME=DT*FLOAT(N)/TCASER(NI)+TBEGIN*(TCON/TCISER(NI))
        ELSE
          TIME=TIMESEC/TCISER(NA)
        ENDIF
        M1=MITLAST(NI)    
  500   CONTINUE
        M2=M1+1
        IF(TIME.GT.TISER(M2,NI))THEN
          M1=M2
          GOTO 500
        ELSE
          MITLAST(NI)=M1
        ENDIF      
        TDIFF=TISER(M2,NI)-TISER(M1,NI)
        WTM1=(TISER(M2,NI)-TIME)/TDIFF
        WTM2=(TIME-TISER(M1,NI))/TDIFF
        RICECOVT(NI)=WTM1*RICECOVS(M1,NI)+WTM2*RICECOVS(M2,NI)
        RICETHKT(NI)=WTM1*RICETHKS(M1,NI)+WTM2*RICETHKS(M2,NI)
	ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      DO L=2,LA
        RICECOVL(L)=0.
	  RICETHKL(L)=0.
	ENDDO

      IF(NISER.EQ.1)THEN
        DO L=2,LA
          RICECOVL(L)=RICECOVT(1)
	    RICETHKL(L)=RICETHKT(1)
	  ENDDO
	ENDIF
C
      IF(NISER.GT.1)THEN
        DO NI=1,NISER
          DO L=2,LA
            RICECOVL(L)=RICECOVL(L)+RICEWHT(L,NI)*RICECOVT(NI)
            RICETHKL(L)=RICETHKL(L)+RICEWHT(L,NI)*RICETHKT(NI)
	    ENDDO
        ENDDO
	ENDIF
C
      IF(NISER.GT.0)THEN
        DO L=2,LA
          ICECOVL(L)=NINT(RICECOVL(L))
	  ENDDO
        DO L=2,LA
          RICECOVL(L)=FLOAT(ICECOVL(L))
	  ENDDO
        DO L=2,LA
          IF(ICECOVL(L).EQ.0) RICETHKL(L)=0.0
	  ENDDO
	ENDIF
C
C**********************************************************************C
C
       IF(ISICE.EQ.1)THEN
        DO L=2,LA
	   TAUICE=-0.001*SQRT(U(L,KC)*U(L,KC)+V(L,KC)*V(L,KC))
         TSX(L)=RICECOVL(L)*TAUICE*U(L,KC)+(1.-RICECOVL(L))*TSX(L)
         TSY(L)=RICECOVL(L)*TAUICE*V(L,KC)+(1.-RICECOVL(L))*TSY(L)
	   WINDST(L)=(1.-RICECOVL(L))*WINDST(L)
        ENDDO
	 ENDIF
C
C**********************************************************************C
C
      RETURN
      END
