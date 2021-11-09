C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALHEATGVC(ISTL)
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
C **  SUBROUTINE CALHEAT CALCULATES SURFACE AND INTERNAL HEAT SOURCES 
C **  AND SINKS IN THE HEAT (TEM) TRANSPORT EQUATION
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION TEMOLD(LCM,KCM)
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF(ISTL.EQ.2)THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      ENDIF
C
C**********************************************************************C
C
      IF(ISTOPT(2).EQ.1)THEN
C
      DO K=1,KCM
      DO L=2,LA
        TEMOLD(L,K)=TEM(L,K)
      ENDDO
      ENDDO
C
      DO L=2,LA
        TVAR3S(L)=0.
      ENDDO
C
C **  ADSORB SW SOLR RAD TO ALL LAYERS AND BED
C
      IF(IASWRAD.EQ.0)THEN
C
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     &              (1.+0.00412*TEM(L,KC))))
C    &    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
CX        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)
     &  -(DELT*DZIC(KC)*GVCSCLPI(L))*( 1.312E-14*((TEM(L,KC)+273.)**4)
     &                  *(0.39-0.05*SQRT(VPA(L)))*(1.-.8*CLOUDT(L))
     &   +5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))
     &   +CCNHTT(L)*0.288E-3*WNDTMP*(TEM(L,KC)-TATMT(L))
     &   +CLEVAP(L)*0.445*WNDTMP*(SVPW-VPA(L))/PATMT(L) )
     &  +(DELT*DZIC(KC)*GVCSCLPI(L))*0.2385E-6*SOLSWRT(L)*(
     &     FSWRATF*EXP(SWRATNF*HDEP*GVCSCLP(L)*(Z(KC)-1.))
     &    +(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLP(L)*(Z(KC)-1.))
     &    -FSWRATF*EXP(SWRATNF*HDEP*GVCSCLP(L)*(Z(KC-1)-1.))
     &    -(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLP(L)*(Z(KC-1)-1.)) )
       ENDDO

c     &   +1.5*0.288E-6*WNDTMP*(TEM(L,KC)-TATMT(L))
c     &   +1.5*0.445E-3*WNDTMP*(SVPW-VPA(L))/PATMT(L) )

C
       IF(KC.GT.1)THEN
        DO K=1,KS
         DO L=2,LA
          HDEP=MAX(HP(L),0.)
          TVAR1S(L,K)=HDEP*TEM(L,K)
     &     +(DELT*DZIC(K)*GVCSCLPI(L))*0.2385E-6*SOLSWRT(L)*(
     &      FSWRATF*EXP(SWRATNF*HDEP*GVCSCLP(L)*(Z(K)-1.))
     &     +(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLP(L)*(Z(K)-1.))
     &     -FSWRATF*EXP(SWRATNF*HDEP*GVCSCLP(L)*(Z(K-1)-1.))
     &     -(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLP(L)*(Z(K-1)-1.)) )
         ENDDO
        ENDDO
       ENDIF
c
       DO K=1,KS
	  DO L=2,LA
	    IF(K.LT.KGVCP(L))TVAR1S(L,K)=0.0
        ENDDO
	 ENDDO
C
c adsorb remaining solar radiation into the bottom layer
c
       IF(DABEDT.GT.0.0)THEN
       DO L=2,LA
	  K=KGVCP(L)
        HDEP=MAX(HP(L),0.)
        TVAR1S(L,K)=TVAR1S(L,K)
     &     +(DELT*DZIC(K)*GVCSCLPI(L))*0.2385E-6*SOLSWRT(L)*(
     &      FSWRATF*EXP(SWRATNF*HDEP*GVCSCLP(L)*(Z(K-1)-1.))
     &     +(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLP(L)*(Z(K-1)-1.)) )
       ENDDO
       ENDIF

C
c relax bottom layer temperature toward deep bed
c
C       IF(HTBED2.GT.0.0)THEN
C       DO L=2,LA
C	  K=KGVCP(L)
C        IF(HP(L).GT.0.)THEN
C          TVAR1S(L,K)=TVAR1S(L,K)
C     &     +(DELT*DZIC(K)*GVCSCLPI(L))*HTBED2*
C     &       (TEMB(L,KBHM)-(TVAR1S(L,K)*HPI(L)))
C        ENDIF
C       ENDDO
C       ENDIF
C
c
c       DO L=2,LA
c	  K=KGVCP(L)
c        UBED=0.5*( U(L,K)+U(L+1,K) )
c        VBED=0.5*( V(L,K)+V(LNC(L),K) )
c        USPD=SQRT( UBED*UBED+VBED*VBED )
c        TMPVAL=(HTBED1*USPD+HTBED2)*(TEM(L,K)-TEMB(L))
c        TVAR1S(L,K)=TVAR1S(L,K)-DELT*DZIC(K)*GVCSCLPI(L)*TMPVAL
c        TEMB(L)=TEMB(L)+DELT*TMPVAL/DABEDT
c     &     +(DELT/DABEDT)*0.2385E-6*SOLSWRT(L)*(
c     &     +FSWRATF*EXP(SWRATNF*HDEP*GVCSCLPI(L)*(Z(K-1)-1.))
c     &     +(1.-FSWRATF)*EXP(SWRATNS*HDEP*GVCSCLPI(L)*(Z(K-1)-1.)) )
c       ENDDO
C
       DO K=1,KC
        DO L=2,LA
         IF(HP(L).GT.0.) TEM(L,K)=HPI(L)*TVAR1S(L,K)
	   TEM(L,K)=MAX(TEM(L,K),0.0)
        ENDDO
       ENDDO
C
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
C
       IF(ISDRY.GT.0)THEN
         DO K=1,KC
           DO L=2,LA
             IF(IMASKDRY(L).EQ.1) TEM(L,K)=TATMT(L)
           ENDDO
         ENDDO
C         IF(HTBED2.EQ.0.0)THEN
C           DO L=2,LA
C             IF(IMASKDRY(L).EQ.1) TEMB(L,KBHM)=TATMT(L)
C           ENDDO
C         ENDIF
       ENDIF
C      
       NTSTBCM1=NTSTBC-1
	 IF(DABEDT.GT.0.0) THEN
         IF(IS2TIM.EQ.1) CALL CALHEATBGVC(ISTL)	
	   IF(IS2TIM.EQ.0) THEN
	     IF(NCTBC.EQ.NTSTBCM1) CALL CALHEATBGVC(ISTL)	
         ENDIF
       ENDIF
C
      ENDIF
C
C
C **  ADSORB SW SOLR RAD TO SURFACE LAYER
C
      IF(IASWRAD.EQ.1)THEN
C
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     &              (1.+0.00412*TEM(L,KC))))
C    &    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
CX        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)
     &  -(DELT*DZIC(KC)*GVCSCLPI(L))*( 1.312E-14*((TEM(L,KC)+273.)**4)
     &         *(0.39-0.05*SQRT(VPA(L)))*(1.-.8*CLOUDT(L))
     &   +5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))
     &   +CCNHTT(L)*0.288E-3*WNDTMP*(TEM(L,KC)-TATMT(L))
     &   +CLEVAP(L)*0.445*WNDTMP*(SVPW-VPA(L))/PATMT(L) )
     &  +(DELT*DZIC(KC)*GVCSCLPI(L))*0.2385E-6*SOLSWRT(L)
       ENDDO
C
c relax bottom layer temperature toward deep bed
c
c       IF(HTBED2.GT.0.0)THEN
c       DO L=2,LA
c	  K=KGVCP(L)
c        IF(HP(L).GT.0.)THEN
c          TVAR1S(L,K)=TVAR1S(L,K)
c     &     +(DELT*DZIC(K)*GVCSCLPI(L))*HTBED2*
c     &       (TEMB(L,KBHM)-(TVAR1S(L,K)*HPI(L)))
c        ENDIF
c       ENDDO
c       ENDIF

       DO L=2,LA
        IF(HP(L).GT.0.) TEM(L,KC)=HPI(L)*TVAR1S(L,KC)
	  TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
C
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
C
       IF(ISDRY.GT.0)THEN
         DO K=1,KC
           DO L=2,LA
             IF(IMASKDRY(L).EQ.1.) TEM(L,K)=TATMT(L)
           ENDDO
         ENDDO
c         IF(HTBED2.EQ.0.0)THEN
c           DO L=2,LA
c             IF(IMASKDRY(L).EQ.1) TEMB(L,KBHM)=TATMT(L)
c           ENDDO
c         ENDIF
       ENDIF
C
       NTSTBCM1=NTSTBC-1
	 IF(DABEDT.GT.0.0) THEN
         IF(IS2TIM.EQ.1) CALL CALHEATBGVC(ISTL)	
	   IF(IS2TIM.EQ.0) THEN
	     IF(NCTBC.EQ.NTSTBCM1) CALL CALHEATBGVC(ISTL)	
         ENDIF
       ENDIF
C
      ENDIF
C
C
C **  ENDIF ISOPT(2) EQ 1
C
      ENDIF
C
  600 FORMAT(4I5,2E12.4)
C
C**********************************************************************C
C
      IF(ISTOPT(2).EQ.2)THEN
C
C ** IMPLEMENT EXTERNALLY SPECIFIED EQUILIBRIUM TEMPERATURE FROMULATION
C
      TMPKC=DELT/DZC(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-TMPKC*CLOUDT(L)*HPI(L)*(TEM(L,KC)-TATMT(L))
	  TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
      ENDDO
C
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
C

      ENDIF
C
C**********************************************************************C
C
      IF(ISTOPT(2).EQ.3)THEN
C
C ** IMPLEMENT CONSTANT COEFFICIENT EQUILIBRIUM TEMPERATURE FROMULATION
C
      DTHEQT=DELT*HEQT*FLOAT(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-DTHEQT*HPI(L)*(TEM(L,KC)-TEMO)
	  TEM(L,KC)=MAX(TEM(L,KC),0.0)
       ENDDO
      ENDDO
C
       IF(ISICE.EQ.1)THEN
       DO K=1,KC
        DO L=2,LA
         TEM(L,K)=RICECOVL(L)*TEMPICE+(1.-RICECOVL(L))*TEM(L,K)
        ENDDO
       ENDDO
	 ENDIF
C
      ENDIF
C
C**********************************************************************C
C
      RETURN
      END
