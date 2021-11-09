C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALUVWGVC (ISTL,IS2TL)
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
C **  CALCULATE THE INTERNAL SOLUTION AT TIME LEVEL (N+1)
C **  THE VALUE OF ISTL INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C
C      IF(ISTL.EQ.3)THEN
C       DELT=DT2
C       DELTD2=DT
C      ELSE
C       DELT=DT
C       DELTD2=0.5*DT
C      ENDIF
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
      IF(ISCDMA.EQ.3)THEN
        DELT=0.5*DELT
        DELTD2=0.5*DELTD2
        DELTI=1./DELT
      ENDIF
C
      IF(ISCDMA.EQ.4)THEN
        DELT=0.5*DELT
        DELTD2=0.5*DELTD2
        DELTI=1./DELT
      ENDIF
C
      IF(KC.EQ.1) GOTO 30
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM FRICTION COEFFICIENT
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
       DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        RCX(L)=AVCON1/H1U(L)+STBX(L)*SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)
     &        +V1U(L)*V1U(L))
        RCY(L)=AVCON1/H1V(L)+STBY(L)*SQRT(U1V(L)*U1V(L)
     &        +V1(L,KBTMPV)*V1(L,KBTMPV))
       ENDDO
c
CGVCDIAG
c      IF(N.LE.4)THEN
c      DO L=2,LA
c	  KBTMPU=KGVCU(L)
c	  WRITE(8,8881)N,IL(L),JL(L),KBTMPU,RCX(L),H1U(L),STBX(L),
c     &               U1(L,KBTMPU),V1U(L)
c      ENDDO
c      ENDIF
c 8881 FORMAT('KBU,RCX,H1U,ST,U,V GVC ',4I5,10E14.6)
CGVCDIAG
C
C      IF(ISVEG.GE.1.AND.KC.GT.1)
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C        RCX(L)=RCX(L)+0.5*FXVEG(L,1)*SQRT(U1(L,1)*U1(L,1)
C     &        +V1U(L)*V1U(L))
C        RCY(L)=RCY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*U1V(L)
C     &        +V1(L,1)*V1(L,1))
C       ENDDO
C      ENDDO
C      ENDIF
C
      ELSE
C
       IF(AVCON1.LT.0.00001)THEN
       DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        Q1=SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,KBTMPU)*U(L,KBTMPU)+VU(L)*VU(L))
        RCX(L)=STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,KBTMPV)*V1(L,KBTMPV))
        Q2=SQRT(UV(L)*UV(L)+V(L,KBTMPV)*V(L,KBTMPV))
        RCY(L)=STBY(L)*SQRT(Q1*Q2)
       ENDDO
	 ELSE
       DO L=2,LA
	  KBTMPU=KGVCU(L)
	  KBTMPV=KGVCV(L)
        Q1=SQRT(U1(L,KBTMPU)*U1(L,KBTMPU)+V1U(L)*V1U(L))
        Q2=SQRT(U(L,KBTMPU)*U(L,KBTMPU)+VU(L)*VU(L))
        RCX(L)=AVCON1/SQRT(H1U(L)*HU(L))+STBX(L)*SQRT(Q1*Q2)
        Q1=SQRT(U1V(L)*U1V(L)+V1(L,KBTMPV)*V1(L,KBTMPV))
        Q2=SQRT(UV(L)*UV(L)+V(L,KBTMPV)*V(L,KBTMPV))
        RCY(L)=AVCON1/SQRT(H1V(L)*HV(L))+STBY(L)*SQRT(Q1*Q2)
       ENDDO
	 ENDIF
C
C      IF(ISVEG.GE.1.AND.KC.GT.1)
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C        Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))
C        Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))
C        RCX(L)=RCX(L)+0.5*FXVEG(L,1)*SQRT(Q1*Q2)
C        Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))
C        Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
C        RCY(L)=RCY(L)+0.5*FYVEG(L,1)*SQRT(Q1*Q2)
C       ENDDO
C      ENDDO
C      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE THE U AND V SHEARS
C
C----------------------------------------------------------------------C
C
c **  working definitions
c
c           CDZM(K)=0.5*DZC(K)*DZC(K+1)
c           CDZU(K)=-DZC(K)/(DZC(K)+DZC(K+1))
c           CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
c
c
c           CDZR(1)=DZC(1)-1.
c           CDZD(1)=DZC(1)
c           DO K=2,KS
c           CDZR(K)=DZC(K)+CDZR(K-1)
c           CDZD(K)=DZC(K)+CDZD(K-1)
c           ENDDO
c
c           DO K=1,KS
c           CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
c           ENDDO
C
c
c **  first interface at k=1, cases for null and kbu/kbv =1
c
CGVCDIAG
c      IF(N.LE.4)THEN
c      DO L=2,LA
c	  WRITE(8,881)N,IL(L),JL(L),RCX(L),UHE(L),(DU(L,K),K=1,KS)
c      ENDDO
c      ENDIF
c  881 FORMAT('R,UHE,DU GVC ',3I5,10E14.6)
CGVCDIAG

C
      RCDZM=CDZM(1)*DELTI
      RCDZU=CDZU(1)
      RCDZL=CDZL(1)
       DO L=2,LA
	   IF(LGVCU(L,1))THEN
           CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,1)
           EU=1./CMU
           CU1(L,1)=RCDZU*EU
           DU(L,1)=(DU(L,1)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU
           UUU(L,1)=EU
         ELSE
           CMU=1.
C           EU=0.
           CU1(L,1)=0.
           DU(L,1)=0.
           UUU(L,1)=0.
         ENDIF
       ENDDO
       DO L=2,LA
	   IF(LGVCV(L,1))THEN
           CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,1)
           EV=1./CMV
           CU2(L,1)=RCDZU*EV
           DV(L,1)=(DV(L,1)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV
           VVV(L,1)=EV
         ELSE
           CMV=1.
C           EV=0.
           CU2(L,1)=0.
           DV(L,1)=0.
           VVV(L,1)=0.
         ENDIF
       ENDDO
C
       DO K=2,KS
        RCDZM=CDZM(K)*DELTI
        RCDZU=CDZU(K)
        RCDZL=CDZL(K)
        DO L=2,LA
	    IF(LGVCU(L,K))THEN
	      IF(K.EQ.KGVCU(L))THEN
              CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,K)
              EU=1./(CMU-RCDZL*CU1(L,K-1))
              CU1(L,K)=RCDZU*EU
              DU(L,K)=(DU(L,K)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU
              UUU(L,K)=EU  
		  ELSE
              CMU=1.+RCDZM*GVCSCLU(L)*HU(L)*AVUI(L,K)
              EU=1./(CMU-RCDZL*CU1(L,K-1))
              CU1(L,K)=RCDZU*EU
              DU(L,K)=(DU(L,K)-RCDZL*DU(L,K-1))*EU
              UUU(L,K)=-RCDZL*UUU(L,K-1)*EU  
		  ENDIF    
          ELSE
            CMU=1.
C            EU=0.
            CU1(L,K)=0.
            DU(L,K)=0.
            UUU(L,K)=0.     
          ENDIF
        ENDDO
        DO L=2,LA
	    IF(LGVCV(L,K))THEN
	      IF(K.EQ.KGVCV(L))THEN
              CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,K)
              EV=1./(CMV-RCDZL*CU2(L,K-1))
              CU2(L,K)=RCDZU*EV
              DV(L,K)=(DV(L,K)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV
              VVV(L,K)=EV      
		  ELSE
              CMV=1.+RCDZM*GVCSCLV(L)*HV(L)*AVVI(L,K)
              EV=1./(CMV-RCDZL*CU2(L,K-1))
              CU2(L,K)=RCDZU*EV
              DV(L,K)=(DV(L,K)-RCDZL*DV(L,K-1))*EV
              VVV(L,K)=-RCDZL*VVV(L,K-1)*EV      
		  ENDIF    
          ELSE
            CMV=1.
C            EV=0.
            CU2(L,K)=0.
            DV(L,K)=0.
            VVV(L,K)=0.     
          ENDIF
        ENDDO
       ENDDO
C
c **  back substitution
c
       DO K=KS-1,1,-1
        DO L=2,LA
         DU(L,K)=DU(L,K)-CU1(L,K)*DU(L,K+1)
         DV(L,K)=DV(L,K)-CU2(L,K)*DV(L,K+1)
         UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
         VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
        ENDDO
       ENDDO
C
CGVCDIAG
c      IF(N.LE.4)THEN
c      DO L=2,LA
c	  WRITE(8,884)N,IL(L),JL(L),(DU(L,K),K=1,KS)
c      ENDDO
c      DO L=2,LA
c	  WRITE(8,885)N,IL(L),JL(L),(UUU(L,K),K=1,KS)
c      ENDDO
c      ENDIF
c  884 FORMAT('BKSUB DU  GVC ',3I5,10E14.6)
c  885 FORMAT('BKSUB UUU GVC ',3I5,10E14.6)
CGVCDIAG
C
c **  set SHERMAN-MORRISON BACK SUBSTITUTION
c
       DO L=2,LA
        AAU(L)=0.
        AAV(L)=0.
        BBU(L)=1.
        BBV(L)=1.
       ENDDO
c
c           CDZR(1)=DZC(1)-1.
c
c           DO K=2,KS
c           CDZR(K)=DZC(K)+CDZR(K-1)
c           ENDDO
c
c           DO K=1,KS
c           CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
c           ENDDO
c
c           CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
c
       DO K=1,KS
        DO L=2,LA
	   IF(LGVCU(L,K))THEN
           CRU=CDZRGVCU(L,K)*RCX(L)*AVUI(L,K)
           AAU(L)=AAU(L)+CRU*DU(L,K)
           BBU(L)=BBU(L)+CRU*UUU(L,K)
	   ENDIF
        ENDDO
       ENDDO
c
       DO K=1,KS
        DO L=2,LA
	   IF(LGVCV(L,K))THEN
           CRV=CDZRGVCV(L,K)*RCY(L)*AVVI(L,K)
           AAV(L)=AAV(L)+CRV*DV(L,K)
           BBV(L)=BBV(L)+CRV*VVV(L,K)
	   ENDIF
        ENDDO
       ENDDO
C
       DO L=2,LA
        AAU(L)=AAU(L)/BBU(L)
        AAV(L)=AAV(L)/BBV(L)
       ENDDO
C
      DO K=1,KS
       RDZG=DZG(K)
       DO L=2,LA
        DU(L,K)=RDZG*HU(L)*AVUI(L,K)*(DU(L,K)-AAU(L)*UUU(L,K)) 
        DV(L,K)=RDZG*HV(L)*AVVI(L,K)*(DV(L,K)-AAV(L)*VVV(L,K)) 
       ENDDO
      ENDDO
C 
CGVCDIAG
c      IF(N.LE.4)THEN
c      DO L=2,LA
c	  WRITE(8,882)N,IL(L),JL(L),(DU(L,K),K=1,KS)
c      ENDDO
c      ENDIF
c  882 FORMAT('SMW DU GVC ',3I5,10E14.6)
CGVCDIAG
C
c     note that du and dv are actually the depth time the velocity 
c     difference
C
C**********************************************************************C
C
C **  CALCULATED U AND V       
C
C **  DUSUM+UHE=UHE, DVSUM+VHE=VHE
C
C----------------------------------------------------------------------C
C
C
C        CDZD(1)=DZC(1)
C         DO K=2,KS
C         CDZD(K)=DZC(K)+CDZD(K-1)
C         ENDDO
C
C	   CDZDGVCU(L,K)=GVCSCLU(L)*(CDZD(K) FIXED FOR GVC)
C	   CDZDGVCU(L,K)=GVCSCLV(L)*(CDZD(K) FIXED FOR GVC)
C
c     solve for the surface layer velocity * depth and temporarily
c     store results in uhe and vhe
c
       DO K=1,KS
C        RCDZD=CDZD(K)
        DO L=2,LA
C	   IF(LGVCU(L,K)) UHE(L)=UHE(L)+RCDZD*GVCSCLU(L)*DU(L,K)
	   IF(LGVCU(L,K)) UHE(L)=UHE(L)+CDZDGVCU(L,K)*DU(L,K)
        ENDDO
       ENDDO
c
       DO K=1,KS
C        RCDZD=CDZD(K)
        DO L=2,LA
C         IF(LGVCV(L,K)) VHE(L)=VHE(L)+RCDZD*GVCSCLV(L)*DV(L,K)
         IF(LGVCV(L,K)) VHE(L)=VHE(L)+CDZDGVCV(L,K)*DV(L,K)
        ENDDO
       ENDDO
C
C      move surface layer velocity times depth to uhdy and vhdx
C
       DO L=2,LA
        UHDY(L,KC)=UHE(L)*SUB3D(L,KC)
        VHDX(L,KC)=VHE(L)*SVB3D(L,KC)
       ENDDO
C
       DO K=KS,1,-1
        DO L=2,LA
         IF(LGVCU(L,K))THEN
	     UHDY(L,K)=SUB3D(L,K)*(UHDY(L,K+1)-DU(L,K))
         ELSE
	     UHDY(L,K)=0.0
         ENDIF
        ENDDO
       ENDDO
C
       DO K=KS,1,-1
        DO L=2,LA
         IF(LGVCV(L,K))THEN
	     VHDX(L,K)=SVB3D(L,K)*(VHDX(L,K+1)-DV(L,K))
         ELSE
	     VHDX(L,K)=0.0
         ENDIF
        ENDDO
       ENDDO
C
      DO K=1,KC
       DO L=2,LA
        U(L,K)=SUB3D(L,K)*UHDY(L,K)*HUI(L)
        V(L,K)=SVB3D(L,K)*VHDX(L,K)*HVI(L)
        UHDY(L,K)=SUB3D(L,K)*UHDY(L,K)*DYU(L)
        VHDX(L,K)=SVB3D(L,K)*VHDX(L,K)*DXV(L)
       ENDDO
      ENDDO
C
CGVCDIAG
c      IF(N.LE.4)THEN
c      DO L=2,LA
c	  WRITE(8,883)N,IL(L),JL(L),(UHDY(L,K),K=1,KC)
c      ENDDO
c      ENDIF
c  883 FORMAT('UHDY GVC ',3I5,10E14.6)
CGVCDIAG
C
C
C **  ADD ADJUSTMENT TO 3D HORIZONTAL TRANSPORT 
C
      DO L=2,LA
        TVAR3E(L)=0.
        TVAR3N(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
       ENDDO
      ENDDO
C
      UERMX=-1.E+12
      UERMN=1.E+12
      VERMX=-1.E+12
      VERMN=1.E+12
      DO L=2,LA
        TVAR3E(L)=TVAR3E(L)-UHDYE(L)
        TVAR3N(L)=TVAR3N(L)-VHDXE(L)
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        IF(LGVCU(L,K)) UHDY(L,K)=UHDY(L,K)-GVCSCLUI(L)*TVAR3E(L)*DZIC(K)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        IF(LGVCV(L,K)) VHDX(L,K)=VHDX(L,K)-GVCSCLVI(L)*TVAR3N(L)*DZIC(K)
       ENDDO
      ENDDO
C

C **  RESET VELOCITIES
C
      DO L=2,LA
       UHE(L)=0.
       VHE(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        UHE(L)=UHE(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
        VHE(L)=VHE(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
        U(L,K)=UHDY(L,K)*HUI(L)
        V(L,K)=VHDX(L,K)*HVI(L)
       ENDDO
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        U(L,K)=U(L,K)*DYIU(L)
        V(L,K)=V(L,K)*DXIV(L)
       ENDDO
      ENDDO
C
      DO L=2,LA
       UHE(L)=UHE(L)*DYIU(L)
       VHE(L)=VHE(L)*DXIV(L)
      ENDDO
C
C **  UNCOMMENT BELOW TO WRITE CONTINUITY DIAGNOSITCS
C
c      DO L=2,LA
c        TVAR3E(L)=0.
c        TVAR3N(L)=0.
c      ENDDO
C
c      DO K=1,KC
c       DO L=2,LA
c        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY(L,K)*DZC(K)
c        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX(L,K)*DZC(K)
c       ENDDO
c      ENDDO
C
c      UERMX=-1.E+12
c      UERMN=1.E+12
c      VERMX=-1.E+12
c      VERMN=1.E+12
c      DO L=2,LA
c        TVAR3E(L)=TVAR3E(L)-UHDYE(L)
c        TVAR3N(L)=TVAR3N(L)-VHDXE(L)
c        IF(TVAR3E(L).GT.UERMX)THEN
c          UERMX=TVAR3E(L)
c          LUMX=L
c        ENDIF
c        IF(TVAR3E(L).LT.UERMN)THEN
c          UERMN=TVAR3E(L)
c          LUMN=L
c        ENDIF
c        IF(TVAR3N(L).GT.VERMX)THEN
c          VERMX=TVAR3N(L)
c          LVMX=L
c        ENDIF
c        IF(TVAR3N(L).LT.VERMN)THEN
c          VERMN=TVAR3N(L)
c          LVMN=L
c        ENDIF
c      ENDDO
C
c      WRITE(8,6661)IL(LUMX),JL(LUMX),UERMX
c      WRITE(8,6662)IL(LUMN),JL(LUMN),UERMN
c      WRITE(8,6663)IL(LVMX),JL(LVMX),VERMX
c      WRITE(8,6664)IL(LVMN),JL(LVMN),VERMN
C
 6660 FORMAT(' INT MOD ERR AT = ',I8)
 6661 FORMAT(' I,J,UHDYERMX = ',2I5,10E12.4)
 6662 FORMAT(' I,J,UHDYERMN = ',2I5,10E12.4)
 6663 FORMAT(' I,J,VHDXERMX = ',2I5,10E12.4)
 6664 FORMAT(' I,J,VHDXERMX = ',2I5,10E12.4)
 6665 FORMAT(' I,J,UMXSTUFF = ',2I5,10E12.4)
C
C**********************************************************************C
C
C **  CALCULATE W
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3)THEN
C
       DO L=2,LA
	  LN=LNC(L)
	  LE=L+1
        TVAR3E(L)=GVCSCLP(L)*UHDYE(LE)
        TVAR3N(L)=GVCSCLP(L)*VHDXE(LN)
        TVAR3W(L)=GVCSCLP(L)*UHDY2E(LE)
        TVAR3S(L)=GVCSCLP(L)*VHDX2E(LN)
       ENDDO
C
c      mode here

       DO K=1,KC
        DO L=2,LA
	   LN=LNC(L)
	   LE=L+1
         TVAR1E(L,K)=GVCSCLU(LE)*UHDY(LE,K)
         TVAR1N(L,K)=GVCSCLV(LN)*VHDX(LN,K)
         TVAR1W(L,K)=GVCSCLU(LE)*UHDY2(LE,K)
         TVAR1S(L,K)=GVCSCLV(LN)*VHDX2(LN,K)
        ENDDO
       ENDDO
C
c020906       DO K=1,KS
c        DO L=2,LA
c	   IF(K.GE.KGVCP(L))THEN
c         LN=LNC(L)
c         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
c     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
c     &          -TVAR3E(L)+GVCSCLU(L)*UHDYE(L)
c     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)
c     &          -TVAR3W(L)+GVCSCLU(L)*UHDY2E(L)
c     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
c     &          -TVAR3N(L)+GVCSCLV(L)*VHDXE(L)
c     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)
c     &          -TVAR3S(L)+GVCSCLV(L)*VHDX2E(L) 
c     &          )*DXYIP(L) )
c     &          +2.*SWB3D(L,K)*( QSUM(L,K)
c     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
c        ENDIF
c        ENDDO
c       ENDDO
C
       DO K=1,KS
        DO L=2,LA
	   IF(K.GE.KGVCP(L))THEN
         LN=LNC(L)
         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
     &          -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)
     &          -TVAR3W(L)+GVCSCLP(L)*UHDY2E(L)
     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
     &          -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)
     &          -TVAR3S(L)+GVCSCLP(L)*VHDX2E(L) 
     &          )*DXYIP(L) )
     &          +2.*SWB3D(L,K)*( QSUM(L,K)
     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
        ENDIF
        ENDDO
       ENDDO
C
C **  UNCOMMENT BELOW FOR CONTINUITY DIAGNOSTICS
C
C      WSFMAX=-1.E+12
C      WSFMIN=1.E+12
C      SURFOUT=0.
C      K=KC
C      DO L=2,LA
C      IF(SWB(L).GT.0.5)THEN
C         WSURF=SWB(L)*( W(L,K-1)-W2(L,K)+W2(L,K-1)
C     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
C     &   -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
C     &       +TVAR1W(L,K)-GVCSCLU(L)*UHDY2(L,K)
C     &   -TVAR3W(L)+GVCSCLP(L)*UHDY2E(L)
C     &       +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
C     &   -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
C     &       +TVAR1S(L,K)-GVCSCLV(L)*VHDX2(L,K)
C     &   -TVAR3S(L)+GVCSCLP(L)*VHDX2E(L) )*DXYIP(L) )
C     &       +2.*SWB(L)*( QSUM(L,K)-
C     &   DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
C      SURFOUT=SURFOUT+WSURF*DXYP(L)
C      IF(WSURF.GT.WSFMAX)THEN
C        WSFMAX=WSURF
C        IMAX=IL(L)
C        JMAX=JL(L)
C      ENDIF
C      IF(WSURF.LT.WSFMIN)THEN
C        WSFMIN=WSURF
C        IMIN=IL(L)
C        JMIN=JL(L)
C      ENDIF
C      ENDIF
C      ENDDO
C
C      L=LIJ(IMAX,JMAX)
C      WSFMAX=WSFMAX*DXYP(L)
C      L=LIJ(IMIN,JMIN)
C      WSFMIN=WSFMIN*DXYP(L)
C      WRITE(8,601)IMAX,JMAX,WSFMAX
C      WRITE(8,602)IMIN,JMIN,WSFMIN
C      WRITE(8,603)SURFOUT
C
      ENDIF
C
      IF(ISTL.EQ.2)THEN
C
       DO L=2,LA
	  LN=LNC(L)
	  LE=L+1
        TVAR3E(L)=GVCSCLP(L)*UHDYE(LE)
        TVAR3N(L)=GVCSCLP(L)*VHDXE(LN)
        TVAR3W(L)=GVCSCLP(L)*UHDY1E(LE)
        TVAR3S(L)=GVCSCLP(L)*VHDX1E(LN)
       ENDDO
C
       DO K=1,KC
        DO L=2,LA
	   LN=LNC(L)
	   LE=L+1
         TVAR1E(L,K)=GVCSCLU(LE)*UHDY(LE,K)
         TVAR1N(L,K)=GVCSCLV(LN)*VHDX(LN,K)
         TVAR1W(L,K)=GVCSCLU(LE)*UHDY1(LE,K)
         TVAR1S(L,K)=GVCSCLV(LN)*VHDX1(LN,K)
        ENDDO
       ENDDO
C
c020906       DO K=1,KS
c        DO L=2,LA
c	   IF(K.GE.KGVCP(L))THEN
c         LN=LNC(L)
c         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
c     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
c     &          -TVAR3E(L)+GVCSCLU(L)*UHDYE(L)
c     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)
c     &          -TVAR3W(L)+GVCSCLU(L)*UHDY1E(L)
c     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
c     &          -TVAR3N(L)+GVCSCLV(L)*VHDXE(L)
c     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)
c     &          -TVAR3S(L)+GVCSCLV(L)*VHDX1E(L)
c     &          )*DXYIP(L) )
c     &          +2.*SWB3D(L,K)*( QSUM(L,K)
c     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
c         ENDIF
c        ENDDO
c       ENDDO
C
       DO K=1,KS
        DO L=2,LA
	   IF(K.GE.KGVCP(L))THEN
         LN=LNC(L)
         W(L,K)=SWB3D(L,K)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
     &          -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
     &          +TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)
     &          -TVAR3W(L)+GVCSCLP(L)*UHDY1E(L)
     &          +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
     &          -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
     &          +TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)
     &          -TVAR3S(L)+GVCSCLP(L)*VHDX1E(L)
     &          )*DXYIP(L) )
     &          +2.*SWB3D(L,K)*( QSUM(L,K)
     &          -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
         ENDIF
        ENDDO
       ENDDO
C
C  UNCOMMENT BELOW FOR SURFACE LAYER CONTINUITY DIAGNOSTIC
C
c      WSFMAX=-1.E+12
c      WSFMIN=1.E+12
c      SURFOUT=0.
c      K=KC
c      DO L=2,LA
c      IF(SWB(L).GT.0.5)THEN
c         WSURF=SWB(L)*( W(L,K-1)-W1(L,K)+W1(L,K-1)
c     &  -DZC(K)*(TVAR1E(L,K)-GVCSCLU(L)*UHDY(L,K)
c     &   -TVAR3E(L)+GVCSCLP(L)*UHDYE(L)
c     &       +TVAR1W(L,K)-GVCSCLU(L)*UHDY1(L,K)
c     &   -TVAR3W(L)+GVCSCLP(L)*UHDY1E(L)
c     &       +TVAR1N(L,K)-GVCSCLV(L)*VHDX(L,K)
c     &   -TVAR3N(L)+GVCSCLP(L)*VHDXE(L)
c     &       +TVAR1S(L,K)-GVCSCLV(L)*VHDX1(L,K)
c     &   -TVAR3S(L)+GVCSCLP(L)*VHDX1E(L) )*DXYIP(L) )
c     &       +2.*SWB(L)*( QSUM(L,K)
c     &   -DZC(K)*GVCSCLP(L)*QSUME(L) )*DXYIP(L)
c      SURFOUT=SURFOUT+WSURF*DXYP(L)
c      IF(WSURF.GT.WSFMAX)THEN
c        WSFMAX=WSURF
c        IMAX=IL(L)
c        JMAX=JL(L)
c      ENDIF
c      IF(WSURF.LT.WSFMIN)THEN
c        WSFMIN=WSURF
c        IMIN=IL(L)
c        JMIN=JL(L)
c      ENDIF
c      ENDIF
c      ENDDO
C
c      L=LIJ(IMAX,JMAX)
c      WSFMAX=WSFMAX*DXYP(L)
c      L=LIJ(IMIN,JMIN)
c      WSFMIN=WSFMIN*DXYP(L)
c      WRITE(8,601)IMAX,JMAX,WSFMAX
c      WRITE(8,602)IMIN,JMIN,WSFMIN
c      WRITE(8,603)SURFOUT
C
      ENDIF
C
  601 FORMAT(' IMAX,JMAX,QWSFMAX = ',2I5,E14.5)      
  602 FORMAT(' IMIN,JMIN,QWSFMIN = ',2I5,E14.5)      
  603 FORMAT(' TOTAL SURF Q ERR = ',E14.5)      
C
C**********************************************************************C
C
C **  CALCULATE U AND V ON OPEN BOUNDARIES
C
   30 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.EQ.0)THEN
C
      DO K=1,KC
      DO LL=1,NCBS
	IF(ISPBS(LL).LE.1)THEN
      L=LCBS(LL)      
      LN=LNC(L)
      LNN=LNC(LN)
      IF(LN.NE.LC)THEN
        VHDX(LN,K)=SVB3D(LN,K)*(VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN))
        V(LN,K)=VHDX(LN,K)/(HV(LN)*DXV(LN))
       ELSE
        VHDX(LN,K)=0.
        V(LN,K)=0.
      ENDIF
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBW
	IF(ISPBW(LL).LE.1)THEN
      L=LCBW(LL)
      LP=L+1
      LPP=L+2
      IF(LP.NE.LC)THEN      
        UHDY(LP,K)=SUB3D(LP,K)*(UHDY(LPP,K)-UHDYE(LPP)+UHDYE(LP))
        U(LP,K)=UHDY(LP,K)/(HU(LP)*DYU(LP))
       ELSE
        UHDY(LP,K)=0.
        U(LP,K)=0.
      ENDIF
      ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBE
	IF(ISPBE(LL).LE.1)THEN
      L=LCBE(LL)      
      UHDY(L,K)=SUB3D(L,K)*(UHDY(L-1,K)-UHDYE(L-1)+UHDYE(L))
      U(L,K)=UHDY(L,K)/(HU(L)*DYU(L))
	ENDIF
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBN
	IF(ISPBN(LL).LE.1)THEN
      L=LCBN(LL)      
      LS=LSC(L)
      VHDX(L,K)=SVB3D(L,K)*(VHDX(LS,K)-VHDXE(LS)+VHDXE(L))
      V(L,K)=VHDX(L,K)/(HV(L)*DXV(L))
	ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C
      IF(IS1DCHAN.GE.1)THEN
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)      
      LN=LNC(L)
      LNN=LNC(LN)
      VHDX(LN,K)=VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN)
      V(LN,K)=VHDX(LN,K)/FADXV(LN)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      LN=LNC(L)
      UHDY(L+1,K)=UHDY(L+2,K)-UHDYE(L+2)+UHDYE(L+1)
      U(L+1,K)=UHDY(L+1,K)/FADYU(L+1)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)      
      UHDY(L,K)=UHDY(L-1,K)-UHDYE(L-1)+UHDYE(L)
      U(L,K)=UHDY(L,K)/FADYU(L)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)      
      LS=LSC(L)
      VHDX(L,K)=VHDX(LS,K)-VHDXE(LS)+VHDXE(L)
      V(L,K)=VHDX(L,K)/FADXV(L)
      ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  CALCULATE AVERAGE CELL FACE TRANSPORTS FOR SALT, TEMPERATURE AND 
C **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
C
C----------------------------------------------------------------------C
C
      IF(ISCDMA.LE.4)THEN
      IF(ISTL.EQ.2)THEN
C
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))
          U2(L,K)=0.5*(U(L,K)+U1(L,K))
          V2(L,K)=0.5*(V(L,K)+V1(L,K))
          W2(L,K)=0.5*(W(L,K)+W1(L,K))
         ENDDO
        ENDDO
C
       ELSE
C
        DO K=1,KC
         DO L=2,LA
C         DU(L,K)=0.25*(UHDY(L,K)-UHDY2(L,K))
C         DV(L,K)=0.25*(VHDX(L,K)-VHDX2(L,K))
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY2(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX2(L,K))
          U2(L,K)=0.5*(U(L,K)+U2(L,K))
          V2(L,K)=0.5*(V(L,K)+V2(L,K))
C         DW(L,K)=0.25*(W(L,K)-W2(L,K))
          W2(L,K)=0.5*(W(L,K)+W2(L,K))
         ENDDO
        ENDDO
C
      ENDIF
      ENDIF
C
C
      IF(ISCDMA.GE.5)THEN
      IF(ISTL.EQ.2)THEN
C
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))
          U2(L,K)=0.5*(U(L,K)+U1(L,K))
          V2(L,K)=0.5*(V(L,K)+V1(L,K))
          W2(L,K)=0.5*(W(L,K)+W1(L,K))
         ENDDO
        ENDDO
C
       ELSE
C
        DO K=1,KC
         DO L=2,LA
C         DU(L,K)=0.25*(UHDY(L,K)-UHDY2(L,K))
C         DV(L,K)=0.25*(VHDX(L,K)-VHDX2(L,K))
          UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY2(L,K))
          VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX2(L,K))
          U2(L,K)=0.5*(U(L,K)+U2(L,K))
          V2(L,K)=0.5*(V(L,K)+V2(L,K))
C         DW(L,K)=0.25*(W(L,K)-W2(L,K))
          W2(L,K)=0.5*(W(L,K)+W2(L,K))
         ENDDO
        ENDDO
C
      ENDIF
      ENDIF
C
C
      IF(ISWVSD.GE.1)THEN
C
        DO K=1,KC
         DO L=2,LA
          UHDY2(L,K)=UHDY2(L,K)+DYU(L)*UVPT(L,K)
          VHDX2(L,K)=VHDX2(L,K)+DXV(L)*VVPT(L,K)
          U2(L,K)=U2(L,K)+UVPT(L,K)/HMU(L)
          V2(L,K)=U2(L,K)+VVPT(L,K)/HMV(L)
          W2(L,K)=W2(L,K)+WVPT(L,K)
         ENDDO
        ENDDO
C
      ENDIF
C
C **  ADDITIONAL 3D CONTINUITY ADJUSTED ADDED BELOW
C
      IF(KC.GT.1)THEN
C
      DO L=2,LA
        TVAR3E(L)=0.
        TVAR3N(L)=0.
      ENDDO
C
      DO K=1,KC
       DO L=2,LA
        TVAR3E(L)=TVAR3E(L)+GVCSCLU(L)*UHDY2(L,K)*DZC(K)
        TVAR3N(L)=TVAR3N(L)+GVCSCLV(L)*VHDX2(L,K)*DZC(K)
       ENDDO
      ENDDO
C
      IF(ISTL.EQ.3)THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*( QSUME(L)
     &             -TVAR3E(L+1)+TVAR3E(L)
     &             -TVAR3N(LN) +TVAR3N(L) )
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP
     &                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
	  ERRHPP=ABS(HP(L)-HPPTMP)
C	  IF(ERRHPP.GT.0.00001.AND.SPB(L).GT.0.5) 
C     &    WRITE(8,6999)N,IL(L),JL(L),HPPTMP,HP(L)
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        HPI(L)=1./HP(L)
        ENDDO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*( QSUME(L)
     &             -TVAR3E(L+1)+TVAR3E(L)
     &             -TVAR3N(LN) +TVAR3N(L) )
        IF(ISGWIE.GE.1) HPPTMP=HPPTMP
     &                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
	  ERRHPP=ABS(HP(L)-HPPTMP)
C	  IF(ERRHPP.GT.0.00001.AND.SPB(L).GT.0.5) 
C     &    WRITE(8,6999)N,IL(L),JL(L),HPPTMP,HP(L)
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        HPI(L)=1./HP(L)
        ENDDO
      ENDIF
C
      IF(MDCHH.GE.1)THEN
        RLAMN=QCHERR
        RLAMO=1.-RLAMN
        DO NMD=1,MDCHH
        LHOST=LMDCHH(NMD)
        LCHNU=LMDCHU(NMD)
        LCHNV=LMDCHV(NMD)
        IF(MDCHTYP(NMD).EQ.1)THEN
          TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUN(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
          HPI(LHOST)=1./HP(LHOST)
          HPI(LCHNU)=1./HP(LCHNU)
        ENDIF            
        IF(MDCHTYP(NMD).EQ.2)THEN
          TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVN(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
          HPI(LHOST)=1./HP(LHOST)
          HPI(LCHNV)=1./HP(LCHNV)
        ENDIF            
        ENDDO
      ENDIF
C
      ENDIF
C
C
CGVCDIAG
c      IF(N.LE.NTSPTC)THEN
c      DO L=2,LA
c	  WRITE(8,891)N,IL(L),JL(L),UHDYE(L),(U(L,K),K=1,KC)
c      ENDDO
c      DO L=2,LA
c	  WRITE(8,892)N,IL(L),JL(L),HP(L),(W(L,K),K=1,KC)
c      ENDDO
c      ENDIF
c  891 FORMAT('FINAL QX,U GVC ',3I5,10E14.6)
c  892 FORMAT('FINAL HP,W GVC ',3I5,10E14.6)
CGVCDIAG
C
 6999 FORMAT('  DEP ADJ ERR ',3I5,2E14.6)
C
C**********************************************************************C
C
C **  ACCUMULTATE MAX COURANT NUMBERS
C
      DO K=1,KC
      DO L=2,LA
       CFLUUUT=DELT*ABS(DXIU(L)*U(L,K))
       CFLUUU(L,K)=MAX(CFLUUUT,CFLUUU(L,K))
       CFLVVVT=DELT*ABS(DYIV(L)*V(L,K))
       CFLVVV(L,K)=MAX(CFLVVVT,CFLVVV(L,K))
       CFLWWWT=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))
       CFLWWW(L,K)=MAX(CFLWWWT,CFLWWW(L,K))
       CFLCACT=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))
       CFLCAC(L,K)=MAX(CFLCACT,CFLCAC(L,K))
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C ** CALCULATE NONHYDROSTATIC PRESSURE
C
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) CALL CALPNHS
C
C**********************************************************************C
C
C **  WRITE TO DIAGNOSTIC FILE CFL.OUT WITH DIAGNOSTICS OF MAXIMUM
C **  TIME STEP 
C **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
C
C----------------------------------------------------------------------C
C
      IF(ISCFL.GE.1.AND.ISTL.EQ.3)THEN
C
        OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')
        IF(ISCFLM.GE.1.AND.N.EQ.1)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          CLOSE(2,STATUS='DELETE')
          DO L=1,LC
          ICFLMP(L)=0
          ENDDO
        ENDIF
C
        DTCFL=1.E+18
C
        K=1
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WBTMP=0.
        WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)
     &         +1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
C
        IF(KC.GT.1)THEN
        K=KC
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WTTMP=0.
        WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)
     &         +1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
        ENDIF
C
        IF(KC.GT.2)THEN
        DO K=2,KS
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))
        WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K  ))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)
     &         +1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        ENDIF
        ENDDO
        ENDDO
        ENDIF
C
        IVAL=MOD(N,ISCFL)
        IDTCFL=NINT(DTCFL)
        IF(ISCFL.EQ.1) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL
        IF(ISCFL.GE.2.AND.IVAL.EQ.0 )WRITE(1,1213)IDTCFL
        IF(ISCFLM.GE.1 )THEN
          LTMP=LIJ(ICFL,JCFL)
          ICFLMP(LTMP)=ICFLMP(LTMP)+1
        ENDIF
C
        IF(ISCFLM.GE.1.AND.N.EQ.NTS)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          TMPVALN=1./FLOAT(NTS)
          DO L=2,LA
           TMPVAL=TMPVALN*FLOAT(ICFLMP(L))
           WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL
          ENDDO
          CLOSE(2)
        ENDIF
C
        CLOSE(1)
      ENDIF
C
C----------------------------------------------------------------------C
C
      IF(ISCFL.GE.1.AND.IS1DCHAN.GE.1)THEN
C
        OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')
        IF(ISCFLM.GE.1.AND.N.EQ.1)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          CLOSE(2,STATUS='DELETE')
          DO L=1,LC
          ICFLMP(L)=0
          ENDDO
        ENDIF
C
        DTCFL=1.E+18
C
        K=1
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))
        UETMP=ABS(DXIU(L+1)*U2(L+1,K))
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)
     &         +1.0E-12
        DTMAX=0.5/DTMAXI
        IF(DTMAX.LT.DTCFL)THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
        ENDIF
        ENDDO
C
        IVAL=MOD(N,ISCFL)
        IDTCFL=NINT(DTCFL)
        IF(ISCFL.EQ.1) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL
        IF(ISCFL.GE.2.AND.IVAL.EQ.0 )WRITE(1,1213)IDTCFL
        IF(ISCFLM.GE.1 )THEN
          LTMP=LIJ(ICFL,JCFL)
          ICFLMP(LTMP)=ICFLMP(LTMP)+1
        ENDIF
C
        IF(ISCFLM.GE.1.AND.N.EQ.NTS)THEN
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')
          TMPVALN=1./FLOAT(NTS)
          DO L=2,LA
           TMPVAL=TMPVALN*FLOAT(ICFLMP(L))
           WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL
          ENDDO
          CLOSE(2)
        ENDIF
C
        CLOSE(1)
      ENDIF
C
 1212 FORMAT(' MAX TIME STEP =',F10.2,' SEC FOR N,I,J,K =',I8,3I5)
 1213 FORMAT(I4)
 1214 FORMAT(2I5,I12,F10.2)
C
C**********************************************************************C
C
      RETURN
      END
