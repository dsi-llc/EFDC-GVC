C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE FILTRAN
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
C **  SUBROUTINE FILTRAN SPATIALLY FILTERS THE MEAN MASS TRANSPORT 
C **  VELOCITY AND DIFFUSIVITY FIELDS AND COMPUTES SUBSCALE DIFFUSION
C **  COEFFICIENTS  AXZ AND AYZ REPLACED BY UUU AND VVV
C **  COEFFICIENTS  UHPF AND VHPF REPLACED BY DU AND DV
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  FILTER HORIZONTAL VELOCITY COMPONENTS IN VERTICAL DIRECTION
C **  AND CALCULATE HORIZONTAL DISPERSION COEFFICIENTS
C     
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LE=L+1
      LW=L-1
      LN=LNC(L)
      LS=LSC(L)
C
C----------------------------------------------------------------------C
C
      ULPF(L,1)=0.25*(3.*UHDY(L,1)+UHDY(L,2))
      VLPF(L,1)=0.25*(3.*VHDX(L,1)+VHDX(L,2))
      ULPF(L,KC)=0.25*(3.*UHDY(L,KC)+UHDY(L,KS))
      VLPF(L,KC)=0.25*(3.*VHDX(L,KC)+VHDX(L,KS))
      DO K=2,KS
      ULPF(L,K)=0.25*(UHDY(L,K-1)+2.*UHDY(L,K)+UHDY(L,K+1))
      VLPF(L,K)=0.25*(VHDX(L,K-1)+2.*VHDX(L,K)+VHDX(L,K+1))
      ENDDO
C     
      DO K=1,KC
      DU(L,K)=UHDY(L,K)-ULPF(L,K)
      DV(L,K)=VHDX(L,K)-VLPF(L,K)
      ENDDO
C
      ISFILAZ=0
      IF(ISFILAZ.EQ.0)GOTO 100
C
C----------------------------------------------------------------------C
C
C     WTB=0.5*DZIS*HPI(L)
C
      ERU=0.
      ERV=0.
      DO K=1,KC
      DU(L,K)=DU(L,K)*DZIC(K)/(DYU(L)*HU(L))
      DV(L,K)=DV(L,K)*DZIC(K)/(DXV(L)*HV(L))
      ERU=ERU+DU(L,K)
      ERV=ERV+DV(L,K)
      ENDDO
C
      ERU=ERU*DZ
      ERV=ERV*DZ
      DO K=1,KC
      DU(L,K)=DU(L,K)-ERU
      DV(L,K)=DV(L,K)-ERV
      ENDDO
C
      FZU(0)=0.
      FZV(0)=0.
C
      DO K=1,KS
      FZU(K)=FZU(K-1)+DU(L,K)
      FZV(K)=FZV(K-1)+DV(L,K)
      ENDDO
C
      CUB(1)=0.
      CLB(1)=0.
C
      DO K=1,KS
      CUB(K+1)=CUB(K)+FZU(K)/(AB(L,K)+AB(L-1,K))
      CLB(K+1)=CLB(K)+FZV(K)/(AB(L,K)+AB(LS,K))
      ENDDO
C
      ERU=0.
      ERV=0.
      DO K=1,KC
      ERU=ERU+CUB(K)
      ERV=ERV+CLB(K)
      ENDDO
C
      ERU=ERU*DZ
      ERV=ERV*DZ
      DO K=1,KC
      CUB(K)=CUB(K)-ERU
      CLB(K)=CLB(K)-ERV
      ENDDO
C
      AXZEX=0.
      AYZEX=0.
      DO K=1,KC
      DU(L,K)=-2.*DZ*DZ*SUB(L)*HU(L)*DU(L,K)*CUB(K)
      DV(L,K)=-2.*DZ*DZ*SVB(L)*HV(L)*DV(L,K)*CLB(K)
      AXZEX=AXZEX+DZ*DU(L,K)
      AYZEX=AYZEX+DZ*DV(L,K)
      ENDDO
C
      UUU(L,1)=0.25*(3.*DU(L,1)+DU(L,1))
      VVV(L,1)=0.25*(3.*DV(L,1)+DV(L,2))
      UUU(L,KC)=0.25*(3.*DU(L,KC)+DU(L,KS))
      VVV(L,KC)=0.25*(3.*DV(L,KC)+DV(L,KS))
      DO K=2,KS
      UUU(L,K)=0.25*(DU(L,K-1)+2.*DU(L,K)+DU(L,K+1))
      VVV(L,K)=0.25*(DV(L,K-1)+2.*DV(L,K)+DV(L,K+1))
      ENDDO
C
      ISAZCON=0
      IF(ISAZCON.EQ.1)THEN
      DO K=1,KC
      UUU(L,K)=AXZEX
      VVV(L,K)=AYZEX
      ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
  100 CONTINUE
C
C----------------------------------------------------------------------C
C
      ISFILUH=0
      IF(ISFILUH.EQ.1)THEN
      DO K=1,KC
      UHDY(L,K)=ULPF(L,K)
      U(L,K)=UHDY(L,K)/(HU(L)*DYU(L))
      U1(L,K)=U(L,K)
      VHDX(L,K)=VLPF(L,K)
      V(L,K)=VHDX(L,K)/(HV(L)*DXV(L))
      V1(L,K)=V(L,K)
      ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
      ENDDO
C
C**********************************************************************C
C

      DO LL=1,NCBE
      L=LCBE(LL)
      LE=L+1
      DO K=1,KC
      UHDY(LE,K)=UHDY(L,K)
      UHDY(LE,K)=UHDY(L,K)
      U(LE,K)=UHDY(L,K)/(HU(L)*DYU(L))
      U1(LE,K)=U(LE,K)
      ENDDO
      ENDDO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LN=LNC(L)
      DO K=1,KC
      VHDX(LN,K)=VHDX(L,K)
      VHDX(LN,K)=VHDX(L,K)
      V(LN,K)=DZI*VHDX(L,K)/(HV(L)*DXV(L))
      V1(LN,K)=V(LN,K)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  PERFORM DIAGNOSTIC ON HORIZONTAL DISPERSION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      AXZMAX=0.
      AXZMIN=100000.
      AYZMAX=0.
      AYZMIN=100000.
C
      DO L=1,LC
      DO K=1,KC
      IF(UUU(L,K).GT.AXZMAX)THEN
       AXZMAX=UUU(L,K)
       IAXZMA=IL(L)
       JAXZMA=JL(L)
       KAXZMA=K
      ENDIF
      IF(UUU(L,K).LT.AXZMIN)THEN
       AXZMIN=UUU(L,K)
       IAXZMI=IL(L)
       JAXZMI=JL(L)
       KAXZMI=K
      ENDIF
      IF(VVV(L,K).GT.AYZMAX)THEN
       AYZMAX=VVV(L,K)
       IAYZMA=IL(L)
       JAYZMA=JL(L)
       KAYZMA=K
      ENDIF
      IF(VVV(L,K).LT.AYZMIN)THEN
       AYZMIN=VVV(L,K)
       IAYZMI=IL(L)
       JAYZMI=JL(L)
       KAYZMI=K
      ENDIF
      ENDDO
      ENDDO
C
      WRITE(6,686)IAXZMA,JAXZMA,KAXZMA,AXZMAX
      WRITE(6,687)IAXZMI,JAXZMI,KAXZMI,AXZMIN
      WRITE(6,688)IAYZMA,JAYZMA,KAYZMA,AYZMAX
      WRITE(6,689)IAYZMI,JAYZMI,KAYZMI,AYZMIN
  686 FORMAT('  IJKAXZMAX=',3I10,3X,E12.4//)
  687 FORMAT('  IJKAXZMIN=',3I10,3X,E12.4//)
  688 FORMAT('  IJKAYZMAX=',3I10,3X,E12.4//)
  689 FORMAT('  IJKAYZMIN=',3I10,3X,E12.4//)
C
C**********************************************************************C
C
      RETURN
      END
