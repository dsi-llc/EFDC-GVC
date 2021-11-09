C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE REDKC
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
C **  SUBROUTINE REDKC REDUCES THE NUMBER OF VERTICAL LAYERS BY 1/2 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      KCD2=KC/2
C
C**********************************************************************C
C
C **  LOOP OVER ALL CELLS
C     
C----------------------------------------------------------------------C
C
      DO L=2,LA
C
      DO K=2,KC,2
      KM=K-1
      KK=K/2
      BH(L,KK)=UHDY(L,K)+UHDY(KM,L)
      ENDDO
      DO K=1,KCD2
      UHDY(L,K)=BH(L,K)
      UHDY(L,K)=UHDY(L,K)
      U(L,K)=0.5*DZI*UHDY(L,K)/(HU(L)*DYU(L))
      U1(L,K)=U(L,K)
      ENDDO
C
      DO K=2,KC,2
      KM=K-1
      KK=K/2
      BH(L,KK)=VHDX(L,K)+VHDX(KM,L)
      ENDDO
      DO K=1,KCD2
      VHDX(L,K)=BH(L,K)
      VHDX(L,K)=VHDX(L,K)
      V(L,K)=0.5*VHDX(L,K)/(HV(L)*DXV(L))
      V1(L,K)=V(L,K)
      ENDDO
C
      DO K=2,KC,2
      KM=K-1
      KK=K/2
      BH(L,KK)=0.5*(UUU(L,K)+UUU(KM,L))
      ENDDO
      DO K=1,KCD2
      UUU(L,K)=BH(L,K)
      ENDDO
C
      DO K=2,KC,2
      KM=K-1
      KK=K/2
      BH(L,KK)=0.5*(VVV(L,K)+VVV(KM,L))
      ENDDO
      DO K=1,KCD2
      VVV(L,K)=BH(L,K)
      ENDDO
C
      DO K=2,KC,2
      KM=K-1
      KK=K/2
      BH(L,KK)=0.5*(SAL(L,K)+SAL(L,K-1))
      ENDDO
      DO K=1,KCD2
      SAL(L,K)=BH(L,K)
      SAL1(L,K)=BH(L,K)
      ENDDO
C
      DO K=2,KC-2,2
      KM=K-1
      KP=K+1
      KK=K/2
      BH(L,KK)=W(L,K)
      ENDDO
      DO K=1,KCD2-1
      W(L,K)=BH(L,K)
      W1(L,K)=W(L,K)
      ENDDO
      K=KCD2
      W(L,K)=0.
      W1(L,K)=W(L,K)
C
      ISFILAB=1 
      IF(ISFILAB.EQ.0)THEN
C
      DO K=2,KC-2,2
      KM=K-1
      KP=K+1
      KK=K/2
      BH(L,KK)=AB(L,K)
      ENDDO
      DO K=1,KCD2-1
      AB(L,K)=BH(L,K)
      ENDDO
      AB(KCD2,L)=ABO
C
      ELSE
C
      DO K=2,KC-2,2
      KM=K-1
      KP=K+1
      KK=K/2
      BH(L,KK)=4./(1./AB(L,K-1)+2./AB(L,K)+1./AB(KP,L))
      ENDDO
      DO K=1,KCD2-1
      AB(L,K)=BH(L,K)
      ENDDO
      AB(KCD2,L)=ABO
C
      ENDIF
C
      ENDDO
C
C**********************************************************************C
C
      KC=KCD2
      KS=KC-1
C
      IF(KS.EQ.0) KS=1
C
C     DZI=FLOAT(KC)
C     DZIS=DZI*DZI
C     DZISD4=DZIS/4.
      DZ=1./DZI
      DZ2=2.*DZ
      DZS=DZ*DZ
      DZDKC=DZ/FLOAT(KC)
C
      DZDDT=DZ/DT
      DZDDT2=0.5*DZ/DT
      DZSDDT=DZ*DZ/DT
C
C**********************************************************************C
C
      RETURN
      END
