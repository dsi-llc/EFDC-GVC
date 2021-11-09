C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T2
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
C 05/01/2002        john hamrick       05/01/2002       john hamrick
C  subroutine added for 2 time-level balances including sed,snd,tox
C----------------------------------------------------------------------C
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT
      ELSE
        DELT=DTDYN
      END IF
C
C**********************************************************************C
C
C **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
C
      VOLOUT=VOLOUT-DELT*VHDX2(LN,K)*DZC(K)
      WVOLOUT=WVOLOUT-DELT*VHDX2(LN,K)*DZC(K)
      SALOUT=SALOUT-DELT*MIN(VHDX2(LN,K),0.)*SAL(LN,K)*DZC(K)
     &             -DELT*MAX(VHDX2(LN,K),0.)*SAL(L,K)*DZC(K)
      DYEOUT=DYEOUT-DELT*MIN(VHDX2(LN,K),0.)*DYE(LN,K)*DZC(K)
     &             -DELT*MAX(VHDX2(LN,K),0.)*DYE(L,K)*DZC(K)
      DO NS=1,NSED
        SEDOUT2T(NS)=SEDOUT2T(NS)
     &              -DELT*MIN(VHDX2(LN,K),0.)*SED(LN,K,NS)*DZC(K)
     &              -DELT*MAX(VHDX2(LN,K),0.)*SED(L,K,NS)*DZC(K)
	ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS)=SNDOUT2T(NS)
     &              -DELT*MIN(VHDX2(LN,K),0.)*SND(LN,K,NS)*DZC(K)
     &              -DELT*MAX(VHDX2(LN,K),0.)*SND(L,K,NS)*DZC(K)
	ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT)=TOXOUT2T(NT)
     &              -DELT*MIN(VHDX2(LN,K),0.)*TOX(LN,K,NT)*DZC(K)
     &              -DELT*MAX(VHDX2(LN,K),0.)*TOX(L,K,NT)*DZC(K)
	ENDDO
      PPEOUT=PPEOUT-DELT*VHDX2(LN,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LN))
     &       +0.125*(HP(L)+H1P(L)+HP(LN)+H1P(LN))*(Z(K)+Z(K-1)) )
      BBEOUT=BBEOUT-DELT*MIN(VHDX2(LN,K),0.)*DZC(K)*GP*( BELV(LN)
     &                  +0.5*HP(LN)*(Z(K)+Z(K-1)) )*B(LN,K)
     &             -DELT*MAX(VHDX2(LN,K),0.)*DZC(K)*GP*( BELV(L)
     &                  +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
C
      ENDDO
      ENDDO
C
C      DO LL=1,NCBS
C      L=LCBS(LL)
C      LN=LNC(L)
C      VOLOUT=VOLOUT-VHDX2E(LN)
C      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)
C
      VOLOUT=VOLOUT-DELT*UHDY2(L+1,K)*DZC(K)
      WVOLOUT=WVOLOUT-DELT*UHDY2(L+1,K)*DZC(K)
      SALOUT=SALOUT-DELT*MIN(UHDY2(L+1,K),0.)*SAL1(L+1,K)*DZC(K)
     &             -DELT*MAX(UHDY2(L+1,K),0.)*SAL1(L,K)*DZC(K)
      DYEOUT=DYEOUT-DELT*MIN(UHDY2(L+1,K),0.)*DYE1(L+1,K)*DZC(K)
     &             -DELT*MAX(UHDY2(L+1,K),0.)*DYE1(L,K)*DZC(K)
      DO NS=1,NSED
        SEDOUT2T(NS)=SEDOUT2T(NS)
     &              -DELT*MIN(UHDY2(L+1,K),0.)*SED(L+1,K,NS)*DZC(K)
     &              -DELT*MAX(UHDY2(L+1,K),0.)*SED(L,K,NS)*DZC(K)
	ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS)=SNDOUT2T(NS)
     &              -DELT*MIN(UHDY2(L+1,K),0.)*SND(L+1,K,NS)*DZC(K)
     &              -DELT*MAX(UHDY2(L+1,K),0.)*SND(L,K,NS)*DZC(K)
	ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT)=TOXOUT2T(NT)
     &              -DELT*MIN(UHDY2(L+1,K),0.)*TOX(L+1,K,NT)*DZC(K)
     &              -DELT*MAX(UHDY2(L+1,K),0.)*TOX(L,K,NT)*DZC(K)
	ENDDO
      PPEOUT=PPEOUT-DELT*UHDY2(L+1,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(L+1))
     &       +0.125*(HP(L)+H2P(L)+HP(L+1)+H2P(L+1))*(Z(K)+Z(K-1)) )
      BBEOUT=BBEOUT-DELT*MIN(UHDY2(L+1,K),0.)*DZC(K)*GP*( BELV(L+1)
     &                  +0.5*HP(L+1)*(Z(K)+Z(K-1)) )*B1(L+1,K)
     &             -DELT*MAX(UHDY2(L+1,K),0.)*DZC(K)*GP*( BELV(L)
     &                  +0.5*HP(L)*(Z(K)+Z(K-1)) )*B1(L,K)
C
      ENDDO
      ENDDO
C
C      DO LL=1,NCBW
C      L=LCBW(LL)
C      VOLOUT=VOLOUT-UHDY2E(L+1)
C      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
C
      VOLOUT=VOLOUT+DELT*UHDY2(L,K)*DZC(K)
      WVOLOUT=WVOLOUT+DELT*UHDY2(L,K)*DZC(K)
      SALOUT=SALOUT+DELT*MIN(UHDY2(L,K),0.)*SAL(L,K)*DZC(K)
     &             +DELT*MAX(UHDY2(L,K),0.)*SAL(L-1,K)*DZC(K)
      DYEOUT=DYEOUT+DELT*MIN(UHDY2(L,K),0.)*DYE(L,K)*DZC(K)
     &             +DELT*MAX(UHDY2(L,K),0.)*DYE(L-1,K)*DZC(K)
      DO NS=1,NSED
        SEDOUT2T(NS)=SEDOUT2T(NS)
     &              +DELT*MIN(UHDY2(L,K),0.)*SED(L,K,NS)*DZC(K)
     &              +DELT*MAX(UHDY2(L,K),0.)*SED(L-1,K,NS)*DZC(K)
	ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS)=SNDOUT2T(NS)
     &              +DELT*MIN(UHDY2(L,K),0.)*SND(L,K,NS)*DZC(K)
     &              +DELT*MAX(UHDY2(L,K),0.)*SND(L-1,K,NS)*DZC(K)
	ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT)=TOXOUT2T(NT)
     &              +DELT*MIN(UHDY2(L,K),0.)*TOX(L,K,NT)*DZC(K)
     &              +DELT*MAX(UHDY2(L,K),0.)*TOX(L-1,K,NT)*DZC(K)
	ENDDO
      PPEOUT=PPEOUT+DELT*UHDY2(L,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(L-1))
     &       +0.125*(HP(L)+H1P(L)+HP(L-1)+H1P(L-1))*(Z(K)+Z(K-1)) )
      BBEOUT=BBEOUT+DELT*MIN(UHDY2(L,K),0.)*DZC(K)*GP*(BELV(L)
     &                  +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
     &             +DELT*MAX(UHDY2(L,K),0.)*DZC(K)*GP*(BELV(L-1)
     &                  +0.5*HP(L-1)*(Z(K)+Z(K-1)) )*B(L-1,K)
C
      ENDDO
      ENDDO
C
C
C      DO LL=1,NCBE
C      L=LCBE(LL)
C      VOLOUT=VOLOUT+UHDY2E(L)
C      ENDDO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
C
      VOLOUT=VOLOUT+DELT*VHDX2(L,K)*DZC(K)
      WVOLOUT=WVOLOUT+DELT*VHDX2(L,K)*DZC(K)
      SALOUT=SALOUT+DELT*MIN(VHDX2(L,K),0.)*SAL(L,K)*DZC(K)
     &             +DELT*MAX(VHDX2(L,K),0.)*SAL(LS,K)*DZC(K)
      DYEOUT=DYEOUT+DELT*MIN(VHDX2(L,K),0.)*DYE(L,K)*DZC(K)
     &             +DELT*MAX(VHDX2(L,K),0.)*DYE(LS,K)*DZC(K)
      DO NS=1,NSED
        SEDOUT2T(NS)=SEDOUT2T(NS)
     &              +DELT*MIN(VHDX2(L,K),0.)*SED(L,K,NS)*DZC(K)
     &              +DELT*MAX(VHDX2(L,K),0.)*SED(LS,K,NS)*DZC(K)
	ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS)=SNDOUT2T(NS)
     &              +DELT*MIN(VHDX2(L,K),0.)*SND(L,K,NS)*DZC(K)
     &              +DELT*MAX(VHDX2(L,K),0.)*SND(LS,K,NS)*DZC(K)
	ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT)=TOXOUT2T(NT)
     &              +DELT*MIN(VHDX2(L,K),0.)*TOX(L,K,NT)*DZC(K)
     &              +DELT*MAX(VHDX2(L,K),0.)*TOX(LS,K,NT)*DZC(K)
	ENDDO
      PPEOUT=PPEOUT+DELT*VHDX2(L,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LS))
     &       +0.125*(HP(L)+H1P(L)+HP(LS)+H1P(LS))*(Z(K)+Z(K-1)) )
      BBEOUT=BBEOUT+DELT*MIN(VHDX2(L,K),0.)*DZC(K)*GP*( BELV(L)
     &                  +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
     &             +DELT*MAX(VHDX2(L,K),0.)*DZC(K)*GP*( BELV(LS)
     &                  +0.5*HP(LS)*(Z(K)+Z(K-1)) )*B(LS,K)
C
      ENDDO
      ENDDO
C
C      DO LL=1,NCBN
C      L=LCBN(LL)
C      LS=LSC(L)
C      VOLOUT=VOLOUT+VHDX2E(L)
C      ENDDO
C
C**********************************************************************C
C
      RETURN
      END
