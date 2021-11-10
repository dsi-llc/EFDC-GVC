C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTIN1
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
C  added code to properly initial restart input for drying and wetting

C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE RESTIN1 READS A RESTART FILE
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION TDUMMY(KCM)
C
C**********************************************************************C
C
      OPEN(1,FILE='RESTART.INP',STATUS='UNKNOWN')
C
C**********************************************************************C
C
      ISBELVC=0
C
      READ(1,908,ERR=1000)NREST
C
      DO L=2,LA
      IF(ISRESTI.EQ.1)THEN
	  IF(IMORPH.EQ.0)THEN
          READ(1,*,ERR=1001)HP(L),H1P(L),HWQ(L),H2WQ(L)
        ELSE
          READ(1,*,ERR=1001)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
        ENDIF
      ELSE
        READ(1,*,ERR=1002)HP(L),H1P(L),HWQ(L),H2WQ(L),BELTMP
        IF(BELTMP.NE.BELV(L))THEN
          ISBELVC=1
          WRITE(6,600)IL(L),JL(L),BELTMP,BELV(L)
          HP(L)=HP(L)+BELTMP-BELV(L)
          H1P(L)=H1P(L)+BELTMP-BELV(L)
          HWQ(L)=HWQ(L)+BELTMP-BELV(L)
          H2WQ(L)=H2WQ(L)+BELTMP-BELV(L)
        ENDIF
      ENDIF
C
      IF(HP(L).LT.0.0.OR.H1P(L).LT.0.0)THEN
C      IF(HP(L).LT.0.0)THEN
        WRITE(6,9696)L,IL(L),JL(L),HP(L),H1P(L)
        STOP
      ENDIF
C
      READ(1,*,ERR=1003)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      READ(1,*,ERR=1004)(U(L,K),K=1,KC)
      READ(1,*,ERR=1005)(U1(L,K),K=1,KC)
      READ(1,*,ERR=1006)(V(L,K),K=1,KC)
      READ(1,*,ERR=1007)(V1(L,K),K=1,KC)
c
C      READ(1,*,ERR=1008)(QQ(L,K),K=0,KS)
C      READ(1,*,ERR=1008) QQ(L,KC)
C      READ(1,*,ERR=1009)(QQ1(L,K),K=0,KS)
C      READ(1,*,ERR=1009) QQ1(L,KC)
C      READ(1,*,ERR=1010)(QQL(L,K),K=0,KS)
C      READ(1,*,ERR=1010) QQL(L,KC)
C      READ(1,*,ERR=1011)(QQL1(L,K),K=0,KS)
C      READ(1,*,ERR=1011) QQL1(L,KC)
C      READ(1,*,ERR=1012)(DML(L,K),K=0,KS)
C      READ(1,*,ERR=1012) DML(L,KC)

      READ(1,*,ERR=1008)(QQ(L,K),K=0,KC)
      READ(1,*,ERR=1009)(QQ1(L,K),K=0,KC)
      READ(1,*,ERR=1010)(QQL(L,K),K=0,KC)
      READ(1,*,ERR=1011)(QQL1(L,K),K=0,KC)
      READ(1,*,ERR=1012)(DML(L,K),K=0,KC)
c
      IF(ISCI(1).EQ.1)THEN
       READ(1,*,ERR=1013)(SAL(L,K),K=1,KC)
       READ(1,*,ERR=1014)(SAL1(L,K),K=1,KC)
      ENDIF
      IF(ISCI(2).GE.1)THEN  ! *** HANDLE QEA TEMPB BYPASS
       READ(1,*,ERR=1015)(TEM(L,K),K=1,KC)
       READ(1,*,ERR=1016)(TEM1(L,K),K=1,KC)
      ENDIF
      IF(ISCI(3).EQ.1)THEN
       READ(1,*,ERR=1017)(DYE(L,K),K=1,KC)
       READ(1,*,ERR=1018)(DYE1(L,K),K=1,KC)
      ENDIF
      IF(ISCI(4).EQ.1)THEN
C       READ(1,*,ERR=1021)(SFL(L,K),K=1,KC)
C       READ(1,*,ERR=1022)(SFL2(L,K),K=1,KC)
       READ(1,*,ERR=1021)SFLSBOT(L),(SFL(L,K),K=1,KC)
       READ(1,*,ERR=1022)SFLSBOT(L),(SFL2(L,K),K=1,KC)
      ENDIF
      IF(ISCI(5).EQ.1)THEN
       DO NT=1,NTOX
       READ(1,*,ERR=1019)(TOXB(L,K,NT),K=1,KB)
       READ(1,*,ERR=1019)(TOX(L,K,NT),K=1,KC)
       READ(1,*,ERR=1020)(TOXB1(L,K,NT),K=1,KB)
       READ(1,*,ERR=1020)(TOX1(L,K,NT),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCI(6).EQ.1)THEN
       DO NS=1,NSED
       READ(1,*,ERR=1019)(SEDB(L,K,NS),K=1,KB)
       READ(1,*,ERR=1019)(SED(L,K,NS),K=1,KC)
       READ(1,*,ERR=1020)(SEDB1(L,K,NS),K=1,KB)
       READ(1,*,ERR=1020)(SED1(L,K,NS),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCI(7).EQ.1)THEN
       DO NS=1,NSND
       READ(1,*,ERR=1019)(SNDB(L,K,NS),K=1,KB)
       READ(1,*,ERR=1019)(SND(L,K,NS),K=1,KC)
       READ(1,*,ERR=1020)(SNDB1(L,K,NS),K=1,KB)
       READ(1,*,ERR=1020)(SND1(L,K,NS),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCI(6).EQ.1.OR.ISCI(7).EQ.1)THEN
       READ(1,*,ERR=1019)(HBED(L,K),K=1,KB)
       READ(1,*,ERR=1019)(HBED1(L,K),K=1,KB)
       READ(1,*,ERR=1019)(VDRBED(L,K),K=1,KB)
       READ(1,*,ERR=1019)(VDRBED1(L,K),K=1,KB)
      ENDIF
      ENDDO  ! *** END OF DO L=2,LA
C
      DO M=1,4
      IF(ISCI(M).EQ.1)THEN
C
       DO LL=1,NCBS
       READ(1,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
       READ(1,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       READ(1,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
       READ(1,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       READ(1,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
       READ(1,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       READ(1,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
       READ(1,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDIF
      ENDDO
C
      IF(ISCI(5).EQ.1)THEN
      DO NT=1,NTOX
      M=MSVTOX(NT)
C
       DO LL=1,NCBS
       READ(1,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
       READ(1,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       READ(1,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
       READ(1,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       READ(1,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
       READ(1,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       READ(1,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
       READ(1,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
      IF(ISCI(6).EQ.1)THEN
      DO NT=1,NSED
      M=MSVSED(NT)
C
       DO LL=1,NCBS
       READ(1,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
       READ(1,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       READ(1,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
       READ(1,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       READ(1,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
       READ(1,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       READ(1,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
       READ(1,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
      IF(ISCI(7).EQ.1)THEN
      DO NT=1,NSND
      M=MSVSND(NT)
C
       DO LL=1,NCBS
       READ(1,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
       READ(1,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       READ(1,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
       READ(1,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       READ(1,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
       READ(1,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       READ(1,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
       READ(1,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
c     IF(ISDRY.GT.0)THEN
        DO L=2,LA
        READ(1,*,ERR=1031)QSUME(L),(QSUM(L,K),K=1,KC)
C       WRITE(6,6666)L,IL(L),JL(L),QSUME(L)
        ENDDO
c      ENDIF
C
      IF(MDCHH.GE.1)THEN
c      IF(MDCHH.EQ.2.AND.MDCHHD.EQ.1)THEN
        DO NMD=1,MDCHH
        READ(1,*,ERR=1032)ITMP1,JTMP1,ITMP2,JTMP2,
     &                    ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
C       WRITE(6,6667)NMD,ITMP1,JTMP1,ITMP2,JTMP2,
C    &                   ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
        ENDDO
       ELSE
        DO NMD=1,MDCHH
         QCHANU(NMD)=0.
         QCHANV(NMD)=0.
         QCHANUN(NMD)=0.
         QCHANVN(NMD)=0.
        ENDDO
c      ENDIF
      ENDIF
C
      IF(ISGWIE.GE.1)THEN
        DO L=2,LA
        READ(1,*,ERR=1033)AGWELV(L),AGWELV1(L)
        ENDDO
      ENDIF
C
      CLOSE(1)
C
 6666 FORMAT(3I10,F12.6)
 6667 FORMAT(7I5,2X,E12.4,2X,E12.4)
C
      DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      SED(1,K,1)=0.
      SND(1,K,1)=0.
      TOX(1,K,1)=0.
      SFL(1,K)=0.
      CWQ(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SED1(1,K,1)=0.
      SND1(1,K,1)=0.
      TOX1(1,K,1)=0.
      SFL2(1,K)=0.
      CWQ2(1,K)=0.
      VHDX1(1,K)=0.
      UHDY1(1,K)=0.
      VHDXWQ(1,K)=0.
      UHDYWQ(1,K)=0.
      SAL(LC,K)=0.
      TEM(LC,K)=0.
      DYE(LC,K)=0.
      SED(LC,K,1)=0.
      SND(LC,K,1)=0.
      TOX(LC,K,1)=0.
      SFL(LC,K)=0.
      CWQ(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SED1(LC,K,1)=0.
      SND1(LC,K,1)=0.
      TOX1(LC,K,1)=0.
      SFL2(LC,K)=0.
      CWQ2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXWQ(LC,K)=0.
      UHDYWQ(LC,K)=0.
      ENDDO
C
      DO L=2,LA
      UHDYE(L)=SUB(L)*UHDYE(L)
      UHDY1E(L)=SUB(L)*UHDY1E(L)
      VHDXE(L)=SVB(L)*VHDXE(L)
      VHDX1E(L)=SVB(L)*VHDX1E(L)
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
      U(L,K)=SUB(L)*U(L,K)
      U1(L,K)=SUB(L)*U1(L,K)
      V(L,K)=SVB(L)*V(L,K)
      V1(L,K)=SVB(L)*V1(L,K)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
      DO L=2,LA
      LS=LSC(L)
C      H1U(L)=0.5*(H1P(L)+H1P(L-1))
C      H1V(L)=0.5*(H1P(L)+H1P(LS))
       H1U(L)=0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(L-1)*DYP(L-1)*H1P(L-1))
     &           /(DXU(L)*DYU(L))
       H1V(L)=0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(LS )*DYP(L-1)*H1P(LS ))
     &           /(DXV(L)*DYV(L))
       P1(L)=G*(H1P(L)+BELV(L))
C      HU(L)=0.5*(HP(L)+HP(L-1))
C      HV(L)=0.5*(HP(L)+HP(LS))
       HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(L-1)*DYP(L-1)*HP(L-1))
     &           /(DXU(L)*DYU(L))
       HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LS )*DYP(L-1)*HP(LS ))
     &           /(DXV(L)*DYV(L))
       P(L)=G*(HP(L)+BELV(L))
       HPI(L)=1./HP(L)
       HUI(L)=1./HU(L)
       HVI(L)=1./HV(L)
       H1UI(L)=1./H1U(L)
       H1VI(L)=1./H1V(L)
      ENDDO
C
      H1U(1)=H1U(2)
      H1V(1)=H1V(2)
      P1(1)=P1(2)
      HU(1)=HU(2)
      HV(1)=HV(2)
      P(1)=P(2)
      HPI(1)=1./HP(2)
      HUI(1)=1./HU(2)
      HVI(1)=1./HV(2)
      H1UI(1)=1./H1U(2)
      H1VI(1)=1./H1V(2)
C
      H1U(LC)=H1U(LA)
      H1V(LC)=H1V(LA)
      P1(LC)=P1(LA)
      HU(LC)=HU(LA)
      HV(LC)=HV(LA)
      P(LC)=P(LA)
      HPI(LC)=1./HP(LA)
      HUI(LC)=1./HU(LA)
      HVI(LC)=1./HV(LA)
      H1UI(LC)=1./H1U(LA)
      H1VI(LC)=1./H1V(LA)
C
      DO K=1,KC
      DO L=2,LA
      UHDY1(L,K)=DYU(L)*H1U(L)*U1(L,K)
      VHDX1(L,K)=DXV(L)*H1V(L)*V1(L,K)
      UHDY(L,K)=DYU(L)*HU(L)*U(L,K)
      VHDX(L,K)=DXV(L)*HV(L)*V(L,K)
      SAL(L,K)=MAX(SAL(L,K),0.)
      SAL1(L,K)=MAX(SAL1(L,K),0.)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CORRECT FOR CHANGED BOTTOM ELEV
C
      IF(ISRESTI.EQ.-1.AND.ISBELVC.EQ.1)THEN
C
        DO L=2,LA
        UHE(L)=0.
        VHE(L)=0.
        ENDDO
C
        DO K=1,KC
        DO L=2,LA
        UHE(L)=UHE(L)+UHDY1(L,K)
        VHE(L)=VHE(L)+VHDX1(L,K)
        ENDDO
        ENDDO
C
        DO L=2,LA
        IF(UHE(L).NE.0.)THEN
          TMPVAL=UHDY1E(L)/UHE(L)
          DO K=1,KC
          U1(L,K)=TMPVAL*U1(L,K)
          ENDDO
        ENDIF
        IF(VHE(L).NE.0.)THEN
          TMPVAL=VHDX1E(L)/VHE(L)
          DO K=1,KC
          V1(L,K)=TMPVAL*V1(L,K)
          ENDDO
        ENDIF
        ENDDO
C      
        DO L=2,LA
        UHE(L)=0.
        VHE(L)=0.
        ENDDO
C
        DO K=1,KC
        DO L=2,LA
        UHE(L)=UHE(L)+UHDY(L,K)
        VHE(L)=VHE(L)+VHDX(L,K)
        ENDDO
        ENDDO
C
        DO L=2,LA
        IF(UHE(L).NE.0.)THEN
          TMPVAL=UHDYE(L)/UHE(L)
          DO K=1,KC
          U(L,K)=TMPVAL*U(L,K)
          ENDDO
        ENDIF
        IF(VHE(L).NE.0.)THEN
          TMPVAL=VHDXE(L)/VHE(L)
          DO K=1,KC
          V(L,K)=TMPVAL*V(L,K)
          ENDDO
        ENDIF
        ENDDO
C      
        DO K=1,KC
        DO L=2,LA
        UHDY1(L,K)=DYU(L)*H1U(L)*U1(L,K)
        VHDX1(L,K)=DXV(L)*H1V(L)*V1(L,K)
        UHDY(L,K)=DYU(L)*HU(L)*U(L,K)
        VHDX(L,K)=DXV(L)*HV(L)*V(L,K)
        ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      N=0
      IF(ISDRY.EQ.0)THEN
        CALL CALTSXY
        CALL CALQVS (2)
      ENDIF
C
      DO K=1,KS
      RDZC=DZC(K)
      DO L=2,LA
      LN=LNC(L)
      W(L,K)=SWB(L)*(W(L,K-1)
     &       -RDZC*(UHDY(L+1,K)-UHDY(L,K)-UHDYE(L+1)+UHDYE(L)
     &       +VHDX(LN,K)-VHDX(L,K)-VHDXE(LN)+VHDXE(L))*DXYIP(L))
     &        +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      W1(L,K)=SWB(L)*(W1(L,K-1)
     &       -RDZC*(UHDY1(L+1,K)-UHDY1(L,K)-UHDY1E(L+1)+UHDY1E(L)
     &       +VHDX1(LN,K)-VHDX1(L,K)-VHDX1E(LN)+VHDX1E(L))*DXYIP(L))
     &        +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      ENDDO       
      ENDDO       
C
      OPEN(1,FILE='NEWSALT.INP',STATUS='UNKNOWN')
      IONE=1
      WRITE(1,101)IONE
      DO L=2,LA
      WRITE(1,102)L,IL(L),JL(L),(SAL(L,K),K=1,KC)
      ENDDO
      CLOSE(1)
C
C ** READ BED TEMPERATURE  HTBED1 HTBED2
C
      IF(ISCI(2).EQ.1)THEN                       ! *** Bypass for AQ Compatibililty ( ISCI(2) = 2 Bypasses )
       OPEN(1,FILE='TEMPB.RST',STATUS='UNKNOWN')
       DO L=2,LA
C        READ(1,*)LDUM,IDUM,JDUM,(TDUMMY(K),K=1,KC),TEMB(L) 
        READ(1,*)LDUM,IDUM,JDUM,(TEMB(L,K),K=1,KBHM)
       ENDDO
       CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  WRITE NEW DXDY FILES TO INCORPORATE RESTART DEPTH
C
      IDUMMY=0
      OPEN(1,FILE='RSTDXDY.INP',STATUS='UNKNOWN')
      DO J=1,JC
      DO I=1,IC
      L=LIJ(I,J)
      IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LT.9)THEN
        WRITE(1,339)IL(L),JL(L),DXP(L),DYP(L),HP(L),BELV(L),
     &              ZBR(L),IDUMMY
      ENDIF
      ENDDO
      ENDDO
      CLOSE(1)
C
  339 FORMAT(2I5,5F14.5,I5)
C
C**********************************************************************C
C
C **  SET DRYING AND WETTING FLAGS
C
C
      IF(ISDRY.EQ.0)THEN
      DO L=2,LA
        IMASKDRY(L)=0
        LMASKDRY(L)=.TRUE.
	ENDDO
	ENDIF
C
      IF(ISDRY.GT.0.AND.ISDRY.LT.97)THEN
C
      DO L=2,LA
        ISCDRY(L)=0
        LS=LSC(L)
        LN=LNC(L)
        IF(HP(L).LE.HDRY)THEN
          ISCDRY(L)=1
          SUB(L)=0.
          SVB(L)=0.
          SUB(L+1)=0.
          SVB(LN)=0.
          SBX(L)=0.
          SBY(L)=0.
          SBX(L+1)=0.
          SBY(LN)=0.
        ENDIF
      ENDDO
C
      ENDIF
C
      IF(ISDRY.EQ.98)THEN
C
        HDRY2=2.*HDRY
        DO L=2,LA
        LS=LSC(L)
        LN=LNC(L)
        IF(HP(L).LE.HDRY)THEN
CZZ	    SUBW=SUB(L)
CZZ	    SUBE=SUB(L+1)
CZZ	    SVBS=SVB(L)
CZZ	    SVBN=SVB(LN)
	    DHPDT=(HP(L)-H1P(L))/DT
	    IF(DHPDT.GT.0.0)THEN
CZZ
	      SUBW=SUB(L)
	      SUBE=SUB(L+1)
	      SVBS=SVB(L)
	      SVBN=SVB(LN)
            SUB(L)=0.0  
            SUB(L+1)=0.0  
            SVB(L)=0.0  
            SVB(LN)=0.0  
            SBX(L)=0.0  
            SBX(L+1)=0.0  
            SBY(L)=0.0  
            SBY(LN)=0.0  
CZZ
	      IF(SUBO(L).GT.0.5)THEN
	        IF(UHDYE(L).GT.0.0.AND.HP(L-1).GT.HDRY2)THEN
	          SUB(L)=1.
                SBX(L)=1.
              ENDIF
            ENDIF
	      IF(SUBO(L+1).GT.0.5)THEN
	        IF(UHDYE(L+1).LT.0.0.AND.HP(L+1).GT.HDRY2)THEN
	          SUB(L+1)=1.
                SBX(L+1)=1.
              ENDIF
            ENDIF
	      IF(SVBO(L).GT.0.5)THEN
	        IF(VHDXE(L).GT.0.0.AND.HP(LS).GT.HDRY2)THEN
	          SVB(L)=1.
                SBY(L)=1.
              ENDIF
            ENDIF
	      IF(SVBO(LN).GT.0.5)THEN
	        IF(VHDXE(LN).LT.0.0.AND.HP(LN).GT.HDRY2)THEN
	          SVB(LN)=1.
                SBY(LN)=1.
              ENDIF
            ENDIF
	      RDRY=SUB(L)+SUB(L+1)+SVB(L)+SVB(LN)
	      IF(RDRY.LT.0.5)THEN
              ISCDRY(L)=1
            ELSE
              ISCDRY(L)=0
            ENDIF
	      IDFLAG=0
	      TMPVAL=ABS(SUB(L)-SUBW)
	      IF(TMPVAL.GT.0.5)ICORDRY=1
	      TMPVAL=ABS(SUB(L+1)-SUBE)
	      IF(TMPVAL.GT.0.5)ICORDRY=1
	      TMPVAL=ABS(SVB(L)-SVBS)
	      IF(TMPVAL.GT.0.5)ICORDRY=1
	      TMPVAL=ABS(SVB(LN)-SVBN)
	      IF(TMPVAL.GT.0.5)ICORDRY=1
CZZ            IF(IDFLAG.GT.0) ICORDRY=1
	    ELSE
	      SUB(L)=0.0
	      SUB(L+1)=0.0
	      SVB(L)=0.0
	      SVB(LN)=0.0
	      SBX(L)=0.0
	      SBX(L+1)=0.0
	      SBY(L)=0.0
	      SBY(LN)=0.0
            IF(ISCDRY(L).EQ.0)THEN
		    ISCDRY(L)=1
            ENDIF
	    ENDIF
	  ENDIF
C
        ENDDO
C
      ENDIF
C
      IF(ISDRY.EQ.99)THEN
       OPEN(1,FILE='RSTWD.INP',STATUS='UNKNOWN')
       OPEN(2,FILE='RSTWD.RCK',STATUS='UNKNOWN')
       CLOSE(2,STATUS='DELETE')
       OPEN(2,FILE='RSTWD.RCK',STATUS='UNKNOWN')
       DO L=2,LA
        READ(1,*)LDUM,IDUM,JDUM,ISCDRY(L),NATDRY(L),
     &           IMASKDRY(L),SUB(L),SVB(L) 
        WRITE(2,913)LDUM,IDUM,JDUM,ISCDRY(L),NATDRY(L),
     &           IMASKDRY(L),SUB(L),SVB(L),SUBO(L),SVBO(L) 
       ENDDO
       CLOSE(1)
       CLOSE(2)
      ENDIF
C
  913 FORMAT(6I5,4F7.3)

C
C      DO L=2,LA
C        IMASKDRY(L)=0
C        LMASKDRY(L)=.TRUE.
C      END DO
C
C
      DO L=2,LA
        IF(IMASKDRY(L).EQ.0) LMASKDRY(L)=.TRUE.
        IF(IMASKDRY(L).GT.0) LMASKDRY(L)=.FALSE.
      END DO
C
C**********************************************************************C
C
  101 FORMAT(I5)
  102 FORMAT(3I5,16F8.2)
C
C**********************************************************************C
C
C **  WRITE READ ERRORS ON RESTART
C
      GOTO 3000
 1000 WRITE(6,2000)
      STOP
 1001 WRITE(6,2001)L
      STOP
 1002 WRITE(6,2002)L
      STOP
 1003 WRITE(6,2003)L
      STOP
 1004 WRITE(6,2004)L
      STOP
 1005 WRITE(6,2005)L
      STOP
 1006 WRITE(6,2006)L
      STOP
 1007 WRITE(6,2007)L
      STOP
 1008 WRITE(6,2008)L
      STOP
 1009 WRITE(6,2009)L
      STOP
 1010 WRITE(6,2010)L
      STOP
 1011 WRITE(6,2011)L
      STOP
 1012 WRITE(6,2012)L
      STOP
 1013 WRITE(6,2013)L
      STOP
 1014 WRITE(6,2014)L
      STOP
 1015 WRITE(6,2015)L
      STOP
 1016 WRITE(6,2016)L
      STOP
 1017 WRITE(6,2017)L
      STOP
 1018 WRITE(6,2018)L
      STOP
 1019 WRITE(6,2019)L
      STOP
 1020 WRITE(6,2020)L
      STOP
 1021 WRITE(6,2021)L
      STOP
 1022 WRITE(6,2022)L
      STOP
 1023 WRITE(6,2023)L
      STOP
 1024 WRITE(6,2024)L
      STOP
 1025 WRITE(6,2025)L
      STOP
 1026 WRITE(6,2026)L
      STOP
 1027 WRITE(6,2027)L
      STOP
 1028 WRITE(6,2028)L
      STOP
 1029 WRITE(6,2029)L
      STOP
 1030 WRITE(6,2030)L
      STOP
 1031 WRITE(6,2031)L
      STOP
 1032 WRITE(6,2032)NMD
      STOP
 1033 WRITE(6,2033)L
      STOP
 3000 CONTINUE
C
C**********************************************************************C
C
  600 FORMAT(2X,'I,J,BELVOLD,BELVNEW',2I5,2F12.2)
  906 FORMAT(5E15.7)
  907 FORMAT(12E12.4)
  908 FORMAT(12I10)
 2000 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1000')
 2001 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1001 L =',I6)
 2002 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1002 L =',I6)
 2003 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1003 L =',I6)
 2004 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1004 L =',I6)
 2005 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1005 L =',I6)
 2006 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1006 L =',I6)
 2007 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1007 L =',I6)
 2008 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1008 L =',I6)
 2009 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1009 L =',I6)
 2010 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1010 L =',I6)
 2011 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1011 L =',I6)
 2012 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1012 L =',I6)
 2013 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1013 L =',I6)
 2014 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1014 L =',I6)
 2015 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1015 L =',I6)
 2016 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1016 L =',I6)
 2017 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1017 L =',I6)
 2018 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1018 L =',I6)
 2019 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1019 L =',I6)
 2020 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1020 L =',I6)
 2021 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1021 L =',I6)
 2022 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1022 L =',I6)
 2023 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1023 L =',I6)
 2024 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1024 L =',I6)
 2025 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1025 L =',I6)
 2026 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1026 L =',I6)
 2027 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1027 L =',I6)
 2028 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1028 L =',I6)
 2029 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1029 L =',I6)
 2030 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1030 L =',I6)
 2031 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1031 L =',I6)
 2032 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1032 NMD =',I6)
 2033 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1033 L =',I6)
 9696 FORMAT('  NEGATIVE DEPTH RESTART, L,I,J,HP,H1P = ',
     &           3I7,2F10.4)
C
C**********************************************************************C
C
      RETURN
      END
