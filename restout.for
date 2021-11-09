C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTOUT(IRSTYP)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C  11/14/2001       john hamric        11/14/2001       john hamric
C   add output of bed load transport QSBDLDX  QSBDLDY
C
C----------------------------------------------------------------------C
C
C **  SUBROUTINE RESTOUT WRITES A RESTART FILE
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION FVOLSED(KBM),FVOLSND(KBM)

C
C**********************************************************************C
C
      IF(IRSTYP.EQ.0)THEN
        OPEN(99,FILE='RESTART.OUT',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='RESTART.OUT',STATUS='UNKNOWN')
      ENDIF
C
      IF(IRSTYP.EQ.1)THEN
        OPEN(99,FILE='CRASHST.OUT',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='CRASHST.OUT',STATUS='UNKNOWN')
      ENDIF
C
C**********************************************************************C
C
      IF(ISRESTO.EQ.-11)THEN
C
      DO L=1,LC
C     HP(L)=HMP(L)
C     H1P(L)=HMP(L)
C     HWQ(L)=HMP(L)
C     H2WQ(L)=HMP(L)
      HP(L)=-BELV(L)
      H1P(L)=-BELV(L)
      HWQ(L)=-BELV(L)
      H2WQ(L)=-BELV(L)
      UHDYE(L)=0.
      UHDY1E(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      ENDDO
C
      DO K=0,KC
      DO L=1,LC
      QQ(L,K)=QQMIN
      QQ1(L,K)=QQMIN
      QQL(L,K)=QQLMIN
      QQL1(L,K)=QQLMIN
      DML(L,K)=DMLMIN
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=1,LC
      U(L,K)=0.
      U1(L,K)=0.
      V(L,K)=0.
      V1(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=MAX(SAL(L,K),0.)
      SAL1(L,K)=MAX(SAL1(L,K),0.)
      ENDDO
      ENDDO
C
      ENDIF
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
      WRITE(99,909)N,TIME
C
      DO L=2,LA
C     WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L)
      WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
      WRITE(99,907)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      WRITE(99,907)(U(L,K),K=1,KC)
      WRITE(99,907)(U1(L,K),K=1,KC)
      WRITE(99,907)(V(L,K),K=1,KC)
      WRITE(99,907)(V1(L,K),K=1,KC)
      WRITE(99,907)(QQ(L,K),K=0,KC)
      WRITE(99,907)(QQ1(L,K),K=0,KC)
      WRITE(99,907)(QQL(L,K),K=0,KC)
      WRITE(99,907)(QQL1(L,K),K=0,KC)
      WRITE(99,907)(DML(L,K),K=0,KC)
      IF(ISCO(1).EQ.1)THEN
       WRITE(99,907)(SAL(L,K),K=1,KC)
       WRITE(99,907)(SAL1(L,K),K=1,KC)
      ENDIF
      IF(ISCO(2).EQ.1)THEN
       WRITE(99,907)(TEM(L,K),K=1,KC)
       WRITE(99,907)(TEM1(L,K),K=1,KC)
      ENDIF
      IF(ISCO(3).EQ.1)THEN
       WRITE(99,907)(DYE(L,K),K=1,KC)
       WRITE(99,907)(DYE1(L,K),K=1,KC)
      ENDIF
      IF(ISCO(4).EQ.1)THEN
       WRITE(99,907)SFLSBOT(L),(SFL(L,K),K=1,KC)
       WRITE(99,907)SFLSBOT(L),(SFL2(L,K),K=1,KC)
      ENDIF
      IF(ISCO(5).EQ.1)THEN
       DO NT=1,NTOX
       WRITE(99,907)(TOXB(L,K,NT),K=1,KB)
       WRITE(99,907)(TOX(L,K,NT),K=1,KC)
       WRITE(99,907)(TOXB1(L,K,NT),K=1,KB)
       WRITE(99,907)(TOX1(L,K,NT),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCO(6).EQ.1)THEN
       DO NS=1,NSED
       WRITE(99,907)(SEDB(L,K,NS),K=1,KB)
       WRITE(99,907)(SED1(L,K,NS),K=1,KC)
       WRITE(99,907)(SEDB1(L,K,NS),K=1,KB)
       WRITE(99,907)(SED1(L,K,NS),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCO(7).EQ.1)THEN
       DO NS=1,NSND
       WRITE(99,907)(SNDB(L,K,NS),K=1,KB)
       WRITE(99,907)(SND(L,K,NS),K=1,KC)
       WRITE(99,907)(SNDB1(L,K,NS),K=1,KB)
       WRITE(99,907)(SND1(L,K,NS),K=1,KC)
       ENDDO
      ENDIF
      IF(ISCO(6).EQ.1.OR.ISCO(7).EQ.1)THEN
       WRITE(99,907)(HBED(L,K),K=1,KB)
       WRITE(99,907)(HBED1(L,K),K=1,KB)
       WRITE(99,907)(VDRBED(L,K),K=1,KB)
       WRITE(99,907)(VDRBED1(L,K),K=1,KB)
      ENDIF
      ENDDO
C
      DO M=1,4
      IF(ISCO(M).EQ.1)THEN
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDIF
      ENDDO
C
      IF(ISCO(5).EQ.1)THEN
      DO NT=1,NTOX
       M=MSVTOX(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
      IF(ISCO(6).EQ.1)THEN
      DO NT=1,NSED
       M=MSVSED(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
      IF(ISCO(7).EQ.1)THEN
      DO NT=1,NSND
       M=MSVSND(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       ENDDO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       ENDDO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       ENDDO
C
      ENDDO
      ENDIF
C
      DO L=2,LA
        WRITE(99,907)QSUME(L),(QSUM(L,K),K=1,KC)
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        WRITE(99,910)IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),
     &               IMDCHV(NMD),JMDCHV(NMD),QCHANU(NMD),QCHANV(NMD)
        ENDDO
      ENDIF
C
      IF(ISGWIE.GE.1)THEN
        DO L=2,LA
        WRITE(99,907)AGWELV(L),AGWELV1(L)
        ENDDO
      ENDIF
C
C
      IF(ISWAVE.GE.1)THEN
       OPEN(1,FILE='WVQWCP.OUT',STATUS='UNKNOWN')
       CLOSE(1, STATUS='DELETE')
       OPEN(1,FILE='WVQWCP.OUT',STATUS='UNKNOWN')
       DO L=2,LA
       WRITE(1,911)IL(L),JL(L),QQWV1(L),QQWV2(L),QQWV3(L),QQWC(L),
     &                         QQWCR(L),QQ(L,0) 
       ENDDO
       CLOSE(1)
      ENDIF     
C
      IF(ISCO(1).GE.1.)THEN
       OPEN(1,FILE='SALT.RST',STATUS='UNKNOWN')
       CLOSE(1, STATUS='DELETE')
       OPEN(1,FILE='SALT.RST',STATUS='UNKNOWN')
       DO L=2,LA
        WRITE(1,912)L,IL(L),JL(L),(SAL(L,K),K=1,KC) 
       ENDDO
       CLOSE(1)
      ENDIF
C
      IF(ISCO(2).GE.1.)THEN
       OPEN(1,FILE='TEMP.RST',STATUS='UNKNOWN')
       OPEN(2,FILE='TEMPB.RST',STATUS='UNKNOWN')
       CLOSE(1, STATUS='DELETE')
       CLOSE(2, STATUS='DELETE')
       OPEN(1,FILE='TEMP.RST',STATUS='UNKNOWN')
       OPEN(2,FILE='TEMPB.RST',STATUS='UNKNOWN')
       DO L=2,LA
        WRITE(1,912)L,IL(L),JL(L),(TEM(L,K),K=1,KC) 
        WRITE(2,912)L,IL(L),JL(L),(TEMB(L,K),K=1,KBHM)
       ENDDO
       CLOSE(1)
       CLOSE(2)
      ENDIF
C
      IF(ISDRY.EQ.99)THEN
       OPEN(1,FILE='RSTWD.OUT',STATUS='UNKNOWN')
       CLOSE(1, STATUS='DELETE')
       OPEN(1,FILE='RSTWD.OUT',STATUS='UNKNOWN')
       DO L=2,LA
        WRITE(1,913)L,IL(L),JL(L),ISCDRY(L),NATDRY(L),
     &           IMASKDRY(L),SUB(L),SVB(L),SUBO(L),SVBO(L) 
       ENDDO
       CLOSE(1)
      ENDIF
C
C**********************************************************************C
C
C **  OUTPUT SALINITY AND TEMPATURE DATA ASSIMILATION
C
      IF(NLCDA.GT.0)THEN
C
      OPEN(1,FILE='STDATAASM.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='STDATAASM.OUT')
      DO J=1,NLCDA
	  DO I=1,NTC
          WRITE(1,5678)J,I,FSALASM(I,J),FVOLASM(I,J),FTEMASM(I,J)
	  ENDDO
	ENDDO
C
      ENDIF
C
C **  OUTPUT CORRECTIVE FLOW FOR WATER SURFACE ELEVATION DATA ASSIMILATION
C
      IF(NLWSEDA.GT.0)THEN
C
      TMP1=1.
      IF(ITIMSOL.EQ.0)THEN
       TMP1=(FLOAT(NTSTBC)+1.)/(FLOAT(NTSTBC))
	ENDIF
	TMP1=TMP1*FLOAT(NTSPTC)
C
      OPEN(1,FILE='WSEDATAASM.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WSEDATAASM.OUT')
      DO J=1,NLWSEDA
	  DO I=1,NTC
	    QWSEASMT=QWSEASM(I,J)/TMP1
          TIMDAY=TBEGIN+FLOAT(I)-0.5
          WRITE(1,5679)J,TIMDAY,QWSEASMT
	  ENDDO
	ENDDO
C
      ENDIF
C
 5678 FORMAT(2I6,3E14.5)
 5679 FORMAT(I6,F10.2,3E14.5)
C
C**********************************************************************C
C
      CLOSE(99)
C
C     IF(ISRESTO.EQ.-2) CALL OUT3D
C
C**********************************************************************C
C
C     WRITE BED CONDITIONS AT RESTART
C
      IWRTBED=0
C
      IF(ISTRAN(6).EQ.1.OR.ISTRAN(7).EQ.1)THEN
      IF(ISDTXBUG.EQ.1.OR.N.EQ.NTS)THEN
         IWRTBED=1
      ENDIF
	ENDIF
C
      IF(ISTRAN(6).EQ.1.OR.ISTRAN(7).EQ.1)THEN
      IF(IRSTYP.EQ.1)THEN
         IWRTBED=1
      ENDIF
	ENDIF
C
      IF(IWRTBED.EQ.1)THEN
C
      OPEN(1,FILE='BEDRST_SED.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_SED.OUT')
      WRITE(1,111)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SEDB(L,K,1),K=1,KB)
        IF(NSED.GT.1) THEN
          DO NX=2,NSED
           WRITE(1,102)(SEDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_SND.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_SND.OUT')
      WRITE(1,112)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SNDB(L,K,1),K=1,KB)
        IF(NSND.GT.1)THEN
          DO NX=2,NSND
           WRITE(1,102)(SNDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_VDR.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_VDR.OUT')
      WRITE(1,113)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(VDRBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_VDRS.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_VDRS.OUT')
      WRITE(1,113)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(VDRBED(L,K),K=1,KB)
        WRITE(1,102)(VDRBEDSED(L,K),K=1,KB)
        WRITE(1,102)(VDRBEDSND(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_VFS.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_VFS.OUT')
      WRITE(1,113)  
      DO L=2,LA
	  DO K=1,KB
          FVOLSED(K)=0.0
	    FVOLSND(K)=0.0
	    IF(K.LE.KBT(L))THEN
	      DO NS=1,NSED
              FVOLSED(K)=FVOLSED(K)+VFRBED(L,K,NS)
	      ENDDO
            DO NX=1,NSND
              NS=NSED+NX
	        FVOLSND(K)=FVOLSND(K)+VFRBED(L,K,NS)
	      ENDDO
          ENDIF
	  ENDDO
        WRITE(1,101)IL(L),JL(L),(FVOLSED(K),K=1,KB)
        WRITE(1,102)(FVOLSND(K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_POR.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_POR.OUT')
      WRITE(1,114)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(PORBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_ZHB.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_ZHB.OUT')
      WRITE(1,115)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),(HBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_BDN.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_BDN.OUT')
      WRITE(1,116)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(BDENBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_ELV.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_ELV.OUT')
      WRITE(1,117)
	RVAL=0.
	TMP1=0.  
	TMP2=0.  
	TMP3=0.  
	TMP4=0.  
      DO L=2,LA
	RVAL=RVAL+1.
	TMP1=TMP1+ZELBEDA(L) 
	TMP2=TMP2+HBEDA(L) 
	TMP3=TMP3+BELV(L) 
	TMP4=TMP4+HP(L)  
        SURF=HP(L)+BELV(L)
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),BELV(L),HP(L),SURF
      ENDDO
	TMP1=TMP1/RVAL
	TMP2=TMP2/RVAL
	TMP3=TMP3/RVAL
	TMP4=TMP4/RVAL
	IDUM=0
	JDUM=0
        WRITE(1,101)IDUM,JDUM,TMP1,TMP2,TMP3,TMP4
      CLOSE(1)
C
      OPEN(1,FILE='WATRST_SED.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WATRST_SED.OUT')
      WRITE(1,118)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SED(L,K,1),K=1,KC)
        IF(NSED.GT.1) THEN
          DO NX=2,NSED
           WRITE(1,102)(SED(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='WATRST_SND.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WATRST_SND.OUT')
      WRITE(1,119)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SND(L,K,1),K=1,KC)
        IF(NSND.GT.1)THEN
          DO NX=2,NSND
           WRITE(1,102)(SND(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_BDL.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_BDL.OUT')
      WRITE(1,120)  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),QSBDLDX(L,1),QSBDLDY(L,1)
        IF(NSND.GT.1)THEN
          DO NX=2,NSND
           WRITE(1,102)QSBDLDX(L,NX),QSBDLDY(L,NX)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
C
      OPEN(1,FILE='BEDRST_TOX.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='BEDRST_TOX.OUT')
	DO NT=1,NTOX
      WRITE(1,121)NT  
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(TOXB(L,K,NT),K=1,KB)
      ENDDO
      ENDDO
      CLOSE(1)
C
      ENDIF
c
c      OPEN(1,FILE='RSTDXDY.INP',STATUS='UNKNOWN')
c      DO J=1,JC
c      DO I=1,IC
c      L=LIJ(I,J)
c      IF(IJCT(I,J).GE.1.AND.IJCT(I,J).LT.9)THEN
c        WRITE(1,339)IL(L),JL(L),DXP(L),DYP(L),HP(L),BELV(L),
c     &              ZBR(L)
C       WRITE(1,339)I,J,DXTMP,DYTMP,H(I,J),BELVIJ(I,J),
C    &              ZBRIJ(I,J)
c      ENDIF
c      ENDDO
c      ENDDO
c      CLOSE(1)
c
C
  339 FORMAT(2I5,6F14.5)
C
  101 FORMAT(2I5,18E13.5)
  102 FORMAT(10X,18E13.5)
  111 FORMAT('   IL   JL    SEDBT(K=1,KB)')
  112 FORMAT('   IL   JL    SNDBT(K=1,KB)')
  113 FORMAT('   IL   JL    VRDBED(K=1,KB)')
  114 FORMAT('   IL   JL    PORBED(K=1,KB)')
  115 FORMAT('   IL   JL    ZBEDB        HBEDT        HBED(K=1,KB)')
  116 FORMAT('   IL   JL    BDENBED(K=1,KB)')
  117 FORMAT('   IL   JL    ZBEDB        HBEDT        BELV',
     &       '        HWCOL        SELV')
  118 FORMAT('   IL   JL    SEDT(K=1,KC)')
  119 FORMAT('   IL   JL    SNDT(K=1,KC)')
  120 FORMAT('   IL   JL    QSBDLDX      QSBDLDY')
  121 FORMAT('   IL   JL    TOXB(K=1,KB,NT)  NT = ',I5)
C
C**********************************************************************C
C
  906 FORMAT(26E17.8)
  907 FORMAT(26E17.8)
CMRM      906 FORMAT(1P,5E17.8)
CMRM      907 FORMAT(1P,13E17.8)
  908 FORMAT(26I12)
  909 FORMAT(I20,4X,F12.4)
  910 FORMAT(6I5,2X,E17.8,2X,E17.8)
  911 FORMAT(2I5,2X,26E13.4)
CMRM      910 FORMAT(6I5,2X,1P,E17.8,2X,E17.8)
CMRM      911 FORMAT(2I5,2X,1P,6E13.4)
  912 FORMAT(3I5,26F7.3)
  913 FORMAT(6I5,26F7.3)
C
C**********************************************************************C
C
      RETURN
      END
