C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE NEGDEPHOUS(QCHANUT,QCHANVT,ISTL)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 02/15/2002        John Hamrick       02/11/2002       John Hamrick
C  added alternate sor equation solver relax2t
C----------------------------------------------------------------------C  
C
C**********************************************************************C
C
C ** SUBROUTINE NEGDEP CHECK EXTERNAL SOLUTION FOR NEGATIVE DEPTHS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION QCHANUT(NCHANM),QCHANVT(NCHANM)
C
C**********************************************************************C
C
C **  CHECK FOR NEGATIVE DEPTHS
C
C----------------------------------------------------------------------C
C
      IF(ISNEGH.GE.1)THEN
      INEGFLG=0
C
      DO L=2,LA
      IF(HP(L).LT.0.)THEN
cjah  INEGFLG=1
      LN=LNC(L)
      WRITE(6,1111)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE(8,1111)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
cjah start
      HP(L)=HDRY
      IMASKDRY(L)=2
CJMH CHANGED TO LOGICAL
C      LMASKDRY(L)=0
	LMASKDRY(L)=.FALSE.
CJMH CHANGED TO LOGICAL
cjah end
      ENDIF
      ENDDO
C
      DO L=2,LA
      IF(HU(L).LT.0.AND.SUBO(L).GT.0.5)THEN
cjah  INEGFLG=1
      LN=LNC(L)
      WRITE(6,1112)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE(8,1112)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
cjah
      HU(L)=HDRY
cjah
      ENDIF
      ENDDO
C
C
      DO L=2,LA
      IF(HV(L).LT.0.AND.SVBO(L).GT.0.5)THEN
cjah  INEGFLG=1
      LN=LNC(L)
      WRITE(6,1113)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE(8,1113)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
cjah
      HV(L)=HDRY
cjah
      ENDIF
      ENDDO
C
C      DO L=2,LA
C      IF(HU(L).LT.0.)THEN
C      INEGFLG=1
C      LN=LNC(L)
C      LS=LSC(L)
C      WRITE (6,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (6,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (6,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (6,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (6,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      WRITE (8,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (8,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (8,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (8,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (8,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      ENDIF
C      ENDDO
C
C
C      DO L=2,LA
C      IF(HV(L).LT.0.)THEN
C      INEGFLG=1
C      LN=LNC(L)
C      LS=LSC(L)
C      WRITE (6,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (6,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (6,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (6,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (6,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      WRITE (8,6060)IL(L),JL(L),ISCDRY(L),HP(L),H1P(L),H2P(L)
C      WRITE (8,6061)IL(L),JL(L),ISCDRY(L-1),HU(L),H1U(L)
C      WRITE (8,6062)IL(L),JL(L),ISCDRY(L+1),HU(L+1),H1U(L+1)
C      WRITE (8,6063)IL(L),JL(L),ISCDRY(LS),HV(L),H1V(L)
C      WRITE (8,6064)IL(L),JL(L),ISCDRY(LN),HV(LN),H1V(LN)
C      ENDIF
C      ENDDO
C
      ENDIF
C
      IF(ISNEGH.EQ.2)THEN
      IF(INEGFLG.EQ.1)THEN
C
C
        IF(MDCHH.GT.0)THEN
        DO NMD=1,MDCHH
          WRITE(8,8000)
          LHOST=LMDCHH(NMD)
          IHOST=IL(LHOST)        
          JHOST=JL(LHOST)        
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
          ICHNU=IL(LCHNU)
          JCHNU=JL(LCHNU)
          SRFCHAN=HP(LCHNU)+BELV(LCHNU)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNU)+BELV(LCHNU)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(8,8001)N,NMD,MDCHTYP(NMD),ICHNU,JCHNU,ISCDRY(LCHNU),
     &                 SRFCHAN,HP(LCHNU),P1(LCHNU),H1P(LCHNU)
          WRITE(8,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),P1(LHOST),H1P(LHOST)
          WRITE(8,8003)QCHANU(NMD),QCHANUT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
          ICHNV=IL(LCHNV)
          JCHNV=JL(LCHNV)
          SRFCHAN=HP(LCHNV)+BELV(LCHNV)
          SRFHOST=HP(LHOST)+BELV(LHOST)
          SRFCHAN1=H1P(LCHNV)+BELV(LCHNV)
          SRFHOST1=H1P(LHOST)+BELV(LHOST)
          WRITE(8,8001)N,NMD,MDCHTYP(NMD),ICHNV,JCHNV,ISCDRY(LCHNV),
     &                 SRFCHAN,HP(LCHNV),SRFCHAN1,H1P(LCHNV)
          WRITE(8,8002)IHOST,JHOST,ISCDRY(LHOST),
     &                 SRFHOST,HP(LHOST),SRFHOST1,H1P(LHOST)
          WRITE(8,8003)QCHANV(NMD),QCHANVT(NMD),CCCCHU(NMD),CCCCHV(NMD)
          ENDIF
          WRITE(8,8004)
        ENDDO
        ENDIF
C
        CALL RESTOUT(1)
C
        OPEN(1,FILE='EQCOEF.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     &               FPTMP(L),SURFTMP
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='EQTERM.OUT',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     &               HRVO(L),HUTMP(L),HVTMP(L)
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='CFLMAX.OUT')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='CFLMAX.OUT')
        DO L=2,LA
         WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
         WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
         WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
         WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
        ENDDO
        CLOSE(1)
C
        OPEN(1,FILE='NEGDEPDIA.OUT')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='NEGDEPDIA.OUT')
        DO L=2,LA
         WRITE(1,1993)IL(L),JL(L),IMASKDRY(L),BELV(L),HP(L),H1P(L),
     &                SUB(L),SUB(L+1),SVB(L),SVB(LNC(L))
        ENDDO
        CLOSE(1)
C
        CALL EE_LINKAGE(0)
        CALL SURFPLT
        STOP
C
      ENDIF
      ENDIF
C
 1001 FORMAT(2I5,10(1X,E12.4))
 1002 FORMAT(3I4,10(1X,E9.2))
 1991 FORMAT(2I5,12F8.3)
 1993 FORMAT(3I5,12F8.3)
 1992 FORMAT(10X,12F8.3)
 1111 FORMAT(' NEG DEPTH AT CELL CENTER')
 1112 FORMAT(' NEG DEPTH AT WEST FACE')
 1113 FORMAT(' NEG DEPTH AT SOUTH FACE')
 6060 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HP,H1P,H2P =',3(2X,E12.4))
 6061 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUW,H1UW =',2(2X,E12.4))
 6062 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HUE,H1UE =',2(2X,E12.4))
 6063 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVS,H1VS =',2(2X,E12.4))
 6064 FORMAT('  NEG DEPTH AT I,J =',2I4,'  HVN,H1VN =',2(2X,E12.4))
 8001 FORMAT(I7,5I5,4E13.4)
 8002 FORMAT(17X,3I5,4E13.4)
 8003 FORMAT(32X,4E13.4)
 8000 FORMAT('    N    NMD  MTYP   I    J  IDRY      P           H',
     &       '           P1           H1')
 8004 FORMAT('                                     QCHANU',
     &       '       QCHANUT      CCCCHU       CCCCHV ')   
C
C**********************************************************************C
C
      RETURN
      END
