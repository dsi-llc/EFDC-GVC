C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALSTEPD
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
C **  SUBROUTINE CALSTEP ESTIMATE THE CURRENT MAXIMUM TIME STEP SIZE
C **  FORM LINEAR STABILITY CRITERIA AND A FACTOR OF SAFETY
C 
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION QSUBOUT(LCM,KCM),QSUBINN(LCM,KCM)
      DIMENSION DTL1(LCM),DTL2(LCM),DTL3(LCM),DTL4(LCM)
C
C**********************************************************************C
C
      ITRNTMP=0
      DO NX=1,7
        ITRNTMP=ITRNTMP+ISTRAN(NX)
      ENDDO
C
      IF(N.LE.0)DTDYN=0.0
C
C **  CLEAN TIME STEP LOG FILE
C
      IF(N.LE.1) THEN
        OPEN(1,FILE='TIMSTP.LOG')
        CLOSE(1,STATUS='DELETE')
      ENDIF
C
      DTMIN=DT
      DTMAX=TIDALP
C
      DO L=2,LA
        DTL1(L)=DTMAX
        DTL2(L)=DTMAX
        DTL3(L)=DTMAX
        DTL4(L)=DTMAX
      ENDDO
C
C **  DETERMINE SOURCE/SINKS FOR SUBGRID SCALE CHANNEL EXCHANGES
C
      DO K=1,KC
        DO L=2,LA
          QSUBOUT(L,K)=0.0
          QSUBINN(L,K)=0.0
        ENDDO
      ENDDO
C
      IF(MDCHH.GE.1)THEN
        DO K=1,KC
        DO NMD=1,MDCHH
          LMDCHHT=LMDCHH(NMD)
          LMDCHUT=LMDCHU(NMD)
          LMDCHVT=LMDCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            QUKTMP=QCHANU(NMD)*DZC(K)
            QVKTMP=0.
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            QVKTMP=QCHANV(NMD)*DZC(K)
            QUKTMP=0.
          ENDIF
          IF(MDCHTYP(NMD).EQ.3)THEN
            QUKTMP=QCHANU(NMD)*DZC(K)
            QVKTMP=QCHANV(NMD)*DZC(K)
          ENDIF
          QSUBOUT(LMDCHHT,K)=QSUBOUT(LMDCHHT,K)
     &               +MIN(QUKTMP,0.)
     &               +MIN(QVKTMP,0.)
          QSUBINN(LMDCHHT,K)=QSUBINN(LMDCHHT,K)
     &               +MAX(QUKTMP,0.)
     &               +MAX(QVKTMP,0.)
          QSUBOUT(LMDCHUT,K)=QSUBOUT(LMDCHUT,K)
     &               -MAX(QUKTMP,0.)
          QSUBINN(LMDCHUT,K)=QSUBINN(LMDCHUT,K)
     &               -MIN(QUKTMP,0.)
          QSUBOUT(LMDCHVT,K)=QSUBOUT(LMDCHVT,K)
     &               -MAX(QVKTMP,0.)
          QSUBINN(LMDCHVT,K)=QSUBINN(LMDCHVT,K)
     &               -MIN(QVKTMP,0.)
        ENDDO
        ENDDO
      ENDIF
C
C **  METHOD 1: UPWIND DIFF IN MOMENTUM EQUATIONS
C
      DO K=1,KC
        DO L=2,LA
          IF(LMASKDRY(L))THEN
          LE=L+1
          LN=LNC(L)
          LS=LSC(L)
          LSE=LSEC(L)
          LNW=LNWC(L)
          KM=K-1
          VATUUU=0.25*(V(L,K)+V(L-1,K)+V(LN,K)+V(LNW,K))
          TMPUUU=ABS(U(L,K)/DXU(L))+ABS(VATUUU/DYU(L))
          DTTMP=DTMAX
          IF(TMPUUU.GT.0.0) DTTMP=1./TMPUUU
          DTL1(L)=MIN(DTL1(L),DTTMP)
c          write(*,*)L,LN,LS
          UATVVV=0.25*(U(L,K)+U(LS,K)+U(L+1,K)+U(LSE,K))
          TMPVVV=ABS(V(L,K)/DYV(L))+ABS(UATVVV/DXV(L))
          DTTMP=DTMAX
          IF(TMPVVV.GT.0.0) DTTMP=1./TMPVVV
          DTL1(L)=MIN(DTL1(L),DTTMP)
          UEAST=ABS(U(L,K))
          UWEST=ABS(U(L+1,K))
          VSOUTH=ABS(V(L,K))
          VNORTH=ABS(V(LN,K))
          TMPVVV=MAX(VSOUTH,VNORTH)
          TMPUUU=MAX(UEAST,UWEST)
          TMPVAL=TMPUUU/DXP(L)+TMPVVV/DYP(L)
          DTTMP=DTMAX
          IF(TMPVAL.GT.0.0)DTTMP=1./TMPVAL
          DTL1(L)=MIN(DTL1(L),DTTMP)
          ENDIF
        ENDDO
      ENDDO
C
C **  METHOD 2: POSITIVITY OF ADVECTED MATERIAL, DTL2
C
      IF(ITRNTMP.GE.1)THEN
      DO K=1,KC
        DO L=2,LA
          IF(LMASKDRY(L))THEN
          LE=L+1
          LN=LNC(L)
          LS=LSC(L)
          KM=K-1
          TOP=DZC(K)*H1P(L)*DXYP(L)
          QXPLUS=UHDY2(LE,K)*DZC(K)
          QXPLUS=MAX(QXPLUS,0.0)
          QYPLUS=VHDX2(LN,K)*DZC(K)
          QYPLUS=MAX(QYPLUS,0.0)
          QZPLUS=W2(L,K)*DXYP(L)
          QZPLUS=MAX(QZPLUS,0.0)
          QXMINS=UHDY2(L,K)*DZC(K)
          QXMINS=-MIN(QXMINS,0.0)
          QYMINS=VHDX2(L,K)*DZC(K)
          QYMINS=-MIN(QYMINS,0.0)
          QZMINS=W2(L,KM)*DXYP(L)
          QZMINS=-MIN(QZMINS,0.0)
          QTOTAL=QSUM(L,K)+QSUBOUT(L,K)+QSUBINN(L,K)
          QSRC=-MIN(QTOTAL,0.0)
          BOT=QXPLUS+QYPLUS+QZPLUS+QXMINS+QYMINS+QZMINS+QSRC
          IF(BOT.GT.0.0)THEN
            DTTMP=TOP/BOT
            DTL2(L)=MIN(DTL2(L),DTTMP)
            IF(DTTMP.LT.0.0)THEN
              WRITE(6,880)IL(L),JL(L),K,TOP,QXPLUS,QYPLUS,QZPLUS,
     &                    QXMINS,QYMINS,QZMINS,QSRC
              WRITE(8,880)IL(L),JL(L),K,TOP,QXPLUS,QYPLUS,QZPLUS,
     &                    QXMINS,QYMINS,QZMINS,QSRC
            ENDIF
          ENDIF
          ENDIF
        ENDDO
      ENDDO
      ENDIF
C
C **  METHOD 3: implicit BOTTOM FRICTION AND ROTATIONAL ACCELERATION DAMPING
C
      DO L=2,LA
        IF(LMASKDRY(L))THEN
          TMPVAL=SUB(L)+SUB(L+1)+SVB(L)+SVB(LNC(L))
          IF(TMPVAL.LT.0.5)THEN
            LN=LNC(L)
            TAUBC=QQ(L,0)/CTURB2
            UCTR=0.5*(U(L,1)+U(L+1,1))
            VCTR=0.5*(V(L,1)+V(LN,1))
            UHMAG=HP(L)*SQRT(UCTR*UCTR+VCTR*VCTR)
            IF(UHMAG.GT.0.0)THEN
              FRIFRE=TAUBC/UHMAG
              FRIFRE2=FRIFRE*FRIFRE
              ACACTMP=(CAC(L,KC)*HPI(L)*DXYIP(L))**2
              IF(ACACTMP.GT.FRIFRE2)THEN
                DTTMP=2.*FRIFRE/(ACACTMP-FRIFRE2)
                DTL3(L)=MIN(DTL3(L),DTTMP)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
C **  METHOD 4: LIMIT RATE OF DEPTH CHANGE
C
      DO L=2,LA
        IF(LMASKDRY(L))THEN
          IF(HP(L).LT.HDRY/10.) HP(L)=HDRY/10.
          IF(H1P(L).LT.HDRY/10.) H1P(L)=HDRY/10.
cjah avoid divide by zero
          TESTTEMP=MAX(ABS(HP(L)-H1P(L)),1.E-06)
          TMPVAL=DTDYN*HP(L)/TESTTEMP
cjah      TMPVAL=DTDYN*HP(L)/ABS(HP(L)-H1P(L)+1.E-6)
            DTL4(L)=DTSSDHDT*TMPVAL
        ENDIF
      ENDDO
C
C **  CHOOSE THE MINIMUM OF THE THREE METHODS
C
      DTL1MN=2.*DTMAX
      DTL2MN=2.*DTMAX
      DTL3MN=2.*DTMAX
      DTL4MN=2.*DTMAX
      DTTMP=2.*DTMAX
      DO L=2,LA
        IF(LMASKDRY(L))THEN
          IF(DTL1MN.GT.DTL1(L))THEN
            DTL1MN=DTL1(L)
            L1LOC=L
          ENDIF
          IF(DTL2MN.GT.DTL2(L))THEN
            DTL2MN=DTL2(L)
            L2LOC=L
          ENDIF
          IF(DTL3MN.GT.DTL3(L))THEN
            DTL3MN=DTL3(L)
            L3LOC=L
          ENDIF
          IF(DTL4MN.GT.DTL4(L))THEN
            DTL4MN=DTL4(L)
            L4LOC=L
          ENDIF
        ENDIF
      ENDDO
      IF(DTTMP.GT.DTL1MN)THEN
        DTTMP=DTL1MN
        LLOC=L1LOC
      ENDIF
      IF(DTTMP.GT.DTL2MN)THEN
        DTTMP=DTL2MN
        LLOC=L2LOC
      ENDIF
      IF(DTTMP.GT.DTL3MN)THEN
        DTTMP=DTL3MN
        LLOC=L3LOC
      ENDIF
      IF(DTTMP.GT.DTL4MN)THEN
        DTTMP=DTL4MN
        LLOC=L4LOC
      ENDIF
      DTRAW=DTTMP
C
C **  CHECK IF CURVATURE INSTABILITY IS CONTROLLED
C
c      CACDTMX=-1000.
c      DO L=2,LA
c        IF(LMASKDRY(L))THEN
c          CACTMP=ABS(DTTMP*CAC(L,KC)*HPI(L)*DXYIP(L))
c          CACDTMX=MAX(CACDTMX,CACTMP)
c        ENDIF
c      ENDDO
c      CACAMP=SQRT(1.+CACDTMX*CACDTMX)
C
C **  APPLY A SAFTY FACTOR
C
      DTTMP=DTSSFAC*DTTMP
C
C **  MAKE A MULTIPLE OF OF DTMIN
C
      TIMEDAY=TIMESEC/86400.
      IF(DTTMP.LE.DTMIN)THEN
        WRITE(8,800)TIMEDAY,DTTMP,DTMIN,IL(LLOC),JL(LLOC)
        WRITE(6,800)TIMEDAY,DTTMP,DTMIN,IL(LLOC),JL(LLOC)
        WRITE(8,801)IL(L1LOC),JL(L1LOC),DTL1MN
        WRITE(6,801)IL(L1LOC),JL(L1LOC),DTL1MN
        WRITE(8,802)IL(L2LOC),JL(L2LOC),DTL2MN
        WRITE(6,802)IL(L2LOC),JL(L2LOC),DTL2MN
        WRITE(8,803)IL(L3LOC),JL(L3LOC),DTL3MN
        WRITE(6,803)IL(L3LOC),JL(L3LOC),DTL3MN
        WRITE(8,804)IL(L4LOC),JL(L4LOC),DTL4MN
        WRITE(6,804)IL(L4LOC),JL(L4LOC),DTL4MN
        DTTMP=DTMIN
      ELSE
        TMPVAL=DTTMP/DTMIN
        ITMPR=NINT(TMPVAL)
        RTMPR=FLOAT(ITMPR)
        IF(RTMPR.LT.TMPVAL)THEN
          DTTMP=RTMPR*DTMIN
        ELSE
          DTTMP=(RTMPR-1.)*DTMIN
        ENDIF
      ENDIF
C
C **  SET TO MINIMUM TIME STEP ON STARTUP
C
      IF(N.EQ.0)DTTMP=DTMIN
C
C **  RESTRICT INCREASE IN TIME STEP TO DTMIN
C
C      DTDYN2=2.*DTDYN
C      IF(DTTMP.GT.DTDYN2)THEN
C        DTTMP=DTDYN2
C      ENDIF
C
      DTDYNP=DTDYN+DTMIN
      IF(DTTMP.GT.DTDYNP)THEN
        DTTMP=DTDYNP
      ENDIF
C
      DTDYN=DTTMP
C
C **  SET INCREMENTAL INCREASE IN OUTPUT COUNTER
C
      NINCRMT=NINT(DTDYN/DTMIN)
C
C **  ADJUST INCREMENT FOR N TO LAND EVENLY ON NTSPTC
C
      RTCTMP=FLOAT(N)/FLOAT(NTSPTC)
      NTCTMP=RTCTMP
      NTMP=(1+NTCTMP)*NTSPTC-N
      IF(NINCRMT.GT.NTMP)THEN
        NINCRMT=NTMP
        DTDYN=FLOAT(NTMP)*DTMIN
      ENDIF
C
C **  WRITE TO TIME STEP LOG FILE
C
c      CACAMP=CACAMP-1.0
C      OPEN(1,FILE='TIMSTP.LOG',POSITION='APPEND')
C      WRITE(1,100)N,NINCRMT,NTCTMP,IL(LLOC),JL(LLOC),DTRAW,DTDYN,
C     &            DTL1MN,DTL2MN,DTL3MN,CACAMP
C      CLOSE(1)
C
  100 FORMAT(5I5,5F12.5,E13.5)     
C  100 FORMAT('  TIME STEP N,TDY,TMN,TMX,CAM = ',2I5,4F12.5,E13.5)     
  101 FORMAT(3I5,E13.5)     
  800 FORMAT('  TIME,DTDYN,DTMIN,I,J = ',F12.5,2E12.4,2I7)     
  801 FORMAT('  MOM  ADV,I,J,DTM = ',2I5,E13.4)     
  802 FORMAT('  MASS ADV,I,J,DTM = ',2I5,E13.4)     
  803 FORMAT('  CURV ACC,I,J,DTM = ',2I5,E13.4) 
  804 FORMAT('  LIM DHDT,I,J,DTM = ',2I5,E13.4) 
  880 FORMAT(3I5,8E13.4) 
 8899 FORMAT(' DT3 ERROR ',2I5,6E13.5)   
C
C**********************************************************************C
C
      RETURN
      END
