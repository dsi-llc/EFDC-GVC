C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SUBCHAN(QCHANUT,QCHANVT,IACTIVE,RLAMN,RLAMO,DELT,
     &IACTALL)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C----------------------------------------------------------------------C  
c
C**********************************************************************C
C
C ** SUBROUTINE SUBCHAN CALCULATES SUBGRID CHANNEL INTERACTIONS AND IS
C ** CALLED FROM CALPUV2TC 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
c
      DIMENSION QCHANUT(NCHANM),QCHANVT(NCHANM)
      DIMENSION IACTIVE(NCHANM)
C
C----------------------------------------------------------------------C
C
      IF(MDCHH.GE.1)THEN
        IACTALL=0
C
        DO NMD=1,MDCHH
          CCCCHU(NMD)=0.0
          CCCCHV(NMD)=0.0
          CCCCHH(NMD)=0.0
          LHOST=LMDCHH(NMD)
          LCHNU=LMDCHU(NMD)
          LCHNV=LMDCHV(NMD)
C         X-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.1)THEN
            IACTIVE(NMD)=0
            SRFHOST=H1P(LHOST)+BELV(LHOST)
            SRFCHAN=H1P(LCHNU)+BELV(LCHNU)
            IF(SRFCHAN.GT.SRFHOST)THEN
c            IF(H1P(LCHNU).GT.HWET)THEN
            IF(H1P(LCHNU).GT.HDRY)THEN
              HCHNCOR=HWET
              IACTIVE(NMD)=1
              IACTALL=IACTALL+1
C              IF(ISCDRY(LCHNU).EQ.1) IACTIVE=1
            ENDIF
            ENDIF
            IF(SRFHOST.GT.SRFCHAN)THEN
            IF(H1P(LHOST).GT.HDRY)THEN
              HCHNCOR=HDRY
              IACTIVE(NMD)=1
              IACTALL=IACTALL+1
C              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            ENDIF
            IF(HP(LHOST).LE.0.0.OR.HP(LCHNU).LE.0.0)THEN
              IF(IACTIVE(NMD).EQ.1)THEN
                IACTALL=IACTALL-1
                IACTIVE(NMD)=0  !PMC
              ENDIF
            ENDIF
            IF(IACTIVE(NMD).EQ.1)THEN
              WCHAN=DXP(LCHNU)
C              RLCHN=0.5*DYP(LCHNU)+0.25*DYP(LHOST)
C              HCHAN=0.5*DYP(LCHNU)*HP(LCHNU)+0.25*DYP(LHOST)*HP(LHOST)
              RLCHN=0.5*DYP(LCHNU)+CHANLEN(NMD)
              HCHAN=0.5*DYP(LCHNU)*H1P(LCHNU)+CHANLEN(NMD)*H1P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBX(LCHNU)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*ABS(QCHANUT(NMD)))
C                HCHANT=HCHAN-HCHNCOR
C                HCHANT=MAX(HCHANT,0.0)
C                CCCCHV(NMD)=DELT*HCHANT*WCHAN/RLCHN
                CCCCHV(NMD)=DELT*HCHAN*WCHAN/RLCHN
              ENDIF
            ENDIF
          ENDIF
C         Y-DIRECTION CHANNEL
          IF(MDCHTYP(NMD).EQ.2)THEN
            IHCHMX=0
            IHCHMN=0
            IACTIVE(NMD)=0
            SRFHOST=H1P(LHOST)+BELV(LHOST)
            SRFCHAN=H1P(LCHNV)+BELV(LCHNV)
            IF(SRFCHAN.GT.SRFHOST)THEN
C            IF(H1P(LCHNV).GT.HWET)THEN
            IF(H1P(LCHNV).GT.HDRY)THEN
C              HCHNMX=H1P(LCHNV)+BELV(LCHNV)-BELV(LHOST)
              HCHNMX=-H1P(LCHNV)
              IHCHMX=1
              IACTIVE(NMD)=1
              IACTALL=IACTALL+1
C              IF(ISCDRY(LCHNV).EQ.0) IACTIVE=1
            ENDIF
            ENDIF
            IF(SRFHOST.GT.SRFCHAN)THEN
            IF(H1P(LHOST).GT.HDRY)THEN
              HCHNMN=H1P(LHOST)
              IHCHMN=1
              IACTIVE(NMD)=1
              IACTALL=IACTALL+1
C              IF(ISCDRY(LHOST).EQ.0) IACTIVE=1
            ENDIF
            ENDIF
            IF(HP(LHOST).LE.0.0.OR.HP(LCHNU).LE.0.0)THEN
              IF(IACTIVE(NMD).EQ.1)THEN
                IACTALL=IACTALL-1
                IACTIVE(NMD)=0  !PMC
              ENDIF
            ENDIF
            IF(IACTIVE(NMD).EQ.1)THEN
              WCHAN=DYP(LCHNV)
C              RLCHN=0.5*DXP(LCHNV)+0.25*DXP(LHOST)
C              HCHAN=0.5*DXP(LCHNV)*HP(LCHNV)+0.25*DXP(LHOST)*HP(LHOST)
              RLCHN=0.5*DXP(LCHNV)+CHANLEN(NMD)
              HCHAN=0.5*DXP(LCHNV)*H1P(LCHNV)+CHANLEN(NMD)*H1P(LHOST)
              HCHAN=HCHAN/RLCHN
              IF(IHCHMX.EQ.1) HCHAN=MAX(HCHAN,HCHNMX)
              IF(IHCHMN.EQ.1) HCHAN=MIN(HCHAN,HCHNMN)
              IF(HCHAN.GT.0.)THEN
C                TMPVAL=STBY(LCHNV)*DELT/(HCHAN*HCHAN*WCHAN)
                TMPVAL=CHANFRIC(NMD)*DELT/(HCHAN*HCHAN*WCHAN)
                CCCCHU(NMD)=1./(1.+TMPVAL*ABS(QCHANVT(NMD)))
C                HCHANT=HCHAN-HCHNCOR
C                HCHANT=MAX(HCHANT,0.0)
C                CCCCHV(NMD)=DELT*HCHANT*WCHAN/RLCHN
                CCCCHV(NMD)=DELT*HCHAN*WCHAN/RLCHN
              ENDIF
            ENDIF
          ENDIF
          CC(LHOST)=CC(LHOST)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
          CCCCHH(NMD)=G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
          IF(MDCHTYP(NMD).EQ.1)THEN
            CC(LCHNU)=CC(LCHNU)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
            TMPVAL=G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANUT(NMD)
     &     -G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNU))
            FP(LHOST)=FP(LHOST)+TMPVAL
            FP(LCHNU)=FP(LCHNU)-TMPVAL
          ENDIF
          IF(MDCHTYP(NMD).EQ.2)THEN
            CC(LCHNV)=CC(LCHNV)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
            TMPVAL=G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANVT(NMD)
     &     -G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNV))
            FP(LHOST)=FP(LHOST)+TMPVAL
            FP(LCHNV)=FP(LCHNV)-TMPVAL
          ENDIF
        ENDDO
C
        WRITE(8,1949)N,IACTALL
C
c        DO NMD=1,MDCHH
c	    WRITE(8,1948)NMD,CCCCHH(NMD),CCCCHU(NMD),CCCCHV(NMD)
c        ENDDO
C
      ENDIF
C
 1949 FORMAT(' N, # ACTIVE 2 GRID FLOWS = ',2I8)
 1948 FORMAT(I5,3E12.4)
C
C**********************************************************************C
C
      RETURN
      END
