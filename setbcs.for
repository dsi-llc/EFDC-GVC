C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SETBCS
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 01/15/2002        John Hamrick       01/11/2002       John Hamrick
C  modified boundary condition flags for type 2 open boundaries
C 05/02/2002        John Hamrick       05/01/2002       John Hamrick
C  added real flags RSSBCE(L),RSSBCW(L),RSSBCN(L),RSSBCS(L)
C  to modified calculation of cell center bed stress (stored as QQ(l,0))
C  and the outputed cell center velocity for cells have source/sinks
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
C **  SUBROUTINE SETBCS SETS BOUNDARY CONDITION SWITCHES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION SUBEW(LCM),SVBNS(LCM)
C
C**********************************************************************C
C
C **  SET LAND-WATER BOUNDARY SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      I=IL(L)
      J=JL(L)
      IF(LCT(L).EQ.1)THEN
       STCUV(L)=0.
       STCAP(L)=0.5
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.2)THEN
       STCUV(L)=0.
       STCAP(L)=0.5
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.3)THEN
       STCUV(L)=0.
       STCAP(L)=0.5
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.4)THEN
       STCUV(L)=0.
       STCAP(L)=0.5
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.5)THEN
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.6)THEN
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.6) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.7) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.6) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.7) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      IF(LCT(L).EQ.7)THEN
       IF(IJCT(I-1,J).EQ.1) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.2) SUB(L)=0.
       IF(IJCT(I-1,J).EQ.3) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.4) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.5) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.6) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.7) SUB(L)=1.
       IF(IJCT(I-1,J).EQ.9) SUB(L)=0.
       IF(IJCT(I,J-1).EQ.1) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.2) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.3) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.4) SVB(L)=0.
       IF(IJCT(I,J-1).EQ.5) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.6) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.7) SVB(L)=1.
       IF(IJCT(I,J-1).EQ.9) SVB(L)=0.
      ENDIF
      ENDDO
C
      SUB(1)=0.
      SVB(1)=0.
      SUB(LC)=0.
      SVB(LC)=0.
C
C**********************************************************************C
C
C **  MODIFY LAND-WATER BNDRY CONDS FOR PERIOD GRID IN N-S DIRECTION
C
      IF(ISPGNS.GE.1)THEN
      DO NPN=1,NPNSBP
C
C     SET SOUTH CELL SVB'S TO 1
C
      LS=LIJ(ISPNS(NPN),JSPNS(NPN))
      SVB(LS)=1.
      SVBO(LS)=1.
C
      ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  SET WATER-WATER (P OR SURFACE ELEVATION) BOUNDARY SWITCHES  
C     
C----------------------------------------------------------------------C
C
      DO LL=1,NPBW
      I=IPBW(LL)
      J=JPBW(LL)
      L=LIJ(I,J)
      LPBW(LL)=L      
      SPB(L)=0.
      SUB(L)=0.
      SVB(L)=0.
      SWB(L)=0.
      IF(ISPBW(LL).LE.1) THEN
        SVB(L+1)=0.
        SWB(L+1)=0.
        SCAX(L+1)=0.
      END IF
      SAAX(L+1)=0.
      SDX(L+1)=0.
      ENDDO
C
      DO LL=1,NPBE
      I=IPBE(LL)
      J=JPBE(LL)
      L=LIJ(I,J)
      LPBE(LL)=L
      SPB(L)=0.
      SVB(L)=0.
      SWB(L)=0.
      IF(ISPBE(LL).LE.1) THEN
        SWB(L-1)=0.
        SVB(L-1)=0.
      SCAX(L)=0.
      END IF
      SAAX(L)=0.
      SDX(L)=0.
      ENDDO
C
      DO LL=1,NPBS
      I=IPBS(LL)
      J=JPBS(LL)
      L=LIJ(I,J)
      LPBS(LL)=L
      LN=LNC(L)
      SPB(L)=0.
      SVB(L)=0.
      SUB(L)=0.
      SWB(L)=0.
      IF(ISPBS(LL).LE.1) THEN
        SUB(LN)=0.
        SWB(LN)=0.
        SCAY(LN)=0.
      END IF
      SAAY(LN)=0.
      SDY(LN)=0.
      ENDDO     
C
      DO LL=1,NPBN
      I=IPBN(LL)
      J=JPBN(LL)
      L=LIJ(I,J)
      LPBN(LL)=L
      LS=LSC(L)
      SPB(L)=0.
      SUB(L)=0.
      SWB(L)=0.
      IF(ISPBN(LL).LE.1) THEN
        SUB(LS)=0.
        SWB(LS)=0.
        SCAY(L)=0.
      END IF
      SAAY(L)=0.
      SDY(L)=0.
      ENDDO
C
C
      DO L=1,LC
      SUBO(L)=SUB(L)
      SVBO(L)=SVB(L)
      SUB1(L)=SUB(L)
      SVB1(L)=SVB(L)
      ENDDO
C
C**********************************************************************C
C
C ** RESET DXU,DYU,DXV,DYV BASED ON BOUNDARY CONDITION SWITCHES
C
      DO L=2,LA
        IF(SUB(L).GT.0.5)THEN
          DXU(L)=0.5*(DXP(L)+DXP(L-1))
          DYU(L)=0.5*(DYP(L)+DYP(L-1))
        ENDIF
        IF(SUB(L).LT.0.5.AND.SUB(L+1).GT.0.5)THEN
          write(8,*)'DDYDDDX',L,IL(L),JL(L),DXP(L),DXP(L+1)
          DXU(L)=DXP(L)
          DDYDDDX=2.*(DYP(L+1)-DYP(L))/(DXP(L)+DXP(L+1))
          DYU(L)=DYP(L)-0.5*DXP(L)*DDYDDDX
        ENDIF
        IF(SUB(L).LT.0.5.AND.SUB(L+1).LT.0.5)THEN
          DXU(L)=DXP(L)
          DYU(L)=DYP(L)
        ENDIF
      ENDDO
C      
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        IF(SVB(L).GT.0.5)THEN
          DXV(L)=0.5*(DXP(L)+DXP(LS))
          DYV(L)=0.5*(DYP(L)+DYP(LS))
        ENDIF
        IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
          write(8,*)'DDXDDDY',L,IL(L),JL(L),DXP(L),DXP(LN)
          DDXDDDY=2.*(DXP(LN)-DXP(L))/(DYP(L)+DYP(LN))
          DXV(L)=DXP(L)-0.5*DYP(L)*DDXDDDY
          DYV(L)=DYP(L)
        ENDIF
        IF(SVB(L).LT.0.5.AND.SVB(LN).LT.0.5)THEN
          DXV(L)=DXP(L)
          DYV(L)=DYP(L)
        ENDIF
      ENDDO
C
C**********************************************************************C
C
C **  SET THIN BARRIERS BY CALLING CELLMASK
C **  CALL MOVED FROM AAEFDC ON 23 JAN 2004
C
      IF(ISMASK.EQ.1) CALL CELLMASK	 
C
C**********************************************************************C
C     
C **  SET VOLUMETRIC & CONCENTRATION SOURCE LOCATIONS AND BED STRESS
C **  AND CELL CENTER BED STRESS AND VELOCITY MODIFERS 
C
C----------------------------------------------------------------------C
C
      DO LL=1,NQSIJ
      I=IQS(LL)
      J=JQS(LL)
      LTMP=LIJ(I,J)
      LQS(LL)=LTMP
      IF(NQSMUL(LL).EQ.0)RQSMUL(LL)=1.
      IF(NQSMUL(LL).EQ.1)RQSMUL(LL)=DYP(LTMP)
      IF(NQSMUL(LL).EQ.2)RQSMUL(LL)=DXP(LTMP)
      IF(NQSMUL(LL).EQ.3)RQSMUL(LL)=DXP(LTMP)+DYP(LTMP)
      ENDDO
C
      DO NCTL=1,NQCTL
      RQDW=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LTMP=LIJ(IU,JU)
      IF(NQCMUL(NCTL).EQ.0)RQCMUL(NCTL)=1.
      IF(NQCMUL(NCTL).EQ.1)RQCMUL(NCTL)=DYP(LTMP)
      IF(NQCMUL(NCTL).EQ.2)RQCMUL(NCTL)=DXP(LTMP)
      IF(NQCMUL(NCTL).EQ.3)RQCMUL(NCTL)=DXP(LTMP)+DYP(LTMP)
      ENDDO
C
      DO L=2,LA
        RSSBCE(L)=1.0
	  RSSBCW(L)=1.0
	  RSSBCN(L)=1.0
	  RSSBCS(L)=1.0
	  SUBEW(L)=SUB(L)+SUB(L+1)
	  SVBNS(L)=SVB(L)+SVB(LNC(L))
	ENDDO
C
      DO LL=1,NQSIJ
        L=LQS(LL)
	  LE=L+1
	  LN=LNC(L)
	  IF(SVBNS(L).LT.0.5)THEN
	    IF(SUB(L).LT.0.5.AND.SUB(LE).GT.0.5)THEN
	      RSSBCW(L)=0.0
            RSSBCE(L)=2.0
          ENDIF
	    IF(SUB(L).GT.0.5.AND.SUB(LE).LT.0.5)THEN
	      RSSBCW(L)=2.0
            RSSBCE(L)=0.0
          ENDIF
        ENDIF
	  IF(SUBEW(L).LT.0.5)THEN
	    IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
	      RSSBCS(L)=0.0
            RSSBCN(L)=2.0
          ENDIF
	    IF(SVB(L).GT.0.5.AND.SVB(LN).LT.0.5)THEN
	      RSSBCS(L)=2.0
            RSSBCN(L)=0.0
          ENDIF
        ENDIF
      ENDDO
C
      DO NCTL=1,NQCTL
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        L=LIJ(IU,JU)
	  LE=L+1
	  LN=LNC(L)
	  IF(SVBNS(L).LT.0.5)THEN
	    IF(SUB(L).LT.0.5.AND.SUB(LE).GT.0.5)THEN
	      RSSBCW(L)=0.0
            RSSBCE(L)=2.0
          ENDIF
	    IF(SUB(L).GT.0.5.AND.SUB(LE).LT.0.5)THEN
	      RSSBCW(L)=2.0
            RSSBCE(L)=0.0
          ENDIF
        ENDIF
	  IF(SUBEW(L).LT.0.5)THEN
	    IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
	      RSSBCS(L)=0.0
            RSSBCN(L)=2.0
          ENDIF
	    IF(SVB(L).GT.0.5.AND.SVB(LN).LT.0.5)THEN
	      RSSBCS(L)=2.0
            RSSBCN(L)=0.0
          ENDIF
        ENDIF
      ENDDO
C
      DO NCTL=1,NQCTL
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF(ID.NE.0.AND.JD.NE.0)THEN
          L=LIJ(ID,JD)
	    LE=L+1
	    LN=LNC(L)
	    IF(SVBNS(L).LT.0.5)THEN
	      IF(SUB(L).LT.0.5.AND.SUB(LE).GT.0.5)THEN
	        RSSBCW(L)=0.0
              RSSBCE(L)=2.0
            ENDIF
	      IF(SUB(L).GT.0.5.AND.SUB(LE).LT.0.5)THEN
	        RSSBCW(L)=2.0
              RSSBCE(L)=0.0
            ENDIF
          ENDIF
	    IF(SUBEW(L).LT.0.5)THEN
	      IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
	        RSSBCS(L)=0.0
              RSSBCN(L)=2.0
            ENDIF
	      IF(SVB(L).GT.0.5.AND.SVB(LN).LT.0.5)THEN
	        RSSBCS(L)=2.0
              RSSBCN(L)=0.0
            ENDIF
          ENDIF
	  END IF
      ENDDO
C
      DO NWR=1,NQWR
        IU=IQWRU(NWR)
        JU=JQWRU(NWR)
        L=LIJ(IU,JU)
	  LE=L+1
	  LN=LNC(L)
	  IF(SVBNS(L).LT.0.5)THEN
	    IF(SUB(L).LT.0.5.AND.SUB(LE).GT.0.5)THEN
	      RSSBCW(L)=0.0
            RSSBCE(L)=2.0
          ENDIF
	    IF(SUB(L).GT.0.5.AND.SUB(LE).LT.0.5)THEN
	      RSSBCW(L)=2.0
            RSSBCE(L)=0.0
          ENDIF
        ENDIF
	  IF(SUBEW(L).LT.0.5)THEN
	    IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
	      RSSBCS(L)=0.0
            RSSBCN(L)=2.0
          ENDIF
	    IF(SVB(L).GT.0.5.AND.SVB(LN).LT.0.5)THEN
	      RSSBCS(L)=2.0
            RSSBCN(L)=0.0
          ENDIF
        ENDIF
      ENDDO
C
      DO NWR=1,NQWR
        ID=IQWRD(NWR)
        JD=JQWRD(NWR)
        L=LIJ(ID,JD)
	  LE=L+1
	  LN=LNC(L)
	  IF(SVBNS(L).LT.0.5)THEN
	    IF(SUB(L).LT.0.5.AND.SUB(LE).GT.0.5)THEN
	      RSSBCW(L)=0.0
            RSSBCE(L)=2.0
          ENDIF
	    IF(SUB(L).GT.0.5.AND.SUB(LE).LT.0.5)THEN
	      RSSBCW(L)=2.0
            RSSBCE(L)=0.0
          ENDIF
        ENDIF
	  IF(SUBEW(L).LT.0.5)THEN
	    IF(SVB(L).LT.0.5.AND.SVB(LN).GT.0.5)THEN
	      RSSBCS(L)=0.0
            RSSBCN(L)=2.0
          ENDIF
	    IF(SVB(L).GT.0.5.AND.SVB(LN).LT.0.5)THEN
	      RSSBCS(L)=2.0
            RSSBCN(L)=0.0
          ENDIF
        ENDIF
      ENDDO
C
C**********************************************************************C
C     
C **  SET OPEN BOUNDARY FLAGS
C
C----------------------------------------------------------------------C
C
      DO LL=1,NCBS
      I=ICBS(LL)
      J=JCBS(LL)
      LCBS(LL)=LIJ(I,J)
      L=LIJ(I,J)
      SCB(L)=0.
      ENDDO
C
      DO LL=1,NCBW
      I=ICBW(LL)
      J=JCBW(LL)
      LCBW(LL)=LIJ(I,J)
      L=LIJ(I,J)
      SCB(L)=0.
      ENDDO
C
      DO LL=1,NCBE
      I=ICBE(LL)
      J=JCBE(LL)
      LCBE(LL)=LIJ(I,J)
      L=LIJ(I,J)
      SCB(L)=0.
      ENDDO
C
      DO LL=1,NCBN
      I=ICBN(LL)
      J=JCBN(LL)
      LCBN(LL)=LIJ(I,J)
      L=LIJ(I,J)
      SCB(L)=0.
      ENDDO
C
C**********************************************************************C
C
C **  SET CHANNEL HOST AND GUEST LOCATION MAPPINGS
C
      IF(MDCHH.GE.1)THEN
        DO NMD=1,MDCHH
        LMDCHH(NMD)=LIJ(IMDCHH(NMD),JMDCHH(NMD))
        IF(IMDCHU(NMD).EQ.1.AND.JMDCHU(NMD).EQ.1)THEN
          LMDCHU(NMD)=1
         ELSE
          LMDCHU(NMD)=LIJ(IMDCHU(NMD),JMDCHU(NMD))
        ENDIF
        IF(IMDCHV(NMD).EQ.1.AND.JMDCHV(NMD).EQ.1)THEN
          LMDCHV(NMD)=1
         ELSE
          LMDCHV(NMD)=LIJ(IMDCHV(NMD),JMDCHV(NMD))
        ENDIF
        ENDDO
      ENDIF
C
C**********************************************************************C
C
C **  SET CELL FACE WET DEPTHS
C
      HUWET(1)=HWET
      HUWET(LC)=HWET
      HVWET(1)=HWET
      HVWET(LC)=HWET
      HUDRY(1)=HDRY
      HUDRY(LC)=HDRY
      HVDRY(1)=HDRY
      HVDRY(LC)=HDRY
C
      DO L=2,LA
      LS=LSC(L)
      HUDRY(L)=HDRY+0.5*ABS(BELV(L)-BELV(L-1))
      HVDRY(L)=HDRY+0.5*ABS(BELV(L)-BELV(LS))
      HUWET(L)=HWET+0.5*ABS(BELV(L)-BELV(L-1))
      HVWET(L)=HWET+0.5*ABS(BELV(L)-BELV(LS))
      ENDDO
C
      IF(ISDRY.GT.0)THEN
      NDRYTMP=MOD(ISDRY,2)
      IF(NDRYTMP.NE.0)THEN
        DO L=2,LA
        HUWET(L)=HWET
        HVWET(L)=HWET
        HUDRY(L)=HDRY
        HVDRY(L)=HDRY
        ENDDO
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  DIAGNOSTIC OUTPUT
C
      OPEN(1,FILE='SETBC.DIA')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='SETBC.DIA')
C
      DO L=2,LA
       WRITE(1,1001)IL(L),JL(L),SUB(L),SUB(L+1),SVB(L),SVB(LNC(L)),
     &              SPB(L)
      ENDDO
C
      CLOSE(1)
C
 1001 FORMAT(2I5,8E13.4)
C
C**********************************************************************C
C
      RETURN
      END
