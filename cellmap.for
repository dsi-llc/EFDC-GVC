C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CELLMAP
C
C **  SUBROUTINE CELLMAP GENERATES CELL MAPPINGS 
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
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  SET 1D CELL INDEX SEQUENCE AND MAPPINGS
C
C----------------------------------------------------------------------C
C
cSO   OPEN(1,FILE='CELL9.OUT',STATUS='UNKNOWN')
cSO   WRITE (1,1616)IC,JC
C
      L=2
      DO J=1,JC
      DO I=1,IC
       IF(IJCT(I,J).GT.0.AND.IJCT(I,J).LT.9)THEN
        IL(L)=I
        JL(L)=J
        LCT(L)=IJCT(I,J)
        LIJ(I,J)=L
        L=L+1
       ENDIF
cSO   IF(IJCT(I,J).EQ.9) WRITE(1,1616)I,J
      ENDDO
      ENDDO
      LA=L-1
      LCTT=L
      IF(LCTT.NE.LC)THEN
        WRITE(6,1617)LCTT,LC
        WRITE(7,1617)LCTT,LC
        WRITE(8,1617)LCTT,LC
        STOP
      ENDIF
      IL(1)=0
      IL(LC)=0
      JL(1)=0
      JL(LC)=0
cSO   WRITE(1,601)LA
C
C      WRITE(6,601)LA
      WRITE(7,601)LA
      WRITE(8,601)LA
C
cSO   CLOSE(1)
  601 FORMAT('  LA=',I10,//)
 1616 FORMAT(2I10)
 1617 FORMAT(' LC =',I5,' DETERMINED IN CELLMAP',
     &        ' INCONSITENT WITH INPUT VALUE =',I5//)
C
C----------------------------------------------------------------------C
C
      IF(ISWASP.LE.5)THEN 
      L=2
      DO J=1,JC
      DO I=1,IC
       IF(IJCTLT(I,J).GT.0.AND.IJCTLT(I,J).LT.9)THEN
        ILLT(L)=I
        JLLT(L)=J
        LCTLT(L)=IJCTLT(I,J)
        LIJLT(I,J)=L
        L=L+1
       ENDIF
      ENDDO
      ENDDO
      LALT=L-1
      LCLT=L
      ENDIF
C
      IF(ISWASP.GE.6)THEN 
      L=2
      DO J=1,JC
      DO I=1,IC
       IF(IJCTLT(I,J).GT.0.AND.IJCTLT(I,J).LT.7)THEN
        ILLT(L)=I
        JLLT(L)=J
        LCTLT(L)=IJCTLT(I,J)
        LIJLT(I,J)=L
        L=L+1
       ENDIF
      ENDDO
      ENDDO
      LALT=L-1
      LCLT=L
      ENDIF
C
      WRITE(7,1616)LALT,LCLT
      WRITE(8,1616)LALT,LCLT
C
C----------------------------------------------------------------------C
C
C **  ASSIGN RED AND BLACK CELL SEQUENCES
C
      IF(IRVEC.NE.9)THEN
C
      LR=1
      LB=1
C
      DO L=2,LA
      K=IL(L)+JL(L)
       IF(MOD(K,2) .EQ. 0)THEN
        LRC(LR)=L
        LLRC(L)=LR
        ISRED(L)=1
        LR=LR+1
       ELSE
        LBC(LB)=L
        LLBC(L)=LB
        ISRED(L)=0
        LB=LB+1
       ENDIF
      ENDDO
C
      NRC=LR-1
      NBC=LB-1
C
      ENDIF
C 
C----------------------------------------------------------------------C
C
C **  SET NORTH AND SOUTH CELL IDENTIFIER ARRAYS
C
      LNC(1)=LC
      LSC(1)=LC
      LNEC(1)=LC
      LNWC(1)=LC
      LSEC(1)=LC
      LSWC(1)=LC
      LNC(LC)=1
      LSC(LC)=1
      LNEC(LC)=1
      LNWC(LC)=1
      LSEC(LC)=1
      LSWC(LC)=1
C
      DO L=2,LA
      I=IL(L)
      J=JL(L)
C
      IF(IJCT(I,J+1).EQ.9)THEN
       LNC(L)=LC
      ELSE
       LNC(L)=LIJ(I,J+1)
	 LNC(L)=MAX(LNC(L),1)
      ENDIF
C
      IF(IJCT(I,J-1).EQ.9)THEN
       LSC(L)=LC
      ELSE
       LSC(L)=LIJ(I,J-1)
	 LSC(L)=MAX(LSC(L),1)
      ENDIF
C
      IF(IJCT(I+1,J+1).EQ.9)THEN
       LNEC(L)=LC
      ELSE
       LNEC(L)=LIJ(I+1,J+1)
	 LNEC(L)=MAX(LNEC(L),1)
      ENDIF
C
      IF(IJCT(I-1,J+1).EQ.9)THEN
       LNWC(L)=LC
      ELSE
       LNWC(L)=LIJ(I-1,J+1)
	 LNWC(L)=MAX(LNWC(L),1)
      ENDIF
C
      IF(IJCT(I+1,J-1).EQ.9)THEN
       LSEC(L)=LC
      ELSE
       LSEC(L)=LIJ(I+1,J-1)
	 LSEC(L)=MAX(LSEC(L),1)
      ENDIF
C
      IF(IJCT(I-1,J-1).EQ.9)THEN
       LSWC(L)=LC
      ELSE
       LSWC(L)=LIJ(I-1,J-1)
	 LSWC(L)=MAX(LSWC(L),1)
      ENDIF
C

      ENDDO
C
C----------------------------------------------------------------------C
C
C **  MODIFY NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
C
      IF(ISPGNS.GE.1)THEN
      DO NPN=1,NPNSBP
C
C     SET NORTH CELL SOUTH OF SOUTH CELL
C
      LS=LIJ(ISPNS(NPN),JSPNS(NPN))
      LSC(LS)=LIJ(INPNS(NPN),JNPNS(NPN))
      IF( IJCT(INPNS(NPN)+1,JNPNS(NPN)).EQ.9)THEN
        LSEC(LS)=LC
       ELSE
        LSEC(LS)=LIJ(INPNS(NPN)+1,JNPNS(NPN))
      ENDIF
      IF( IJCT(INPNS(NPN)-1,JNPNS(NPN)).EQ.9)THEN
        LSWC(LS)=LC
       ELSE
        LSWC(LS)=LIJ(INPNS(NPN)-1,JNPNS(NPN)) 
      ENDIF
C
C     SET SOUTH CELL NORTH OF NORTH CELL
C
      LN=LIJ(INPNS(NPN),JNPNS(NPN))
      LNC(LN)=LIJ(ISPNS(NPN),JSPNS(NPN))
      IF( IJCT(ISPNS(NPN)+1,JSPNS(NPN)).EQ.9)THEN
        LNEC(LN)=LC
       ELSE
        LNEC(LN)=LIJ(ISPNS(NPN)+1,JSPNS(NPN))
      ENDIF
      IF( IJCT(ISPNS(NPN)-1,JSPNS(NPN)).EQ.9)THEN
        LNWC(LN)=LC
       ELSE
        LNWC(LN)=LIJ(ISPNS(NPN)-1,JSPNS(NPN)) 
      ENDIF
C
      ENDDO
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET LT NORTH AND SOUTH CELL IDENTIFIER ARRAYS
C
      LNCLT(1)=LCLT
      LSCLT(1)=LCLT
      LNCLT(LC)=1
      LSCLT(LC)=1
C
      DO L=2,LALT
      I=ILLT(L)
      J=JLLT(L)
C
      IF(IJCTLT(I,J+1).EQ.9)THEN
       LNCLT(L)=LCLT
      ELSE
       LNCLT(L)=LIJLT(I,J+1)
      ENDIF
C
      IF(IJCTLT(I,J-1).EQ.9)THEN
       LSCLT(L)=LCLT
      ELSE
       LSCLT(L)=LIJLT(I,J-1)
      ENDIF
C
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  SET NORTH, SOUTH, EAST AND WEST RED CELL IDENTIFIER ARRAYS
C
      IF(IRVEC.NE.9)THEN
C
      DO LR=1,NRC
      L=LRC(LR)
C
      LN=LNC(L)
      IF(LN.EQ.LC)THEN
       LBNRC(LR)=NBC+1
      ELSE
       LBNRC(LR)=LLBC(LN)
      ENDIF
C
      LS=LSC(L)
      IF(LS.EQ.LC)THEN
       LBSRC(LR)=NBC+1
      ELSE
       LBSRC(LR)=LLBC(LS)
      ENDIF
C
      LE=L+1
      IF(LE.EQ.LC)THEN
       LBERC(LR)=NBC+1
      ELSE
       IF(ISRED(LE).EQ.1)THEN
        LBERC(LR)=NBC+1
       ELSE
        LBERC(LR)=LLBC(LE)
       ENDIF
      ENDIF
C
      LW=L-1
      IF(LW.EQ.1)THEN
       LBWRC(LR)=NBC+1
      ELSE
       IF(ISRED(LW).EQ.1)THEN
        LBWRC(LR)=NBC+1
       ELSE
        LBWRC(LR)=LLBC(LW)
       ENDIF
      ENDIF
C
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  SET NORTH, SOUTH, EAST AND WEST BLACK CELL IDENTIFIER ARRAYS
C
      IF(IRVEC.NE.9)THEN
C
      DO LB=1,NBC
      L=LBC(LB)
C
      LN=LNC(L)
      IF(LN.EQ.LC)THEN
       LRNBC(LB)=NRC+1
      ELSE
       LRNBC(LB)=LLRC(LN)
      ENDIF
C
      LS=LSC(L)
      IF(LS.EQ.LC)THEN
       LRSBC(LB)=NRC+1
      ELSE
       LRSBC(LB)=LLRC(LS)
      ENDIF
C
      LE=L+1
      IF(LE.EQ.LC)THEN
       LREBC(LB)=NRC+1
      ELSE
       IF(ISRED(LE).EQ.0)THEN
        LREBC(LB)=NRC+1
       ELSE
        LREBC(LB)=LLRC(LE)
       ENDIF
      ENDIF
C
      LW=L-1
      IF(LW.EQ.1)THEN
       LRWBC(LB)=NRC+1
      ELSE
       IF(ISRED(LW).EQ.0)THEN
        LRWBC(LB)=NRC+1
       ELSE
        LRWBC(LB)=LLRC(LW)
       ENDIF
      ENDIF
C
      ENDDO
C
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  DIAGNOSE OF RED-BLACK CELL MAPPING
C
      IF(IRVEC.EQ.2.OR.IRVEC.EQ.9) GOTO 220
C
      OPEN(1,FILE='RBCM.DIA',STATUS='UNKNOWN')
C
C **  RED CELL LOOP
C
      DO LR=1,NRC
      LTMP=LRC(LR)
      LBNRTMP=LBNRC(LR)
      LBSRTMP=LBSRC(LR)
      LBERTMP=LBERC(LR)
      LBWRTMP=LBWRC(LR)
      LNTMP=LBC(LBNRTMP)
      LSTMP=LBC(LBSRTMP)
      LETMP=LBC(LBERTMP)
      LWTMP=LBC(LBWRTMP)
      IF(LTMP.EQ.1) WRITE(1,101)LR,LTMP
      IF(LTMP.EQ.LC)WRITE(1,102)LR,LTMP
      IF(LTMP.GT.1.AND.LTMP.LT.LC)THEN
        ITMP=IL(LTMP)
        JTMP=JL(LTMP)
        IF(LNTMP.EQ.1) WRITE(1,103)LR,LTMP,ITMP,JTMP
        IF(LNTMP.EQ.LC)WRITE(1,104)LR,LTMP,ITMP,JTMP
        IF(LNTMP.GT.1.AND.LNTMP.LT.LC)THEN
          INTMP=IL(LNTMP)
          JNTMP=JL(LNTMP)
          JNTMPM1=JNTMP-1
          IF(ITMP.NE.INTMP)  WRITE(1,105)LR,LTMP,ITMP,JTMP,INTMP,JNTMP
          IF(JTMP.NE.JNTMPM1)WRITE(1,106)LR,LTMP,ITMP,JTMP,INTMP,JNTMP
        ENDIF
        IF(LSTMP.EQ.1) WRITE(1,107)LR,LTMP,ITMP,JTMP
        IF(LSTMP.EQ.LC)WRITE(1,108)LR,LTMP,ITMP,JTMP
        IF(LSTMP.GT.1.AND.LSTMP.LT.LC)THEN
          ISTMP=IL(LSTMP)
          JSTMP=JL(LSTMP)
          JSTMPP1=JSTMP+1
          IF(ITMP.NE.ISTMP)  WRITE(1,109)LR,LTMP,ITMP,JTMP,ISTMP,JSTMP
          IF(JTMP.NE.JSTMPP1)WRITE(1,110)LR,LTMP,ITMP,JTMP,ISTMP,JSTMP
        ENDIF
        IF(LETMP.EQ.1) WRITE(1,111)LR,LTMP,ITMP,JTMP
        IF(LETMP.EQ.LC)WRITE(1,112)LR,LTMP,ITMP,JTMP
        IF(LETMP.GT.1.AND.LETMP.LT.LC)THEN
          IETMP=IL(LETMP)
          IETMPM1=IETMP-1
          JETMP=JL(LETMP)
          IF(ITMP.NE.IETMPM1.AND.SUB(LETMP).NE.0.0)THEN
            WRITE(1,113)LR,LTMP,ITMP,JTMP,IETMP,JETMP
            WRITE(1,119)SUB(LETMP)
          ENDIF
          IF(JTMP.NE.JETMP.AND.SUB(LETMP).NE.0.0)   THEN
            WRITE(1,114)LR,LTMP,ITMP,JTMP,IETMP,JETMP
            WRITE(1,119)SUB(LETMP)
          ENDIF
        ENDIF
        IF(LWTMP.EQ.1) WRITE(1,115)LR,LTMP,ITMP,JTMP
        IF(LWTMP.EQ.LC)WRITE(1,116)LR,LTMP,ITMP,JTMP
        IF(LWTMP.GT.1.AND.LWTMP.LT.LC)THEN
          IWTMP=IL(LWTMP)
          IWTMPP1=IWTMP+1
          JWTMP=JL(LWTMP)
          IF(ITMP.NE.IWTMPP1.AND.SUB(LTMP).NE.0.0)THEN
            WRITE(1,117)LR,LTMP,ITMP,JTMP,IWTMP,JWTMP
            WRITE(1,120)SUB(LTMP)
          ENDIF
          IF(JTMP.NE.JWTMP.AND.SUB(LTMP).NE.0.0)   THEN
            WRITE(1,118)LR,LTMP,ITMP,JTMP,IWTMP,JWTMP
            WRITE(1,120)SUB(LTMP)
          ENDIF
        ENDIF
      ENDIF
      ENDDO  
C
C **  BLACK CELL LOOP
C
      DO LB=1,NBC
      LTMP=LBC(LB)
      LRNBTMP=LRNBC(LB)
      LRSBTMP=LRSBC(LB)
      LREBTMP=LREBC(LB)
      LRWBTMP=LRWBC(LB)
      LNTMP=LRC(LRNBTMP)
      LSTMP=LRC(LRSBTMP)
      LETMP=LRC(LREBTMP)
      LWTMP=LRC(LRWBTMP)
      IF(LTMP.EQ.1) WRITE(1,101)LR,LTMP
      IF(LTMP.EQ.LC)WRITE(1,102)LR,LTMP
      IF(LTMP.GT.1.AND.LTMP.LT.LC)THEN
        ITMP=IL(LTMP)
        JTMP=JL(LTMP)
        IF(LNTMP.EQ.1) WRITE(1,103)LB,LTMP,ITMP,JTMP
        IF(LNTMP.EQ.LC)WRITE(1,104)LB,LTMP,ITMP,JTMP
        IF(LNTMP.GT.1.AND.LNTMP.LT.LC)THEN
          INTMP=IL(LNTMP)
          JNTMP=JL(LNTMP)
          JNTMPM1=JNTMP-1
          IF(ITMP.NE.INTMP)  WRITE(1,105)LB,LTMP,ITMP,JTMP,INTMP,JNTMP
          IF(JTMP.NE.JNTMPM1)WRITE(1,106)LB,LTMP,ITMP,JTMP,INTMP,JNTMP
        ENDIF
        IF(LSTMP.EQ.1) WRITE(1,107)LB,LTMP,ITMP,JTMP
        IF(LSTMP.EQ.LC)WRITE(1,108)LB,LTMP,ITMP,JTMP
        IF(LSTMP.GT.1.AND.LSTMP.LT.LC)THEN
          ISTMP=IL(LSTMP)
          JSTMP=JL(LSTMP)
          JSTMPP1=JSTMP+1
          IF(ITMP.NE.ISTMP)  WRITE(1,109)LB,LTMP,ITMP,JTMP,ISTMP,JSTMP
          IF(JTMP.NE.JSTMPP1)WRITE(1,110)LB,LTMP,ITMP,JTMP,ISTMP,JSTMP
        ENDIF
        IF(LETMP.EQ.1) WRITE(1,111)LB,LTMP,ITMP,JTMP
        IF(LETMP.EQ.LC)WRITE(1,112)LB,LTMP,ITMP,JTMP
        IF(LETMP.GT.1.AND.LETMP.LT.LC)THEN
          IETMP=IL(LETMP)
          IETMPM1=IETMP-1
          JETMP=JL(LETMP)
          IF(ITMP.NE.IETMPM1.AND.SUB(LETMP).NE.0.0)THEN
            WRITE(1,113)LB,LTMP,ITMP,JTMP,IETMP,JETMP
            WRITE(1,119)SUB(LETMP)
          ENDIF
          IF(JTMP.NE.JETMP.AND.SUB(LETMP).NE.0.0)   THEN
            WRITE(1,114)LB,LTMP,ITMP,JTMP,IETMP,JETMP
            WRITE(1,119)SUB(LETMP)
          ENDIF
        ENDIF
        IF(LWTMP.EQ.1) WRITE(1,115)LB,LTMP,ITMP,JTMP
        IF(LWTMP.EQ.LC)WRITE(1,116)LB,LTMP,ITMP,JTMP
        IF(LWTMP.GT.1.AND.LWTMP.LT.LC)THEN
          IWTMP=IL(LWTMP)
          IWTMPP1=IWTMP+1
          JWTMP=JL(LWTMP)
          IF(ITMP.NE.IWTMPP1.AND.SUB(LTMP).NE.0.0)THEN
            WRITE(1,117)LB,LTMP,ITMP,JTMP,IWTMP,JWTMP
            WRITE(1,120)SUB(LTMP)
          ENDIF
          IF(JTMP.NE.JWTMP.AND.SUB(LTMP).NE.0.0)   THEN
            WRITE(1,118)LB,LTMP,ITMP,JTMP,IWTMP,JWTMP
            WRITE(1,120)SUB(LTMP)
          ENDIF
        ENDIF
      ENDIF
      ENDDO
C
      CLOSE(1)
C
  220 CONTINUE
C
  101 FORMAT(' LR,LTMP = ',2I6/)
  102 FORMAT(' LR,LTMP = ',2I6/)
  103 FORMAT(' LN= 1, LR,LTMP,ITMP,JTMP = ',4I6/)
  104 FORMAT(' LN=LC, LR,LTMP,ITMP,JTMP = ',4I6/)
  105 FORMAT(' NERR, LR,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)
  106 FORMAT(' NERR, LR,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)
  107 FORMAT(' LS= 1, LR,LTMP,ITMP,JTMP = ',4I6/)
  108 FORMAT(' LS=LC, LR,LTMP,ITMP,JTMP = ',4I6/)
  109 FORMAT(' SERR, LR,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)
  110 FORMAT(' SERR, LR,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)
  111 FORMAT(' LE= 1, LR,LTMP,ITMP,JTMP = ',4I6/)
  112 FORMAT(' LE=LC, LR,LTMP,ITMP,JTMP = ',4I6/)
  113 FORMAT(' EERR, LR,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)
  114 FORMAT(' EERR, LR,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)
  115 FORMAT(' LW= 1, LR,LTMP,ITMP,JTMP = ',4I6/)
  116 FORMAT(' LW=LC, LR,LTMP,ITMP,JTMP = ',4I6/)
  117 FORMAT(' WERR, LR,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)
  118 FORMAT(' WERR, LR,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)
  201 FORMAT(' LB,LTMP = ',2I6/)
  202 FORMAT(' LB,LTMP = ',2I6/)
  203 FORMAT(' LN= 1, LB,LTMP,ITMP,JTMP = ',4I6/)
  204 FORMAT(' LN=LC, LB,LTMP,ITMP,JTMP = ',4I6/)
  205 FORMAT(' NERR, LB,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)
  206 FORMAT(' NERR, LB,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)
  207 FORMAT(' LS= 1, LB,LTMP,ITMP,JTMP = ',4I6/)
  208 FORMAT(' LS=LC, LB,LTMP,ITMP,JTMP = ',4I6/)
  209 FORMAT(' SERR, LB,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)
  210 FORMAT(' SERR, LB,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)
  211 FORMAT(' LE= 1, LB,LTMP,ITMP,JTMP = ',4I6/)
  212 FORMAT(' LE=LC, LB,LTMP,ITMP,JTMP = ',4I6/)
  213 FORMAT(' EERR, LB,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)
  214 FORMAT(' EERR, LB,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)
  215 FORMAT(' LW= 1, LB,LTMP,ITMP,JTMP = ',4I6/)
  216 FORMAT(' LW=LC, LB,LTMP,ITMP,JTMP = ',4I6/)
  217 FORMAT(' WERR, LB,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)
  218 FORMAT(' WERR, LB,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)
  119 FORMAT('       SUB(LETMP) = ',F10.2/)
  120 FORMAT('       SUB(LTMP)  = ',F10.2/)
C 
C**********************************************************************C
C
C **  DEFINE MAPPING TO 3D GRAPHICS GRID
C
      IF(ISCLO.EQ.0.OR.NWGG.EQ.0)THEN
        IG=IC
        JG=JC
       ELSE
        OPEN(1,FILE='GCELLMP.INP',STATUS='UNKNOWN')
        READ(1,1111)
        READ(1,1111)
        READ(1,1111)
        READ(1,1111)
        READ(1,*)IG,JG
        DO NW=1,NWGG
        READ(1,*)IGTMP,JGTMP,ICOMP,JCOMP
        LTMP=LIJ(ICOMP,JCOMP)
        IWGG(NW)=IGTMP
        JWGG(NW)=JGTMP
        LWGG(NW)=LTMP        
        ENDDO
        CLOSE(1)
      ENDIF
C
 1111 FORMAT(80X)
C 
C**********************************************************************C
C
      RETURN
      END
