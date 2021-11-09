C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
C      SUBROUTINE WSMTS(TINDAY)
      SUBROUTINE WSMTS
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
C
C**********************************************************************C
C
C WRITE TIME-SERIES OUTPUT
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      TINDAY=TINDAY
C
      OPEN(1,FILE='WQSDTS1.OUT',STATUS='UNKNOWN',POSITION='APPEND')
      OPEN(2,FILE='WQSDTS2.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN
        TIMTMP=TIMTMP/TCTMSR  
      ELSE
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF
C
      DO M=1,ISMTS
        LL=LSMTS(M)
        TSMPON = SMPON(LL,1)+SMPON(LL,2)+SMPON(LL,3)
        TSMPOP = SMPOP(LL,1)+SMPOP(LL,2)+SMPOP(LL,3)
        TSMPOC = SMPOC(LL,1)+SMPOC(LL,2)+SMPOC(LL,3)
        WRITE(1,71)IL(LL),JL(LL),TIMTMP,SM1NH4(LL),SM2NH4(LL),
     *   WQBFNH4(LL),SM1NO3(LL),SM2NO3(LL),WQBFNO3(LL),SM1PO4(LL),
     *   SM2PO4(LL),WQBFPO4D(LL),SM1H2S(LL),SM2H2S(LL),WQBFO2(LL),
     *   WQBFCOD(LL),SM1SI(LL),SM2SI(LL),WQBFSAD(LL),SMT(LL),SMBST(LL),
     *   TSMPON,TSMPOP,TSMPOC
        WRITE(2,71)IL(LL),JL(LL),TIMTMP,SMCSOD(LL),SMNSOD(LL),
     *   SMD1PO4(LL),SMD1SI(LL),SMSS(LL),SMJNIT(LL),SMJDEN(LL),
     *   SMJAQH2S(LL),SMJGCH4(LL),SMDGFN(LL),SMDGFP(LL),SMDGFC(LL),
     *   SMDFN(LL,1),SMDFN(LL,2),SMDFN(LL,3),SMDFP(LL,1),SMDFP(LL,2),
     *   SMDFP(LL,3),SMDFC(LL,1),SMDFC(LL,2),SMDFC(LL,3)
      ENDDO
C
CMRM   71 FORMAT(2I5, F11.5, 23E12.4)
   71 FORMAT(2I5, F11.5, 1P, 23E11.3)
C
      CLOSE(1)
      CLOSE(2)
C
      RETURN
      END
