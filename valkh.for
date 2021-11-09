C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION VALKH(HFFDG)
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      IF(HFFDG.LE.0.02)THEN
        VALKH=HFFDG*HFFDG
        RETURN
      ENDIF
C
      IF(HFFDG.GE.10.)THEN
        VALKH=HFFDG
        RETURN
      ENDIF
C
      DO NTAB=2,1001
      FTMPM1=FUNKH(NTAB-1)
      FTMP  =FUNKH(NTAB  )
      IF(FTMPM1.LE.HFFDG.AND.HFFDG.LT.FTMP)THEN
        VALKH=RKHTAB(NTAB)
     &   -(RKHTAB(NTAB)-RKHTAB(NTAB-1))*(FTMP-HFFDG)/(FTMP-FTMPM1)
        RETURN
      ENDIF
      ENDDO
C
      IF(NTAB.EQ.1001)THEN
        WRITE(6,600) RKHTAB(1001)
        WRITE(8,600) RKHTAB(1001)
        STOP
      ENDIF
C
C **  INITIALIZE WAVE DISPERSION RELATION TABLE
C
C      DO N=1,501
C       RKH(N)=0.02*FLOAT(N-1)
C       FRKH(N)=RKH(N)*TANH(RKH(N))
C      ENDDO
C
C      WFRKH=WVFRQ*WVFRQ*DEP(N)/9.8
C      IF(WFRKH.LE.0.02) WRKH=WFRKH*WFRKH
C      IF(WFRKH.GE.10.) WRKH=WFRKH
C      IF(WFRKH.GT.0.02.AND.WFRKH.LT.10.)THEN
C        DO M=1,500
C         IF(WFRKH.GT.FRKH(M).AND.WFRKH.LT.FRKH(M+1))THEN
C           DRDF=(RKH(M+1)-RKH(M))/(FRKH(M+1)-FRKH(M))
C           WRKH=DRDF*(WFRKH-FRKH(M))+RKH(M)
C           GOTO 200 
C         ENDIF
C        ENDDO
C      ENDIF
C  200 CONTINUE
C
C
C
  600 FORMAT(' WAVE DISPERSION TABLE OUT OF BOUNDS KH = ',E12.4)
C 
      RETURN
      END
