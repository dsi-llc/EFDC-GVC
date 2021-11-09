C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION SETSTVEL(D,SSG)
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
C
C **  NONCOHEASIVE SEDIMENT SETTLING AND SHIELDS CRITERIA
C **  USING VAN RIJN'S EQUATIONS
C
      VISC=1.E-6
      GP=(SSG-1.)*9.82
      GPD=GP*D
      SQGPD=SQRT(GPD)
      RD=SQGPD*D/VISC
C
C **  SETTLING VELOCITY  
C
      IF(D.LT.1.0E-4)THEN
        WSET=SQGPD*RD/18.
      ENDIF
C
      IF(D.GE.1.0E-4.AND.D.LT.1.E-3)THEN
        TMP=SQRT(1.+0.01*RD*RD)-1.
        WSET=10.0*SQGPD*TMP/RD
      ENDIF
C
      IF(D.GE.1.E-3)THEN
        WSET=1.1*SQGPD
      ENDIF
C
      SETSTVEL=WSET
C
      RETURN
      END
