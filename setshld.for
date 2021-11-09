C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SETSHLD(TSC,THETA,D,SSG,DSR,USC)
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
      TMP=GP/(VISC*VISC)
      DSR=D*(TMP**0.333333)
      GPD=GP*D
C
C **  SHIELDS 
C
      IF(DSR.LE.4.0)THEN
        THETA=0.24/DSR
      ENDIF
C
      IF(DSR.GT.4.0.AND.DSR.LE.10.0)THEN
        THETA=0.14/(DSR**0.64)
      ENDIF
C
      IF(DSR.GT.10.0.AND.DSR.LE.20.0)THEN
        THETA=0.04/(DSR**0.1)
      ENDIF
C
      IF(DSR.GT.20.0.AND.DSR.LE.150.0)THEN
        THETA=0.013*(DSR**0.29)
      ENDIF
C
      IF(DSR.GT.150.0)THEN
        THETA=0.055
      ENDIF
C
      TSC=GPD*THETA
      USC=SQRT(TSC)
C
      RETURN
      END
