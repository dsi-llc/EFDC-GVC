C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FUNDEN(SAL,SED,TEM)
C
C**********************************************************************C
C
C **  FUNDEN CALCULATED DENSITY AS A FUNTION OF SAL,TEM,AND SED
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
C      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'EFDC.PAR'
C
      SSG=2.5
C
      SDEN=1./2500000.
C
C **  DENSITY AT GIVEN VALUES OF SAL AND TEM
C
      SSTMP=SAL
      TTMP=TEM
      RHTMP=999.842594+6.793952E-2*TTMP-9.095290E-3*TTMP*TTMP
     &    +1.001685E-4*TTMP*TTMP*TTMP-1.120083E-6*TTMP*TTMP*TTMP*TTMP
     &    +6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
      RHO=RHTMP+SSTMP*(0.824493-4.0899E-3*TTMP+7.6438E-5*TTMP*TTMP
     &   -8.2467E-7*TTMP*TTMP*TTMP+5.3875E-9*TTMP*TTMP*TTMP*TTMP)
     &   +SQRT(SSTMP)*SSTMP*(-5.72466E-3+1.0227E-4*TTMP
     &   -1.6546E-6*TTMP*TTMP)+4.8314E-4*SSTMP*SSTMP
C
C **  CORRECTION FOR SEDIMENT
C **  REVISED 6/18/07 BY JMH
C
C       RHO=RHO*( (1.-SDEN*SED)+(SSG-1.)*SDEN*SED )
       RHO=RHO*( (1.-SDEN*SED)+SSG*SDEN*SED )
C
C **  RETURN DENSITY
C
      FUNDEN=RHO
C
C**********************************************************************C
C
      RETURN
      END
