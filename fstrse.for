C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FSTRSE(VOID,BMECH1,BMECH2,BMECH3)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C 11/15/2001        john hamrick       11/15/2001       john hamrick
C  added standard exponential form constitutive relationship
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
C **  FSTRSE IS WATER SPECIFIC WEIGHT NORMALIZED EFFECTIVE STRESS 
C
C     UNITS ARE METERS
C
      IF(BMECH1.GT.0.0)THEN
C 
c      st=sto*exp(-(void-voidso)/voidss)
c      BMECH1=sto
c      BMECH2=voidso
c      BMECH3=voidss
C     
        TMP=-(VOID-BMECH2)/BMECH3
        FSTRSE=BMECH1*EXP(TMP)
C
      ELSE
C 
c        FSTRSELOG=-0.0147351*(VOID**3)+0.333957*(VOID**2)-3.28661*VOID
c     &          +8.90864
        FSTRSELOG=-0.0147351*(VOID**3)+0.311854*(VOID**2)-2.96371*VOID
     &          +7.34698
        FSTRSE=EXP(FSTRSELOG)
C
      END IF
C
      RETURN
      END
