C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FDSTRSE(VOID,BMECH1,BMECH2,BMECH3)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 11/15/2001        john hamrick       11/15/2001       john hamrick
C  added standard exponential form constitutive relationship
C
C----------------------------------------------------------------------C
C
      INCLUDE 'EFDC.PAR'
C
C **  FDSTRSE IS COMPRESSION LENGTH SCALE 
C
C     =-D(EFFECTIVE STRESS/WATER SPECIFIC WEIGHT)/D(VOID RATIO)
C
C     IE DERIVATIVE OF SPECIFIC WEIGHT NORMALIZED EFFECTIVE 
C        STRESS WITH RESPECT TO VOID RATIO
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
        TMP=-(VOID-BMECH2)/BMECH3
        FDSTRSE=(BMECH1/BMECH3)*EXP(TMP)
C
      ELSE
C
c        FDSTRSEL=-0.00835084*(VOID**3)+0.207061*(VOID**2)
c     &          -2.57907*VOID+ 8.79215
        FSTRSEL=-0.0147351*(VOID**3)+0.311854*(VOID**2)
     &          -2.96371*VOID+7.34698
        DFSTRSEL=-0.0442053*(VOID**2)+0.623708*VOID
     &          -2.96371
        FDSTRSE=DFSTRSEL*EXP(FSTRSEL)
C
      END IF
C
      RETURN
      END
