C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FHYDCN(VOID,BMECH4,BMECH5,BMECH6,IBMECHK)
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
C
      INCLUDE 'EFDC.PAR'
C
C **  FHYDCN IS HYDRAULIC CONDUCTIVITY DIVIDED BY 1+VOID RATIO
C
C     UNITS ARE METERS/SECOND
C
      IF(BMECH4.GT.0.0)THEN
C
c      k=ko*exp((void-voidko)/voidkk)
c      BMECH4=ko
c      BMECH5=voidko
c      BMECH6=voidkk
C
        TMP=(VOID-BMECH5)/BMECH6
        IF(IBMECHK.EQ.0) THEN
          FHYDCN=BMECH4*EXP(TMP)/(1.+VOID)
        ELSE
          FHYDCN=BMECH4*EXP(TMP)
        END IF
C
      ELSE
C
        FHYDCNLOG=0.00816448*(VOID**3)-0.232453*(VOID**2)+2.5759*VOID
     &         -28.581
        FHYDCN=EXP(FHYDCNLOG)
C
      END IF
C
      RETURN
      END
