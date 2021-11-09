C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SMLCFG(NOTUSED, XSOL, FVAL, GRAD, HAVEG)
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      LOGICAL HAVEG
      DIMENSION XSOL(LCM-2),GRAD(LCM-2),PVAL(LCM)
C
      HAVEG=.TRUE.
C
C  LOAD XSOL INTO PVAL
C
      PVAL(1)=0.
      PVAL(LC)=0.
C
      DO L=2,LA
        PVAL(L)=XSOL(L-1)
      ENDDO 
C
C  CALCULATE GRADIENT AND FUNCTION
C
      FVAL=0.
      DO L=2,LA
        LN=LNC(L)
        LS=LSC(L)
        GVAL=CCC(L)*PVAL(L)+CCS(L)*PVAL(LS)+CCW(L)*PVAL(L-1)
     &        +CCE(L)*PVAL(L+1)+CCN(L)*PVAL(LN)
        FVAL=FVAL+PVAL(L)*( 0.5*GVAL-FPTMP(L) )
        GRAD(L-1)=GVAL-FPTMP(L)
      ENDDO
C
      RETURN
      END
