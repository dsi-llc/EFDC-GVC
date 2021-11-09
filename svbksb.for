C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
C
C **  FORM NUMERICAL RECIPES
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
C
      PARAMETER (NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
      S=0.
      IF(W(J).NE.0.)THEN
       DO 11 I=1,M
       S=S+U(I,J)*B(I)
   11  CONTINUE
       S=S/W(J)
      ENDIF
      TMP(J)=S
   12 CONTINUE
      DO 14 J=1,N
      S=0.
      DO 13 JJ=1,N
      S=S+V(J,JJ)*TMP(JJ)
   13 CONTINUE
      X(J)=S
   14 CONTINUE
      RETURN
      END
