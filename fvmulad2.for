C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE FVMULAD2(N,LC,A,B,C,D)
      DIMENSION  A(N),B(N),C(N),D(N)
C
C **  VECTOR MULTIPLY-ADD IN FORTRAN
C
      LM3=LC-3
      DO L=1,LM3,4
        A(L  )=B(L  )+C(L  )*D(L  )
        A(L+1)=B(L+1)+C(L+1)*D(L+1)
        A(L+2)=B(L+2)+C(L+2)*D(L+2)
        A(L+3)=B(L+3)+C(L+3)*D(L+3)
      END DO
C
      RETURN
      END
