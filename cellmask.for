C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CELLMASK
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
C **  SUBROUTINE CELLMASK CONVERTS LAND CELLS TO WATER CELLS BY 
C **  MASKING VARIABLES.  DEPTHS IN THE MASKED CELLS SHOULD BE INPUT AT
C **  THE END OF THE DXDY.INP FILE.
C 
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      OPEN(1,FILE='MASK.INP',STATUS='UNKNOWN')
C
      DO NS=1,6
      READ(1,1) 
      ENDDO
    1 FORMAT(80X)
      READ(1,*,ERR=1000) MMASK
C
C**********************************************************************C
C
      DO M=1,MMASK
      READ(1,*,ERR=1000) I,J,MTYPE
      L=LIJ(I,J)
	IF(IJCT(I,J).EQ.0.OR.IJCT(I,J).EQ.9)THEN
	  WRITE(6,*)'INVALID MASK LOCATION AT ',M,I,J
	ENDIF
	IF(L.LT.2.OR.L.GT.LA)THEN
	  WRITE(6,*)'INVALID MASK LOCATION AT ',M,I,J
	ENDIF
C
      IF(MTYPE.EQ.1)THEN
      SUB(L)=0.
      SUBO(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      DO K=1,KC
      U(L,K)=0.     
      U1(L,K)=0.     
      U2(L,K)=0.     
      UHDY(L,K)=0.     
      UHDY1(L,K)=0.     
      UHDY2(L,K)=0.     
      ENDDO
      ENDIF
C
      IF(MTYPE.EQ.2)THEN
      SVB(L)=0.
      SVBO(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      DO K=1,KC
      V(L,K)=0.     
      V1(L,K)=0.     
      V2(L,K)=0.     
      VHDX(L,K)=0.
      VHDX1(L,K)=0.
      VHDX2(L,K)=0.
      ENDDO
      ENDIF
C
      IF(MTYPE.EQ.3)THEN
      LN=LNC(L)
      SUB(L)=0.
      SUBO(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      SUB(L+1)=0.
      SUBO(L+1)=0.
      UHDYE(L+1)=0.
      UHDY1E(L+1)=0.
      UHDY2E(L+1)=0.
      SVB(L)=0.
      SVBO(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      SVB(LN)=0.
      SVBO(LN)=0.
      VHDXE(LN)=0.
      VHDX1E(LN)=0.
      VHDX2E(LN)=0.
      P(L)=0.
      P1(L)=0.
      DO K=1,KC
      B(L,K)=0.
      B1(L,K)=0.
      SAL(L,K)=0.
      SAL1(L,K)=0.
      TEM(L,K)=TEMO
      TEM1(L,K)=TEMO
      DYE(L,K)=0.
      DYE1(L,K)=0.
      SED(L,K,1)=0.
      SED1(L,K,1)=0.
      QQ(L,K)=0.
      QQ1(L,K)=0.
      QQL(L,K)=0.
      QQL1(L,K)=0.
      U(L,K)=0.     
      U1(L,K)=0.     
      U2(L,K)=0.     
      UHDY(L,K)=0.     
      UHDY1(L,K)=0.     
      UHDY2(L,K)=0.     
      U(L+1,K)=0.     
      U1(L+1,K)=0.     
      U2(L+1,K)=0.     
      UHDY(L+1,K)=0.     
      UHDY1(L+1,K)=0.     
      UHDY2(L+1,K)=0.     
      V(L,K)=0.     
      V1(L,K)=0.     
      V2(L,K)=0.     
      VHDX(L,K)=0.
      VHDX1(L,K)=0.
      VHDX2(L,K)=0.
      V(LN,K)=0.     
      V1(LN,K)=0.     
      V2(LN,K)=0.     
      VHDX(LN,K)=0.
      VHDX1(LN,K)=0.
      VHDX2(LN,K)=0.
      ENDDO
      ENDIF
C
      ENDDO
C
C**********************************************************************C
C
      CLOSE(1)
C
C**********************************************************************C
C
C **  WRITE READ ERRORS ON CELLMASK
C
      GOTO 1002
 1000 WRITE(6,1001)
 1001 FORMAT('  READ ERROR ON FILE MASK.INP ')
      STOP
 1002 CONTINUE
C
C**********************************************************************C
C
      RETURN
      END
