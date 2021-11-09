C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE OUTPUT2
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
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  OUTPUT RESULTS OF RELAXATION SOLUTION
C
      WRITE (7,40) RP
   40 FORMAT (1H1,' RESULTS OF RELAX SOLUTION - RP=',F5.2,//)
C
      WRITE (7,41)
   41 FORMAT (' GLOBAL SQUARED ERROR',//)
C
      WRITE(7,43)ERRMAX,ERRMIN
   43 FORMAT('ERRMAX =',3X,E12.4,5X,'ERRMIN =',3X,E12.4)
C
C     DO N=1,NTS,10
C     WRITE (7,20) N, (ERR(NN),NN=N,N+9)
C     ENDDO
   20 FORMAT (1X,I5,3X,10E12.4)
C
      WRITE(7,40)RP
      WRITE (7,42)
   42 FORMAT (' ITERATIONS TO CONVERGENCE',//)
C
      WRITE(7,44)ITRMAX,ITRMIN
   44 FORMAT('ITRMAX =',I5,5X,'ITRMIN =',I5)
C
C     DO N=1,NTS,10
C     WRITE (7,21) N, (ITR(NN),NN=N,N+9)
C     ENDDO
   21 FORMAT (1X,I5,5X,10I10)
C
   30 FORMAT (10E12.4)
C
C**********************************************************************C
C
C **  OUTPUT HARMONIC ANALYSIS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      PAM(L)=(AMCP(L)*AMCP(L)+AMSP(L)*AMSP(L))**.5
      IF(AMSP(L).EQ.0.0.AND.AMCP(L).EQ.0.0)THEN
       PPH(L)=999999.
      ELSE
       PPH(L)=ATAN2(AMSP(L),AMCP(L))
      ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      PAM(L)=PAM(L)*GI
      ENDDO
      WRITE (7,55)
      CALL PPLOT (1)
C
   55 FORMAT (1H1,'TIDAL SURFACE DISPLACEMENT AMPLITUDE IN METERS',//)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      PAM(L)=0.5*TIDALP*PPH(L)/PI
      ENDDO
      WRITE(7,588)
      CALL PPLOT (1)
C
  588 FORMAT (1H1,'TIDAL SURFACE DISPLACEMENT PHASE IN SEC',//)
C
C----------------------------------------------------------------------C
C
C     DO KK=1,KC,KS
C
C     DO L=2,LA
C     PAM(L)=SQRT(AMCU(L,KK)*AMCU(L,KK)
C    &               +AMSU(L,KK)*AMSU(L,KK))
C     ENDDO
C     WRITE (7,56) KK
C     CALL PPLOT (2)
C
C     ENDDO
C
C  56 FORMAT (1H1,' X TIDAL VELOCITY AMPLITUDE, M/S, LAYER',I5,//)
C
C----------------------------------------------------------------------C
C
C     DO KK=1,KC,KS
C
C     DO L=2,LA
C     LN=LNC(L)
C     PAM(L)=SQRT(AMCV(L,KK)*AMCV(L,KK)
C    &               +AMSV(L,KK)*AMSV(L,KK))
C     ENDDO
C     WRITE (7,57) KK
C     CALL PPLOT (2)
C
C     ENDDO
C
C  57 FORMAT (1H1,' Y TIDAL VELOCITY AMPLITUDE, M/S, LAYER',I5,//)
C
C**********************************************************************C
C
C **  PRINTED OUTPUT OF P,U,AND V AMPLITUDES
C
C     WRITE(7,71)
C  71 FORMAT(1H1,'   L   IL   JL',9X,'PA',12X,'PP',13X,'U',14X,'V',//)
C
C     DO L=2,LA
C     WRITE(7,72)L,IL(L),JL(L),PAM(L),PPH(L),U(L),V(L)
C     ENDDO
   72 FORMAT(3I5,4(3X,E12.4))
C
C**********************************************************************C
C
C **  OUTPUT VECTOR POTENTIAL TRANSPORT VELOCITY
C
C     IF(ISVPTHA.NE.1) GOTO 100
C
C----------------------------------------------------------------------C
C
C     DO KK=1,KC,KS
C
C     DO L=2,LA
C     PAM(L)=0.5*(UVPT(L,KK)+UVPT(L+1,KK))/HMU(L)
C     ENDDO
C     WRITE (7,1458) KK
C     CALL PPLOT (2)
C
C     ENDDO
C
 1458 FORMAT(1H1,' X VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
C
C----------------------------------------------------------------------C
C
C     DO KK=1,KC,KS
C
C     DO L=2,LA
C     LN=LNC(L)
C     PAM(L)=0.5*(VVPT(L,KK)+VVPT(LN,KK))/HMV(L)
C     ENDDO
C     WRITE (7,1459) KK
C     CALL PPLOT (2)
C
C     ENDDO
C
 1459 FORMAT(1H1,' Y VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
C
C**********************************************************************C
C
  100 CONTINUE
C
      RETURN
      END
