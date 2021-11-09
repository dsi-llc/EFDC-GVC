C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE OUTPUT1
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
C **  PLOT SURFACE ELEVATION
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      PAM(L)=P(L)*GI
      ENDDO
      WRITE (7,50) N  
      CALL PPLOT (1)
C
   50 FORMAT(1H1,' SURFACE ELEVATION IN METERS AT TIMESTEP  ',I5,//)
C
C**********************************************************************C
C
C **  PLOT SURFACE AND BOTTOM SALINITY
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      PAM(L)=SAL(L,KK)
      ENDDO
      WRITE (7,509) KK,N  
      CALL PPLOT (1)
C
      ENDDO
C
  509 FORMAT(1H1,' SALINITY, PPT, LAYER',I5,3X,'AT TIME STEP  ',I5,//)
C
C**********************************************************************C
C
C **  PLOT SALINITY STRATIFICATION
C
C----------------------------------------------------------------------C
C
      IF(KC.GT.1)THEN
C
      DO L=2,LA
      PAM(L)=SAL(L,1)-SAL(L,KC)
      ENDDO
      WRITE (7,549) N  
      CALL PPLOT (1)
C
      ENDIF
C
  549 FORMAT(1H1,' SALINITY STRATIFICATION, PPT, AT TIME STEP  ',I5,//)
C
C**********************************************************************C
C
C **  OUTPUT RESIDUAL TOTAL DEPTH
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      PAM(L)=HLPF(L)
      ENDDO
      WRITE (7,585)
      WRITE (7,757)N
      CALL PPLOT (2)
C
  585 FORMAT (1H1,'RESIDUAL TOTAL DEPTH IN METERS',//)
  757 FORMAT (1H1,'AVERAGED OVER TWO TIDAL CYCLES ENDING AT N=',I5,//)
C
C**********************************************************************C
C
C **  OUTPUT EULERIAN RESIDUAL TRANSPORT VELOCITY
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      PAM(L)=0.5*(UHLPF(L,KK)+UHLPF(L+1,KK))/HMP(L)
      ENDDO
      WRITE (7,58) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      LN=LNC(L)
      PAM(L)=0.5*(VHLPF(L,KK)+VHLPF(LN,KK))/HMP(L)
      ENDDO
      WRITE (7,59) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
   58 FORMAT(1H1,' X EULERIAN RESID TRANSPORT VEL, M/S, LAYER',I5,//)
   59 FORMAT(1H1,' Y EULERIAN RESID TRANSPORT VEL, M/S, LAYER',I5,//)
C
C**********************************************************************C
C
C **  OUTPUT VECTOR POTENTIAL TRANSPORT VELOCITY
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      PAM(L)=0.5*(UVPT(L,KK)+UVPT(L+1,KK))/HMP(L)
      ENDDO
      WRITE (7,1458) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      LN=LNC(L)
      PAM(L)=0.5*(VVPT(L,KK)+VVPT(LN,KK))/HMP(L)
      ENDDO
      WRITE (7,1459) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
 1458 FORMAT(1H1,' X VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
 1459 FORMAT(1H1,' Y VECTOR POTENTIAL TRANSPORT VEL, M/S, LAYER',I5,//)
C
C**********************************************************************C
C
C **  OUTPUT LAGRANGIAN RESIDUAL TRANSPORT VELOCITIES
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      PAM(L)=0.5*(UHLPF(L,KK)+UHLPF(L+1,KK)+UVPT(L,KK)
     &       +UVPT(L+1,KK))/HMP(L)
      ENDDO
      WRITE (7,1558) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      LN=LNC(L)
      PAM(L)=0.5*(VHLPF(L,KK)+VHLPF(LN,KK)+VVPT(L,KK)
     &       +VVPT(LN,KK))/HMP(L)
      ENDDO
      WRITE (7,1559) KK
      WRITE (7,757)N
      CALL PPLOT (2)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
 1558 FORMAT(1H1,' X LAGRANGIAN RESID TRANSPORT VEL, M/S, LAYER',I5,//)
 1559 FORMAT(1H1,' Y LAGRANGIAN RESID TRANSPORT VEL, M/S, LAYER',I5,//)
C
C**********************************************************************C
C
C **  OUTPUT RESIDUAL VOLUMETRIC FLOW ACROSS OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
      WRITE (7,60)
      WRITE (7,757)N
      WRITE (7,61) QXW,QXE
      WRITE (7,62) QYS,QYN
      WRITE (7,1661) QXWVP,QXEVP
      WRITE (7,1662) QYSVP,QYNVP
C
C----------------------------------------------------------------------C
C
   60 FORMAT (1H1,' RESIDUAL FLOW ACROSS OPEN BOUNDARIES, M3/S',//)
   61 FORMAT (5X,' QXW=',5X,E12.4,10X,' QXE=',5X,E12.4,//)
   62 FORMAT (5X,' QYS=',5X,E12.4,10X,' QYN=',5X,E12.4,//)
 1661 FORMAT (5X,' QXWVP=',5X,E12.4,10X,' QXEVP=',5X,E12.4,//)
 1662 FORMAT (5X,' QYSVP=',5X,E12.4,10X,' QYNVP=',5X,E12.4,//)
C
C**********************************************************************C
C
C **  OUTPUT RESIDUAL BUOYANCY AND STRATIFICATION
C
C----------------------------------------------------------------------C
C
      DO KK=1,KC,KS
C
      DO L=2,LA
      PAM(L)=SALLPF(L,KK)
      ENDDO
      WRITE (7,569) KK
      WRITE (7,757) N
      CALL PPLOT (1)
C
      ENDDO
C
C----------------------------------------------------------------------C
C
      IF(KC.GT.1)THEN
C
      DO L=2,LA
      LN=LNC(L)
      PAM(L)=SALLPF(L,1)-SALLPF(L,KC)
      ENDDO
      WRITE (7,568) 
      WRITE (7,757) N
      CALL PPLOT (1)
C
      ENDIF
C
C----------------------------------------------------------------------C
C
  568 FORMAT(1H1,' RESIDUAL SALINITY STRATIFICATION, PPT ',//)
  569 FORMAT(1H1,' RESIDUAL SALINITY, PPT, LAYER',I5,//)
C
C**********************************************************************C
C
      RETURN
      END
