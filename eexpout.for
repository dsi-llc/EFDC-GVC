C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE EEXPOUT(JSEXPLORER)
C
C**********************************************************************C
C
C      SUBROUTINE BEDTOP  
C  
C **  SUBROUTINE BEDTOP OUTPUTS THE SEDIMENT AND TOXIC VARIABLES  
C **  FOR THE TOP LAYER OF THE SEDIMENT BED  
C  
      INCLUDE'EFDC.PAR'  
      INCLUDE'EFDC.CMN'  
C
      INTEGER*4 IVER     
C
C**********************************************************************C
C
C **  INITIAL CALL
C
      IF(JSEXPLORER.EQ.1)THEN
        ! *** SETUP SEDIMENT AND WATER COLUMN FILES

        ! *** SEDIMENT TOP LAYER AND WATER COLUMN
        OPEN(95,FILE='BED_TOP.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     .       FORM='BINARY')
        CLOSE(95,STATUS='DELETE')
        OPEN(95,FILE='BED_TOP.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     .       FORM='BINARY')
        IVER=103
        WRITE(95)IVER
        WRITE(95)ISTRAN(1),ISTRAN(2),ISTRAN(3),ISTRAN(4)
        WRITE(95)ISTRAN(5),ISTRAN(6),ISTRAN(7)
        WRITE(95)NSED,NSND,KB,KC,NTOX
        NSXD=NSED+NSND  
        DO NS=1,NSXD
         WRITE(95)SEDDIA(NS)
        ENDDO 

        ! *** BED SEDIMENT LAYERS
        IF(ISBEXP.GE.1)THEN
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN

C**********************************************************************C
C *** HQI BEGIN BLOCK
C     HQI 08/25/03, SO and RM
C     Change to output model computed critical shear stress based on
C     option selected for IWRSP
            OPEN(1001,FILE='TAU_CRIT_COH.OUT',FORM='BINARY',
     &                ACCESS='SEQUENTIAL')
            CLOSE(1001,STATUS='DELETE')
            OPEN(1001,FILE='TAU_CRIT_COH.OUT',FORM='BINARY',
     &                ACCESS='SEQUENTIAL')
C *** HQI END BLOCK
C**********************************************************************C

            OPEN(96,FILE='BED_LAY.OUT',STATUS='UNKNOWN',
     .           ACCESS='SEQUENTIAL',FORM='BINARY')
              CLOSE(96,STATUS='DELETE')
            OPEN(96,FILE='BED_LAY.OUT',STATUS='UNKNOWN',
     .      ACCESS='SEQUENTIAL',FORM='BINARY')
            WRITE(96)IVER
            WRITE(96)ISTRAN(1),ISTRAN(2),ISTRAN(3),ISTRAN(4)
            WRITE(96)ISTRAN(5),ISTRAN(6),ISTRAN(7)
            WRITE(96)NSED,NSND,KB,KC,NTOX  
              NSXD=NSED+NSND  
            DO NS=1,NSXD
             WRITE(96)SEDDIA(NS)
            ENDDO 
          ENDIF 
        ENDIF 
      ENDIF 
C
C**********************************************************************C
C
C **  SUBSEQUENT CALLS
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
      ELSE
        TIME=TIMESEC
      ENDIF
      IF(JSEXPLORER.EQ.1)TIME=TCON*TBEGIN

      ! *** TIME IN SECONDS CONIVERTED TO DAYS (EE TIME UNITS)
      TIME=TIME/86400.

      IF(ISSPH(8).GE.1)THEN
C
      WRITE(95)TIME,LA-1

      DO L=2,LA  
        N1=KBT(L)
        IF(ISTRAN(6).GT.0.OR.ISTRAN(7).GT.0)THEN  
          IF(ISBEDSTR.GE.1)THEN
            WRITE(95)TAUBSED(L)
            IF(ISBEDSTR.EQ.1)THEN
              WRITE(95)TAUBSND(L)  
            ENDIF
          ELSE
            WRITE(95)TAUB(L)  
          ENDIF
        ELSE  
          SHEAR=MAX(QQ(L,0),QQMIN)/CTURB2  
          WRITE(95)SHEAR  
        ENDIF  
        IF(ISTRAN(1).EQ.1) WRITE(95)(SAL(L,K),K=1,KC)  
        IF(ISTRAN(2).EQ.1) WRITE(95)(TEM(L,K),K=1,KC)  
        IF(ISTRAN(3).EQ.1) WRITE(95)(DYE(L,K),K=1,KC)  
        IF(ISTRAN(4).EQ.1) WRITE(95)(SFL(L,K),K=1,KC)  

        IF(ISTRAN(5).EQ.1)THEN
           WRITE(95)(TOXB(L,N1,NT),NT=1,NTOX)
           WRITE(95)((TOX(L,K,NT),K=1,KC),NT=1,NTOX)
        ENDIF

        IF(ISTRAN(6).EQ.1.OR.ISTRAN(7).GE.1)THEN
          WRITE(95) N1,BELV(L),HBED(L,N1),BDENBED(L,N1),PORBED(L,N1)
          IF(ISTRAN(6).EQ.1)THEN 
             WRITE(95)(SEDB(L,N1,NS),VFRBED(L,N1,NS),NS=1,NSED)
             WRITE(95)((SED(L,K,NS),K=1,KC),NS=1,NSED)
          ENDIF  
          IF(ISTRAN(7).EQ.1)THEN
            WRITE(95)(SNDB(L,N1,NX),VFRBED(L,N1,NX+NSED),NX=1,NSND)
            WRITE(95)((SND(L,K,NX),K=1,KC),NX=1,NSND)
            IF(ISBDLDBC.GT.0)THEN
              WRITE(95)(CQBEDLOADX(L,NX),CQBEDLOADY(L,NX),NX=1,NSND)
            ENDIF
          ENDIF
        ENDIF

      ENDDO  
      ENDIF
C  
C *** NOW OUTPUT ALL THE BEDINFO TO A SINGLE FILE  
C
      IF(ISBEXP.GE.1)THEN
        IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1.AND.KB.GT.1)THEN

          ! *** WRITE TIME STEP RESULTS
          WRITE(96)TIME,LA-1  
  
          DO L=2,LA  
            WRITE(96)KBT(L)  
          ENDDO

          DO L=2,LA  
C102            WRITE(96)KBT(L)  
            DO K=1,KBT(L)  
              WRITE(96)HBED(L,K),BDENBED(L,K),PORBED(L,K)  
              IF(ISTRAN(6).GE.1)THEN  
                DO NS=1,NSED  
                  WRITE(96)SEDB(L,K,NS),VFRBED(L,K,NS)  
                ENDDO  
              ENDIF  
              IF(ISTRAN(7).GE.1)THEN  
                DO NX=1,NSND  
                  NS=NSED+NX  
                  WRITE(96)SNDB(L,K,NX),VFRBED(L,K,NS)  
                ENDDO  
              ENDIF  
              IF(ISTRAN(5).GE.1)THEN  
                DO NT=1,NTOX  
                  WRITE(96)TOXB(L,K,NT)  
                ENDDO  
              ENDIF  
            ENDDO  
          ENDDO
        ENDIF  
      ENDIF 
C  
C**********************************************************************C
C
c moved to bedplth.for      IF(JSBPHA.EQ.1)THEN
c moved to bedplth.for        OPEN(1,FILE='BEDARD.OUT')
c moved to bedplth.for        CLOSE(1,STATUS='DELETE')
c moved to bedplth.for        OPEN(1,FILE='BEDARD.OUT')
c moved to bedplth.for        WRITE(1,131)
c moved to bedplth.for        CLOSE(1)
c moved to bedplth.for       JSBPHA=0
c moved to bedplth.for      ENDIF
C
c moved to bedplth.for      IF(ISBARD.GE.1)THEN
c moved to bedplth.for        OPEN(1,FILE='BEDARD.OUT',POSITION='APPEND')
c moved to bedplth.for        WRITE(1,122)TIME
c moved to bedplth.for        DO L=2,LA
c moved to bedplth.for          WRITE(1,101)IL(L),JL(L),
c moved to bedplth.for     &               (SEDFDTAP(L,NS),SEDFDTAN(L,NS),NS=1,NSED),
c moved to bedplth.for     &               (SNDFDTAP(L,NS),SNDFDTAN(L,NS),NS=1,NSND)
c moved to bedplth.for        ENDDO
c moved to bedplth.for        CLOSE(1)
c moved to bedplth.for      ENDIF
C
c moved to bedplth.for  101 FORMAT(2I5,18E13.5)
c moved to bedplth.for  131 FORMAT(' IL JL (SEDFDTAP SEDFDTAN)(1,NSED)',
c moved to bedplth.for     &             ' (SNDFDTAP SNDFDTAN)(1,NSND)')
c moved to bedplth.for  122 FORMAT(F12.5,' TIME OF OUTPUT')
C
C**********************************************************************C
C
      RETURN
      END
