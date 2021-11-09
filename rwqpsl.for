C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
CTT      OLD SUBROUTINE RWQPSL
C
C**********************************************************************C
C
C READ IN TEMPORALLY VARYING POINT SOURCE INPUT (UNIT INWQPSL).
C: IN (KG/D) EXCEPT XPSQ(M^3/S), XPO2(G/M^3),
C                  XPTAM(KMOL/D), XPFCB(MPN/100ML).
C: IN GES, LOAD IS IN (G/D) EXCEPT TAM IN (MOL/D),
C                  FCB IN (MPN/100ML)*(M^3/S).
C: TO CONVERT KG/D TO G/D, (XW KG/D)*(10^3 G/KG) = (CONV1*XW G/D) WITH
C   CONV1=1.0E3.
C: FOR O2, (XPSQ M^3/S)*(XPO2 G/M^3)*(86400 S/D)=(CONV2*XPSQ*XPO2 G/D)
C   WITH CONV2=8.64E4.
C: FOR TAM, (XPTAM KMOL/D)*(10^3 MOL/KMOL) = (CONV1*XPTAM MOL/D).
C: FOR FCB, (XPFCB MPN/100ML)*(XPSQ M^3/S)*(86400 S/D) =
C   (CONV2*XPSQ*XPFCB (MPN/100ML)*M^3/D).
C
C**********************************************************************C
C
CTT     INCLUDE 'EFDC.PAR'
CTT      INCLUDE 'EFDC.CMN'
C
CTT      PARAMETER (CONV1=1.0E3,CONV2=8.64E4)
CTT      CHARACTER TITLE(3)*79, PSLCONT*3
C
CTT      OPEN(1,FILE=PSLFN,STATUS='UNKNOWN')
CTT      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
CTT      IF(IWQTPSL.EQ.0)THEN
CTT        READ(1,50) (TITLE(M),M=1,3)
CTT        WRITE(2,999)
CTT        WRITE(2,50) (TITLE(M),M=1,3)
CTT        READ(1,999)
CTT        READ(1,94) IWQPS
CTT        WRITE(2,83)'* NUMBER OF CELLS FOR POINT SOURCE INPUT      = ',
CTT     *    IWQPS
CTT      ENDIF
CTT      WRITE(2,60)'* POINT SOURCE INPUT AT ', IWQTPSL,
CTT     *  ' TH TC FROM MODEL START'
C
CTT      READ(1,999)
CTT      READ(1,50) (TITLE(M),M=1,3)
CTT      WRITE(2,50) (TITLE(M),M=1,3)
CTT      DO M=1,IWQPS
CTT        READ(1,94) I,J,K,WQPSQ(M),(WQWPSL(M,NW),NW=1,NWQV)
CTT        WRITE(2,94) I,J,K,WQPSQ(M),(WQWPSL(M,NW),NW=1,NWQV)
CTT        IF(IJCT(I,J).LT.1 .OR. IJCT(I,J).GT.8)THEN
CTT          PRINT*, 'I, J, K, M, WQPSQ = ', I,J,K,M,WQPSQ(M)
CTT          STOP 'ERROR!! INVALID (I,J) IN FILE 1'
CTT        ENDIF
CTT        L=LIJ(I,J)
CTT        IWQPSC(L,K)=M
CTT        DO NW=1,18
CTT          WQWPSL(M,NW) = WQWPSL(M,NW) * CONV1
CTT        ENDDO
CTT        WQTT = WQPSQ(M)*CONV2
CTT        WQWPSL(M,19) = WQWPSL(M,19) * WQTT
CTT        WQWPSL(M,20) = WQWPSL(M,20) * CONV1
CTT        WQWPSL(M,NWQV) = WQWPSL(M,NWQV) * WQTT
CTT      ENDDO
C
CTT      READ(1,52) IWQTPSL, PSLCONT
CTT      WRITE(2,52) IWQTPSL, PSLCONT
CTT      IF(PSLCONT.EQ.'END')THEN
CTT        CLOSE(1)
CTT        IWQPSL = 0
CTT      ENDIF
C
CTT      CLOSE(2)
C
CTT  999 FORMAT(1X)
CTT   50 FORMAT(A79)
CTT   52 FORMAT(I7, 1X, A3)
CTT   60 FORMAT(/, A24, I5, A24)
CTT   83 FORMAT(A48, I5)
CTT   94 FORMAT(3I5,1X,7F8.3, /, 16X,7F8.3, /, 16X, 8F8.3)
C
CTT      RETURN
CTT      END
C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQPSL
C
C**********************************************************************C
C
C **  NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 8 AUGUST 2001
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
C
C**********************************************************************C
C
C **  ALL LOADS ARE IN KG/DAY EXPECT COLIFORM IN MPN/DAY
C **  INTERNAL CONVERSION TO GM/DAY FOR FIRST 20 STATE VARIABLES
C **  FIX FOR TOTAL ACTIVE METAL AND COLIFORM IN FUTURE VERSIONS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
       DIMENSION RLDTMP(21)
C
C**********************************************************************C
C
      IF(ITNWQ.GT.0) GOTO 1000
C
C**********************************************************************C
C
C **  READ IN LOADING SERIES FROM FILE 'WQPSL.INP'
C
C----------------------------------------------------------------------C
C
      IF( NPSTMSR.GE.1)THEN
        OPEN(1,FILE='WQPSL.INP',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,13
          READ(1,1)
        ENDDO
C
        DO NS=1,NPSTMSR
         MWQPTLT(NS)=1
         READ(1,*,IOSTAT=ISO)MWQPSR(NS),TCWQPSR(NS),
     &                   TAWQPSR(NS),RMULADJ,ADDADJ
         IF(ISO.GT.0) GOTO 900
         RMULADJ=1000.*RMULADJ
         ADDADJ=ADDADJ
         DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(NW),NW=1,7)
          IF(ISO.GT.0) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=8,14)
          IF(ISO.GT.0) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=15,21)
          IF(ISO.GT.0) GOTO 900
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
C  FOR FECAL COLIFORM BACTERIA:
          WQPSSER(M,21,NS)=RMULADJ*RLDTMP(21)
          DO NW=1,20
            WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)
          ENDDO
          DO NW=1,21
            WQPSSER(M,NW,NS)=MAX(WQPSSER(M,NW,NS),0.0)
          ENDDO
         ENDDO
        ENDDO
C
        CLOSE(1)
      ENDIF
C
C     WRITE(6,602)
C
      GOTO 901
C
  900 CONTINUE
      WRITE(6,601)NS,M
      STOP
C
  901 CONTINUE
C
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQPS TIME SERIES, NSER,MDATA = ',2I5)
  602 FORMAT(' READ OF FILE WQPSL.INP SUCCESSFUL'/)
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  INITIALIZE NULL SERIES LOADING TO ZERO
C
      DO NW=1,NWQV
       WQPSSRT(NW,0)=0.
      ENDDO
C
C**********************************************************************C
C
C **  LOADING SERIES INTERPOLTATION
C
      DO NS=1,NPSTMSR
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCWQPSR(NS)
      ELSE
        TIME=TIMESEC/TCWQPSR(NS)
      ENDIF
C
       M1=MWQPTLT(NS)
  100  CONTINUE
       M2=M1+1
       IF(TIME.GT.TWQPSER(M2,NS))THEN
         M1=M2
         GOTO 100
        ELSE
         MWQPTLT(NS)=M1
       ENDIF
C
       TDIFF=TWQPSER(M2,NS)-TWQPSER(M1,NS)
       WTM1=(TWQPSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TWQPSER(M1,NS))/TDIFF
        DO NW=1,NWQV
         WQPSSRT(NW,NS)=WTM1*WQPSSER(M1,NW,NS)
     &                   +WTM2*WQPSSER(M2,NW,NS)
        ENDDO
C      WRITE(6,6000)N,CSERTWQ(1,NS,NW),CSERTWQ(KC,NS,NW)
C
      ENDDO
C
C
      IF(ITNWQ.EQ.0)THEN
C
       OPEN(1,FILE='WQPSLT.DIA',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WQPSLT.DIA',STATUS='UNKNOWN')
C
       WRITE(1,112)N,TIME
C
       DO NS=1,NPSTMSR
            WRITE(1,111)NS,(WQPSSRT(NW,NS),NW=1,NWQV)
       ENDDO
C
      CLOSE(1)
C
      ENDIF
C
C**********************************************************************C
C
C **  COMBINE CONSTANT AND TIME VARIABLE PS LOADS
C
C M.R. MORTON 02/20/1999
C MODIFIED SO MULTIPLE POINT SOURCES CAN BE ADDED TO ANY GRID CELL
C AND ANY LAYER (HAD TO CHANGE WQWPSL ARRAY FROM 2D TO 3D).
CMRM      DO NW=1,NWQV
CMRM       DO K=1,KC
CMRM        DO L=2,LA
CMRM         WQWPSL(IWQPSC(L,K),NW)=WQWPSLC(IWQPSC(L,K),NW)
CMRM     &                        +WQPSSRT(NW,IWQPSV(L,K))
CMRM        ENDDO
CMRM       ENDDO
CMRM      ENDDO
C
CMRM  INITIALIZE COMBINED PSL ARRAY TO ZERO:
C
      DO NW=1,NWQV
       DO K=1,KC
        DO L=2,LA
         WQWPSL(L,K,NW) = 0.0
        ENDDO
       ENDDO
      ENDDO
C
CMRM  ADD CONSTANT AND VARIABLE PSLS TO APPROPRIATE GRID CELLS:
C
      IF(ITNWQ.EQ.0)THEN
C
       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
       WRITE(1,112)N,TIME
C 
      ENDIF
C
      DO NS=1,IWQPS
        L = LIJ(ICPSL(NS), JCPSL(NS))
        K = KCPSL(NS)
        ITMP = MVPSL(NS)
        IF(ITNWQ.EQ.0) WRITE(1,121)NS,L,ICPSL(NS),JCPSL(NS),K,ITMP
        IF(K.GE.1)THEN
          DO NW=1,NWQV
           WQWPSL(L,K,NW) = WQWPSL(L,K,NW)
     +         + WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP)
          ENDDO
         ELSE

C original even distribution of load was for sigma KC is the total layers
C Sen Bai modified to add even distribution of load when k=0 for GVC code
        if (IGRIDV .eq. 0) then  
          TMPVAL=1./FLOAT(KC)
          DO KK=1,KC
          DO NW=1,NWQV
           WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW)
     +         + TMPVAL*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
          ENDDO
          ENDDO
C when    IGRIDV=1 
        else 
          TMPVAL=1.0/ (FLOAT(KC)/GVCSCLP(L))
		TMPVAL=1./ (GVCSCLPI(L)*FLOAT(KC))
		DO KK=KC,KC-int(GVCSCLPI(L)*KC)+1,-1
          DO NW=1,NWQV
           WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW)
     +         + TMPVAL*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
          ENDDO
          ENDDO
        end if 


        ENDIF
      ENDDO
C
      IF(ITNWQ.EQ.0)THEN
C
C       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
C       CLOSE(1,STATUS='DELETE')
C       OPEN(1,FILE='WQPSL.DIA',STATUS='UNKNOWN')
C
C       WRITE(1,112)N,TIME
C
       DO L=2,LA
        ITMP=IWQPSC(L,1)
        IF(ITMP.GT.0)THEN
          DO K=1,KC
            WRITE(1,110)ITMP,IL(L),JL(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
          ENDDO
        ENDIF
       ENDDO
C
C       DO K=1,KC
C        DO L=2,LA
C         ITMP=IWQPSC(L,K)
C         IF(ITMP.GT.0)THEN
C            WRITE(1,110)IL(L),JL(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
C         ENDIF
C        ENDDO
C       ENDDO
C
      CLOSE(1)
C
      ENDIF
C
  110 FORMAT(1X,4I4,2X,7E12.4,/,19X,7E12.4,/,19X,7E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,7E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
  121 FORMAT(' NS,L,I,J,K,ITMP = ', 6I5/)
C
C**********************************************************************C
C
      RETURN
      END
