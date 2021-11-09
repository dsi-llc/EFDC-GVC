C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQCSR
C
C**********************************************************************C
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
C     READ AND UPDATE WATER COLUMN WATER QUALITY VARIABLE TIME SERIES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
       CHARACTER*11 FNWQSR(21)
C
C**********************************************************************C
C
      IF(ITNWQ.GT.0) GOTO 1000
C
C**********************************************************************C
C
       FNWQSR( 1)='CWQSR01.INP'
       FNWQSR( 2)='CWQSR02.INP'
       FNWQSR( 3)='CWQSR03.INP'
       FNWQSR( 4)='CWQSR04.INP'
       FNWQSR( 5)='CWQSR05.INP'
       FNWQSR( 6)='CWQSR06.INP'
       FNWQSR( 7)='CWQSR07.INP'
       FNWQSR( 8)='CWQSR08.INP'
       FNWQSR( 9)='CWQSR09.INP'
       FNWQSR(10)='CWQSR10.INP'
       FNWQSR(11)='CWQSR11.INP'
       FNWQSR(12)='CWQSR12.INP'
       FNWQSR(13)='CWQSR13.INP'
       FNWQSR(14)='CWQSR14.INP'
       FNWQSR(15)='CWQSR15.INP'
       FNWQSR(16)='CWQSR16.INP'
       FNWQSR(17)='CWQSR17.INP'
       FNWQSR(18)='CWQSR18.INP'
       FNWQSR(19)='CWQSR19.INP'
       FNWQSR(20)='CWQSR20.INP'
       FNWQSR(21)='CWQSR21.INP'
C
C**********************************************************************C
C
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
C **  TIME SERIES FROM THE FILES CWQSRNN.INP
C
C----------------------------------------------------------------------C
C
      DO NW=1,NWQV
C
      IF(NWQCSR(NW).GE.1)THEN
        OPEN(1,FILE=FNWQSR(NW),STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      ENDDO
C
        DO NS=1,NWQCSR(NW)
        MWQCTLT(NS,NW)=1
        READ(1,*,IOSTAT=ISO)ISTYP,MWQCSR(NS,NW),TCWQCSR(NS,NW),
     &                   TAWQCSR(NS,NW),RMULADJ,ADDADJ
        IF(ISO.GT.0) GOTO 900
        IF(ISTYP.EQ.1)THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GOTO 900
           DO M=1,MWQCSR(NS,NW)
           READ(1,*,IOSTAT=ISO)TWQCSER(M,NS,NW),CSERTMP
           IF(ISO.GT.0) GOTO 900
           TWQCSER(M,NS,NW)=TWQCSER(M,NS,NW)+TAWQCSR(NS,NW)
            DO K=1,KC
            WQCSER(M,K,NS,NW)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
           ENDDO
         ELSE
          DO M=1,MWQCSR(NS,NW)
          READ(1,*,IOSTAT=ISO)TWQCSER(M,NS,NW),
     &                       (WQCSER(M,K,NS,NW), K=1,KC)
          IF(ISO.GT.0) GOTO 900
          TWQCSER(M,NS,NW)=TWQCSER(M,NS,NW)+TAWQCSR(NS,NW)
           DO K=1,KC
           WQCSER(M,K,NS,NW)=RMULADJ*(WQCSER(M,K,NS,NW)+ADDADJ)
           ENDDO
          ENDDO
        ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
      ENDDO
C
C     WRITE(6,602)
      GOTO 901
C
  900 CONTINUE
      WRITE(6,601)NW,NS,M
      STOP
C
  901 CONTINUE
C
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)
  602 FORMAT(' READ OF FILES CWQSRNN.INP SUCCESSFUL'/)
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  INITIALIZE NULL SERIES CONCENTRATIONS
C
      DO NW=1,NWQV
       DO K=1,KC
        CSERTWQ(K,0,NW)=0.
       ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CONCENTRATION SERIES INTERPOLTATION FOR WATER QUALITY VARIABLES
C
      DO NW=1,NWQV
C
       DO NS=1,NWQCSR(NW)
       IF(ISDYNSTP.EQ.0)THEN
         TIME=DT*FLOAT(N)+TCON*TBEGIN
         TIME=TIME/TCWQCSR(NS,NW)
       ELSE
         TIME=TIMESEC/TCWQCSR(NS,NW)
       ENDIF
C
       M1=MWQCTLT(NS,NW)
  100  CONTINUE
       M2=M1+1
       IF(TIME.GT.TWQCSER(M2,NS,NW))THEN
         M1=M2
         GOTO 100
        ELSE
         MWQCTLT(NS,NW)=M1
       ENDIF
C
       TDIFF=TWQCSER(M2,NS,NW)-TWQCSER(M1,NS,NW)
       WTM1=(TWQCSER(M2,NS,NW)-TIME)/TDIFF
       WTM2=(TIME-TWQCSER(M1,NS,NW))/TDIFF
        DO K=1,KC
         CSERTWQ(K,NS,NW)=WTM1*WQCSER(M1,K,NS,NW)
     &                   +WTM2*WQCSER(M2,K,NS,NW)
        ENDDO
C      WRITE(6,6000)N,CSERTWQ(1,NS,NW),CSERTWQ(KC,NS,NW)
       ENDDO
C
      ENDDO
C
C**********************************************************************C
C
C **  ON FIRST CALL, INITIALIZE OUTFLOW CONCENTRATIONS
C
C----------------------------------------------------------------------C
C
      IF(ITNWQ.EQ.0)THEN
C
      DO NW=1,NWQV
      DO K=1,KC
       DO LL=1,NWQOBS
        NSID=IWQOBS(LL,NW)
        L=LIJ( IWQCBS(LL),JWQCBS(LL) )
        CWQLOS(LL,K,NW)=WTCI(K,1)*WQOBCS(LL,1,NW)
     &    +WTCI(K,2)*WQOBCS(LL,2,NW)+CSERTWQ(K,NSID,NW)
        NWQLOS(LL,K,NW)=0
       ENDDO
       DO LL=1,NWQOBW
        NSID=IWQOBW(LL,NW)
        L=LIJ( IWQCBW(LL),JWQCBW(LL) )
        CWQLOW(LL,K,NW)=WTCI(K,1)*WQOBCW(LL,1,NW)
     &    +WTCI(K,2)*WQOBCW(LL,2,NW)+CSERTWQ(K,NSID,NW)
        NWQLOW(LL,K,NW)=0
       ENDDO
       DO LL=1,NWQOBE
        NSID=IWQOBE(LL,NW)
        L=LIJ( IWQCBE(LL),JWQCBE(LL) )
        CWQLOE(LL,K,NW)=WTCI(K,1)*WQOBCE(LL,1,NW)
     &    +WTCI(K,2)*WQOBCE(LL,2,NW)+CSERTWQ(K,NSID,NW)
        NWQLOE(LL,K,NW)=0
       ENDDO
       DO LL=1,NWQOBN
        NSID=IWQOBN(LL,NW)
        L=LIJ( IWQCBN(LL),JWQCBN(LL) )
        CWQLON(LL,K,NW)=WTCI(K,1)*WQOBCN(LL,1,NW)
     &    +WTCI(K,2)*WQOBCN(LL,2,NW)+CSERTWQ(K,NSID,NW)
        NWQLON(LL,K,NW)=0
       ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      RETURN
      END
