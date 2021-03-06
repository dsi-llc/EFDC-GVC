C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQ3DINP
C
C**********************************************************************C
C
C  READ WATER QUALITY SUBMODEL INPUT FILES
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER*3 CWQHDR(NWQVM)
C
      CHARACTER*11  HHMMSS
C
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
      IWQTICI=IWQTICI
      IWQTAGR=IWQTAGR
      IWQTSTL=IWQTSTL
      IWQTSUN=IWQTSUN
      IWQTBEN=IWQTBEN
      IWQTPSL=IWQTPSL
      IWQTNPL=IWQTNPL
      ISMTICI=ISMTICI
C
CTT      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
CTT      CLOSE(1,STATUS='DELETE')
C
      OPEN(1,FILE='WQ3D.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
C
C **  HARDWIRE BY PASS OF RATE COEFFICIENT MAPS
C
C MRM      ISWQCMAP=0
C MRM      ISWQSMAP=0
C
      NWQKCNT=0
      NWQKDPT=1
C
      UHEQ(1)=0.0
      UHEQ(LC)=0.0
      DO ND=1,NDMWQ
       LF=2+(ND-1)*LDMWQ
       LL=LF+LDM-1
       DO L=LF,LL
        UHEQ(L)=1.0
       ENDDO
      ENDDO
CH-
      TINDAY = 0.0
      TINDAY=TINDAY
      ITNWQ = 0
C
      RKCWQ = 1.0/REAL(KC)
      DO K=1,KC
        WQHT(K)=REAL(KC-K)*RKCWQ
      ENDDO
C
CXH      INWQRST=21
CXH      INWQICI=21
CXH      INWQAGR=22
CXH      INWQSTL=23
CXH      INWQSUN=24
CXH      INWQBEN=25
CXH      INWQPSL=26
CXH      INWQNPL=27
CXH      IWQONC=30
CXH      IWQORST=31
C
C      WQTSNAME(1)  = 'CHL'
C      WQTSNAME(2)  = 'TPC'
C      WQTSNAME(3)  = 'DOC'
C      WQTSNAME(4)  = 'TPP'
C      WQTSNAME(5)  = 'DOP'
C      WQTSNAME(6)  = 'P4T'
C      WQTSNAME(7)  = 'P4D'
C      WQTSNAME(8)  = 'APC'
C      WQTSNAME(9)  = 'TNN'
C      WQTSNAME(10) = 'DON'
C      WQTSNAME(11) = 'NH4'
C      WQTSNAME(12) = 'NO3'
C      WQTSNAME(13) = 'TSI'
C      WQTSNAME(14) = 'SUU'
C      WQTSNAME(15) = 'SAA'
C      WQTSNAME(16) = 'SAD'
C      WQTSNAME(17) = 'COD'
C      WQTSNAME(18) = 'DOO'
C      WQTSNAME(19) = 'TAM'
C      WQTSNAME(20) = 'TMP'
C      WQTSNAME(21) = 'FCB'
C
      WQTSNAME(1)  = 'CHC'
      WQTSNAME(2)  = 'CHG'
      WQTSNAME(3)  = 'CHD'
      WQTSNAME(4)  = 'ROC'
      WQTSNAME(5)  = 'LOC'
      WQTSNAME(6)  = 'DOC'
      WQTSNAME(7)  = 'ROP'
      WQTSNAME(8)  = 'LOP'
      WQTSNAME(9)  = 'DOP'
      WQTSNAME(10) = 'P04'
      WQTSNAME(11) = 'RON'
      WQTSNAME(12) = 'LON'
      WQTSNAME(13) = 'DON'
      WQTSNAME(14) = 'NHX'
      WQTSNAME(15) = 'NOX'
      WQTSNAME(16) = 'SUU'
      WQTSNAME(17) = 'SAA'
      WQTSNAME(18) = 'COD'
      WQTSNAME(19) = 'DOX'
      WQTSNAME(20) = 'TAM'
      WQTSNAME(21) = 'FCB'
C Sen Bai 
      WQTSNAME(22) = 'MAC'
C End of Sen Bai

      DO M=0,NWQPS
       WQPSQ(M)=0.0
       WQPSQC(M)=0.0
       DO J=1,NWQV
        WQWPSLC(M,J)=0.0
       ENDDO
      ENDDO
C
      DO K=1,KC
       IWQPSC(1,K)=0
       WQDSQ(1,K)=0.0
       IWQPSC(LC,K)=0
       WQDSQ(LC,K)=0.0
      ENDDO
C
      DO ND=1,NDMWQ
       LF=2+(ND-1)*LDMWQ
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         IWQPSC(L,K)=0
         IWQPSV(L,K)=0
         WQDSQ(L,K)=0.0
        ENDDO
       ENDDO
      ENDDO
C
      DO J=1,NWQV
       DO K=1,KC
        DO L=1,LC
         WQWDSL(L,K,J)=0.0
         WQWPSL(L,K,J)=0.0
        ENDDO
       ENDDO
      ENDDO
C
      CALL RWQC1(IWQDT)
C      CALL RWQC2
C      CALL RWQMAP
C
C
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN')
C
      NWQVOUT=0
      DO NW=1,NWQV
       IF(ISTRWQ(NW).EQ.1)THEN
         NWQVOUT=NWQVOUT+1
         CWQHDR(NWQVOUT)=WQTSNAME(NW)
       ENDIF
      ENDDO
C
C Sen Bai
C Add macroalgae output      
	NWQVOUT=NWQVOUT+1
      CWQHDR(NWQVOUT)=WQTSNAME(NWQV+1)
C End of Sen Bai
      WRITE(1,1969)(CWQHDR(NW),NW=1,NWQVOUT)
C
 1969 FORMAT('C   I    J    K    TIME',7X,A3,8X,A3,8X,A3,
     &  8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,    
     &  8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,    
     &  8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3)
C    
C 1969 FORMAT('C   I    J    K    TIME       CHL        TPC',
C     &       '        DOC        TPP        DOP        P4T',
C     &       '        P4D        APC        TNN        DON',
C     &       '        NH4        NO3        TSI        SUU',
C     &       '        SAA        SAD        COD        DOO',
C     &       '        TAM        TMP        FCB        MALG')
C
      CLOSE(1)
C
C **  INITIALIZE DIURNAL DO ANALYSIS
C
      IF(NDDOAVG.GE.1)THEN
        OPEN(1,FILE='DIURNDO.OUT')
        CLOSE(1,STATUS='DELETE')
        DO K=1,KC
         DO L=2,LA
          DDOMAX(L,K)=-1.E6
          DDOMIN(L,K)=1.E6
         ENDDO
        ENDDO
      ENDIF
C
C **  INITIALIZE LIGHT EXTINCTION ANALYSIS
C
      IF(NDLTAVG.GE.1)THEN
        OPEN(1,FILE='LIGHT.OUT')
        CLOSE(1,STATUS='DELETE')
        NDLTCNT=0
        DO K=1,KC
         DO L=2,LA
          RLIGHTT(L,K)=0.
          RLIGHTC(L,K)=0.
         ENDDO
        ENDDO
      ENDIF
C
C ** INITIALIZE WATER QUALITY AVERAGING SUMMATION ARRAYS:
C
      CALL WQZERO
C
      IF(IWQICI.EQ.2) CALL RWQRST
C
      IF(IWQBEN.EQ.1)THEN
        CALL SMINIT
        CALL SMRIN1(IWQDT)
C        CALL SMRIN2
C        CALL RSMMAP
        IF(ISMICI.EQ.2) CALL RSMRST
      ENDIF
C
C      OPEN(1,FILE='TFILE',STATUS='UNKNOWN')
C      CALL TIME(HHMMSS)
C      CLOSE(1,STATUS='DELETE')
C      OPEN(1,FILE='TFILE',STATUS='UNKNOWN')
C      WRITE(1,100) HHMMSS
C      CLOSE(1)
C
  100 FORMAT('  TIME = ',A11,' HH.MM.SS.HH')
CMRM   100 FORMAT('  TIME = ',A11,' HH:MM:SS.SS')
C
      RETURN
      END
