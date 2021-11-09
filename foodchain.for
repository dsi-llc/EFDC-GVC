C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE FOODCHAIN(IFINISH)
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
C **  SUBROUTINES OUTPUT SPACE AND TIME AVERAGE TOXICS CONCENTRATIONS
C **  FOR FOOD CHAIN MODEL
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION KBFC(LCM),VALPOCW(LCM,KCM),TMPVOLW(LCM,KCM)
      DIMENSION WTBED(LCM,KBM),VALPOCB(LCM,KBM),TMPVOLB(LCM,KBM)
C
      DIMENSION TMPTXWF(NFDCHZM,NTXM),TMPTXWC(NFDCHZM,NTXM),
     &          TMPTXWP(NFDCHZM,NTXM),TMPTXBF(NFDCHZM,NTXM),
     &          TMPTXBC(NFDCHZM,NTXM),TMPTXBP(NFDCHZM,NTXM),
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
     &          TMPTXBPD(NFDCHZM,NTXM)
      DIMENSION VALBCONC(LCM,KBM)
C####################################################################################
      DIMENSION TMPDOCW(NFDCHZM),TMPPOCW(NFDCHZM),
     &          TMPDOCB(NFDCHZM),TMPPOCB(NFDCHZM),
C####################################################################################
C RM 07/17/04
C Change to average data-based foc in WC and sediment in mg/Kg instead of mg/L
     &          VALWCSS(LCM,KCM),VALBEDSS(LCM,KBM),
     &          TMPWCFC(NFDCHZM),TMPBEDFC(NFDCHZM),
     &          VALWCFC(LCM,KCM),VALBEDFC(LCM,KBM)
C####################################################################################
C 
      DIMENSION VOLFCW(NFDCHZM),VOLFCB(NFDCHZM)
C
      LOGICAL LMASKFC(LCM) 
C
C**********************************************************************C
C
C
C**** HAMRICK ADDED INITIALIZATION OF TEMPORARY SPATIAL VARIABLES
C
      DO L=1,LC
        KBFC(L)=0
      ENDDO
C
      DO K=1,KC
      DO L=1,LC
        VALPOCW(L,K)=0.
        TMPVOLW(L,K)=0.
        VALWCSS(L,K)=0.
        VALWCFC(L,K)=0.
      ENDDO
      ENDDO
C
      DO K=1,KB
      DO L=1,LC
        VALBEDSS(L,K)=0.
        VALBEDFC(L,K)=0.
        WTBED(L,K)=0.
        VALPOCB(L,K)=0.
        TMPVOLB(L,K)=0.
      ENDDO
      ENDDO
C
C**** END HAMRICK ADDED INITIALIZATION OF TEMPORARY SPATIAL VARIABLES
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON 
      ENDIF
C
      IF(IFINISH.EQ.1) GO TO 2000
      IF(JSFDCH.EQ.0)GO TO 1000
C
C      WRITE(8,*)' FIRST ENTRY TO FOODCHAIN.FOR '
C
      OPEN(1,FILE='FOODCHAIN.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='FOODCHAIN.OUT')
      WRITE(1,121)
      WRITE(1,122)
      WRITE(1,123)
      CLOSE(1)
C
C     JSFDCH=0
C
      DO M=1,NFDCHZ
        FDCHDOCW(M)=0. 
        FDCHPOCW(M)=0. 
        FDCHDOCB(M)=0.
        FDCHPOCB(M)=0.
      ENDDO
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        FDCHTXWF(M,NT)=0.
        FDCHTXWC(M,NT)=0.
        FDCHTXWP(M,NT)=0.  
        FDCHTXBF(M,NT)=0.
        FDCHTXBC(M,NT)=0. 
        FDCHTXBP(M,NT)=0.
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        FDCHTXBD(M,NT)=0.        
C####################################################################################
      ENDDO
      ENDDO
C
      TIMFDCH=0.0
C
C**********************************************************************C
C
 1000 CONTINUE
C
      TIMFDCH=TIMFDCH+DTSED
C
C      TIMETMP=TIMESEC-TCON*TBEGIN
C      NEQUIVAL=TIMETMP/DT
C      WRITE(8,*)'FOODCHAIN ETSEC,TIMFDCH,DT,N,NE ',TIMETMP,TIMFDCH,
C     &           DTSED,N,NEQUIVAL
C
C **  INITIALIZE VOLUMES AND VOLUME AVERAGES
C
      DO M=1,NFDCHZ
        VOLFCW(M)=0.
        VOLFCB(M)=0.
      ENDDO
C
      DO M=1,NFDCHZ
        TMPDOCW(M)=0.  
        TMPPOCW(M)=0.  
        TMPDOCB(M)=0. 
        TMPPOCB(M)=0.  
C####################################################################################
C RM 07/17/04
C Change to average foc in WC and sediment in mg/Kg instead of mg/L
        TMPWCFC(M)=0.
        TMPBEDFC(M)=0.
C####################################################################################
      ENDDO
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        TMPTXWF(M,NT)=0.
        TMPTXWC(M,NT)=0.  
        TMPTXWP(M,NT)=0.  
        TMPTXBF(M,NT)=0.  
        TMPTXBC(M,NT)=0.  
        TMPTXBP(M,NT)=0. 
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        TMPTXBPD(M,NT)=0.  
C#################################################################################### 
      ENDDO
      ENDDO
C
C **  INITIALIZE MASK
C
      DO L=2,LA
        LMASKFC(L)=.FALSE.
      ENDDO
C
      DO L=2,LA
        IF(LMASKDRY(L))THEN
          IF(MFDCHZ(L).GT.0)LMASKFC(L)=.TRUE.
        ENDIF
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  VOLUME WEIGHTED AVERAGE OVER WATER COLUMN ZONES
C
C     STDOCW(L,K) HAS UNITS: MG/L OR GM/M**3
C     STPOCW(L,K) AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
C
      IF(ISTPOCW.LE.1)THEN
        DO K=1,KC
          DO L=2,LA
            IF(LMASKFC(L)) VALPOCW(L,K)=STPOCW(L,K)
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTPOCW.GE.2)THEN
        DO K=1,KC
          DO L=2,LA
            IF(LMASKFC(L)) THEN
              VALPOCW(L,K)=0.
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
              VALWCSS(L,K)=0.
              VALWCFC(L,K)=0.
            ENDIF
C####################################################################################
          ENDDO
        ENDDO
        DO NS=1,NSED
          DO K=1,KC
            DO L=2,LA
              IF(LMASKFC(L)) THEN
                VALPOCW(L,K)=VALPOCW(L,K)+SED(L,K,NS)*STFPOCW(L,K,NS)
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
                VALWCSS(L,K)=VALWCSS(L,K)+SED(L,K,NS)
C####################################################################################
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DO NX=1,NSND
          NS=NX+NSED
          DO K=1,KC
            DO L=2,LA
              IF(LMASKFC(L)) THEN
                VALPOCW(L,K)=VALPOCW(L,K)+SND(L,K,NX)*STFPOCW(L,K,NS)
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
                VALWCSS(L,K)=VALWCSS(L,K)+SND(L,K,NX)
C####################################################################################
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
      DO K=1,KC
         DO L=2,LA
cjah start
cjah        VALWCFC(L,K)=VALPOCW(L,K)*1000000./VALWCSS(L,K)
            IF(VALWCSS(L,K).GT.0.) THEN
              VALWCFC(L,K)=VALPOCW(L,K)*1000000./VALWCSS(L,K)
            ELSE
              VALWCFC(L,K)=0.
            ENDIF
cjah end
         ENDDO
      ENDDO
C####################################################################################
C
c     changed to areal weighting from voloume weighting
      DO K=1,KC
      DO L=2,LA
c        IF(LMASKFC(L)) TMPVOLW(L,K)=DXYP(L)*HP(L)*DZC(K)
        IF(LMASKFC(L)) TMPVOLW(L,K)=DXYP(L)*DZC(K)
      ENDDO
      ENDDO
C
      DO K=1,KC
      DO L=2,LA
        IF(LMASKFC(L))THEN
          M=MFDCHZ(L)
          VOLFCW(M)=VOLFCW(M)+TMPVOLW(L,K)
        ENDIF  
      ENDDO
      ENDDO
C
C     TMPTXWF(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
C     TMPTXWC(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
C     TMPTXWP(M,NT)/TMPVOLW HAS UNITS: UG/MG
C     TMPDOCW(M,NT)/TMPVOLW AND STDOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
C     TMPPOCW(M,NT)/TMPVOLW AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
C
      DO K=1,KC
      DO L=2,LA
        IF(LMASKFC(L))THEN
          M=MFDCHZ(L)
          TMPDOCW(M)=TMPDOCW(M)+TMPVOLW(L,K)*STDOCW(L,K)
          TMPPOCW(M)=TMPPOCW(M)+TMPVOLW(L,K)*VALPOCW(L,K)
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
          TMPWCFC(M)=TMPWCFC(M)+TMPVOLW(L,K)*VALWCFC(L,K)
C####################################################################################
        ENDIF  
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=1,KC
      DO L=2,LA
        IF(LMASKFC(L))THEN
        M=MFDCHZ(L)
        TMPTXWF(M,NT)=TMPTXWF(M,NT)
     &               +TMPVOLW(L,K)*TOXFDFW(L,K,NT)*TOX(L,K,NT)
        TMPTXWC(M,NT)=TMPTXWC(M,NT)
     &               +TMPVOLW(L,K)*TOXCDFW(L,K,NT)*TOX(L,K,NT) 
        IF(VALPOCW(L,K).GT.0.) TMPTXWP(M,NT)=TMPTXWP(M,NT)
     &               +TMPVOLW(L,K)*TOXPFTW(L,K,NT)*TOX(L,K,NT)
     &               /VALPOCW(L,K) 
        ENDIF  
      ENDDO
      ENDDO
      ENDDO
C
      DO M=1,NFDCHZ
        IF(VOLFCW(M).GT.0.0)THEN
          TMPDOCW(M)=TMPDOCW(M)/VOLFCW(M)
          TMPPOCW(M)=TMPPOCW(M)/VOLFCW(M)
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
          TMPWCFC(M)=TMPWCFC(M)/VOLFCW(M)
C####################################################################################
        ENDIF
      ENDDO
C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
      DO M=1,NFDCHZ
         TMPPOCW(M)=TMPWCFC(M)
      ENDDO
C####################################################################################
C

      DO NT=1,NTOX
      DO M=1,NFDCHZ
        IF(VOLFCW(M).GT.0.0)THEN
          TMPTXWF(M,NT)=TMPTXWF(M,NT)/VOLFCW(M)
          TMPTXWC(M,NT)=TMPTXWC(M,NT)/VOLFCW(M)
          TMPTXWP(M,NT)=TMPTXWP(M,NT)/VOLFCW(M)
        ENDIF
      ENDDO
      ENDDO
C
C     CONVERT PARTICULATE FROM UG/MG TO UG/GM
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        TMPTXWP(M,NT)=1000.*TMPTXWP(M,NT)
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  VOLUME WEIGHTED AVERAGE OVER BED ZONES
C
C     STDOCB(L,K) HAS UNITS: MG/L OR GM/M**3 (MASS PER VOLUME OF PORE WATER)
C     STPOCB(L,K) AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)
C
      IF(ISTPOCB.LE.1)THEN
        DO K=1,KB
          DO L=2,LA
            IF(LMASKFC(L)) VALPOCB(L,K)=STPOCB(L,K)
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISTPOCB.GE.2)THEN
        DO K=1,KB
          DO L=2,LA
            VALPOCB(L,K)=0.
C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
            VALBEDSS(L,K)=0.
            VALBEDFC(L,K)=0.
C####################################################################################
          ENDDO
        ENDDO
        DO NS=1,NSED
           DO K=1,KB
              DO L=2,LA
                 IF(LMASKFC(L))THEN
                   IF(K.LE.KBT(L)) THEN
                     VALPOCB(L,K)=VALPOCB(L,K)
C####################################################################################
C RM 05/14/04
C Change to average using data-based foc rather than partitioning foc
C Sorry, go back to partitioning based foc
     &                     +SEDB(L,K,NS)*STFPOCB(L,K,NS)/HBED(L,K)  
c     &                    +SEDB(L,K,NS)*FPOCB(L,K)/HBED(L,K)  
C####################################################################################
C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
                     VALBEDFC(L,K)=VALBEDFC(L,K)+SEDB(L,K,NS)*FPOCB(L,K)/HBED(L,K)  
                     VALBEDSS(L,K)=VALBEDSS(L,K)+SEDB(L,K,NS)/HBED(L,K)  
C####################################################################################
                   ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO NX=1,NSND
           NS=NX+NSED
           DO K=1,KB
              DO L=2,LA
                 IF(LMASKFC(L))THEN
                   IF(K.LE.KBT(L)) THEN
                     VALPOCB(L,K)=VALPOCB(L,K)
C####################################################################################
C RM 05/14/04
C Change to average using data-based foc rather than partitioning foc
C Sorry, go back to partitioning based foc
     &                     +SNDB(L,K,NX)*STFPOCB(L,K,NS)/HBED(L,K)
c     &                    +SNDB(L,K,NX)*FPOCB(L,K)/HBED(L,K)
C####################################################################################
C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
                     VALBEDFC(L,K)=VALBEDFC(L,K)+SNDB(L,K,NX)*FPOCB(L,K)/HBED(L,K)  
                     VALBEDSS(L,K)=VALBEDSS(L,K)+SNDB(L,K,NX)/HBED(L,K)  
C####################################################################################
                   ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
      ENDIF

C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
      DO K=1,KB
         DO L=2,LA
cjah start
cjah        VALBEDFC(L,K)=VALBEDFC(L,K)*1000000./VALBEDSS(L,K)
            IF(VALBEDSS(L,K).GT.0.) THEN
              VALBEDFC(L,K)=VALBEDFC(L,K)*1000000./VALBEDSS(L,K)
            ELSE
              VALBEDFC(L,K)=0.
            ENDIF
cjah end
         ENDDO
      ENDDO
C####################################################################################
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCB
      DO K=1,KB
        DO L=2,LA
           VALBCONC(L,K)=0.
        ENDDO
      ENDDO
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            IF(LMASKFC(L))THEN
              IF(K.LE.KBT(L)) VALBCONC(L,K)=VALBCONC(L,K) + SEDB(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO NX=1,NSND
        NS=NX+NSED
        DO K=1,KB
          DO L=2,LA
            IF(LMASKFC(L))THEN
              IF(K.LE.KBT(L)) VALBCONC(L,K)=VALBCONC(L,K) + SNDB(L,K,NX)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C####################################################################################
        
C
      DO K=1,KB
      DO L=2,LA
        WTBED(L,K)=0.0
      ENDDO
      ENDDO
C
      DO L=2,LA
        IF(LMASKFC(L))THEN
          KBFC(L)=KBT(L)
          HBEDTMP=0.0
          KSTOP=0
          DO K=KBT(L),1,-1
            HBEDTMP=HBEDTMP+HBED(L,K)
C####################################################################################
C HQI change to input spatially foodchain averaging depths
C RM 09/16/05
            if(MFDCHZ(L).EQ.1) then
C     Reach 5A
              HBFDCH = FCMHBDUS
            else
              if(MFDCHZ(L).EQ.5) then
C     Reach 6
                HBFDCH = FCMHBDWP
              else
C     Reach 5B, 5C and 5D
                HBFDCH = FCMHBDDS
              endif
            endif
C####################################################################################
            IF(HBEDTMP.GT.HBFDCH.AND.KSTOP.EQ.0)THEN
              KBFC(L)=K
              KSTOP=1
              HBSTOP=HBED(L,K)-HBEDTMP+HBFDCH
C####################################################################################
C RM 05/12/07
C Correction to weightage factor
C              WTBED(L,K)=HBSTOP/HBFDCH
              WTBED(L,K)=HBSTOP/HBED(L,K)
C####################################################################################
            ENDIF
          ENDDO
          KTMP=KBFC(L)+1
          DO K=KTMP,KBT(L)
C####################################################################################
C RM 05/14/04
C Weightages greater than 1 could occur with this method of depth-weighting.
C When the thickness of the top layer is greater than HBFDCH (0.1524
C meters), the weightage assigned to this layer could become greater
C than 1. Need to confirm this with JH.
C####################################################################################
C####################################################################################
C RM 05/12/07
C Correction to weightage factor
C            WTBED(L,K)=HBED(L,K)/HBFDCH
             WTBED(L,K)=1.0
C####################################################################################
          ENDDO
        ENDIF
      ENDDO
C
      IF(JSFDCH.EQ.1)THEN
        OPEN(1,FILE='FOODCHAIN.DIA')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='FOODCHAIN.DIA')
        DO L=2,LA
          IF(LMASKFC(L))THEN
            WRITE(1,111)IL(L),JL(L),KBFC(L),KBT(L),
     &               (WTBED(L,K),K=KBFC(L),KBT(L))
            WRITE(1,112)(HBED(L,K),K=KBFC(L),KBT(L))
          ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
C
      DO K=1,KB
      DO L=2,LA
        IF(LMASKFC(L))THEN
          IF(K.GE.KBFC(L).AND.K.LE.KBT(L))THEN
            TMPVOLB(L,K)=DXYP(L)*WTBED(L,K)*HBED(L,K)
          ENDIF  
        ENDIF  
      ENDDO
      ENDDO
C
      DO K=1,KB
      DO L=2,LA
        IF(LMASKFC(L))THEN
          M=MFDCHZ(L)
          IF(K.GE.KBFC(L).AND.K.LE.KBT(L))THEN
            VOLFCB(M)=VOLFCB(M)+TMPVOLB(L,K)
          ENDIF  
        ENDIF  
      ENDDO
      ENDDO
C
C     TMPTXBF(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
C     TMPTXBC(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
C     TMPTXWP(M,NT)/TMPVOLB HAS UNITS: UG/MG
C     TMPDOCW(M,NT)/TMPVOLB AND STDOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER VOLUME PORE WATER)
C     TMPPOCW(M,NT)/TMPVOLB AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)
C
      DO K=KB,1,-1
      DO L=2,LA
        IF(LMASKFC(L))THEN
          M=MFDCHZ(L)
          IF(K.GE.KBFC(L).AND.K.LE.KBT(L))THEN
            TMPDOCB(M)=TMPDOCB(M)+TMPVOLB(L,K)*STDOCB(L,K)
            TMPPOCB(M)=TMPPOCB(M)+TMPVOLB(L,K)*VALPOCB(L,K)
C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
            TMPBEDFC(M)=TMPBEDFC(M)+TMPVOLB(L,K)*VALBEDFC(L,K)
C####################################################################################
          ENDIF  
        ENDIF  
      ENDDO
      ENDDO
C
      DO NT=1,NTOX
      DO K=KB,1,-1
      DO L=2,LA
        IF(LMASKFC(L))THEN
          M=MFDCHZ(L)
          IF(K.GE.KBFC(L).AND.K.LE.KBT(L))THEN
            TMPVAL=HBED(L,K)*VALPOCB(L,K)
            PORHINV=1.0/(HBED(L,K)*PORBED(L,K))
            TMPTXBF(M,NT)=TMPTXBF(M,NT)
     &                +TMPVOLB(L,K)*PORHINV*TOXFDFB(L,K,NT)*TOXB(L,K,NT)
            TMPTXBC(M,NT)=TMPTXBC(M,NT)
     &                +TMPVOLB(L,K)*PORHINV*TOXCDFB(L,K,NT)*TOXB(L,K,NT)
            IF(TMPVAL.GT.0.) TMPTXBP(M,NT)=TMPTXBP(M,NT)
     &                +TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)
     &                /TMPVAL
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
            TMPTXBPD(M,NT)=TMPTXBPD(M,NT) + 
     &                TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)
     &                /VALBCONC(L,K)
C####################################################################################
          ENDIF  
        ENDIF  
      ENDDO
      ENDDO
      ENDDO
C
      DO M=1,NFDCHZ
         IF(VOLFCB(M).GT.0.0)THEN
          TMPDOCB(M)=TMPDOCB(M)/VOLFCB(M)
          TMPPOCB(M)=TMPPOCB(M)/VOLFCB(M)
C####################################################################################
C RM 07/17/04
C Change to average foc in BED in mg/Kg instead of mg/L
          TMPBEDFC(M)=TMPBEDFC(M)/VOLFCB(M)
C####################################################################################
        ENDIF
      ENDDO

C####################################################################################
C RM 07/17/04
C Change to average foc in WC in mg/Kg instead of mg/L
      DO M=1,NFDCHZ
         TMPPOCB(M)=TMPBEDFC(M)
      ENDDO
C####################################################################################
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        IF(VOLFCB(M).GT.0.0)THEN
          TMPTXBF(M,NT)=TMPTXBF(M,NT)/VOLFCB(M)
          TMPTXBC(M,NT)=TMPTXBC(M,NT)/VOLFCB(M)
          TMPTXBP(M,NT)=TMPTXBP(M,NT)/VOLFCB(M)
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
          TMPTXBPD(M,NT)=TMPTXBPD(M,NT)/VOLFCB(M)
C####################################################################################
        ENDIF
      ENDDO
      ENDDO
C
C     CONVERT PARTICULATE FROM UG/MG TO UG/GM
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        TMPTXBP(M,NT)=1000.*TMPTXBP(M,NT)
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        TMPTXBPD(M,NT)=1000.*TMPTXBPD(M,NT)
C####################################################################################
      ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C
C **  ACCUMULATE THE TIME AVERAGE
C
      DO M=1,NFDCHZ
        FDCHDOCW(M)=FDCHDOCW(M)+DTSED*TMPDOCW(M)  
        FDCHPOCW(M)=FDCHPOCW(M)+DTSED*TMPPOCW(M)  
        FDCHDOCB(M)=FDCHDOCB(M)+DTSED*TMPDOCB(M) 
        FDCHPOCB(M)=FDCHPOCB(M)+DTSED*TMPPOCB(M)   
      ENDDO
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        FDCHTXWF(M,NT)=FDCHTXWF(M,NT)+DTSED*TMPTXWF(M,NT)
        FDCHTXWC(M,NT)=FDCHTXWC(M,NT)+DTSED*TMPTXWC(M,NT)  
        FDCHTXWP(M,NT)=FDCHTXWP(M,NT)+DTSED*TMPTXWP(M,NT)  
        FDCHTXBF(M,NT)=FDCHTXBF(M,NT)+DTSED*TMPTXBF(M,NT)
        FDCHTXBC(M,NT)=FDCHTXBC(M,NT)+DTSED*TMPTXBC(M,NT)  
        FDCHTXBP(M,NT)=FDCHTXBP(M,NT)+DTSED*TMPTXBP(M,NT)  
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        FDCHTXBD(M,NT)=FDCHTXBD(M,NT)+DTSED*TMPTXBPD(M,NT)  
C####################################################################################
        ENDDO
      ENDDO
C
      JSFDCH=0
C
      IF(TIMFDCH.LT.TFCAVG) RETURN
C
C**********************************************************************C
C
C **  COMPLETE AVERAGING AND OUTPUT RESULTS
C
 2000 CONTINUE
C
C      WRITE(8,*)'ENTRY TO FOODCHAIN OUTPUT'
C      TIMETMP=TIMESEC-TCON*TBEGIN
C      NEQUIVAL=TIMETMP/DT
C      WRITE(8,*)'FOODCHAIN ETSEC,TIMFDCH,DT,N,NE ',TIMETMP,TIMFDCH,
C     &           DTSED,N,NEQUIVAL
C
      FDCHVAL=1./TIMFDCH
      DO M=1,NFDCHZ
        FDCHDOCW(M)=FDCHVAL*FDCHDOCW(M)
        FDCHPOCW(M)=FDCHVAL*FDCHPOCW(M)
        FDCHDOCB(M)=FDCHVAL*FDCHDOCB(M)
        FDCHPOCB(M)=FDCHVAL*FDCHPOCB(M)
      ENDDO
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        FDCHTXWF(M,NT)=FDCHVAL*FDCHTXWF(M,NT)
        FDCHTXWC(M,NT)=FDCHVAL*FDCHTXWC(M,NT)
        FDCHTXWP(M,NT)=FDCHVAL*FDCHTXWP(M,NT)
        FDCHTXBF(M,NT)=FDCHVAL*FDCHTXBF(M,NT)
        FDCHTXBC(M,NT)=FDCHVAL*FDCHTXBC(M,NT)
        FDCHTXBP(M,NT)=FDCHVAL*FDCHTXBP(M,NT)
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        FDCHTXBD(M,NT)=FDCHVAL*FDCHTXBD(M,NT)
C####################################################################################
      ENDDO
      ENDDO
C
      OPEN(1,FILE='FOODCHAIN.OUT',POSITION='APPEND')
C
      WRITE(1,101)TIME,NTOX,NFDCHZ,TIMFDCH
C     WRITE(1,103)
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        WRITE(1,102)NT,M,FDCHTXWF(M,NT),FDCHTXWC(M,NT),FDCHTXWP(M,NT),
     &                   FDCHDOCW(M),FDCHPOCW(M),FDCHTXBF(M,NT),
     &                   FDCHTXBC(M,NT),FDCHTXBP(M,NT),FDCHDOCB(M),
     &                   FDCHPOCB(M),FDCHTXBD(M,NT)
      ENDDO
      ENDDO
C
      CLOSE(1)
C
C**********************************************************************C
C
C **  INITIALIZE FOR NEXT AVERAGING PERIOD
C
      DO M=1,NFDCHZ
        FDCHDOCW(M)=0. 
        FDCHPOCW(M)=0. 
        FDCHDOCB(M)=0.
        FDCHPOCB(M)=0.
      ENDDO
C
      DO NT=1,NTOX
      DO M=1,NFDCHZ
        FDCHTXWF(M,NT)=0.
        FDCHTXWC(M,NT)=0.
        FDCHTXWP(M,NT)=0.  
        FDCHTXBF(M,NT)=0.
        FDCHTXBC(M,NT)=0. 
        FDCHTXBP(M,NT)=0. 
C####################################################################################
C RM 05/14/04
C Change to average dry weight PCBs
        FDCHTXBD(M,NT)=0.        
C####################################################################################
      ENDDO
      ENDDO

      TIMFDCH=0.0
C
C**********************************************************************C
C
  111 FORMAT(4I5,10F10.4)
  112 FORMAT(20X,10F10.4)
  101 FORMAT(F12.4,2I7,F12.3)
  102 FORMAT(1X,2I6,11E13.5)
  103 FORMAT('              TXWF         TXWC         TXWP',
     &       '         DOCW         POCW         TXBF         TXBC',
     &       '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
  121 FORMAT('DATA: OUTPUT TIME (DAYS), NTOX, NZONES, ',
     &       'AERAGING PERIOD (SECS)')
  122 FORMAT('DATA: NT    NZ   TXWF         TXWC         TXWP',
     &     '         DOCW         POCW         TXBF         TXBC',
     &       '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
  123 FORMAT('DATA:            UG/L         UG/L         UG/GM OC',
     &       '      MG/L         MG/KG         UG/L         UG/L',
     &       '         UG/GM OC     MG/L         MG/KG       UG/GM Dry')
C
C**********************************************************************C
C
      RETURN
      END

