C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WASP7
C
C     MIKE MORTONS WASP5MRM
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
C ==========
C REVISIONS:
C ==========
C  M. MORTON 06/06/94: THIS VERSION WRITES DISPERSION TO THE WASPDH.OUT
C     HYDRODYNAMIC FILE INSTEAD OF WASPB.OUT.  A MODIFIED VERSION OF
C     WASP5 IS NECESSARY TO UTILIZE THESE DISPERSIONS.
C  M. MORTON 06/07/94: WRITES HYDRODYNAMIC INFORMATION AND DISPERSION TO
C     AN UNFORMATTED BINARY FILE, WASPDHU.OUT.  THE CORRECT FILES
C     FOR WASP DATA GROUPS B, C, AND D ARE:
C        DATA GROUP B USE WASPB.MRM  (DO NOT USE WASPB.OUT)
C        DATA GROUP C USE WASPC.OUT
C        DATA GROUP D USE WASPD.MRM  (DO NOT USE WASPD.OUT)
C===========
C
C **  SUBROUTINE WASP5 WRITES OUTPUT FILES PROVIDING ADVECTIVE AND
C **  DIFFUSIVE TRANSPORT FIELDS FOR THE WASP5 WATER QUALITY MODEL
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION LUTMP((KCM+1)*LCM),LDTMP((KCM+1)*LCM),QTMP((KCM+1)*LCM)
      CHARACTER*50 TITLEB,TITLEC
      CHARACTER*12 HYDFIL
C
      TITLEB='DATA GROUP B: EXCHANGE COEFFICIENTS'
      TITLEC='DATA GROUP C: VOLUMES'
C
C**********************************************************************C
C
C **  WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C **  THE VALUE OF X IN THE F10.X FORMATS MAY NEED TO BE CHANGED
C **  FROM PROBLEM TO PROBLEM.  A PRELIMINARY RUN USING E10.3
C **  CAN BE USED TO SPEED THE ADJUSTMENT  
C   
C**********************************************************************C
C
C **  READ CONTROL DATA FOR WRITING TO WASP COMPATIBLE FILES     
C
C----------------------------------------------------------------------C
C
      SVPT=1. 
      IF(NTSMMT.LT.NTSPTC)SVPT=0.
C
      IF(JSWASP.EQ.1)THEN
      OPEN(1,FILE='EFDC.WSP',STATUS='UNKNOWN')
C
C1**  READ CELL VOLUME PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
C
C2**  READ DIFFUSION PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) NRFLD,SCALR,CONVR,ISNKH
C
C3**  READ ADVECTION PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD
C
C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
C
      READ(1,1)
      READ(1,1)
      READ(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3                      !MRM
C
      CLOSE(1)
      ENDIF
C
    1 FORMAT (80X)
C
C**********************************************************************C
C
C **  WRITE HORIZONTAL POSITION AND LAYER FILE WASPP.OUT 
C **  WRITE INITIAL VOLUME FILE WASPC.OUT    
C
C **  FILE WASPC.OUT IS CONSISTENT WITH DATA GROUP C SPECIFICATIONS 
C **  ON PAGE 11 OF THE WASP5.1 MANUAL PART B, SEPT 1993
C   
C **  FILE WASPP.OUT DEFINES THE LAYER (1 IS SURFACE WATER LAYER, WITH
C **  LAYER NUMBERING INCREASING WITH DEPTH IN WATER COLUMN) AND
C **  HORIZONTAL POSITIONS IN LON,LAT OR UTME, UTMN OF THE WATER
C **  QUALITY (LONG TERM TRANSPORT) CELLS OR SEGEMENTS
C
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1)THEN
        OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
        OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
        CLOSE(90,STATUS='DELETE')
        CLOSE(93,STATUS='DELETE')
        OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
        OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
C       IVOPT=2
C       IBEDV=0
        WRITE(93,1031)IVOPT,IBEDV,TDINTS,TITLEC
C       SCALV=1.
C       CONVV=1.
        WRITE(93,1032)SCALV,CONVV
C       VMULT=0.
C       VEXP=0.
C       DMULT=0.
C       DEXP=0.
        LCLTM2=LCLT-2
        LWASP=0
        IF(KC.GT.1)THEN
          LTYPE=1
          KWASP=1
           DO LT=2,LALT
           LWASP=LWASP+1
           LBELOW=LWASP+LCLTM2
           I=ILLT(LT)
           J=JLLT(LT)
           L=LIJ(I,J)
           DMULT=HLPF(L)*DZC(KC)
           VOLUME=DXYP(L)*HLPF(L)*DZC(KC)
           IF(NTSMMT.LT.NTSPTC)THEN
             DMULT=HMP(L)*DZC(KC)
             VOLUME=DXYP(L)*HMP(L)*DZC(KC)
           ENDIF
           WRITE(90,1001)LWASP,KWASP,I,J,L,KC
           WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     &                DMULT,DEXP,I,J,L,KC
           ENDDO
          LTYPE=2
           DO K=KS,2,-1
           KWASP=KC-K+1
            DO LT=2,LALT
            LWASP=LWASP+1
            LBELOW=LWASP+LCLTM2
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            DMULT=HLPF(L)*DZC(K)                                             !MRM
            VOLUME=DXYP(L)*HLPF(L)*DZC(K)
            IF(NTSMMT.LT.NTSPTC)THEN
              DMULT=HMP(L)*DZC(KC)
              VOLUME=DXYP(L)*HMP(L)*DZC(KC)
            ENDIF
            WRITE(90,1001)LWASP,KWASP,I,J,L,K
            WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     &                DMULT,DEXP,I,J,L,KC
            ENDDO
           ENDDO
        ENDIF 
        LTYPE=2
        IF(KC.EQ.1) LTYPE=1
        KWASP=KC
        DO LT=2,LALT
         LWASP=LWASP+1
C        LBELOW=0
         LBELOW=LWASP+LCLTM2
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         DMULT=HLPF(L)*DZC(1)                                               !MRM
         VOLUME=DXYP(L)*HLPF(L)*DZC(1)
         IF(NTSMMT.LT.NTSPTC)THEN
              DMULT=HMP(L)*DZC(KC)
              VOLUME=DXYP(L)*HMP(L)*DZC(KC)
         ENDIF
         IONE=1
         WRITE(90,1001)LWASP,KWASP,I,J,L,IONE
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     &                DMULT,DEXP,I,J,L,IONE
        ENDDO
        LTYPE=3
        KWASP=KC+1
         DXYSUM=0.
         LWSPTMP=LWASP+1
         DO LT=2,LALT
         LWSPTMP=LWSPTMP+1
         ENDDO
C THE FOLLOWING THE LOWER BENTHIC LAYER.  ALL UPPER BENTHIC LAYER SEGMENTS   !MRM
C HAVE THIS LAYER IMMEDIATELY BELOW THEM:                                    !MRM
         DO LT=2,LALT
         LWASP=LWASP+1
         LBELOW=LWSPTMP
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         DXYSUM=DXYSUM+DXYP(L)
         VOLUME=DXYP(L)*DEPSED
         IZERO=0
         WRITE(90,1001)LWASP,KWASP,I,J,L,IZERO
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     &                DEPSED,DEXP,I,J,L,IZERO                                            !MRM
C     &                DMULT,DEXP                                            !MRM
         ENDDO
C NEXT DO THE LOWER BENTHIC LAYER:                                           !MRM
        LTYPE=4
        KWASP=KC+2
         LWASP=LWASP+1
         LBELOW=0
         DMULT=DEPSED                                                       !MRM
         VOLUME=DXYSUM*DEPSED
         IM1=-1
         WRITE(90,1001)LWASP,KWASP,I,J,L,IM1
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     &                DMULT,DEXP,I,J,L,IM1
        CLOSE(90)
        CLOSE(93)
      ENDIF
C
 1001 FORMAT(6I5,2F10.4)
 1031 FORMAT(2I5,F10.4,10X,A50)
 1032 FORMAT(2F10.4)
C FORMAT 1033 AS COMMENTED OUT IS TROUBLESOME ... BETTER CHANGE SHOWN
C1033 FORMAT(3I10,5F10.2)
 1033 FORMAT(3I10,F10.1,4F10.3,'   !',4I5)
C
C**********************************************************************C
C
C **  WRITE DIFFUSIVE AND DISPERSIVE TRANSPORT FILE WASPB.OUT 
C
C **  FILE WASPB.OUT IS CONSISTENT WITH DATA GROUP B SPECIFICATIONS          !MRM
C **  ON PAGE 8 OF THE WASP5.1 MANUAL PART B, SEPT 1993
C   
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1)THEN
        OPEN(91,FILE='WASPB.OUT',STATUS='UNKNOWN')
        CLOSE(91,STATUS='DELETE')
        OPEN(91,FILE='WASPB.OUT',STATUS='UNKNOWN')
C       NRFLD=1
        WRITE(91,1011)NRFLD,TITLEB
        NTEX=NTS/NTSMMT
C       SCALR=1.
C       CONVR=1.
        WRITE(91,1012)NTEX,SCALR,CONVR
        CLOSE(91)
C      ENDIF                                                                !MRM
C
      OPEN(91,FILE='WASPB.OUT',POSITION='APPEND',STATUS='UNKNOWN')
C
      LCLTM2=LCLT-2
      NORSH=0
      NORSV=0  
       DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         NORSH=NORSH+INT(SUBO(L))+INT(SVBO(L))
         NORSV=NORSV+INT(SPB(L))
       ENDDO
      NORS=ISNKH*KC*NORSH+KS*NORSV
      WRITE(91,1013)NORS
C
      IF(ISNKH.EQ.1)THEN
        UNITY=1.
        DO K=KC,1,-1
          KMUL=KC-K
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF(SUB(L).EQ.1.)THEN
              LWASP=LT-1+KMUL*LCLTM2
              LWASPW=LWASP-1       
              LW=L-1
              ADDLW=DYU(L)*AHULPF(L,K)*DZC(K)*0.5*(HLPF(L)
     &         +HLPF(LW))*DXIU(L)
              WRITE(91,1014) ADDLW,UNITY,LWASPW,LWASP
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      IF(ISNKH.EQ.1)THEN
        UNITY=1.
        DO K=KC,1,-1
          KMUL=KC-K
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF(SVB(L).EQ.1.)THEN
              LWASP=LT-1+KMUL*LCLTM2
              LSLT=LSCLT(LT)
              LWASPS=LSLT-1+KMUL*LCLTM2
              LS=LSC(L)
              ADDLS=DXV(L)*AHVLPF(L,K)*DZC(K)*0.5*(HLPF(L)
     &         +HLPF(LS))*DYIV(L)
              WRITE(91,1014) ADDLS,UNITY,LWASPS,LWASP
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      IF(KC.GT.1)THEN
        UNITY=1.
        DO K=KS,1,-1
          KMUL1=KS-K
          KMUL2=KMUL1+1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF(SPB(L).EQ.1.)THEN
              LWASP=LT-1+KMUL1*LCLTM2
              LBELOW=LT-1+KMUL2*LCLTM2
              ADDL=DXYP(L)*ABLPF(L,K)*DZIG(K)
              WRITE(91,1014) ADDL,UNITY,LWASP,LBELOW
            ENDIF 
          ENDDO
        ENDDO
      ENDIF 
C
      NBRK=6
      WRITE(91,1015)NBRK
C
C
      IF(ISDYNSTP.EQ.0)THEN
        TSTOP=DT*FLOAT(N)+TCON*TBEGIN
        TSTART=TSTOP-DT*FLOAT(NTSMMT)
      ELSE
        TSTOP=TENDRNSEC
        TSTART=TSTOP-DT*FLOAT(NTSMMT)
      ENDIF
C
CTIME      TSTOP=(DT*FLOAT(N)+TBEGIN*TCON)
CTIME      TSTART=TSTOP-DT*FLOAT(NTSMMT)
      TSTOP=TSTOP/86400.
      TSTART=TSTART/86400.
C
      TSMALL=1.E-5
      D1=0.
      T1=0.-2*TSMALL
      D2=0.
      T2=TSTART-TSMALL
      D3=1.
      T3=TSTART+TSMALL
      D4=1.
      T4=TSTOP-TSMALL
      D5=0.
      T5=TSTOP+TSMALL
      D6=0.
      T6=2*TSMALL+(DT*FLOAT(NTS)+TBEGIN*TCON)/86400.
      WRITE(91,1016)D1,T1,D2,T2,D3,T3,D4,T4
      WRITE(91,1016)D5,T5,D6,T6
C
      CLOSE(91)
C
C **  ADD PORE WATER EXCHANGE FIELD ON LAST CALL
C
C      IF(N.GE.NTS)THEN                                                    !MRM
        OPEN(91,FILE='WASPB.OUT',POSITION='APPEND',STATUS='UNKNOWN')
C
        NTEX=1
        SCALR=1.
        CONVR=1.
        WRITE(91,1012)NTEX,SCALR,CONVR
        NORSV=0  
        DO LT=2,LALT
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        NORSV=NORSV+INT(SPB(L))
        ENDDO
        WRITE(91,1013)NORSV
        IF(KC.GE.1)THEN
          KMUL2=KC+1
          UNITY=1.
          DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          IF(SPB(L).EQ.1.)THEN
            LWASP=LT-1+KC*LCLTM2
            LBELOW=LT-1+KMUL2*LCLTM2
            ADDL=2.*DXYP(L)*SEDIFF/DEPSED
            WRITE(91,1014) ADDL,UNITY,LWASP,LBELOW
          ENDIF 
          ENDDO
        ENDIF
C
        NBRK=6
        WRITE(91,1015)NBRK
C
C
      IF(ISDYNSTP.EQ.0)THEN
        TSTOP=DT*FLOAT(N)+TCON*TBEGIN
        TSTART=TSTOP-DT*FLOAT(NTSMMT)
      ELSE
        TSTOP=TENDRNSEC
        TSTART=TSTOP-DT*FLOAT(NTSMMT)
      ENDIF
C
CTIME      TSTOP=(DT*FLOAT(N)+TBEGIN*TCON)
CTIME      TSTART=TSTOP-DT*FLOAT(NTSMMT)
      TSTOP=TSTOP/86400.
      TSTART=TSTART/86400.
C
        TSMALL=1.E-5
        D1=0.
        T1=0.-2*TSMALL
        D2=0.
        T2=TSTART-TSMALL
        D3=1.
        T3=TSTART+TSMALL
        D4=1.
        T4=TSTOP-TSMALL
        D5=0.
        T5=TSTOP+TSMALL
        D6=0.
        T6=2*TSMALL+(DT*FLOAT(NTS)+TBEGIN*TCON)/86400.
        WRITE(91,1016)D1,T1,D2,T2,D3,T3,D4,T4
        WRITE(91,1016)D5,T5,D6,T6
C
        IBPTMP=0
        WRITE(91,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     &                IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     &                IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     &                IBPTMP,IBPTMP,IBPTMP,IBPTMP
C        
        CLOSE(91)
C      ENDIF                                                                !MRM
      ENDIF                                                                 !MRM
C
 1011 FORMAT(I5,10X,A50)
 1012 FORMAT(I5,2F10.4)
 1013 FORMAT(I5)
C1014 FORMAT(2F10.0,2I5,'   !',3I5,3X,A3)
 1014 FORMAT(2E10.3,2I5,F10.3,'   !',3I5,3X,A3)
 1015 FORMAT(I5)
 1016 FORMAT(4(E10.3,F10.5))
 1017 FORMAT(16I5)
C
C**********************************************************************C
C
C **  WRITE ADVECTIVE TRANSPORT FILE WASPD.OUT    
C
C **  FILE WASPD.OUT IS CONSISTENT WITH DATA GROUP D.1 SPECIFICATIONS 
C **  ON PAGE 13 OF THE WASP5.1 MANUAL PART B, SEPT 1993
C **  THIS FILE IS WRITTEN ONLY IF ISWASPD=1
C   
C----------------------------------------------------------------------C
C
C!!!!!!!!!!CHANGES ON NEXT 2 LINES
      IF(ISWASPD.EQ.1)THEN
C
       IF(JSWASP.EQ.1)THEN
        OPEN(92,FILE='WASPD.OUT',STATUS='UNKNOWN')
        CLOSE(92,STATUS='DELETE')
        OPEN(92,FILE='WASPD.OUT',STATUS='UNKNOWN')
C       IQOPT=1
C       NFIELD=1
        WRITE(92,1021)IQOPT,NFIELD,HYDFIL
        NINQ=NTS/NTSMMT
C       SCALQ=1
C       CONVQ=1
        WRITE(92,1022)NINQ,SCALQ,CONVQ
        CLOSE(92)
       ENDIF
C
      OPEN(92,FILE='WASPD.OUT',POSITION='APPEND',STATUS='UNKNOWN')
      LCLTM2=LCLT-2
      NOQSH=0
      NOQSV=0
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C!!!!!!!!!!!!!!CHANGES ON NEXT 3 LINES
       NOQSH=NOQSH+INT(SUBO(L))+INT(SVBO(L))
CJH       IF(IJCTLT(I+1,J).EQ.6) NOQSH=NOQSH+1
CJH       IF(IJCTLT(I,J+1).EQ.6) NOQSH=NOQSH+1
       IF(IJCTLT(I+1,J).EQ.8) NOQSH=NOQSH+1
       IF(IJCTLT(I,J+1).EQ.8) NOQSH=NOQSH+1
       NOQSV=NOQSV+INT(SWB(L))
       ENDDO
      NOQS=KC*NOQSH+KS*NOQSV
      WRITE(92,1023)NOQS
C
      LL=0
C
      DO K=KC,1,-1
      KMUL=KC-K
C
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C!!!!!!!!!!!!!CHANGES ON NEXT 15 LINES
       IF(SUBO(L).EQ.1.)THEN
         LL=LL+1
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LDTMP(LL)-1
CJH         IF(IJCTLT(I-1,J).EQ.6) LUTMP(LL)=0
         IF(IJCTLT(I-1,J).EQ.8) LUTMP(LL)=0
         QTMP(LL)=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
       ENDIF
CJH       IF(IJCTLT(I+1,J).EQ.6)THEN
       IF(IJCTLT(I+1,J).EQ.8)THEN
         IF(SUBO(L+1).EQ.1.)THEN
           LL=LL+1
           LDTMP(LL)=0
           LUTMP(LL)=LT-1+KMUL*LCLTM2
           QTMP(LL)=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
         ENDIF
       ENDIF
       ENDDO
C
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C!!!!!!!!!!!!!CHANGES ON NEXT 16 LINES
       IF(SVBO(L).EQ.1.)THEN
         LL=LL+1
         LSLT=LSCLT(LT)
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LSLT-1+KMUL*LCLTM2
CJH         IF(IJCTLT(I,J-1).EQ.6) LUTMP(LL)=0
         IF(IJCTLT(I,J-1).EQ.8) LUTMP(LL)=0
         QTMP(LL)=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
       ENDIF
CJH       IF(IJCTLT(I,J+1).EQ.6)THEN
       IF(IJCTLT(I,J+1).EQ.8)THEN
         LN=LNC(L)
         IF(SVBO(LN).EQ.1)THEN
           LL=LL+1
           LDTMP(LL)=0
           LUTMP(LL)=LT-1+KMUL*LCLTM2
           QTMP(LL)=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
         ENDIF
       ENDIF
       ENDDO
C
      ENDDO
C
      IF(KC.GT.1)THEN
        DO K=KS,1,-1
        KMUL1=KS-K
        KMUL2=KMUL1+1
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         IF(SWB(L).EQ.1.)THEN
           LL=LL+1
           LUTMP(LL)=LT-1+KMUL1*LCLTM2
           LDTMP(LL)=LT-1+KMUL2*LCLTM2
           QTMP(LL)=-DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
         ENDIF 
         ENDDO
        ENDDO
      ENDIF
C
      DO L=1,LL,4
      WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),
     &               QTMP(L+1),LUTMP(L+1),LDTMP(L+1),
     &               QTMP(L+2),LUTMP(L+2),LDTMP(L+2),
     &               QTMP(L+3),LUTMP(L+3),LDTMP(L+3)
      ENDDO
C
      NBRKQ=6
      WRITE(92,1025)NBRKQ
      WRITE(92,1026)D1,T1,D2,T2,D3,T3,D4,T4
      WRITE(92,1026)D5,T5,D6,T6
C
      CLOSE(92)
C!!!!!!!!!!CHANGES ON NEXT 2 LINES
C
      ENDIF
C
 1021 FORMAT(2I5,A12)
 1022 FORMAT(I5,2F10.4)
 1023 FORMAT(I5)
C1024 FORMAT(4(F10.0,2I5))
 1024 FORMAT(1P,4(E10.3,2I5))
 1025 FORMAT(I5)
 1026 FORMAT(4(2F10.5))

C**********************************************************************C     !MRM
C M.R. MORTON'S VERSION OF WASP DATA GROUP D                                 !MRM
C **  WRITE ADVECTIVE TRANSPORT FILE WASPD.MRM                               !MRM
C----------------------------------------------------------------------C     !MRM
      IF(JSWASP .EQ. 1)THEN                                                !MRM
        OPEN(92,FILE='WASPD.MRM',STATUS='UNKNOWN')                           !MRM
        WRITE(92,2020) IQOPT,NFIELD,HYDFIL                                   !MRM
        LL=0                                                                 !MRM
        NINQ=0                                                               !MRM
        SCALQ=1.0                                                            !MRM
        CONVQ=1.0/86400.0                                                    !MRM
C DATA BLOCK D.1 (ADVECTIVE FLOWS) IS NOT NEEDED SINCE HYD FILE IS USED:     !MRM
C        WRITE(92,2021) NINQ,SCALQ,CONVQ                                     !MRM
C DATA BLOCK D.2 (PORE WATER FLOWS) NOT NEEDED:                              !MRM
        WRITE(92,2022) NINQ,SCALQ,CONVQ                                      !MRM
C DATA BLOCK D.3 (SEDIMENT #1 TRANSPORT FIELD):                              !MRM
        NINQ=1                                                               !MRM
        WRITE(92,2023) NINQ,SCALQ,CONVQ                                      !MRM
        IF(KC.GT.1)THEN                                                    !MRM
          DO K=KS,0,-1                                                       !MRM
            KMUL1=KS-K                                                       !MRM
            KMUL2=KMUL1+1                                                    !MRM
            DO LT=2,LALT                                                     !MRM
C              WRITE(6,6999)K,KMUL1,KMUL2,LT
C              CALL F_FLUSHNOW(6) 
              I=ILLT(LT)                                                     !MRM
C              WRITE(6,6999)LT,I 
C              CALL F_FLUSHNOW(6) 
              J=JLLT(LT)                                                     !MRM
C              WRITE(6,6999)LT,J
C              CALL F_FLUSHNOW(6) 
              L=LIJ(I,J)                                                     !MRM
C              WRITE(6,6999)K,KMUL1,KMUL2,LT,I,J,L,IL(L),JL(L),SWB(L) 
C              CALL F_FLUSHNOW(6) 
              IF(SWB(L).EQ.1.)THEN                                         !MRM
                LL=LL+1                                                      !MRM
                LUTMP(LL)=LT-1+KMUL1*LCLTM2                                  !MRM
                LDTMP(LL)=LT-1+KMUL2*LCLTM2                                  !MRM
C QTMP ARRAY WILL HOLD THE PLAN VIEW AREA OF EACH CELL:                      !MRM
                QTMP(LL)= DXYP(L)                                            !MRM
              ENDIF                                                         !MRM
            ENDDO                                                           !MRM
          ENDDO                                                             !MRM
        ENDIF                                                               !MRM
C                                                                            !MRM                                               
 6999 FORMAT(9I5,F5.1)                                                 
 6996 FORMAT(9I5,F5.1)                                                 
C                                                               
        WRITE(92,2030) LL                                                    !MRM
        DO L=1,LL,4                                                          !MRM
          WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),                      !MRM
     &                 QTMP(L+1),LUTMP(L+1),LDTMP(L+1),                      !MRM
     &                 QTMP(L+2),LUTMP(L+2),LDTMP(L+2),                      !MRM
     &                 QTMP(L+3),LUTMP(L+3),LDTMP(L+3)                       !MRM
        ENDDO                                                               !MRM
        NBRKQ=2                                                              !MRM
        T1=1.0                                                               !MRM
        T2=366.0                                                             !MRM
        WRITE(92,2030) NBRKQ                                                 !MRM
        WRITE(92,2031) WSS1,T1,WSS1,T2                                       !MRM
C DATA BLOCK D.4 (SEDIMENT #2 TRANSPORT FIELD):                              !MRM
        NINQ=1                                                               !MRM
        WRITE(92,2024) NINQ,SCALQ,CONVQ                                      !MRM
        WRITE(92,2030) LL                                                    !MRM
        DO L=1,LL,4                                                          !MRM
          WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),                      !MRM
     &                 QTMP(L+1),LUTMP(L+1),LDTMP(L+1),                      !MRM
     &                 QTMP(L+2),LUTMP(L+2),LDTMP(L+2),                      !MRM
     &                 QTMP(L+3),LUTMP(L+3),LDTMP(L+3)                       !MRM
        ENDDO                                                               !MRM
        NBRKQ=2                                                              !MRM
        T1=1.0                                                               !MRM
        T2=366.0                                                             !MRM
        WRITE(92,2030) NBRKQ                                                 !MRM
        WRITE(92,2031) WSS2,T1,WSS2,T2                                       !MRM
C DATA BLOCK D.5 (SEDIMENT #3 TRANSPORT FIELD):                              !MRM
        NINQ=1                                                               !MRM
        WRITE(92,2025) NINQ,SCALQ,CONVQ                                      !MRM
        WRITE(92,2030) LL                                                    !MRM
        DO L=1,LL,4                                                          !MRM
          WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),                      !MRM
     &                 QTMP(L+1),LUTMP(L+1),LDTMP(L+1),                      !MRM
     &                 QTMP(L+2),LUTMP(L+2),LDTMP(L+2),                      !MRM
     &                 QTMP(L+3),LUTMP(L+3),LDTMP(L+3)                       !MRM
        ENDDO                                                               !MRM
        NBRKQ=2                                                              !MRM
        T1=1.0                                                               !MRM
        T2=366.0                                                             !MRM
        WRITE(92,2030) NBRKQ                                                 !MRM
        WRITE(92,2031) WSS3,T1,WSS3,T2                                       !MRM
C ADD SYSTEM BYPASS ARRAY TO BOTTOM OF DATA GROUP D:                         !MRM
        WRITE(92,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP,                           !MRM
     +                IBPTMP,IBPTMP,IBPTMP,IBPTMP,                           !MRM
     +                IBPTMP,IBPTMP,IBPTMP,IBPTMP,                           !MRM
     +                IBPTMP,IBPTMP,IBPTMP,IBPTMP                            !MRM
        CLOSE(92)                                                            !MRM
      ENDIF                                                                 !MRM
 2020 FORMAT(2I5,A12,'    DATA GROUP D: FLOWS')                              !MRM
 2021 FORMAT(1P,I5,2E10.3,'    DATA BLOCK D.1 ADVECTIVE FLOWS')              !MRM
 2022 FORMAT(1P,I5,2E10.3,'    DATA BLOCK D.2 PORE WATER FLOWS')             !MRM
 2023 FORMAT(1P,I5,2E10.3,'    DATA BLOCK D.3 SED. #1 TRANSPORT FIELD')      !MRM
 2024 FORMAT(1P,I5,2E10.3,'    DATA BLOCK D.4 SED. #2 TRANSPORT FIELD')      !MRM
 2025 FORMAT(1P,I5,2E10.3,'    DATA BLOCK D.5 SED. #3 TRANSPORT FIELD')      !MRM
 2030 FORMAT(I5)                                                             !MRM
 2031 FORMAT(2(E10.3,F10.5))                                                 !MRM

C
C**********************************************************************C
C
C **  WRITE TO EXTERNAL HYDRO FILE WASPDH.OUT AND DIAGNOSTIC VERSION
C **  OF SAME FILE WASPDHD.OUT 
C
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1)THEN
        OPEN(90,FILE='WASPDHD.OUT',STATUS='UNKNOWN')
        IF(IQOPT.EQ.3) OPEN(94,FILE='WASPDH.OUT',STATUS='UNKNOWN')           !MRM
        IF(IQOPT.EQ.4) OPEN(95,FILE='WASPDHU.OUT',STATUS='UNKNOWN',          !MRM
     &                      FORM='UNFORMATTED')
        OPEN(96,FILE='WASPB.MRM',STATUS='UNKNOWN')                           !MRM
        CLOSE(90,STATUS='DELETE')
        IF(IQOPT.EQ.3) CLOSE(94,STATUS='DELETE')                             !MRM
        IF(IQOPT.EQ.4) CLOSE(95,STATUS='DELETE')                             !MRM
        CLOSE(96,STATUS='DELETE')                                            !MRM
        OPEN(90,FILE='WASPDHD.OUT',STATUS='UNKNOWN')
        IF(IQOPT.EQ.3) OPEN(94,FILE='WASPDH.OUT',STATUS='UNKNOWN')           !MRM
        IF(IQOPT.EQ.4) OPEN(95,FILE='WASPDHU.OUT',STATUS='UNKNOWN',          !MRM
     &                      FORM='UNFORMATTED')
        OPEN(96,FILE='WASPB.MRM',STATUS='UNKNOWN')                           !MRM
        WRITE(96,1011) NRFLD,TITLEB                                          !MRM
        NTEXX=1                                                              !MRM
        WRITE(96,1012) NTEXX,SCALR,CONVR                                     !MRM

C WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 1, DATA OPTIONS:                 !MRM
C  NJUN = NUMBER OF SEGMENTS CONNECTED BY FLOWS FROM THE HYD. FILE           !MRM
C  NCHN = NUMBER OF INTERFACIAL FLOW PAIRS FROM THE HYD. FILE                !MRM
C  DTWASP = WASP5 TIME STEP (SECONDS)                                        !MRM
C  TZERO = BEGIN TIME STEP FOR HYD. FILE (SECONDS)                           !MRM
C  TENDHYD = END TIME STEP FOR HYD. FILE (SECONDS)                           !MRM
C  ISTMP = CONTROL SWITCH, 0=TIME VARIABLE SEGMENT DEPTHS AND VELOCITIES     !MRM
C          ARE READ; 1=TIME VARIABLE SEGMENT DEPTHS AND VELOCITIES ARE NOT   !MRM
C          READ.                                                             !MRM
C
C       KCLC=KC*LCLT
C       LCLTM2=LCLT-2
C        DO KL=1,KCLC 
C        NCHNC(KL)=0
C        ENDDO
C        DO M=1,10
C         DO KL=1,KCLC 
C         LCHNC(KL,M)=0
C         ENDDO
C        ENDDO
        NJUN=KC*(LCLT-2)
        NCHNH=0
        NCHNV=0
C!!!!!!!!CHANGES NEXT 13 LINES
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          NCHNH=NCHNH+INT(SUBO(L))
CJH          IF(IJCTLT(I+1,J).EQ.6)THEN
          IF(IJCTLT(I+1,J).EQ.8)THEN
            IF(SUBO(L+1).EQ.1.) NCHNH=NCHNH+1
          ENDIF 
          NCHNH=NCHNH+INT(SVBO(L))
CJH          IF(IJCTLT(I,J+1).EQ.6)THEN
          IF(IJCTLT(I,J+1).EQ.8)THEN
            IF(SVBO(LNC(L)).EQ.1.) NCHNH=NCHNH+1
          ENDIF 
          NCHNV=NCHNV+INT(SWB(L))
        ENDDO
        NCHN=KC*NCHNH+KS*NCHNV
        ISTMP=0
        NODYN=NFLTMT
        NODYN=NODYN                                                          !MRM
        DTWASP = DT * FLOAT(NTSMMT)                                          !MRM
        TZERO=TBEGIN*TCON
        TENDHYD=TZERO+NTS*DT
        WRITE(90,901)NJUN,NCHN
        IF(IQOPT.EQ.3)THEN                                                  !MRM
          WRITE(94,941) NJUN,NCHN, DTWASP, TZERO,TENDHYD,ISTMP               !MRM
C        WRITE(94,941) NJUN,NCHN,DT,TZERO,TENDHYD,ISTMP
        ENDIF                                                               !MRM
        IF(IQOPT.EQ.4)THEN                                                  !MRM
          WRITE(95) NJUN,NCHN, DTWASP, TZERO,TENDHYD,ISTMP                   !MRM
        ENDIF
        WRITE(96,1013) NCHN                                                  !MRM
C
C **  CHANNEL DATA
C
C WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 2, SEGMENT INTERFACE PAIRS:      !MRM
C   WASP EXPECTS TO SEE BOUNDARY SEGMENTS DESIGNATED AS "0".                 !MRM
C
        RMNDUM=0.
        LCHN=0
        DO K=KC,1,-1
          KMUL=KC-K
C!!!!!!!!!!!!!!CHANGES ON NEXT 38 LINES
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF(SUBO(L).EQ.1.)THEN
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LDTM-1
CJH              IF(IJCTLT(I-1,J).EQ.6) LUTM=0
              IF(IJCTLT(I-1,J).EQ.8) LUTM=0
              RLENTH=DXU(L)
              WIDTH=DYU(L)
              LCHN=LCHN+1
C             NCHNC(LDTM)=NCHNC(LDTM)+1
C             NCHNC(LUTM)=NCHNC(LUTM)+1
C             LCHNC(LDTM,NCHNC(LDTM))=LCHN
C             LCHNC(LUTM,NCHNC(LUTM))=LCHN
              IF(ISDHD. EQ. 1) WRITE(90,902)LCHN,RLENTH,WIDTH,
     +          RMNDUM,LUTM,LDTM
              IF(ISDHD .EQ. 2) WRITE(90,'(2I5)') LUTM,LDTM
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              IF(IQOPT.EQ.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'U 0'                         !MRM
              IF(IQOPT.EQ.4) WRITE(95) LUTM,LDTM                             !MRM
              WRITE(96,1014) UNITY,UNITY,LUTM,LDTM,UNITY,I,J,K,'U 0'                     !MRM
            ENDIF
CJH            IF(IJCTLT(I+1,J).EQ.6)THEN
            IF(IJCTLT(I+1,J).EQ.8)THEN
              IF(SUBO(L+1).EQ.1.)THEN
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                RLENTH=DXU(L+1)
                WIDTH=DYU(L+1)
                LCHN=LCHN+1
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                IF(ISDHD .EQ. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +            RMNDUM,LUTM,LDTM
              IF(ISDHD .EQ. 2) WRITE(90,'(2I5)') LUTM,LDTM
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                IF(IQOPT.EQ.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'U+1'                       !MRM
                IF(IQOPT.EQ.4) WRITE(95) LUTM,LDTM                           !MRM
                UNITY=1.0
                WRITE(96,1014) UNITY,UNITY,LUTM,LDTM,UNITY,I,J,K,'U+1'                   !MRM
              ENDIF
            ENDIF
          ENDDO
C        ENDDO
C        DO K=KC,1,-1
C        KMUL=KC-K
C!!!!!!!!CHANGES NEXT 41 LINES
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF(SVBO(L).EQ.1.)THEN
              LSLT=LSCLT(LT)
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LSLT-1+KMUL*LCLTM2
CJH              IF(IJCTLT(I,J-1).EQ.6) LUTM=0
              IF(IJCTLT(I,J-1).EQ.8) LUTM=0
              RLENTH=DYV(L)
              WIDTH=DXV(L)
              LCHN=LCHN+1
C             NCHNC(LDTM)=NCHNC(LDTM)+1
C             NCHNC(LUTM)=NCHNC(LUTM)+1
C             LCHNC(LDTM,NCHNC(LDTM))=LCHN
C             LCHNC(LUTM,NCHNC(LUTM))=LCHN
              IF(ISDHD .EQ. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +          RMNDUM,LUTM,LDTM
              IF(ISDHD .EQ. 2) WRITE(90,'(2I5)') LUTM,LDTM
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              IF(IQOPT.EQ.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'V 0'                         !MRM
              IF(IQOPT.EQ.4) WRITE(95) LUTM,LDTM                             !MRM
              WRITE(96,1014) UNITY,UNITY,LUTM,LDTM,UNITY,I,J,K,'V 0'                     !MRM
            ENDIF
C            IF(IJCTLT(I,J+1).EQ.6)THEN
            IF(IJCTLT(I,J+1).EQ.8)THEN
              LN=LNC(L)
              IF(SVBO(LN).EQ.1.)THEN
                LSLT=LSCLT(LT)
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                RLENTH=DYV(LN)
                WIDTH=DXV(LN)
                LCHN=LCHN+1
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                IF(ISDHD .EQ. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +            RMNDUM,LUTM,LDTM
                IF(ISDHD .EQ. 2) WRITE(90,'(2I5)') LUTM,LDTM
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                IF(IQOPT.EQ.3) WRITE(94,9941) LUTM,LDTM,UNITY,I,
     +                                       J,K,'V+1'                       !MRM
                IF(IQOPT.EQ.4) WRITE(95) LUTM,LDTM                           !MRM
                WRITE(96,1014) UNITY,UNITY,LUTM,LDTM,UNITY,I,J,K,'V+1'                   !MRM
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(KC.GT.1)THEN
          DO K=KS,1,-1
            KMUL1=KS-K
            KMUL2=KMUL1+1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF(SWB(L).EQ.1.)THEN
                LUTM=LT-1+KMUL1*LCLTM2
                LDTM=LT-1+KMUL2*LCLTM2
                RLENTH=HLPF(L)*DZG(K)
                WIDTH=SQRT(DXYP(L))
                LCHN=LCHN+1
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                WRITE(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                IF(IQOPT.EQ.3) WRITE(94,9941)LUTM,LDTM,I,J,K,'W 0'                        !MRM
                IF(IQOPT.EQ.4) WRITE(95) LUTM,LDTM                           !MRM
                WRITE(96,1014) UNITY,UNITY,LUTM,LDTM,UNITY,I,J,K,'W 0'                   !MRM
              ENDIF
            ENDDO
          ENDDO
C
C WRITE OUT TIME SERIES OF ZERO DISPERSION COEFFICIENTS:
C
          D1=0.0                                                             !MRM
          T1=TZERO/TCON                                                      !MRM
          D2=0.0                                                             !MRM
          T2=TENDHYD/TCON                                                    !MRM
          NBRKQ=2                                                            !MRM
          WRITE(96,905) NBRKQ                                                !MRM
          WRITE(96,1016) D1,T1, D2,T2                                        !MRM
C
C FOR EXCHANGE BETWEEN THE LOWER WATER SURFACE LAYER AND THE UPPER           !MRM
C BENTHIC LAYER, DO THE FOLLOWING:                                           !MRM
C
          WRITE(96,1012) NTEXX,SCALR,CONVR                                   !MRM
          NTEXX=0
          DO K=1,1                                                           !MRM
            DO LT=2,LALT                                                     !MRM
              I=ILLT(LT)                                                     !MRM
              J=JLLT(LT)                                                     !MRM
              L=LIJ(I,J)                                                     !MRM
              IF(SWB(L).EQ.1.)THEN                                         !MRM
                NTEXX=NTEXX+1                                                !MRM
              ENDIF
            ENDDO
          ENDDO
          WRITE(96,1013) NTEXX                                               !MRM
          DO K=1,1                                                           !MRM
            KMUL1=KS-K                                                       !MRM
            KMUL2=KMUL1+1                                                    !MRM
            KMUL3=KMUL2+1                                                    !MRM
            DO LT=2,LALT                                                     !MRM
              I=ILLT(LT)                                                     !MRM
              J=JLLT(LT)                                                     !MRM
              L=LIJ(I,J)                                                     !MRM
              IF(SWB(L).EQ.1.)THEN                                         !MRM
                LUTM=LT-1+KMUL2*LCLTM2                                       !MRM
                LDTM=LT-1+KMUL3*LCLTM2                                       !MRM
                WRITE(96,1014) DXYP(L),DEPSED,LUTM,LDTM                      !MRM
              ENDIF
            ENDDO
          ENDDO
C
C WRITE OUT TIME SERIES OF WATER-BENTHIC EXCHANGE DISPERSION COEFFICIENTS:   !MRM
C
          D1=SEDIFF                                                          !MRM
          T1=TZERO/TCON                                                      !MRM
          D2=SEDIFF                                                          !MRM
          T2=TENDHYD/TCON                                                    !MRM
          NBRKQ=2                                                            !MRM
          WRITE(96,905) NBRKQ                                                !MRM
          WRITE(96,1016) D1,T1, D2,T2                                        !MRM
          IBPTMP=0                                                           !MRM
          WRITE(96,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP,                         !MRM
     +                  IBPTMP,IBPTMP,IBPTMP,IBPTMP,                         !MRM
     +                  IBPTMP,IBPTMP,IBPTMP,IBPTMP,                         !MRM
     +                  IBPTMP,IBPTMP,IBPTMP,IBPTMP                          !MRM
        ENDIF
C
C **  JUNCTION DATA WITH INITIAL CONDITIONS
C
C WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 3, INITIAL SEGMENT PROPERTIES:   !MRM
C
        VELTMP=0.
        DUMVOL=0.
        DO K=KC,1,-1
        KMUL=KC-K
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         LCELL=LT-1+KMUL*LCLTM2
         DEPTMP=HLPF(L)*DZC(K)
         VOLTMP=DEPTMP*DXYP(L)
         IF(NTSMMT.LT.NTSPTC)THEN
           DEPTMP=HMP(L)*DZC(K)
           VOLTMP=DEPTMP*DXYP(L)
         ENDIF
         IF(ISDHD .EQ. 1) WRITE(90,904) LCELL,VOLTMP,I,J
         IF(IQOPT.EQ.3) WRITE(94,9440) VOLTMP,DEPTMP,VELTMP           !MRM
         IF(IQOPT.EQ.4) WRITE(95) VOLTMP,DEPTMP,VELTMP                       !MRM
         ENDDO
        ENDDO
C
        CLOSE(90)
        IF(IQOPT.EQ.3) CLOSE(94)                                             !MRM
        IF(IQOPT.EQ.4) CLOSE(95)                                             !MRM
        CLOSE(96)                                                            !MRM
      ENDIF
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP, VOLUME AND FLOW DATA
C
      OPEN(90,FILE='WASPDHD.OUT',POSITION='APPEND',STATUS='UNKNOWN')
      IF(IQOPT.EQ.3)THEN                                                    !MRM
        OPEN(94,FILE='WASPDH.OUT',POSITION='APPEND',STATUS='UNKNOWN')
      ENDIF                                                                 !MRM
      IF(IQOPT.EQ.4)THEN                                                    !MRM
        OPEN(95,FILE='WASPDHU.OUT',POSITION='APPEND',STATUS='UNKNOWN',
     &                      FORM='UNFORMATTED')
      ENDIF                                                                 !MRM
      LCLTM2=LCLT-2
      IZERO=0
      RZERO=0
      IZERO=IZERO                                                            !MRM
      RZERO=RZERO                                                            !MRM
C
C     NSTEP=N-NTSMMT
C     WRITE(94,945)NSTEP
C
C WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 4, BQ(J) FLOW IN INTERFACE       !MRM
C PAIR "J":                                                                  !MRM
C
C ADVECTION AND DISPERSION IN THE X-DIRECTION:
C
      LCHNUM=0
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:        !MRM
       ADDLW=0.0                                                             !MRM
       IF(SUB(L).EQ.1.)THEN                                                !MRM
         LW=L-1                                                              !MRM
         ADDLW=DYU(L)*AHULPF(L,K)*DZC(K)*0.5*(HLPF(L)                        !MRM
     &         +HLPF(LW))*DXIU(L)                                            !MRM
       ENDIF                                                                !MRM
C +++++ ABOVE ADDED BY M. MORTON                                             !MRM
C!!!!!!!!CHANGES NEXT 12 LINES
       IF(SUBO(L).EQ.1.)THEN
         TMPVAL=UHLPF(L,K)+SVPT*UVPT(L,K)
         FLOWX=DYU(L)*TMPVAL*DZC(K)
         UDDXTMP=2.*TMPVAL*DXIU(L)/(HLPF(L)+HLPF(L-1))
         IMTMP=I-1
         LCHNUM=LCHNUM+1
         IDRTMP=1
         IF(ISDHD .EQ. 1) WRITE(90,944) FLOWX,IMTMP,I,J,K
         IF(IQOPT.EQ.3) WRITE(94,9946) FLOWX,UDDXTMP,ADDLW,IDRTMP                    !MRM/JMH
C         WRITE(94,946) FLOWX                                                !MRM
         IF(IQOPT.EQ.4) WRITE(95) FLOWX,UDDXTMP,ADDLW,IDRTMP                        !MRM/JMH
       ENDIF
CJH       IF(IJCTLT(I+1,J).EQ.6)THEN
       IF(IJCTLT(I+1,J).EQ.8)THEN
         IF(SUBO(L+1).EQ.1.)THEN
           TMPVAL=UHLPF(L+1,K)+SVPT*UVPT(L+1,K)
           FLOWX=DYU(L+1)*TMPVAL*DZC(K)
           UDDXTMP=2.*TMPVAL*DXIU(L+1)/(HLPF(L+1)+HLPF(L))
           IPTMP=I+1
           LCHNUM=LCHNUM+1
           IDRTMP=1
           IF(ISDHD .EQ. 1) WRITE(90,944) LCHNUM,FLOWX,I,IPTMP,J,K
           IF(IQOPT.EQ.3) WRITE(94,9946) FLOWX,UDDXTMP,ADDLW,IDRTMP                  !MRM/JMH
C           WRITE(94,946) FLOWX                                              !MRM
           IF(IQOPT.EQ.4) WRITE(95) FLOWX,UDDXTMP,ADDLW,IDRTMP                      !MRM/JMH
         ENDIF
       ENDIF  
       ENDDO
C     ENDDO
C
C ADVECTION AND DISPERSION IN THE Y-DIRECTION:
C
C     DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:        !MRM
       ADDLS=0.0                                                             !MRM
       IF(SVB(L).EQ.1.)THEN                                                !MRM
         LS=LSC(L)                                                           !MRM
         ADDLS=DXV(L)*AHVLPF(L,K)*DZC(K)*0.5*(HLPF(L)                        !MRM
     &         +HLPF(LS))*DYIV(L)                                            !MRM
       ENDIF                                                                !MRM
C +++++ ABOVE ADDED BY M. MORTON                                             !MRM
C!!!!!!!CHANGES NEXT 13 LINES
       IF(SVBO(L).EQ.1.)THEN
         TMPVAL=VHLPF(L,K)+SVPT*VVPT(L,K)
         FLOWY=DXV(L)*TMPVAL*DZC(K)
         VDDYTMP=2.*TMPVAL*DYIV(L)/(HLPF(L)+HLPF(LSC(L)))
         JMTMP=J-1
         LCHNUM=LCHNUM+1
         IDRTMP=2
         IF(ISDHD .EQ. 1) WRITE(90,944) LCHNUM,FLOWY,I,JMTMP,J,K
         IF(IQOPT.EQ.3) WRITE(94,9946) FLOWY,VDDYTMP,ADDLS,IDRTMP                          !MRM
C         WRITE(94,946) FLOWY                                               !MRM
         IF(IQOPT.EQ.4) WRITE(95) FLOWY,VDDYTMP,ADDLS,IDRTMP                              !MRM
       ENDIF
CJH       IF(IJCTLT(I,J+1).EQ.6)THEN
       IF(IJCTLT(I,J+1).EQ.8)THEN
         LN=LNC(L)
         IF(SVBO(LN).EQ.1.)THEN
           TMPVAL=VHLPF(LN,K)+SVPT*VVPT(LN,K)
           FLOWY=DXV(LN)*TMPVAL*DZC(K)
           VDDYTMP=2.*TMPVAL*DYIV(LN)/(HLPF(LN)+HLPF(L))
           JPTMP=J+1
           LCHNUM=LCHNUM+1
           IDRTMP=2
           IF(ISDHD .EQ. 1) WRITE(90,944) LCHNUM,FLOWY,I,J,JPTMP,K
           IF(IQOPT.EQ.3) WRITE(94,9946) FLOWY,VDDYTMP,ADDLS,IDRTMP                          !MRM
C           WRITE(94,946) FLOWY                                              !MRM
           IF(IQOPT.EQ.4) WRITE(95) FLOWY,VDDYTMP,ADDLS,IDRTMP                              !MRM
         ENDIF
       ENDIF 
       ENDDO
      ENDDO
C
C ADVECTION AND DISPERSION IN THE Z-DIRECTION:
C
      IF(KC.GT.1)THEN
      DO K=KS,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:        !MRM
       ADDL=0.0                                                              !MRM
       IF(SPB(L).EQ.1.)THEN                                                !MRM
         ADDL=DXYP(L)*ABLPF(L,K)*DZIG(K)                                     !MRM
         IF(ISDHD .EQ. 2) WRITE(90, '(4I5,E13.4)') I,J,K,L,ABLPF(L,K)
       ENDIF                                                                !MRM
C +++++ ABOVE ADDED BY M. MORTON                                             !MRM
       IF(SWB(L).EQ.1)THEN
         TMPVAL=WLPF(L,K)+SVPT*WVPT(L,K)
         FLOWZ=-DXYP(L)*TMPVAL
         WDDZTMP=TMPVAL*DZIG(K)/HLPF(L)
         KPTMP=K+1
         IDRTMP=3
         LCHNUM=LCHNUM+1
         IF(ISDHD .EQ. 1) WRITE(90,944) LCHNUM,FLOWZ,I,J,K,KPTMP
         IF(IQOPT.EQ.3) WRITE(94,9946) FLOWZ,WDDZTMP,ADDL,IDRTMP
C         WRITE(94,946) FLOWZ
         IF(IQOPT.EQ.4) WRITE(95) FLOWZ,WDDZTMP,ADDL,IDRTMP                                !MRM
       ENDIF  
       ENDDO
      ENDDO
      ENDIF
C
C WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 5, SEGMENT PROPERTIES:           !MRM
C
      QQSUM=0.
      LCELTMP=0
      DO K=KC,1,-1
       DO LT=2,LALT
       LCELTMP=LCELTMP+1
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       LN=LNC(L)
       VOLUM=DXYP(L)*HLPF(L)*DZC(K)
       IF(NTSMMT.LT.NTSPTC) VOLUM=DXYP(L)*HP(L)*DZC(K)
C      QIN=QSUMELPF(L)*DZC(K)
C      FLOWXI=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
C      FLOWYI=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
C      FLOWZI=DXYP(L)*(WLPF(L,K-1)+SVPT*WVPT(L,K-1))
C      FLOWXO=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
C      FLOWYO=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
C      FLOWZO=DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
C      QQSUM=QIN+FLOWXI+FLOWYI+FLOWZI-FLOWXO-FLOWYO-FLOWZO
       DEPTH=HLPF(L)*DZC(K)
C       IF(NTSMMT.LT.NTSPTC) DEPTH=HP(L)*DZC(K)
       VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K)
     &         +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
       VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K)
     &         +VHLPF(LN,K)+SVPT*VVPT(LN,K))/HLPF(L)
       VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1)
     &         +WLPF(L,K)+SVPT*WVPT(L,K))
       VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
       IF(ISDHD .EQ. 1) WRITE(90,902) LCELTMP,VOLUM,I,J,K
C       IF(IQOPT.EQ.3) WRITE(94,946) VOLUM,QQSUM,DEPTH,VELMAG                 !MRM
       IF(IQOPT.EQ.3) WRITE(94,946) VOLUM,DEPTH,VELMAG                        !MRM
       IF(IQOPT.EQ.4) WRITE(95) VOLUM, DEPTH, VELMAG                         !MRM
       ENDDO
      ENDDO
C
      CLOSE(90)
      IF(IQOPT.EQ.3) CLOSE(94)                                               !MRM
      IF(IQOPT.EQ.4) CLOSE(95)                                               !MRM
C
C----------------------------------------------------------------------C
C
  901 FORMAT(2I5,E12.5,4I5,E12.5)
  902 FORMAT(I5,2X,3F20.8,3I5)
  903 FORMAT(3E12.5,2I5)
  904 FORMAT(I5,2X,F20.8,10I5)
  905 FORMAT(I5)
  906 FORMAT(5E12.5)
  941 FORMAT(2I5,3F20.8,I5)
  942 FORMAT(3E12.5,2I5)
  943 FORMAT(3E12.5,2I5)
  944 FORMAT(I5,2X,F20.8,10I5)
 9440 FORMAT(4F20.8)
  945 FORMAT(I5)
  946 FORMAT(4E17.9)
 9946 FORMAT(3E17.9,I5)
 9941 FORMAT(2I5,'    !',3I5,3X,A3)
C
C**********************************************************************C
C
      JSWASP=0
C
      RETURN
      END
