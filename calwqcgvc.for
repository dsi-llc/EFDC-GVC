C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALWQCGVC (ISTL)
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
C **  SUBROUTINE CALWQC CALCULATES THE CONCENTRATION OF DISSOLVED AND 
C **  SUSPENDED WATER QUALITY CONSTITUTENTS AT TIME LEVEL (N+1). 
C **  CALLED ONLY ON ODD THREE TIME LEVEL STEPS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C     DIMENSION HWQI(LCM)
CDHP  DIMENSION CUBTMP(LCM),CMBTMP(LCM),CLBTMP(LCM),EB(LCM),
CDHP &          VTMP(LCM),ABHWQI(LCM,KCM)
C
C**********************************************************************C
C
      DELT=DT2
c
      IF(IS2TIM.GE.1) THEN
C	  ISTL=2
        IF(ISDYNSTP.EQ.0)THEN
          DELT=DT
          S3TL=0.0
          S2TL=1.0
          ISUD=0
         ELSE
          DELT=DTDYN
          DELTA=DTDYN
          S3TL=0.0
          S2TL=1.0
          ISUD=0
        END IF
	ENDIF
c
C**********************************************************************C
C
C **  UPDATED TIME SERIES CONCENTRATION BOUNDARY CONDITIONS
C
C     CALL CALWQS(ISTL)
C      
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION 
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
      ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
      NVQVT=NWQV
      IF(IWQFCB.EQ.0) NVQVT=NWQV-1
C
      IF(ISTRAN(8).EQ.1)THEN 
        DO NW=1,NVQVT
         IF(ISTRWQ(NW).EQ.1)THEN
           DO K=1,KC
           DO L=2,LA
            CWQ(L,K)=WQV(L,K,NW)
            CWQ2(L,K)=WQV(L,K,NW)
           ENDDO
           ENDDO
           CALL CALTRWQGVC (8,NW,CWQ,CWQ2)
           DO K=1,KC
           DO L=2,LA
            WQV(L,K,NW)=CWQ(L,K)
           ENDDO
           ENDDO
         ENDIF       
        ENDDO
      ENDIF       
C
      IF(ISCRAY.EQ.0)THEN
        TWQADV=TWQADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TWQADV=TWQADV+T2TMP-T1TMP
        WTWQADV=WTWQADV+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
C     IF(ISTRAN(8).EQ.2) CALL CALTRWQI (ISTL,8,CWQ)
C
C**********************************************************************C
C
C **  CALLS TO SOURCE-SINK CALCULATIONS
C
C**********************************************************************C
C
C **  BYPASS OR INITIALIZE VERTICAL DIFFUSION CALCULATION
C
C----------------------------------------------------------------------C
C
      IF(KC.EQ.1) GOTO 2000
C
      DO L=2,LA
      HWQI(L)=1./HWQ(L)
      ENDDO
C
      IF(ISCRAY.EQ.0)THEN
        TTMP=SECNDS(0.0)
      ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      ENDIF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION LEVEL 1
C
C----------------------------------------------------------------------C
C
      IF(ISWQLVL.EQ.1)THEN
C 
      RCDZKK=-DELT*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,1)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB
C       CWQ(L,1)=CWQ(L,1)*EEB
        WQV(L,1, 1)=WQV(L,1, 1)*EEB
        WQV(L,1, 2)=WQV(L,1, 2)*EEB
        WQV(L,1, 3)=WQV(L,1, 3)*EEB
        WQV(L,1, 4)=WQV(L,1, 4)*EEB
        WQV(L,1, 5)=WQV(L,1, 5)*EEB
        WQV(L,1, 6)=WQV(L,1, 6)*EEB
        WQV(L,1, 7)=WQV(L,1, 7)*EEB
        WQV(L,1, 8)=WQV(L,1, 8)*EEB
        WQV(L,1, 9)=WQV(L,1, 9)*EEB
        WQV(L,1,10)=WQV(L,1,10)*EEB
        WQV(L,1,11)=WQV(L,1,11)*EEB
        WQV(L,1,12)=WQV(L,1,12)*EEB
        WQV(L,1,13)=WQV(L,1,13)*EEB
        WQV(L,1,14)=WQV(L,1,14)*EEB
        WQV(L,1,15)=WQV(L,1,15)*EEB
        WQV(L,1,16)=WQV(L,1,16)*EEB
        WQV(L,1,17)=WQV(L,1,17)*EEB
        WQV(L,1,18)=WQV(L,1,18)*EEB
        WQV(L,1,19)=WQV(L,1,19)*EEB
        WQV(L,1,20)=WQV(L,1,20)*EEB
        WQV(L,1,21)=WQV(L,1,21)*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELT*CDZKMK(K)
        RCDZKK=-DELT*CDZKK(K)
        DO L=LF,LL
         CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
         CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K)*AB(L,K)
         CCMBTMP=1.-CCLBTMP-CCUBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         CU1(L,K)=CCUBTMP*EEB
C        CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
         WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
         WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
         WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
         WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
         WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
         WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
         WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
         WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
         WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
         WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
         WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
         WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
         WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
         WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
         WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
         WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
         WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
         WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
         WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
         WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
         WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
        ENDDO
       ENDDO
      ENDDO
C
      K=KC
      RCDZKMK=-DELT*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
        CCMBTMP=1.-CCLBTMP
        EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
C       CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
        WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
        WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
        WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
        WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
        WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
        WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
        WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
        WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
        WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
        WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
        WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
        WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
        WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
        WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
        WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
        WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
        WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
        WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
        WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
        WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
        WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        DO L=LF,LL
C        CWQ(L,K)=CWQ(L,K)-CU1(L,K)*CWQ(L,K+1)
         WQV(L,K, 1)=WQV(L,K, 1)-CU1(L,K)*WQV(L,K+1, 1)
         WQV(L,K, 2)=WQV(L,K, 2)-CU1(L,K)*WQV(L,K+1, 2)
         WQV(L,K, 3)=WQV(L,K, 3)-CU1(L,K)*WQV(L,K+1, 3)
         WQV(L,K, 4)=WQV(L,K, 4)-CU1(L,K)*WQV(L,K+1, 4)
         WQV(L,K, 5)=WQV(L,K, 5)-CU1(L,K)*WQV(L,K+1, 5)
         WQV(L,K, 6)=WQV(L,K, 6)-CU1(L,K)*WQV(L,K+1, 6)
         WQV(L,K, 7)=WQV(L,K, 7)-CU1(L,K)*WQV(L,K+1, 7)
         WQV(L,K, 8)=WQV(L,K, 8)-CU1(L,K)*WQV(L,K+1, 8)
         WQV(L,K, 9)=WQV(L,K, 9)-CU1(L,K)*WQV(L,K+1, 9)
         WQV(L,K,10)=WQV(L,K,10)-CU1(L,K)*WQV(L,K+1,10)
         WQV(L,K,11)=WQV(L,K,11)-CU1(L,K)*WQV(L,K+1,11)
         WQV(L,K,12)=WQV(L,K,12)-CU1(L,K)*WQV(L,K+1,12)
         WQV(L,K,13)=WQV(L,K,13)-CU1(L,K)*WQV(L,K+1,13)
         WQV(L,K,14)=WQV(L,K,14)-CU1(L,K)*WQV(L,K+1,14)
         WQV(L,K,15)=WQV(L,K,15)-CU1(L,K)*WQV(L,K+1,15)
         WQV(L,K,16)=WQV(L,K,16)-CU1(L,K)*WQV(L,K+1,16)
         WQV(L,K,17)=WQV(L,K,17)-CU1(L,K)*WQV(L,K+1,17)
         WQV(L,K,18)=WQV(L,K,18)-CU1(L,K)*WQV(L,K+1,18)
         WQV(L,K,19)=WQV(L,K,19)-CU1(L,K)*WQV(L,K+1,19)
         WQV(L,K,20)=WQV(L,K,20)-CU1(L,K)*WQV(L,K+1,20)
         WQV(L,K,21)=WQV(L,K,21)-CU1(L,K)*WQV(L,K+1,21)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION LEVEL 2
C
C----------------------------------------------------------------------C
C
      IF(ISWQLVL.EQ.2)THEN
C 
      RCDZKK=-DELT*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,1)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB
C       CWQ(L,1)=CWQ(L,1)*EEB
        WQV(L,1, 1)=WQV(L,1, 1)*EEB
        WQV(L,1, 2)=WQV(L,1, 2)*EEB
        WQV(L,1, 3)=WQV(L,1, 3)*EEB
        WQV(L,1, 4)=WQV(L,1, 4)*EEB
        WQV(L,1, 5)=WQV(L,1, 5)*EEB
        WQV(L,1, 6)=WQV(L,1, 6)*EEB
        WQV(L,1, 7)=WQV(L,1, 7)*EEB
        WQV(L,1, 8)=WQV(L,1, 8)*EEB
        WQV(L,1, 9)=WQV(L,1, 9)*EEB
        WQV(L,1,10)=WQV(L,1,10)*EEB
        WQV(L,1,11)=WQV(L,1,11)*EEB
        WQV(L,1,12)=WQV(L,1,12)*EEB
        WQV(L,1,13)=WQV(L,1,13)*EEB
        WQV(L,1,14)=WQV(L,1,14)*EEB
        WQV(L,1,15)=WQV(L,1,15)*EEB
        WQV(L,1,16)=WQV(L,1,16)*EEB
        WQV(L,1,17)=WQV(L,1,17)*EEB
        WQV(L,1,18)=WQV(L,1,18)*EEB
        WQV(L,1,19)=WQV(L,1,19)*EEB
        WQV(L,1,20)=WQV(L,1,20)*EEB
        WQV(L,1,21)=WQV(L,1,21)*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELT*CDZKMK(K)
        RCDZKK=-DELT*CDZKK(K)
        DO L=LF,LL
         CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
         CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K)*AB(L,K)
         CCMBTMP=1.-CCLBTMP-CCUBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         CU1(L,K)=CCUBTMP*EEB
C        CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
         WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
         WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
         WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
         WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
         WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
         WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
         WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
         WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
         WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
         WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
         WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
         WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
         WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
         WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
         WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
         WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
         WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
         WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
         WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
         WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
         WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
        ENDDO
       ENDDO
      ENDDO
C
      K=KC
      RCDZKMK=-DELT*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
        CCMBTMP=1.-CCLBTMP
        EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
C       CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
        WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
        WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
        WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
        WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
        WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
        WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
        WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
        WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
        WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
        WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
        WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
        WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
        WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
        WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
        WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
        WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
        WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
        WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
        WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
        WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
        WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        DO L=LF,LL
C        CWQ(L,K)=CWQ(L,K)-CU1(L,K)*CWQ(L,K+1)
         WQV(L,K, 1)=WQV(L,K, 1)-CU1(L,K)*WQV(L,K+1, 1)
         WQV(L,K, 2)=WQV(L,K, 2)-CU1(L,K)*WQV(L,K+1, 2)
         WQV(L,K, 3)=WQV(L,K, 3)-CU1(L,K)*WQV(L,K+1, 3)
         WQV(L,K, 4)=WQV(L,K, 4)-CU1(L,K)*WQV(L,K+1, 4)
         WQV(L,K, 5)=WQV(L,K, 5)-CU1(L,K)*WQV(L,K+1, 5)
         WQV(L,K, 6)=WQV(L,K, 6)-CU1(L,K)*WQV(L,K+1, 6)
         WQV(L,K, 7)=WQV(L,K, 7)-CU1(L,K)*WQV(L,K+1, 7)
         WQV(L,K, 8)=WQV(L,K, 8)-CU1(L,K)*WQV(L,K+1, 8)
         WQV(L,K, 9)=WQV(L,K, 9)-CU1(L,K)*WQV(L,K+1, 9)
         WQV(L,K,10)=WQV(L,K,10)-CU1(L,K)*WQV(L,K+1,10)
         WQV(L,K,11)=WQV(L,K,11)-CU1(L,K)*WQV(L,K+1,11)
         WQV(L,K,12)=WQV(L,K,12)-CU1(L,K)*WQV(L,K+1,12)
         WQV(L,K,13)=WQV(L,K,13)-CU1(L,K)*WQV(L,K+1,13)
         WQV(L,K,14)=WQV(L,K,14)-CU1(L,K)*WQV(L,K+1,14)
         WQV(L,K,15)=WQV(L,K,15)-CU1(L,K)*WQV(L,K+1,15)
         WQV(L,K,16)=WQV(L,K,16)-CU1(L,K)*WQV(L,K+1,16)
         WQV(L,K,17)=WQV(L,K,17)-CU1(L,K)*WQV(L,K+1,17)
         WQV(L,K,18)=WQV(L,K,18)-CU1(L,K)*WQV(L,K+1,18)
         WQV(L,K,19)=WQV(L,K,19)-CU1(L,K)*WQV(L,K+1,19)
         WQV(L,K,20)=WQV(L,K,20)-CU1(L,K)*WQV(L,K+1,20)
         WQV(L,K,21)=WQV(L,K,21)-CU1(L,K)*WQV(L,K+1,21)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION LEVEL 3
C
C----------------------------------------------------------------------C
C
      IF(ISWQLVL.EQ.3)THEN
C 
      RCDZKK=-DELT*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,1)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB
C       CWQ(L,1)=CWQ(L,1)*EEB
        WQV(L,1, 1)=WQV(L,1, 1)*EEB
        WQV(L,1, 2)=WQV(L,1, 2)*EEB
        WQV(L,1, 3)=WQV(L,1, 3)*EEB
        WQV(L,1, 4)=WQV(L,1, 4)*EEB
        WQV(L,1, 5)=WQV(L,1, 5)*EEB
        WQV(L,1, 6)=WQV(L,1, 6)*EEB
        WQV(L,1, 7)=WQV(L,1, 7)*EEB
        WQV(L,1, 8)=WQV(L,1, 8)*EEB
        WQV(L,1, 9)=WQV(L,1, 9)*EEB
        WQV(L,1,10)=WQV(L,1,10)*EEB
        WQV(L,1,11)=WQV(L,1,11)*EEB
        WQV(L,1,12)=WQV(L,1,12)*EEB
        WQV(L,1,13)=WQV(L,1,13)*EEB
        WQV(L,1,14)=WQV(L,1,14)*EEB
        WQV(L,1,15)=WQV(L,1,15)*EEB
        WQV(L,1,16)=WQV(L,1,16)*EEB
        WQV(L,1,17)=WQV(L,1,17)*EEB
        WQV(L,1,18)=WQV(L,1,18)*EEB
        WQV(L,1,19)=WQV(L,1,19)*EEB
        WQV(L,1,20)=WQV(L,1,20)*EEB
        WQV(L,1,21)=WQV(L,1,21)*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELT*CDZKMK(K)
        RCDZKK=-DELT*CDZKK(K)
        DO L=LF,LL
         CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
         CCUBTMP=RCDZKK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K)*AB(L,K)
         CCMBTMP=1.-CCLBTMP-CCUBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         CU1(L,K)=CCUBTMP*EEB
C        CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
         WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
         WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
         WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
         WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
         WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
         WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
         WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
         WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
         WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
         WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
         WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
         WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
         WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
         WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
         WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
         WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
         WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
         WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
         WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
         WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
         WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
        ENDDO
       ENDDO
      ENDDO
C
      K=KC
      RCDZKMK=-DELT*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCLBTMP=RCDZKMK*GVCSCLPI(L)*HWQI(L)*SWB3D(L,K-1)*AB(L,K-1)
        CCMBTMP=1.-CCLBTMP
        EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
C       CWQ(L,K)=(CWQ(L,K)-CCLBTMP*CWQ(L,K-1))*EEB
        WQV(L,K, 1)=(WQV(L,K, 1)-CCLBTMP*WQV(L,K-1, 1))*EEB
        WQV(L,K, 2)=(WQV(L,K, 2)-CCLBTMP*WQV(L,K-1, 2))*EEB
        WQV(L,K, 3)=(WQV(L,K, 3)-CCLBTMP*WQV(L,K-1, 3))*EEB
        WQV(L,K, 4)=(WQV(L,K, 4)-CCLBTMP*WQV(L,K-1, 4))*EEB
        WQV(L,K, 5)=(WQV(L,K, 5)-CCLBTMP*WQV(L,K-1, 5))*EEB
        WQV(L,K, 6)=(WQV(L,K, 6)-CCLBTMP*WQV(L,K-1, 6))*EEB
        WQV(L,K, 7)=(WQV(L,K, 7)-CCLBTMP*WQV(L,K-1, 7))*EEB
        WQV(L,K, 8)=(WQV(L,K, 8)-CCLBTMP*WQV(L,K-1, 8))*EEB
        WQV(L,K, 9)=(WQV(L,K, 9)-CCLBTMP*WQV(L,K-1, 9))*EEB
        WQV(L,K,10)=(WQV(L,K,10)-CCLBTMP*WQV(L,K-1,10))*EEB
        WQV(L,K,11)=(WQV(L,K,11)-CCLBTMP*WQV(L,K-1,11))*EEB
        WQV(L,K,12)=(WQV(L,K,12)-CCLBTMP*WQV(L,K-1,12))*EEB
        WQV(L,K,13)=(WQV(L,K,13)-CCLBTMP*WQV(L,K-1,13))*EEB
        WQV(L,K,14)=(WQV(L,K,14)-CCLBTMP*WQV(L,K-1,14))*EEB
        WQV(L,K,15)=(WQV(L,K,15)-CCLBTMP*WQV(L,K-1,15))*EEB
        WQV(L,K,16)=(WQV(L,K,16)-CCLBTMP*WQV(L,K-1,16))*EEB
        WQV(L,K,17)=(WQV(L,K,17)-CCLBTMP*WQV(L,K-1,17))*EEB
        WQV(L,K,18)=(WQV(L,K,18)-CCLBTMP*WQV(L,K-1,18))*EEB
        WQV(L,K,19)=(WQV(L,K,19)-CCLBTMP*WQV(L,K-1,19))*EEB
        WQV(L,K,20)=(WQV(L,K,20)-CCLBTMP*WQV(L,K-1,20))*EEB
        WQV(L,K,21)=(WQV(L,K,21)-CCLBTMP*WQV(L,K-1,21))*EEB
       ENDDO
      ENDDO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        DO L=LF,LL
C        CWQ(L,K)=CWQ(L,K)-CU1(L,K)*CWQ(L,K+1)
         WQV(L,K, 1)=WQV(L,K, 1)-CU1(L,K)*WQV(L,K+1, 1)
         WQV(L,K, 2)=WQV(L,K, 2)-CU1(L,K)*WQV(L,K+1, 2)
         WQV(L,K, 3)=WQV(L,K, 3)-CU1(L,K)*WQV(L,K+1, 3)
         WQV(L,K, 4)=WQV(L,K, 4)-CU1(L,K)*WQV(L,K+1, 4)
         WQV(L,K, 5)=WQV(L,K, 5)-CU1(L,K)*WQV(L,K+1, 5)
         WQV(L,K, 6)=WQV(L,K, 6)-CU1(L,K)*WQV(L,K+1, 6)
         WQV(L,K, 7)=WQV(L,K, 7)-CU1(L,K)*WQV(L,K+1, 7)
         WQV(L,K, 8)=WQV(L,K, 8)-CU1(L,K)*WQV(L,K+1, 8)
         WQV(L,K, 9)=WQV(L,K, 9)-CU1(L,K)*WQV(L,K+1, 9)
         WQV(L,K,10)=WQV(L,K,10)-CU1(L,K)*WQV(L,K+1,10)
         WQV(L,K,11)=WQV(L,K,11)-CU1(L,K)*WQV(L,K+1,11)
         WQV(L,K,12)=WQV(L,K,12)-CU1(L,K)*WQV(L,K+1,12)
         WQV(L,K,13)=WQV(L,K,13)-CU1(L,K)*WQV(L,K+1,13)
         WQV(L,K,14)=WQV(L,K,14)-CU1(L,K)*WQV(L,K+1,14)
         WQV(L,K,15)=WQV(L,K,15)-CU1(L,K)*WQV(L,K+1,15)
         WQV(L,K,16)=WQV(L,K,16)-CU1(L,K)*WQV(L,K+1,16)
         WQV(L,K,17)=WQV(L,K,17)-CU1(L,K)*WQV(L,K+1,17)
         WQV(L,K,18)=WQV(L,K,18)-CU1(L,K)*WQV(L,K+1,18)
         WQV(L,K,19)=WQV(L,K,19)-CU1(L,K)*WQV(L,K+1,19)
         WQV(L,K,20)=WQV(L,K,20)-CU1(L,K)*WQV(L,K+1,20)
         WQV(L,K,21)=WQV(L,K,21)-CU1(L,K)*WQV(L,K+1,21)
        ENDDO
       ENDDO
      ENDDO
C
      ENDIF
C
C**********************************************************************C
C
      IF(ISCRAY.EQ.0)THEN
        TWQDIF=TWQDIF+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TWQDIF=TWQDIF+T2TMP-T1TMP
        WTWQDIF=WTWQDIF+(WT2TMP-WT1TMP)*0.001
      ENDIF
C
 2000 CONTINUE
C
C**********************************************************************C
C
      RETURN
      END
