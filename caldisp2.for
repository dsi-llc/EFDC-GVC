C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALDISP2  
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
      DIMENSION UP(KCM),VP(KCM),CDISP(MGM,MGM),CDISPT(KCM,KCM),
     &       CDISPI(KCM,KCM),CCTMP(KCM,KCM),CCUTMP(KCM),CCVTMP(KCM),
     &    CRHS(MGM),CSOL(MGM),WTMP(MGM),VVTMP(MGM,MGM),SVAL(KCM,LCM)
C
C**********************************************************************C
C
      DELT=DT
C
C**********************************************************************C
C
C **  INITIALIZE ON FIRST CALL
C
      IF(N.EQ.NDISP)THEN
C 
      DO L=1,LC
      DO KK=1,KC
      DO K=1,KC
      BDISP(K,KK,L)=0.
      ENDDO
      ENDDO
      ENDDO
C
      DO L=1,LC
      DO K=1,KC
      BDISP(K,K,L)=1.
      FUDISP(K,L)=0.
      FVDISP(K,L)=0.
      ENDDO
      ENDDO
C
      DO L=1,LC
      DXXTCA(L)=0.
      DXYTCA(L)=0.
      DYXTCA(L)=0.
      DYYTCA(L)=0.
      ENDDO
C
      ENDIF
C      
C**********************************************************************C
C
C **  CALCULATE VERTICAL DIFFUSION MATRIX AND INVERSE AND ACCUMULATE
C **  DATA FOR DISPERSION CALCULATION
C
      DO L=2,LA
      IF(LCT(L).EQ.5.AND.SPB(L).NE.0.)THEN
      LN=LNC(L)
C
      DO K=1,KC
      UP(K)=0.5*(U(L,K)+U(L+1,K))
      VP(K)=0.5*(V(L,K)+V(LN,K))
      ENDDO
C      
      UAVG=0.
      VAVG=0.
      DO K=1,KC
      UAVG=UAVG+DZC(K)*UP(K)
      VAVG=VAVG+DZC(K)*VP(K)
      ENDDO
C
      DO K=1,KC
      UP(K)=UP(K)-UAVG
      VP(K)=VP(K)-VAVG
      ENDDO
C      
      DO KK=1,KC
      DO K=1,KC
      CDISP(K,KK)=0.
      ENDDO
      ENDDO
C
      CUTMP=-DELT*CDZKK(1)*AB(L,1)*HPI(L)
      CMTMP=1.-CUTMP
      CDISP(1,1)=CMTMP*DZC(1)
      CDISP(1,2)=CUTMP*DZC(1)
C
      DO K=2,KS
      CLTMP=-DELT*CDZKMK(K)*AB(L,K-1)*HPI(L)
      CUTMP=-DELT*CDZKK(K)*AB(L,K)*HPI(L)
      CMTMP=1.-CLTMP-CUTMP
      CDISP(K,K-1)=CLTMP*DZC(K)
      CDISP(K,K)=CMTMP*DZC(K)
      CDISP(K,K+1)=CUTMP*DZC(K)
      ENDDO
C
      CLTMP=-DELT*CDZKMK(KC)*AB(L,KS)*HPI(L)
      CMTMP=1.-CLTMP
      CDISP(KC,KS)=CLTMP*DZC(KC)
      CDISP(KC,KC)=CMTMP*DZC(KC)
C
      CALL SVDCMP(CDISP,KC,KC,MGM,MGM,WTMP,VVTMP)
C
      DO KK=1,KC
      DO K=1,KC
      CDISPT(K,KK)=CDISP(KK,K)
      ENDDO
      ENDDO
C
      DO KK=1,KC
      DO K=1,KC
      CDISPT(K,KK)=CDISPT(K,KK)/WTMP(K)
      ENDDO
      ENDDO
C
      DO KK=1,KC
      DO K=1,KC
      CTMP=0.
      DO KT=1,KC
      CTMP=CTMP+VVTMP(K,KT)*CDISPT(KT,KK)
      ENDDO
      CDISPI(K,KK)=CTMP
      ENDDO
      ENDDO
C
      DO KK=1,KC
      DO K=1,KC
      CDISPI(K,KK)=CDISPI(K,KK)*DZC(K)
      ENDDO
      ENDDO
C
      DO KK=1,KC
      DO K=1,KC
      CTMP=0.
      DO KT=1,KC
      CTMP=CTMP+CDISPI(K,KT)*BDISP(KT,KK,L)
      ENDDO
      CCTMP(K,KK)=CTMP
      ENDDO
      ENDDO
C
      DO KK=1,KC
      DO K=1,KC
      BDISP(K,KK,L)=CCTMP(K,KK)
      ENDDO
      ENDDO            
C
      DO K=1,KC
      CCUTMP(K)=FUDISP(K,L)-DT*UP(K)/HMIN
      CCVTMP(K)=FVDISP(K,L)-DT*VP(K)/HMIN
      ENDDO
C
      DO K=1,KC
      CCUU=0.
      CCVV=0.
      DO KK=1,KC
      CCUU=CCUU+CDISPI(K,KK)*CCUTMP(KK)
      CCVV=CCVV+CDISPI(K,KK)*CCVTMP(KK)
      ENDDO
      FUDISP(K,L)=CCUU
      FVDISP(K,L)=CCVV
      ENDDO
C
      CCUU=0.
      CCVV=0.
      CCUV=0.
      CCVU=0.
      DO K=1,KC
      CCUU=CCUU+DZC(K)*UP(K)*FUDISP(K,L)
      CCUV=CCUV+DZC(K)*UP(K)*FVDISP(K,L)
      CCVU=CCVU+DZC(K)*VP(K)*FUDISP(K,L)
      CCVV=CCVV+DZC(K)*VP(K)*FVDISP(K,L)
      ENDDO
      DXXTCA(L)=DXXTCA(L)+CCUU*HP(L)
      DXYTCA(L)=DXYTCA(L)+CCUV*HP(L)
      DYXTCA(L)=DYXTCA(L)+CCVU*HP(L)
      DYYTCA(L)=DYYTCA(L)+CCVV*HP(L)
C
      DO K=1,KC
      CCUU=0.
      CCVV=0.
      DO KK=1,KC
      CCUU=CCUU+DZC(KK)*UP(KK)*BDISP(KK,K,L)
      CCVV=CCVV+DZC(KK)*VP(KK)*BDISP(KK,K,L)
      ENDDO
      CUDISPT(K,L)=CCUU*HP(L)
      CVDISPT(K,L)=CCVV*HP(L)
      ENDDO
C
      ENDIF
      ENDDO
C
      IF(N.LT.NTS) RETURN
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF DISPERSION COEFFICIENTS
C
      DO L=2,LA
      IF(LCT(L).EQ.5.AND.SPB(L).NE.0.)THEN
C
      DO KK=1,KC
      DO K=1,KC
      CDISP(K,KK)=-BDISP(K,KK,L)
      ENDDO
      ENDDO
C
      DO K=1,KC
      CDISP(K,K)=1.+CDISP(K,K)
      ENDDO
C
      CALL SVDCMP(CDISP,KC,KC,MGM,MGM,WTMP,VVTMP)
C
      DO K=1,KC
      SVAL(K,L)=WTMP(K)
      ENDDO
C
C     DO K=1,KC
C     AWTMP=ABS(WTMP(K))
C     IF(AWTMP.LT.0.001) WTMP(K)=0.
C     ENDDO
C
      DO K=1,KC
      CRHS(K)=FUDISP(K,L)
      ENDDO
      CALL SVBKSB(CDISP,WTMP,VVTMP,KC,KC,MGM,MGM,CRHS,CSOL)
      CCUU=0.
      CCVU=0.
      DO K=1,KC
      CCUU=CCUU+CUDISPT(K,L)*CSOL(K)
      CCVU=CCVU+CVDISPT(K,L)*CSOL(K)
      ENDDO
      DXXTCA(L)=-(DXXTCA(L)+CCUU)*HMIN/TPN
      DYXTCA(L)=-(DYXTCA(L)+CCVU)*HMIN/TPN
C
      DO K=1,KC
      CRHS(K)=FVDISP(K,L)
      ENDDO
      CALL SVBKSB(CDISP,WTMP,VVTMP,KC,KC,MGM,MGM,CRHS,CSOL)
      CCVV=0.
      CCUV=0.
      DO K=1,KC
      CCVV=CCVV+CVDISPT(K,L)*CSOL(K)
      CCUV=CCUV+CUDISPT(K,L)*CSOL(K)
      ENDDO
      DYYTCA(L)=-(DYYTCA(L)+CCVV)*HMIN/TPN
      DXYTCA(L)=-(DXYTCA(L)+CCUV)*HMIN/TPN
C
      ENDIF
      ENDDO
C
      DO L=2,LA
      DXXTCA(L)=DXXTCA(L)/HLPF(L)
      DXYTCA(L)=DXYTCA(L)/HLPF(L)
      DYXTCA(L)=DYXTCA(L)/HLPF(L)
      DYYTCA(L)=DYYTCA(L)/HLPF(L)
      ENDDO
C
C**********************************************************************C
C
C **  WRITE OUTPUT FILES
C
      OPEN(88,FILE='DISTEN.OUT',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='DISTEN.OUT',STATUS='UNKNOWN')
      WRITE(88,881)
C
      DO L=2,LA
      WRITE(88,2011)IL(L),JL(L),DLON(L),DLAT(L),DXXTCA(L),
     &                    DXYTCA(L),DYXTCA(L),DYYTCA(L)
      ENDDO
C
      CLOSE(88)
C
      OPEN(88,FILE='UVTSC.OUT',STATUS='UNKNOWN')
      WRITE(88,882)
C
      DO L=2,LA
      AMCPT=AMCP(L)*GI
      AMSPT=AMSP(L)*GI
      WRITE(88,2012)IL(L),JL(L),DLON(L),DLAT(L),AMCPT,AMSPT,
     &                    AMCUE(L),AMSUE(L),AMCVE(L),AMSVE(L)
      ENDDO
C
      CLOSE(88)
C
      OPEN(88,FILE='UVERV.OUT',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='UVERV.OUT',STATUS='UNKNOWN')
      WRITE(88,883)
C
      DO L=2,LA
      WRITE(88,2012)IL(L),JL(L),DLON(L),DLAT(L),HLPF(L),UELPF(L),
     &                    VELPF(L),SALLPF(L,1),SALLPF(L,KC)
      ENDDO
C
      CLOSE(88)
C
      OPEN(88,FILE='SINVAL.OUT',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='SINVAL.OUT',STATUS='UNKNOWN')
C
      DO L=2,LA
      WRITE(88,2013)IL(L),JL(L),(SVAL(K,L),K=1,KC)
      ENDDO
C
      CLOSE(88)
C
  881 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'DXX',10X,'DXY',10X,
     &          'DYX',10X,'DYY')
  882 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'AMCPT',8X,'AMSPT',8X,
     &          'AMCUE',8X,'AMSUE',8X,'AMCVE',8X,'AMSVE')
  883 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'HLPF',9X,'UELPF',8X,
     &          'VELPF',8X,'SALLPFBOT',4X,'SALLPFSURF')
 2011 FORMAT(2I4,2X,F10.6,2X,F10.6,4(2X,E12.4))
 2012 FORMAT(2I4,2X,F10.6,2X,F10.6,6(2X,E12.4))
 2013 FORMAT(2I4,8(2X,E12.4))
C
C**********************************************************************C
C
      RETURN
      END
