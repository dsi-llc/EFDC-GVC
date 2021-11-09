C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LAGRES
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
C **  SUBROUTINE LAGRES CALCULATES THE MEAN OR RESIDUAL LAGRANGIAN
C **  VELOCITY FOR VARIOUS PARTICLE RELEASE TIMES OVER THE REFERENCE
C **  TIME PERIOD TIDALP BY INTEGRATING THE THREE DIMENSIONAL
C **  TRAJECTORIES OF NEUTRALLY BUOYANT PARTICLES.  A TRILINEAR
C **  INTERPOLATION WITH AN EXPLICIT FORWARD OR IMPLICIT CENTERED 
C **  INTERGRATION SCHEME IS USED. 
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      TITLE='LAGRANGIAN RESIDUAL DRIFTER DATA'
      IF(JSLRPD.NE.1) GOTO 5
C
      OPEN(80,FILE='LRPD.LOG',STATUS='UNKNOWN')
      CLOSE(80,STATUS='DELETE')
      OPEN(80,FILE='LRPD.LOG',STATUS='UNKNOWN')
      CLOSE(80)
C
      DO M=1,MLRPDRT
      MCNTLR(M)=0
      IACTLR(M)=0
      ENDDO
C
      DO L=1,LC
      NLRPDL(L)=LC
      ILRPD(L)=IC
      JLRPD(L)=JC
      ENDDO
C
      NLRPD=0
C
      DO J=JLRPD1,JLRPD2
      DO I=ILRPD1,ILRPD2
      IF(IJCT(I,J).EQ.5)THEN
        NLRPD=NLRPD+1
        ILRPD(NLRPD)=I
        JLRPD(NLRPD)=J
        L=LIJ(I,J)
        NLRPDL(L)=NLRPD
        DO M=1,MLRPDRT
        DO K=1,KC
        XLRPD(NLRPD,K,M)=FLOAT(I)
        YLRPD(NLRPD,K,M)=FLOAT(J)
        ZLRPD(NLRPD,K,M)=FLOAT(K)
        ENDDO
        ENDDO
      ENDIF 
      ENDDO
      ENDDO
C
      JSLRPD=0
C
C----------------------------------------------------------------------C
C
    5 CONTINUE
      NTSPTCP=NTSPTC+1
      OPEN(80,FILE='LRPD.LOG',POSITION='APPEND',STATUS='UNKNOWN')
C
      DO M=1,MLRPDRT
      IF(N.EQ.NLRPDRT(M))THEN
        IACTLR(M)=1
        WRITE(80,801)M,N
      ENDIF
      IF(N.GT.NLRPDRT(M))THEN
        IF(MCNTLR(M).EQ.NTSPTCP)THEN
          MCNTLR(M)=0
          IACTLR(M)=0
          WRITE(80,802)M,N-1
          DO NLR=1,NLRPD
          L=LIJ(ILRPD(NLR),JLRPD(NLR))
          HLRPD(NLR,M)=H1P(L)
          ENDDO
        ENDIF
      ENDIF
      ENDDO
C      
      DO M=1,MLRPDRT
      IF(IACTLR(M).GE.1)THEN
        MCNTLR(M)=MCNTLR(M)+1
C
        IF(IACTLR(M).EQ.1)THEN
          DO KLR=1,KC
          DO NLR=1,NLRPD
          I=ILRPD(NLR)
          J=JLRPD(NLR)
          L=LIJ(I,J)
          ZILRPD(NLR,KLR,M)=ZZ(KLR)*HP(L)+BELV(L)
          ENDDO 
          ENDDO
          IACTLR(M)=2
        ENDIF 
C
        DO KLR=1,KC
        DO NLR=1,NLRPD
        XTMP=XLRPD(NLR,KLR,M)
        YTMP=YLRPD(NLR,KLR,M)
        ZTMP=ZLRPD(NLR,KLR,M)
        I=NINT(XTMP)
        J=NINT(YTMP)
        K=NINT(ZTMP)
        L=LIJ(I,J)
        LP=L+1
        LM=L-1
        LN=LNC(L)
        LS=LSC(L)
        LNP=LNC(LP)
        LSP=LSC(LP)
        LNM=LNC(LM)
        LSM=LSC(LM)
        APLUS=XTMP-FLOAT(I)+0.5
        AMINU=1.-APLUS
        BPLUS=YTMP-FLOAT(J)+0.5
        BMINU=1.-BPLUS
        CPLUS=ZTMP-FLOAT(K)+0.5
        CMINU=1.-CPLUS
        ABCPPP=APLUS*BPLUS*CPLUS
        ABCPPM=APLUS*BPLUS*CMINU
        ABCPMP=APLUS*BMINU*CPLUS
        ABCPMM=APLUS*BMINU*CMINU
        ABCMPP=AMINU*BPLUS*CPLUS
        ABCMPM=AMINU*BPLUS*CMINU
        ABCMMP=AMINU*BMINU*CPLUS
        ABCMMM=AMINU*BMINU*CMINU
        IF(ISLRPD.EQ.2)THEN
          BCPPP=BPLUS*CPLUS
          BCPPM=BPLUS*CMINU
          BCPMP=BMINU*CPLUS
          BCPMM=BMINU*CMINU
          BCMPP=-BPLUS*CPLUS
          BCMPM=-BPLUS*CMINU
          BCMMP=-BMINU*CPLUS
          BCMMM=-BMINU*CMINU
          ACPPP=APLUS*CPLUS
          ACPPM=APLUS*CMINU
          ACPMP=-APLUS*CPLUS
          ACPMM=-APLUS*CMINU
          ACMPP=AMINU*CPLUS
          ACMPM=AMINU*CMINU
          ACMMP=-AMINU*CPLUS
          ACMMM=-AMINU*CMINU
          ABPPP=APLUS*BPLUS
          ABPPM=-APLUS*BPLUS
          ABPMP=APLUS*BMINU
          ABPMM=-APLUS*BMINU
          ABMPP=AMINU*BPLUS
          ABMPM=-AMINU*BPLUS
          ABMMP=AMINU*BMINU
          ABMMM=-AMINU*BMINU
        ENDIF 
C
        IF(KC.EQ.1)THEN
          ULK   =U(L  ,K )*DXIU(L  )
          ULKP  =ULK
          ULKM  =ULK
          ULPK  =U(LP ,K )*DXIU(LP )
          ULPKP =ULPK
          ULPKM =ULPK
          ULNK  =U(LN ,K )*DXIU(LN )
          ULNKP =ULNK
          ULNKM =ULNK
          ULNPK =U(LNP,K )*DXIU(LNP)
          ULNPKP=ULNPK
          ULNPKM=ULNPK
          ULSK  =U(LS ,K )*DXIU(LS )
          ULSKP =ULSK
          ULSKM =ULSK
          ULSPK =U(LSP,K )*DXIU(LSP)
          ULSPKP=ULSPK
          ULSPKM=ULSPK
          VLK   =V(L  ,K )*DYIV(L  )
          VLKP  =VLK
          VLKM  =VLK
          VLPK  =V(LP ,K )*DYIV(LP )
          VLPKP =VLPK
          VLPKM =VLPK
          VLMK  =V(LM ,K )*DYIV(LM )
          VLMKP =VLMK
          VLMKM =VLMK
          VLNK  =V(LN ,K )*DYIV(LN )
          VLNKP =VLNK
          VLNKM =VLNK
          VLNPK =V(LNP,K )*DYIV(LNP)
          VLNPKP=VLNPK
          VLNPKM=VLNPK
          VLNMK =V(LNM,K )*DYIV(LNM)
          VLNMKP=VLNMK
          VLNMKM=VLNMK
          UPPP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(LN)*(ULNPK+ULNPKP))
     &                                  /(1.+SVB(LN))                     
          UPPM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(LN)*(ULNPK+ULNPKM))
     &                                  /(1.+SVB(LN))
          UPMP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(L )*(ULSPK+ULSPKP))
     &                                  /(1.+SVB(L ))
          UPMM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(L )*(ULSPK+ULSPKM))
     &                                  /(1.+SVB(L ))
          UMPP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(LN)*(ULNK +ULNKP ))
     &                                  /(1.+SVB(LN))
          UMPM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(LN)*(ULNK +ULNKM ))
     &                                  /(1.+SVB(LN))
          UMMP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(L )*(ULSK +ULSKP ))
     &                                  /(1.+SVB(L ))
          UMMM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(L )*(ULSK +ULSKM ))
     &                                  /(1.+SVB(L ))
          VPPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(LP)*(VLNPK+VLNPKP))
     &                                  /(1.+SUB(LP))
          VPPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(LP)*(VLNPK+VLNPKM))
     &                                  /(1.+SUB(LP))
          VPMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(LP)*(VLPK +VLPKP ))
     &                                  /(1.+SUB(LP))
          VPMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(LP)*(VLPK +VLPKM ))
     &                                  /(1.+SUB(LP))
          VMPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(L )*(VLNMK+VLNMKP))
     &                                  /(1.+SUB(L ))
          VMPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(L )*(VLNMK+VLNMKM))
     &                                  /(1.+SUB(L ))
          VMMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(L )*(VLMK +VLMKP ))
     &                                  /(1.+SUB(L ))
          VMMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(L )*(VLMK +VLMKM ))
     &                                  /(1.+SUB(L ))
          WPPP=0.
          WPPM=0.
          WPMP=0.
          WPMM=0.
          WMPP=0.
          WMPM=0.
          WMMP=0.
          WMMM=0.
          GOTO 999
        ENDIF
        IF(K.EQ.1)THEN
          KP=K+1
          ULK   =U(L  ,K )*DXIU(L  )
          ULKP  =U(L  ,KP)*DXIU(L  )
          ULKM  =0.
          ULPK  =U(LP ,K )*DXIU(LP )
          ULPKP =U(LP ,KP)*DXIU(LP )
          ULPKM =0.
          ULNK  =U(LN ,K )*DXIU(LN )
          ULNKP =U(LN ,KP)*DXIU(LN )
          ULNKM =0.
          ULNPK =U(LNP,K )*DXIU(LNP)
          ULNPKP=U(LNP,KP)*DXIU(LNP)
          ULNPKM=0.
          ULSK  =U(LS ,K )*DXIU(LS )
          ULSKP =U(LS ,KP)*DXIU(LS )
          ULSKM =0.
          ULSPK =U(LSP,K )*DXIU(LSP)
          ULSPKP=U(LSP,KP)*DXIU(LSP)
          ULSPKM=0.
          VLK   =V(L  ,K )*DYIV(L  )
          VLKP  =V(L  ,KP)*DYIV(L  )
          VLKM  =0.
          VLPK  =V(LP ,K )*DYIV(LP )
          VLPKP =V(LP ,KP)*DYIV(LP )
          VLPKM =0.
          VLMK  =V(LM ,K )*DYIV(LM )
          VLMKP =V(LM ,KP)*DYIV(LM )
          VLMKM =0.
          VLNK  =V(LN ,K )*DYIV(LN )
          VLNKP =V(LN ,KP)*DYIV(LN )
          VLNKM =0.
          VLNPK =V(LNP,K )*DYIV(LNP)
          VLNPKP=V(LNP,KP)*DYIV(LNP)
          VLNPKM=0.
          VLNMK =V(LNM,K )*DYIV(LNM)
          VLNMKP=V(LNM,KP)*DYIV(LNM)
          VLNMKM=0.
          WLK   =W(L  ,K )*HPI(L  )*DZIG(K )
          WLPK  =W(LP ,K )*HPI(LP )*DZIG(K )
          WLMK  =W(LM ,K )*HPI(LM )*DZIG(K )
          WLNK  =W(LN ,K )*HPI(LN )*DZIG(K )
          WLNPK =W(LNP,K )*HPI(LNP)*DZIG(K )
          WLNMK =W(LNM,K )*HPI(LNM)*DZIG(K )
          WLSK  =W(LS ,K )*HPI(LS )*DZIG(K ) 
          WLSPK =W(LSP,K )*HPI(LSP)*DZIG(K )
          WLSMK =W(LSM,K )*HPI(LSM)*DZIG(K )
          UPPP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(LN)*(ULNPK+ULNPKP))
     &                                  /(1.+SVB(LN))                     
          UPPM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(LN)*(ULNPK+ULNPKM))
     &                                  /(1.+SVB(LN))
          UPMP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(L )*(ULSPK+ULSPKP))
     &                                  /(1.+SVB(L ))
          UPMM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(L )*(ULSPK+ULSPKM))
     &                                  /(1.+SVB(L ))
          UMPP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(LN)*(ULNK +ULNKP ))
     &                                  /(1.+SVB(LN))
          UMPM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(LN)*(ULNK +ULNKM ))
     &                                  /(1.+SVB(LN))
          UMMP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(L )*(ULSK +ULSKP ))
     &                                  /(1.+SVB(L ))
          UMMM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(L )*(ULSK +ULSKM ))
     &                                  /(1.+SVB(L ))
          VPPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(LP)*(VLNPK+VLNPKP))
     &                                  /(1.+SUB(LP))
          VPPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(LP)*(VLNPK+VLNPKM))
     &                                  /(1.+SUB(LP))
          VPMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(LP)*(VLPK +VLPKP ))
     &                                  /(1.+SUB(LP))
          VPMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(LP)*(VLPK +VLPKM ))
     &                                  /(1.+SUB(LP))
          VMPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(L )*(VLNMK+VLNMKP))
     &                                  /(1.+SUB(L ))
          VMPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(L )*(VLNMK+VLNMKM))
     &                                  /(1.+SUB(L ))
          VMMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(L )*(VLMK +VLMKP ))
     &                                  /(1.+SUB(L ))
          VMMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(L )*(VLMK +VLMKM ))
     &                                  /(1.+SUB(L ))
          WPPP=(WLK +SUB(LP)*WLPK +SVB(LN)*WLNK +SUB(LP)*SVB(LN)*WLNPK )
     &        /(1.  +SUB(LP)      +SVB(LN)      +SUB(LP)*SVB(LN)       )
          WPPM=0.
          WPMP=(WLK +SUB(LP)*WLPK +SVB(L )*WLSK +SUB(LP)*SVB(L )*WLSPK )
     &        /(1.  +SUB(LP)      +SVB(L )      +SUB(LP)*SVB(L )       )
          WPMM=0.
          WMPP=(WLK +SUB(L )*WLMK +SVB(LN)*WLNK +SUB(L )*SVB(LN)*WLNMK )
     &        /(1.  +SUB(L )      +SVB(LN)      +SUB(L )*SVB(LN)       )
          WMPM=0.
          WMMP=(WLK +SUB(L )*WLMK +SVB(L )*WLSK +SUB(L )*SVB(L )*WLSMK )
     &        /(1.  +SUB(L )      +SVB(L )      +SUB(L )*SVB(L )*WLSMK )
          WMMM=0.
        ENDIF
        IF(K.GT.1.AND.K.LT.KC)THEN
          KP=K+1
          KM=K-1
          ULK   =U(L  ,K )*DXIU(L  )
          ULKP  =U(L  ,KP)*DXIU(L  )
          ULKM  =U(L  ,KM)*DXIU(L  )
          ULPK  =U(LP ,K )*DXIU(LP )
          ULPKP =U(LP ,KP)*DXIU(LP )
          ULPKM =U(LP ,KM)*DXIU(LP )
          ULNK  =U(LN ,K )*DXIU(LN )
          ULNKP =U(LN ,KP)*DXIU(LN )
          ULNKM =U(LN ,KM)*DXIU(LN )
          ULNPK =U(LNP,K )*DXIU(LNP)
          ULNPKP=U(LNP,KP)*DXIU(LNP)
          ULNPKM=U(LNP,KM)*DXIU(LNP)
          ULSK  =U(LS ,K )*DXIU(LS )
          ULSKP =U(LS ,KP)*DXIU(LS )
          ULSKM =U(LS ,KM)*DXIU(LS )
          ULSPK =U(LSP,K )*DXIU(LSP)
          ULSPKP=U(LSP,KP)*DXIU(LSP)
          ULSPKM=U(LSP,KM)*DXIU(LSP)
          VLK   =V(L  ,K )*DYIV(L  )
          VLKP  =V(L  ,KP)*DYIV(L  )
          VLKM  =V(L  ,KM)*DYIV(L  )
          VLPK  =V(LP ,K )*DYIV(LP )
          VLPKP =V(LP ,KP)*DYIV(LP )
          VLPKM =V(LP ,KM)*DYIV(LP )
          VLMK  =V(LM ,K )*DYIV(LM )
          VLMKP =V(LM ,KP)*DYIV(LM )
          VLMKM =V(LM ,KM)*DYIV(LM )
          VLNK  =V(LN ,K )*DYIV(LN )
          VLNKP =V(LN ,KP)*DYIV(LN )
          VLNKM =V(LN ,KM)*DYIV(LN )
          VLNPK =V(LNP,K )*DYIV(LNP)
          VLNPKP=V(LNP,KP)*DYIV(LNP)
          VLNPKM=V(LNP,KM)*DYIV(LNP)
          VLNMK =V(LNM,K )*DYIV(LNM)
          VLNMKP=V(LNM,KP)*DYIV(LNM)
          VLNMKM=V(LNM,KM)*DYIV(LNM)
          WLK   =W(L  ,K )*HPI(L  )*DZIG(K )
          WLKM  =W(L  ,KM)*HPI(L  )*DZIG(KM)
          WLPK  =W(LP ,K )*HPI(LP )*DZIG(K )
          WLPKM =W(LP ,KM)*HPI(LP )*DZIG(KM)
          WLMK  =W(LM ,K )*HPI(LM )*DZIG(K )
          WLMKM =W(LM ,KM)*HPI(LM )*DZIG(KM)
          WLNK  =W(LN ,K )*HPI(LN )*DZIG(K )
          WLNKM =W(LN ,KM)*HPI(LN )*DZIG(KM)
          WLNPK =W(LNP,K )*HPI(LNP)*DZIG(K )
          WLNPKM=W(LNP,KM)*HPI(LNP)*DZIG(KM)
          WLNMK =W(LNM,K )*HPI(LNM)*DZIG(K )
          WLNMKM=W(LNM,KM)*HPI(LNM)*DZIG(KM)
          WLSK  =W(LS ,K )*HPI(LS )*DZIG(K ) 
          WLSKM =W(LS ,KM)*HPI(LS )*DZIG(KM)
          WLSPK =W(LSP,K )*HPI(LSP)*DZIG(K )
          WLSPKM=W(LSP,KM)*HPI(LSP)*DZIG(KM)
          WLSMK =W(LSM,K )*HPI(LSM)*DZIG(K )
          WLSMKM=W(LSM,KM)*HPI(LSM)*DZIG(KM)
          UPPP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(LN)*(ULNPK+ULNPKP))
     &                                  /(1.+SVB(LN))                     
          UPPM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(LN)*(ULNPK+ULNPKM))
     &                                  /(1.+SVB(LN))
          UPMP=SUB(LP)*(0.5*(ULPK+ULPKP)+0.5*SVB(L )*(ULSPK+ULSPKP))
     &                                  /(1.+SVB(L ))
          UPMM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(L )*(ULSPK+ULSPKM))
     &                                  /(1.+SVB(L ))
          UMPP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(LN)*(ULNK +ULNKP ))
     &                                  /(1.+SVB(LN))
          UMPM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(LN)*(ULNK +ULNKM ))
     &                                  /(1.+SVB(LN))
          UMMP=SUB(L )*(0.5*(ULK +ULKP )+0.5*SVB(L )*(ULSK +ULSKP ))
     &                                  /(1.+SVB(L ))
          UMMM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(L )*(ULSK +ULSKM ))
     &                                  /(1.+SVB(L ))
          VPPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(LP)*(VLNPK+VLNPKP))
     &                                  /(1.+SUB(LP))
          VPPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(LP)*(VLNPK+VLNPKM))
     &                                  /(1.+SUB(LP))
          VPMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(LP)*(VLPK +VLPKP ))
     &                                  /(1.+SUB(LP))
          VPMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(LP)*(VLPK +VLPKM ))
     &                                  /(1.+SUB(LP))
          VMPP=SVB(LN)*(0.5*(VLNK+VLNKP)+0.5*SUB(L )*(VLNMK+VLNMKP))
     &                                  /(1.+SUB(L ))
          VMPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(L )*(VLNMK+VLNMKM))
     &                                  /(1.+SUB(L ))
          VMMP=SVB(L )*(0.5*(VLK +VLKP )+0.5*SUB(L )*(VLMK +VLMKP ))
     &                                  /(1.+SUB(L ))
          VMMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(L )*(VLMK +VLMKM ))
     &                                  /(1.+SUB(L ))
          WPPP=(WLK +SUB(LP)*WLPK +SVB(LN)*WLNK +SUB(LP)*SVB(LN)*WLNPK )
     &        /(1.  +SUB(LP)      +SVB(LN)      +SUB(LP)*SVB(LN)       )
          WPPM=(WLKM+SUB(LP)*WLPKM+SVB(LN)*WLNKM+SUB(LP)*SVB(LN)*WLNPKM)
     &        /(1.  +SUB(LP)      +SVB(LN)      +SUB(LP)*SVB(LN)       )
          WPMP=(WLK +SUB(LP)*WLPK +SVB(L )*WLSK +SUB(LP)*SVB(L )*WLSPK )
     &        /(1.  +SUB(LP)      +SVB(L )      +SUB(LP)*SVB(L )       )
          WPMM=(WLKM+SUB(LP)*WLPKM+SVB(L )*WLSKM+SUB(LP)*SVB(L )*WLSPKM)
     &        /(1.  +SUB(LP)      +SVB(L )      +SUB(LP)*SVB(L )       )
          WMPP=(WLK +SUB(L )*WLMK +SVB(LN)*WLNK +SUB(L )*SVB(LN)*WLNMK )
     &        /(1.  +SUB(L )      +SVB(LN)      +SUB(L )*SVB(LN)       )
          WMPM=(WLKM+SUB(L )*WLMKM+SVB(LN)*WLNKM+SUB(L )*SVB(LN)*WLNMKM)
     &        /(1.  +SUB(L )      +SVB(LN)      +SUB(L )*SVB(LN)       )
          WMMP=(WLK +SUB(L )*WLMK +SVB(L )*WLSK +SUB(L )*SVB(L )*WLSMK )
     &        /(1.  +SUB(L )      +SVB(L )      +SUB(L )*SVB(L )*WLSMK )
          WMMM=(WLKM+SUB(L )*WLMKM+SVB(L )*WLSKM+SUB(L )*SVB(L )*WLSMKM)
     &        /(1.  +SUB(L )      +SVB(L )      +SUB(L )*SVB(L )       )
        ENDIF
        IF(K.EQ.KC)THEN
          KM=K-1
          ULK   =U(L  ,K )*DXIU(L  )
          ULKM  =U(L  ,KM)*DXIU(L  )
          ULPK  =U(LP ,K )*DXIU(LP )
          ULPKM =U(LP ,KM)*DXIU(LP )
          ULNK  =U(LN ,K )*DXIU(LN )
          ULNKM =U(LN ,KM)*DXIU(LN )
          ULNPK =U(LNP,K )*DXIU(LNP)
          ULNPKM=U(LNP,KM)*DXIU(LNP)
          ULSK  =U(LS ,K )*DXIU(LS )
          ULSKM =U(LS ,KM)*DXIU(LS )
          ULSPK =U(LSP,K )*DXIU(LSP)
          ULSPKM=U(LSP,KM)*DXIU(LSP)
          VLK   =V(L  ,K )*DYIV(L  )
          VLKM  =V(L  ,KM)*DYIV(L  )
          VLPK  =V(LP ,K )*DYIV(LP )
          VLPKM =V(LP ,KM)*DYIV(LP )
          VLMK  =V(LM ,K )*DYIV(LM )
          VLMKM =V(LM ,KM)*DYIV(LM )
          VLNK  =V(LN ,K )*DYIV(LN )
          VLNKM =V(LN ,KM)*DYIV(LN )
          VLNPK =V(LNP,K )*DYIV(LNP)
          VLNPKM=V(LNP,KM)*DYIV(LNP)
          VLNMK =V(LNM,K )*DYIV(LNM)
          VLNMKM=V(LNM,KM)*DYIV(LNM)
          WLKM  =W(L  ,KM)*HPI(L  )*DZIG(KM)
          WLPKM =W(LP ,KM)*HPI(LP )*DZIG(KM)
          WLMKM =W(LM ,KM)*HPI(LM )*DZIG(KM)
          WLNKM =W(LN ,KM)*HPI(LN )*DZIG(KM)
          WLNPKM=W(LNP,KM)*HPI(LNP)*DZIG(KM)
          WLNMKM=W(LNM,KM)*HPI(LNM)*DZIG(KM)
          WLSKM =W(LS ,KM)*HPI(LS )*DZIG(KM)
          WLSPKM=W(LSP,KM)*HPI(LSP)*DZIG(KM)
          WLSMKM=W(LSM,KM)*HPI(LSM)*DZIG(KM)
          UPPP=SUB(LP)*( 0.5*(3.*ULPK -ULPKM )
     &          +0.5*SVB(LN)*(3.*ULNPK-ULNPKM) )
     &          /(1.+SVB(LN))
          UPMP=SUB(LP)*( 0.5*(3.*ULPK -ULPKM )
     &          +0.5*SVB(L )*(3.*ULSPK-ULSPKM) )
     &          /(1.+SVB(L ))
          UMPP=SUB(L )*( 0.5*(3.*ULK  -ULKM  )
     &          +0.5*SVB(LN)*(3.*ULNK -ULNKM ) )
     &          /(1.+SVB(LN))
          UMMP=SUB(L )*( 0.5*(3.*ULK  -ULKM  )
     &          +0.5*SVB(L )*(3.*ULSK -ULSKM ) )
     &          /(1.+SVB(L ))
          VPPP=SVB(LN)*( 0.5*(3.*VLNK -VLNKM )
     &          +0.5*SUB(LP)*(3.*VLNPK-VLNPKM) )
     &          /(1.+SUB(LP))
          VPMP=SVB(L )*( 0.5*(3.*VLK  -VLKM  )
     &          +0.5*SUB(LP)*(3.*VLPK -VLPKM ) )
     &          /(1.+SUB(LP))
          VMPP=SVB(LN)*( 0.5*(3.*VLNK -VLNKM )
     &          +0.5*SUB(L )*(3.*VLNMK-VLNMKM) )
     &          /(1.+SUB(L ))
          VMMP=SVB(L )*( 0.5*(3.*VLK  -VLKM  )
     &          +0.5*SUB(L )*(3.*VLMK -VLMKM ) )
     &          /(1.+SUB(L ))
          UPPM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(LN)*(ULNPK+ULNPKM))
     &                                  /(1.+SVB(LN))
          UPMM=SUB(LP)*(0.5*(ULPK+ULPKM)+0.5*SVB(L )*(ULSPK+ULSPKM))
     &                                  /(1.+SVB(L ))
          UMPM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(LN)*(ULNK +ULNKM ))
     &                                  /(1.+SVB(LN))
          UMMM=SUB(L )*(0.5*(ULK +ULKM )+0.5*SVB(L )*(ULSK +ULSKM ))
     &                                  /(1.+SVB(L ))
          VPPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(LP)*(VLNPK+VLNPKM))
     &                                  /(1.+SUB(LP))
          VPMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(LP)*(VLPK +VLPKM ))
     &                                  /(1.+SUB(LP))
          VMPM=SVB(LN)*(0.5*(VLNK+VLNKM)+0.5*SUB(L )*(VLNMK+VLNMKM))
     &                                  /(1.+SUB(L ))
          VMMM=SVB(L )*(0.5*(VLK +VLKM )+0.5*SUB(L )*(VLMK +VLMKM ))
     &                                  /(1.+SUB(L ))
          WPPP=0.
          WPPM=(WLKM+SUB(LP)*WLPKM+SVB(LN)*WLNKM+SUB(LP)*SVB(LN)*WLNPKM)
     &        /(1.  +SUB(LP)      +SVB(LN)      +SUB(LP)*SVB(LN)       )
          WPMP=0.
          WPMM=(WLKM+SUB(LP)*WLPKM+SVB(L )*WLSKM+SUB(LP)*SVB(L )*WLSPKM)
     &        /(1.  +SUB(LP)      +SVB(L )      +SUB(LP)*SVB(L )       )
          WMPP=0.
          WMPM=(WLKM+SUB(L )*WLMKM+SVB(LN)*WLNKM+SUB(L )*SVB(LN)*WLNMKM)
     &        /(1.  +SUB(L )      +SVB(LN)      +SUB(L )*SVB(LN)       )
          WMMP=0.
          WMMM=(WLKM+SUB(L )*WLMKM+SVB(L )*WLSKM+SUB(L )*SVB(L )*WLSMKM)
     &        /(1.  +SUB(L )      +SVB(L )      +SUB(L )*SVB(L )       )
        ENDIF
C
  999 CONTINUE
C
        UXYZ=ABCPPP*UPPP+ABCPPM*UPPM+ABCPMP*UPMP+ABCPMM*UPMM
     &      +ABCMPP*UMPP+ABCMPM*UMPM+ABCMMP*UMMP+ABCMMM*UMMM
        VXYZ=ABCPPP*VPPP+ABCPPM*VPPM+ABCPMP*VPMP+ABCPMM*VPMM
     &      +ABCMPP*VMPP+ABCMPM*VMPM+ABCMMP*VMMP+ABCMMM*VMMM
        WXYZ=ABCPPP*WPPP+ABCPPM*WPPM+ABCPMP*WPMP+ABCPMM*WPMM
     &      +ABCMPP*WMPP+ABCMPM*WMPM+ABCMMP*WMMP+ABCMMM*WMMM
        IF(ISLRPD.EQ.2)THEN
          DXUXYZ=BCPPP*UPPP+BCPPM*UPPM+BCPMP*UPMP+BCPMM*UPMM
     &          +BCMPP*UMPP+BCMPM*UMPM+BCMMP*UMMP+BCMMM*UMMM
          DYUXYZ=ACPPP*UPPP+ACPPM*UPPM+ACPMP*UPMP+ACPMM*UPMM
     &          +ACMPP*UMPP+ACMPM*UMPM+ACMMP*UMMP+ACMMM*UMMM
          DZUXYZ=ABPPP*UPPP+ABPPM*UPPM+ABPMP*UPMP+ABPMM*UPMM
     &          +ABMPP*UMPP+ABMPM*UMPM+ABMMP*UMMP+ABMMM*UMMM
          DXVXYZ=BCPPP*VPPP+BCPPM*VPPM+BCPMP*VPMP+BCPMM*VPMM
     &          +BCMPP*VMPP+BCMPM*VMPM+BCMMP*VMMP+BCMMM*VMMM
          DYVXYZ=ACPPP*VPPP+ACPPM*VPPM+ACPMP*VPMP+ACPMM*VPMM
     &          +ACMPP*VMPP+ACMPM*VMPM+ACMMP*VMMP+ACMMM*VMMM
          DZVXYZ=ABPPP*VPPP+ABPPM*VPPM+ABPMP*VPMP+ABPMM*VPMM
     &          +ABMPP*VMPP+ABMPM*VMPM+ABMMP*VMMP+ABMMM*VMMM
          DXWXYZ=BCPPP*WPPP+BCPPM*WPPM+BCPMP*WPMP+BCPMM*WPMM
     &          +BCMPP*WMPP+BCMPM*WMPM+BCMMP*WMMP+BCMMM*WMMM
          DYWXYZ=ACPPP*WPPP+ACPPM*WPPM+ACPMP*WPMP+ACPMM*WPMM
     &          +ACMPP*WMPP+ACMPM*WMPM+ACMMP*WMMP+ACMMM*WMMM
          DZWXYZ=ABPPP*WPPP+ABPPM*WPPM+ABPMP*WPMP+ABPMM*WPMM
     &          +ABMPP*WMPP+ABMPM*WMPM+ABMMP*WMMP+ABMMM*WMMM
        ENDIF
C
        IF(ISLRPD.EQ.1)THEN
          XLRPD(NLR,KLR,M)=XTMP+DT*UXYZ
          YLRPD(NLR,KLR,M)=YTMP+DT*VXYZ
          ZLRPD(NLR,KLR,M)=ZTMP+DT*WXYZ
        ENDIF
C
        IF(ISLRPD.EQ.2)THEN
          DTD2=0.5*DT
          AAXX=1.-DTD2*DXUXYZ
          AAXY=-DTD2*DYUXYZ
          AAXZ=-DTD2*DZUXYZ
          FFXX=DT*UXYZ
          AAYX=-DTD2*DXVXYZ
          AAYY=1.-DTD2*DYVXYZ
          AAYZ=-DTD2*DZVXYZ
          FFYY=DT*VXYZ
          AAZX=-DTD2*DXWXYZ
          AAZY=-DTD2*DYWXYZ
          AAZZ=1.-DTD2*DZWXYZ
          FFZZ=DT*WXYZ
          DDET=AAXX*AAYY*AAZZ+AAXY*AAYZ*AAZX+AAXZ*AAYX*AAZY
     &        -AAZX*AAYY*AAXZ-AAZY*AAYZ*AAXX-AAZZ*AAYX*AAXY
          DDETI=1./DDET
          DDXX=FFXX*AAYY*AAZZ+AAXY*AAYZ*FFZZ+AAXZ*FFYY*AAZY
     &        -FFZZ*AAYY*AAXZ-AAZY*AAYZ*FFXX-AAZZ*FFYY*AAXY
          DDYY=AAXX*FFYY*AAZZ+FFXX*AAYZ*AAZX+AAXZ*AAYX*FFZZ
     &        -AAZX*FFYY*AAXZ-FFZZ*AAYZ*AAXX-AAZZ*AAYX*FFXX
          DDZZ=AAXX*AAYY*FFZZ+AAXY*FFYY*AAZX+FFXX*AAYX*AAZY
     &        -AAZX*AAYY*FFXX-AAZY*FFYY*AAXX-FFZZ*AAYX*AAXY
          DDXX=DDXX*DDETI
          DDYY=DDYY*DDETI
          DDZZ=DDZZ*DDETI
          XLRPD(NLR,KLR,M)=XTMP+DDXX
          YLRPD(NLR,KLR,M)=YTMP+DDYY
          ZLRPD(NLR,KLR,M)=ZTMP+DDZZ           
        ENDIF
C
        ENDDO
        ENDDO
      ENDIF
      ENDDO
C
      CLOSE(80)
      IF(N.NE.NTS) RETURN
C
C----------------------------------------------------------------------C
C
      MLRAVG=MLRPDRT+1
      DO KLR=1,KC
      DO NLR=1,NLRPD
      XLRPD(NLR,KLR,MLRAVG)=0.
      YLRPD(NLR,KLR,MLRAVG)=0.
      ZLRPD(NLR,KLR,MLRAVG)=0.
      ENDDO 
      ENDDO
C
      DO M=1,MLRPDRT
      DO KLR=1,KC
      DO NLR=1,NLRPD
      I=ILRPD(NLR)
      J=JLRPD(NLR)
      L=LIJ(I,J)
C     LN=LNC(L)
C     LS=LSC(L)
C     XINIT=DLON(L)
C     YINIT=DLAT(L)
      XTMP=XLRPD(NLR,KLR,M)
      YTMP=YLRPD(NLR,KLR,M)
C     IXTMP=NINT(XTMP)
C     JYTMP=NINT(YTMP)
C     LTMP=LIJ(IXTMP,JYTMP)
C     LNRMP=LNC(LTMP)
C     LSTMP=LSC(LTMP)
C     XDIF=XTMP-FLOAT(IXTMP)
C     YDIF=YTMP-FLOAT(JYTMP)
C     IF(XDIF.GT.0.)THEN
C        XTMP=DLON(LTMP)+XDIF*(DLON(LTMP+1)-DLON(LTMP))
C       ELSE
C        XTMP=DLON(LTMP)+XDIF*(DLON(LTMP)-DLON(LTMP-1))
C     ENDIF
C     IF(YDIF.GT.0.)THEN
C        YTMP=DLAT(LTMP)+YDIF*(DLAT(LNTMP)-DLAT(LTMP))
C       ELSE
C        XTMP=DLAT(LTMP)+YDIF*(DLAT(LTMP)-DLON(LSTMP))
C     ENDIF
C     XLRPD(NLR,KLR,M)=(XTMP-XINIT)/TIDALP
C     YLRPD(NLR,KLR,M)=(YTMP-YINIT)/TIDALP
      XLRPD(NLR,KLR,M)=DXP(L)*(XTMP-FLOAT(I))/TIDALP
      YLRPD(NLR,KLR,M)=DYP(L)*(YTMP-FLOAT(J))/TIDALP
      ZTMP=ZLRPD(NLR,KLR,M)
      K=NINT(ZTMP)
      RKTMP=FLOAT(K)
      DKTMP=ZTMP-RKTMP
      ZTMP=( ZZ(K)+DZC(K)*DKTMP )*HLRPD(NLR,M)+BELV(L)
      ZLRPD(NLR,KLR,M)=(ZTMP-ZILRPD(NLR,KLR,M))/TIDALP
      XLRPD(NLR,KLR,MLRAVG)=XLRPD(NLR,KLR,MLRAVG)+XLRPD(NLR,KLR,M)
      YLRPD(NLR,KLR,MLRAVG)=YLRPD(NLR,KLR,MLRAVG)+YLRPD(NLR,KLR,M)
      ZLRPD(NLR,KLR,MLRAVG)=ZLRPD(NLR,KLR,MLRAVG)+ZLRPD(NLR,KLR,M)
      ENDDO
      ENDDO
      ENDDO
C
      TMPVAL=1./FLOAT(MLRPDRT)
      DO KLR=1,KC
      DO NLR=1,NLRPD
      XLRPD(NLR,KLR,MLRAVG)=TMPVAL*XLRPD(NLR,KLR,MLRAVG)
      YLRPD(NLR,KLR,MLRAVG)=TMPVAL*YLRPD(NLR,KLR,MLRAVG)
      ZLRPD(NLR,KLR,MLRAVG)=TMPVAL*ZLRPD(NLR,KLR,MLRAVG)
      ENDDO
      ENDDO
C
      CALL LVELPLTH
      CALL LVELPLTV
C
C**********************************************************************C
C
  801 FORMAT('  LRPDRT ',I5,'  ACTIVATED AT TIME STEP ',I5)
  802 FORMAT('  LRPDRT ',I5,'  DEACTIVATED AT TIME STEP ',I5)
C
C**********************************************************************C
C
      RETURN 
      END
