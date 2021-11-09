C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BAL2T5
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 05/01/2002        john hamrick       05/01/2002       john hamrick
C  subroutine added for 2 time-level balances including sed,snd,tox
C----------------------------------------------------------------------C
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
C **  INCREMENT COUNTER
C
      IF(ISDYNSTP.EQ.0)THEN
        NBAL=NBAL+1
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        NBAL=NBAL+NINCRMT
        TIME=TIMESEC/TCON 
      ENDIF
C
C **  CHECK FOR END OF BALANCE PERIOD
C
      IF(NBAL.GE.NTSMMT)THEN
C
CJH      WRITE(6,6666)N,NBAL
CJH 6666 FORMAT(' ACTIVE CALL TO CALBAL5, N,NBUD = ',2I5)
C
C**********************************************************************C
C
C **  CALCULATE ENDING VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC 
C **  ENERGY AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES
C
C----------------------------------------------------------------------C
C
      RUMERDO=0.
      RVMERDO=0.
      VOLEND=0.
      VOLEND2T=0.
      BVOLEND2T=0.
      WVOLEND2T=0.
      SALEND=0.
      DYEEND=0.
      UMOEND=0.
      VMOEND=0.
      UUEEND=0.
      VVEEND=0.
      PPEEND=0.
      BBEEND=0.

C
      DYEEND2T=0.0
      DO NS=1,NSED
        SEDEND2T(NS)=0.0
        SEDEND2TW(NS)=0.0
        SEDEND2TB(NS)=0.0
      ENDDO
      DO NS=1,NSND
        SNDEND2T(NS)=0.0
        SNDEND2TW(NS)=0.0
        SNDEND2TB(NS)=0.0
      ENDDO
      DO NT=1,NTOX
        TOXEND2T(NT)=0.0
        TOXEND2TW(NT)=0.0
        TOXEND2TB(NT)=0.0
      ENDDO
C
      DO L=2,LA
      LN=LNC(L)
      VOLEND=VOLEND+SPB(L)*DXYP(L)*HP(L)
      VOLEND2T=VOLEND2T+SPB(L)*DXYP(L)*HP(L)
      WVOLEND2T=WVOLEND2T+SPB(L)*DXYP(L)*HP(L)
      UMOEND=UMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)
     &                                 +DYIU(L+1)*HUI(L+1)*UHDYE(L+1))
      VMOEND=VMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)
     &                                 +DXIV(LN)*HVI(LN)*VHDXE(LN))
      PPEEND=PPEEND+SPB(L)*0.5*DXYP(L)
     &             *(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
      ENDDO
C
      AMOEND=SQRT(UMOEND*UMOEND+VMOEND*VMOEND)
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      SALEND=SALEND+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(K)
      DYEEND=DYEEND+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
      DYEEND2T=DYEEND2T+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
C     UUEEND=UUEEND+SPB(L)*0.25*(DXYU(L)*HU(L)*U(L,K)*U(L,K)
C    &      +DXYU(L+1)*HU(L+1)*U(L+1,K)*U(L+1,K))*DZC(K)
C     VVEEND=VVEEND+SPB(L)*0.25*(DXYV(L)*HV(L)*V(L,K)*V(L,K)
C    &      +DXYV(LN)*HV(LN)*V(LN,K)*V(LN,K))*DZC(K)
      UUEEND=UUEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     &      *( (U(L,K)+U(L+1,K))*(U(L,K)+U(L+1,K)) )
      VVEEND=VVEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     &      *( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEEND=BBEEND+SPB(L)*GP*DXYP(L)*HP(L)*DZC(K)*( BELV(L) 
     &      +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
      ENDDO
      ENDDO
C
C
      DO NS=1,NSED
        DO K=1,KC
        DO L=2,LA
          SEDEND2T(NS)=SEDEND2T(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SED(L,K,NS)
          SEDEND2TW(NS)=SEDEND2TW(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SED(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NS=1,NSND
        DO K=1,KC
        DO L=2,LA
          SNDEND2T(NS)=SNDEND2T(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SND(L,K,NS)
          SNDEND2TW(NS)=SNDEND2TW(NS)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*SND(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NT=1,NTOX
        DO K=1,KC
        DO L=2,LA
          TOXEND2T(NT)=TOXEND2T(NT)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*TOX(L,K,NT)
          TOXEND2TW(NT)=TOXEND2TW(NT)
     &                +SCB(L)*DXYP(L)*HP(L)*DZC(K)*TOX(L,K,NT)
        ENDDO
        ENDDO
      ENDDO
C
      DO NS=1,NSED
        DO L=2,LA
        DO K=1,KBT(L)
          SEDEND2T(NS)=SEDEND2T(NS)+SCB(L)*DXYP(L)*SEDB(L,K,NS)
          SEDEND2TB(NS)=SEDEND2TB(NS)+SCB(L)*DXYP(L)*SEDB(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NS=1,NSND
        DO L=2,LA
        DO K=1,KBT(L)
          SNDEND2T(NS)=SNDEND2T(NS)+SCB(L)*DXYP(L)*SNDB(L,K,NS)
          SNDEND2TB(NS)=SNDEND2TB(NS)+SCB(L)*DXYP(L)*SNDB(L,K,NS)
        ENDDO
        ENDDO
      ENDDO
      DO NT=1,NTOX
        DO L=2,LA
        DO K=1,KBT(L)
          TOXEND2T(NT)=TOXEND2T(NT)+SCB(L)*DXYP(L)*TOXB(L,K,NT)
          TOXEND2TB(NT)=TOXEND2TB(NT)+SCB(L)*DXYP(L)*TOXB(L,K,NT)
        ENDDO
        ENDDO
      ENDDO
C
      DO L=2,LA
      DO K=1,KBT(L)
      VOLEND2T=VOLEND2T+SPB(L)*DXYP(L)*HBED(L,K)
      BVOLEND2T=BVOLEND2T+SPB(L)*DXYP(L)*HBED(L,K)
cjmh013103      WVOLEND2T=WVOLEND2T+SPB(L)*DXYP(L)*PORBED(L,K)*HBED(L,K)
      ENDDO
      ENDDO
C
      ENEBEG=UUEBEG+VVEBEG+PPEBEG+BBEBEG
      ENEEND=UUEEND+VVEEND+PPEEND+BBEEND
      ENEOUT=UUEOUT+VVEOUT+PPEOUT+BBEOUT
C
      VOLBMO=VOLBEG-VOLOUT
      VOLBMO2T=VOLBEG2T-VOLOUT
      BVOLBMO2T=BVOLBEG2T-BVOLOUT
      WVOLBMO2T=WVOLBEG2T-WVOLOUT
      SALBMO=SALBEG-SALOUT
      DYEBMO=DYEBEG-DYEOUT
      DYEBMO2T=DYEBEG2T-DYEOUT2T
      UMOBMO=UMOBEG-DYEOUT
      VMOBMO=VMOBEG-DYEOUT
      ENEBMO=ENEBEG-ENEOUT
C
C  TOXOUT2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
c    TO WATER COLUMN VOLUME SOURCES AND SINKS
C  TOXBLB2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
C    DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
C  TOXFLUXW2T(NT) IS WATER COLUMN SIDE TOXIC FLUX DUE TO SUSPENDED LOAD
C    (POSITIVE INTO WATER COLUMN)
C  TOXFLUXB2T(NT) IS BED SIDE TOXIC FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)
C  TADFLUX2T(NT) IS PORE WATER ADVECTION+DIFFUSION FLUX (POSITIVE INTO WATER COLUMN)
C  TOXFBL2T(NT) IS NET TOXIC FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT
C    (SHOULD EQUAL TOXBLB2T(NT)
C
      DO NT=1,NTOX
        TOXOUT2TW(NT)=TOXOUT2T(NT)-TOXFLUXW2T(NT)-TOXFLUXB2T(NT)
     &               -TADFLUX2T(NT)
        TOXOUT2TB(NT)=TOXFLUXW2T(NT)+TOXFLUXB2T(NT)
     &               +TADFLUX2T(NT)+TOXFBL2T(NT)
C MODIFY TOXOUT2T TO INCLUDE BED LOAD BOUNDARY OUT
         TOXOUT2T(NT)=TOXOUT2T(NT)+TOXBLB2T(NT)
        TOXBMO2T(NT)=TOXBEG2T(NT)-TOXOUT2T(NT)
        TOXBMO2TW(NT)=TOXBEG2TW(NT)-TOXOUT2TW(NT)
        TOXBMO2TB(NT)=TOXBEG2TB(NT)-TOXOUT2TB(NT)
      ENDDO
C
C SEDOUT2T(NS) IS IS NET COHESIVE MASS GOING OUT OF DOMAIN DUE
C   TO WATER COLUMN VOLUME SOURCES AND SINKS
C SEDFLUX2T(NS) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED
C   TO WATER COLUMN
C
      DO NS=1,NSED
        SEDOUT2TW(NS)=SEDOUT2T(NS)-SEDFLUX2T(NS)
        SEDOUT2TB(NS)=SEDFLUX2T(NS)
        SEDBMO2T(NS)=SEDBEG2T(NS)-SEDOUT2T(NS)
        SEDBMO2TW(NS)=SEDBEG2TW(NS)-SEDOUT2TW(NS)
        SEDBMO2TB(NS)=SEDBEG2TB(NS)-SEDOUT2TB(NS)
      ENDDO
C
C  SNDOUT2T(NS) IS NET NONCOHESIVE MASS GOING OUT OF DOMAIN DUE
c    TO WATER COLUMN VOLUME SOURCES AND SINKS
C  SBLOUT2T(NS) IS NET NONCOHESIVE SEDIMENT MASS GOING OUT OF DOMAIN DUE
c    DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
C  SNDFLUX2T(NS) IS NET NONCOHESIVE SEDIMENT FLUX DUE TO SUSPENDED LOAD 
C    (POSITIVE INTO WATER COLUMN)
C  SNDFBL2T(NS) IS NET NONCOHESIVE SEDIMENT FLUX FROM BED ASSOCIATED WITH 
C    BED LOAD TRANSPORT (SHOULD EQUAL SBLOUT2T(NSX))
C
      DO NS=1,NSND
        SNDOUT2TW(NS)=SNDOUT2T(NS)-SNDFLUX2T(NS)
        SNDOUT2TB(NS)=SNDFLUX2T(NS)+SNDFBL2T(NS)
C  MODIFY SNDOUT2T TO INCLUDE BED LOAD TRANSPORT OUT OF DOMAIN
        SNDOUT2T(NS)=SNDOUT2T(NS)+SBLOUT2T(NS)
        SNDBMO2T(NS)=SNDBEG2T(NS)-SNDOUT2T(NS)
        SNDBMO2TW(NS)=SNDBEG2TW(NS)-SNDOUT2TW(NS)
        SNDBMO2TB(NS)=SNDBEG2TB(NS)-SNDOUT2TB(NS)
      ENDDO
C
      VOLERR=VOLEND-VOLBMO
      VOLERR2T=VOLEND2T-VOLBMO2T
      BVOLERR2T=BVOLEND2T-BVOLBMO2T
      WVOLERR2T=WVOLEND2T-WVOLBMO2T
      SALERR=SALEND-SALBMO
      DYEERR=DYEEND-DYEBMO
      DYEERR2T=DYEEND2T-DYEBMO2T
      UMOERR=UMOEND-UMOBMO
      VMOERR=VMOEND-VMOBMO
      ENEERR=ENEEND-ENEBMO
C
      DO NS=1,NSED
        SEDERR2T(NS)=SEDEND2T(NS)-SEDBMO2T(NS)
        SEDERR2TW(NS)=SEDEND2TW(NS)-SEDBMO2TW(NS)
        SEDERR2TB(NS)=SEDEND2TB(NS)-SEDBMO2TB(NS)
      ENDDO
      DO NS=1,NSND
        SNDERR2T(NS)=SNDEND2T(NS)-SNDBMO2T(NS)
        SNDERR2TW(NS)=SNDEND2TW(NS)-SNDBMO2TW(NS)
        SNDERR2TB(NS)=SNDEND2TB(NS)-SNDBMO2TB(NS)
      ENDDO
      DO NT=1,NTOX
        TOXERR2T(NT)=TOXEND2T(NT)-TOXBMO2T(NT)
        TOXERR2TW(NT)=TOXEND2TW(NT)-TOXBMO2TW(NT)
        TOXERR2TB(NT)=TOXEND2TB(NT)-TOXBMO2TB(NT)
      ENDDO
C
      RVERDE=-9999.
      RVERDE2T=-9999.
      RWVERDE2T=-9999.
      RSERDE=-9999.
      RDERDE=-9999.
      RDERDE=-9999.
      RUERDE2T=-9999.
      RVERDE=-9999.
      REERDE=-9999.
C
      DO NS=1,NSED
        RSEDERE2T(NS)=-9999.
        RSEDERE2TW(NS)=-9999.
        RSEDERE2TB(NS)=-9999.
      ENDDO
      DO NS=1,NSND
        RSNDERE2T(NS)=-9999.
        RSNDERE2TW(NS)=-9999.
        RSNDERE2TB(NS)=-9999.
      ENDDO
      DO NT=1,NTOX
        RTOXERE2T(NT)=-9999.
        RTOXERE2TW(NT)=-9999.
        RTOXERE2TB(NT)=-9999.
      ENDDO
C
      RVERDO=-9999.
      RVERDO2T=-9999.
      RWVERDO2T=-9999.
      RSERDO=-9999.
      RDERDO=-9999.
      RDERDO2T=-9999.
      RUERDO=-9999.
      RVERDO=-9999.
      REERDO=-9999.
C
      DO NS=1,NSED
        RSEDERO2T(NS)=-9999.
        RSEDERO2TW(NS)=-9999.
        RSEDERO2TB(NS)=-9999.
      ENDDO
      DO NS=1,NSND
        RSNDERO2T(NS)=-9999.
        RSNDERO2TW(NS)=-9999.
        RSNDERO2TB(NS)=-9999.
      ENDDO
      DO NT=1,NTOX
        RTOXERO2T(NT)=-9999.
        RTOXERO2TW(NT)=-9999.
        RTOXERO2TB(NT)=-9999.
      ENDDO
C
      RVERDE=-9999.
      RVERDE2T=-9999.
      RBVERDE2T=-9999.
      RWVERDE2T=-9999.
      RSERDE=-9999.
      RDERDE=-9999.
      RDERDE2T=-9999.
      RUMERDE=-9999.
      RVMERDE=-9999.
      REERDE=-9999.
C
      IF(VOLEND.NE.0.) RVERDE=VOLERR/VOLEND
      IF(VOLEND2T.NE.0.) RVERDE2T=VOLERR2T/VOLEND2T
      IF(BVOLEND2T.NE.0.) RBVERDE2T=BVOLERR2T/BVOLEND2T
      IF(WVOLEND2T.NE.0.) RWVERDE2T=WVOLERR2T/WVOLEND2T
      IF(SALEND.NE.0.) RSERDE=SALERR/SALEND
      IF(DYEEND.NE.0.) RDERDE=DYEERR/DYEEND
      IF(DYEEND2T.NE.0.) RDERDE2T=DYEERR2T/DYEEND2T
      IF(UMOEND.NE.0.) RUMERDE=UMOERR/UMOEND
      IF(VMOEND.NE.0.) RVMERDE=VMOERR/VMOEND
      IF(ENEEND.NE.0.) REERDE=ENEERR/ENEEND
C
c      DO NS=1,NSED
c        IF(SEDEND2T(NS).NE.0.)
c     &    RSEDERE2T(NS)=SEDERR2T(NS)/SEDEND2T(NS)
c        IF(SEDEND2TW(NS).NE.0.)
c     &    RSEDERE2TW(NS)=SEDERR2TW(NS)/SEDEND2TW(NS)
c        IF(SEDEND2TB(NS).NE.0.)
c     &    RSEDERE2TB(NS)=SEDERR2TB(NS)/SEDEND2TB(NS)
c      ENDDO
c      DO NS=1,NSND
c        IF(SNDEND2T(NS).NE.0.)
c     &    RSNDERE2T(NS)=SNDERR2T(NS)/SNDEND2T(NS)
c        IF(SNDEND2TW(NS).NE.0.)
c     &    RSNDERE2TW(NS)=SNDERR2TW(NS)/SNDEND2TW(NS)
c        IF(SNDEND2TB(NS).NE.0.)
c     &    RSNDERE2TB(NS)=SNDERR2TB(NS)/SNDEND2TB(NS)
c      ENDDO
c      DO NT=1,NTOX
c        IF(TOXEND2T(NT).NE.0.)
c     &    RTOXERE2T(NT)=TOXERR2T(NT)/TOXEND2T(NT)
c        IF(TOXEND2TW(NT).NE.0.)
c     &    RTOXERE2TW(NT)=TOXERR2TW(NT)/TOXEND2TW(NT)
c        IF(TOXEND2TB(NT).NE.0.)
c     &    RTOXERE2TB(NT)=TOXERR2TB(NT)/TOXEND2TB(NT)
c      ENDDO
C
C
      DO NS=1,NSED
        TMPVAL=0.5*(SEDEND2T(NS)+SEDBEG2T(NS))
        IF(TMPVAL.NE.0.)
     &    RSEDERE2T(NS)=SEDERR2T(NS)/TMPVAL
        TMPVAL=0.5*(SEDEND2TW(NS)+SEDBEG2TW(NS))
        IF(TMPVAL.NE.0.)
     &    RSEDERE2TW(NS)=SEDERR2TW(NS)/TMPVAL
        TMPVAL=0.5*(SEDEND2TB(NS)+SEDBEG2TB(NS))
        IF(TMPVAL.NE.0.)
     &    RSEDERE2TB(NS)=SEDERR2TB(NS)/TMPVAL
      ENDDO
      DO NS=1,NSND
        TMPVAL=0.5*(SNDEND2T(NS)+SNDBEG2T(NS))
        IF(TMPVAL.NE.0.)
     &    RSNDERE2T(NS)=SNDERR2T(NS)/TMPVAL
        TMPVAL=0.5*(SNDEND2TW(NS)+SNDBEG2TW(NS))
        IF(TMPVAL.NE.0.)
     &    RSNDERE2TW(NS)=SNDERR2TW(NS)/TMPVAL
        TMPVAL=0.5*(SNDEND2TB(NS)+SNDBEG2TB(NS))
        IF(TMPVAL.NE.0.)
     &    RSNDERE2TB(NS)=SNDERR2TB(NS)/TMPVAL
      ENDDO
      DO NT=1,NTOX
        TMPVAL=0.5*(TOXEND2T(NT)+TOXBEG2T(NT))
        IF(TMPVAL.NE.0.)
     &    RTOXERE2T(NT)=TOXERR2T(NT)/TMPVAL
        TMPVAL=0.5*(TOXEND2TW(NT)+TOXBEG2TW(NT))
        IF(TMPVAL.NE.0.)
     &    RTOXERE2TW(NT)=TOXERR2TW(NT)/TMPVAL
        TMPVAL=0.5*(TOXEND2TB(NT)+TOXBEG2TB(NT))
        IF(TMPVAL.NE.0.)
     &    RTOXERE2TB(NT)=TOXERR2TB(NT)/TMPVAL
      ENDDO
C
      IF(VOLOUT.NE.0.) RVERDO=VOLERR/VOLOUT
      IF(VOLOUT.NE.0.) RVERDO2T=VOLERR2T/VOLOUT
      IF(VOLOUT.NE.0.) RWVERDO2T=WVOLERR2T/VOLOUT
      IF(SALOUT.NE.0.) RSERDO=SALERR/SALOUT
      IF(DYEOUT.NE.0.) RDERDO=DYEERR/DYEOUT
      IF(DYEOUT2T.NE.0.) RDERDO2T=DYEERR2T/DYEOUT2T
      IF(UMOOUT.NE.0.) RUMERDO=UMOERR/UMOOUT
      IF(VMOOUT.NE.0.) RVMERDO=VMOERR/VMOOUT
      IF(ENEOUT.NE.0.) REERDO=ENEERR/ENEOUT
C
      DO NS=1,NSED
        IF(SEDOUT2T(NS).NE.0.)
     &    RSEDERO2T(NS)=SEDERR2T(NS)/SEDOUT2T(NS)
        IF(SEDOUT2TW(NS).NE.0.)
     &    RSEDERO2TW(NS)=SEDERR2TW(NS)/SEDOUT2TW(NS)
        IF(SEDOUT2TB(NS).NE.0.)
     &    RSEDERO2TB(NS)=SEDERR2TB(NS)/SEDOUT2TB(NS)
      ENDDO
      DO NS=1,NSND
        IF(SNDOUT2T(NS).NE.0.)
     &    RSNDERO2T(NS)=SNDERR2T(NS)/SNDOUT2T(NS)
        IF(SNDOUT2TW(NS).NE.0.)
     &    RSNDERO2TW(NS)=SNDERR2TW(NS)/SNDOUT2TW(NS)
        IF(SNDOUT2TB(NS).NE.0.)
     &    RSNDERO2TB(NS)=SNDERR2TB(NS)/SNDOUT2TB(NS)
      ENDDO
      DO NT=1,NTOX
        IF(TOXOUT2T(NT).NE.0.)
     &    RTOXERO2T(NT)=TOXERR2T(NT)/TOXOUT2T(NT)
        IF(TOXOUT2TW(NT).NE.0.)
     &    RTOXERO2TW(NT)=TOXERR2TW(NT)/TOXOUT2TW(NT)
        IF(TOXOUT2TB(NT).NE.0.)
     &    RTOXERO2TB(NT)=TOXERR2TB(NT)/TOXOUT2TB(NT)
      ENDDO
C
C**********************************************************************C
C
C **  OUTPUT BALANCE RESULTS TO FILE BAL2T.OUT
C
C----------------------------------------------------------------------C
C
      IF(JSBAL.EQ.1)THEN
        OPEN(89,FILE='BAL2T.OUT')
        CLOSE(89,STATUS='DELETE')
        OPEN(89,FILE='BAL2T.OUT')
        OPEN(81,FILE='BAL2TERSTT.OUT')
        CLOSE(81,STATUS='DELETE')
        OPEN(81,FILE='BAL2TERSTT.OUT')
        OPEN(82,FILE='BAL2TERSTW.OUT')
        CLOSE(82,STATUS='DELETE')
        OPEN(82,FILE='BAL2TERSTW.OUT')
        OPEN(83,FILE='BAL2TERSTB.OUT')
        CLOSE(83,STATUS='DELETE')
        OPEN(83,FILE='BAL2TERSTB.OUT')
        OPEN(84,FILE='BAL2TERVWT.OUT')
        CLOSE(84,STATUS='DELETE')
        OPEN(84,FILE='BAL2TERVWT.OUT')
        JSBAL=0
       ELSE
        OPEN(89,FILE='BAL2T.OUT',POSITION='APPEND')
        OPEN(81,FILE='BAL2TERSTT.OUT',POSITION='APPEND')
        OPEN(82,FILE='BAL2TERSTW.OUT',POSITION='APPEND')
        OPEN(83,FILE='BAL2TERSTB.OUT',POSITION='APPEND')
        OPEN(84,FILE='BAL2TERVWT.OUT',POSITION='APPEND')
      ENDIF
C
      WRITE(89,899)
      WRITE(89,890)NTSMMT,TIME
      WRITE(89,891)
      WRITE(89,892)VOLBEG,SALBEG,DYEBEG,ENEBEG,UMOBEG,VMOBEG,AMOBEG
      WRITE(89,900)
      WRITE(89,893)
      WRITE(89,892)VOLOUT,SALOUT,DYEOUT,ENEOUT,UMOOUT,VMOOUT
      WRITE(89,900)
      WRITE(89,894)
      WRITE(89,892)VOLBMO,SALBMO,DYEBMO,ENEBMO,UMOBMO,VMOBMO
      WRITE(89,900)
      WRITE(89,895)
      WRITE(89,892)VOLEND,SALEND,DYEEND,ENEEND,UMOEND,VMOEND,AMOEND
      WRITE(89,900)
      WRITE(89,896)
      WRITE(89,892)VOLERR,SALERR,DYEERR,ENEERR,UMOERR,VMOERR
      WRITE(89,900)
      WRITE(89,897)
      WRITE(89,892)RVERDE,RSERDE,RDERDE,REERDE,RUMERDE,RVMERDE
      WRITE(89,900)
      WRITE(89,898)
      WRITE(89,892)RVERDO,RSERDO,RDERDO,REERDO,RUMERDO,RVMERDO
      WRITE(89,899)
      UUEBMO=UUEBEG-UUEOUT
      VVEBMO=VVEBEG-VVEOUT
      PPEBMO=PPEBEG-PPEOUT
      BBEBMO=BBEBEG-BBEOUT
      WRITE(89,901)UUEBEG
      WRITE(89,902)UUEOUT
      WRITE(89,903)UUEBMO
      WRITE(89,904)UUEEND
      WRITE(89,900)
      WRITE(89,905)VVEBEG
      WRITE(89,906)VVEOUT
      WRITE(89,907)VVEBMO
      WRITE(89,908)VVEEND
      WRITE(89,900)
      WRITE(89,909)PPEBEG
      WRITE(89,910)PPEOUT
      WRITE(89,911)PPEBMO
      WRITE(89,912)PPEEND
      WRITE(89,900)
      WRITE(89,913)BBEBEG
      WRITE(89,914)BBEOUT
      WRITE(89,915)BBEBMO
      WRITE(89,916)BBEEND
      WRITE(89,900)
      WRITE(89,899)
C
      WRITE(81,8888)TIME,(RSEDERE2T(NS),NS=1,NSED),
     &      (RSNDERE2T(NS),NS=1,NSND),(RTOXERE2T(NT),NT=1,NTOX)
C
      WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE W_COL+BED'
      WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
      WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
      WRITE(89,*)' _BMO_ = BEGINNING - OUT'
      WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
      WRITE(89,*)' _ERR_ = END - BMO'
      WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
      WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
      WRITE(89,900)
      WRITE(89,949)
      NS=0
        WRITE(89,950)NS,DYEBEG2T,DYEOUT2T,DYEBMO2T,
     &             DYEEND2T,DYEERR2T,RDERDE2T,RDERDO2T
      WRITE(89,900)
      WRITE(89,951)
      DO NS=1,NSED
        WRITE(89,950)NS,SEDBEG2T(NS),SEDOUT2T(NS),SEDBMO2T(NS),
     &  SEDEND2T(NS),SEDERR2T(NS),RSEDERE2T(NS),RSEDERO2T(NS),
     &  SEDFLUX2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,952)
      DO NS=1,NSND
        WRITE(89,950)NS,SNDBEG2T(NS),SNDOUT2T(NS),SNDBMO2T(NS),
     &  SNDEND2T(NS),SNDERR2T(NS),RSNDERE2T(NS),RSNDERO2T(NS),
     &  SNDFLUX2T(NS),SBLOUT2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,953)
      DO NT=1,NTOX
        WRITE(89,950)NT,TOXBEG2T(NT),TOXOUT2T(NT),TOXBMO2T(NT),
     &  TOXEND2T(NT),TOXERR2T(NT),RTOXERE2T(NT),RTOXERO2T(NT),
     &  TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
      ENDDO
      WRITE(89,900)
C
      WRITE(82,8888)TIME,(RSEDERE2TW(NS),NS=1,NSED),
     &      (RSNDERE2TW(NS),NS=1,NSND),(RTOXERE2TW(NT),NT=1,NTOX)
C
      WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE WATER COL'
      WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
      WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
      WRITE(89,*)' _BMO_ = BEGINNING - OUT'
      WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
      WRITE(89,*)' _ERR_ = END - BMO'
      WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
      WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
      WRITE(89,900)
      WRITE(89,951)
      DO NS=1,NSED
        WRITE(89,950)NS,SEDBEG2TW(NS),SEDOUT2TW(NS),SEDBMO2TW(NS),
     &  SEDEND2TW(NS),SEDERR2TW(NS),RSEDERE2TW(NS),RSEDERO2TW(NS),
     &  SEDFLUX2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,952)
      SBLOUT2TT=0.0
      DO NS=1,NSND
        DSEDGMM=1./(1.E6*SSG(NS+NSED))
        WRITE(89,950)NS,SNDBEG2TW(NS),SNDOUT2TW(NS),SNDBMO2TW(NS),
     &  SNDEND2TW(NS),SNDERR2TW(NS),RSNDERE2TW(NS),RSNDERO2TW(NS),
     &  SNDFLUX2T(NS),SBLOUT2T(NS)
        SBLOUT2TT=SBLOUT2TT+DSEDGMM*SBLOUT2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,953)
      DO NT=1,NTOX
        WRITE(89,950)NT,TOXBEG2TW(NT),TOXOUT2TW(NT),TOXBMO2TW(NT),
     &  TOXEND2TW(NT),TOXERR2TW(NT),RTOXERE2TW(NT),RTOXERO2TW(NT),
     &  TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
      ENDDO
      WRITE(89,900)
C
      WRITE(83,8888)TIME,(RSEDERE2TB(NS),NS=1,NSED),
     &      (RSNDERE2TB(NS),NS=1,NSND),(RTOXERE2TB(NT),NT=1,NTOX)
C
      WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE  BED'
      WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
      WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
      WRITE(89,*)' _BMO_ = BEGINNING - OUT'
      WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
      WRITE(89,*)' _ERR_ = END - BMO'
      WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
      WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
      WRITE(89,900)
      WRITE(89,951)
      DO NS=1,NSED
        WRITE(89,950)NS,SEDBEG2TB(NS),SEDOUT2TB(NS),SEDBMO2TB(NS),
     &  SEDEND2TB(NS),SEDERR2TB(NS),RSEDERE2TB(NS),RSEDERO2TB(NS),
     &  SEDFLUX2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,952)
      DO NS=1,NSND
        WRITE(89,950)NS,SNDBEG2TB(NS),SNDOUT2TB(NS),SNDBMO2TB(NS),
     &  SNDEND2TB(NS),SNDERR2TB(NS),RSNDERE2TB(NS),RSNDERO2TB(NS),
     &  SNDFLUX2T(NS),SBLOUT2T(NS)
      ENDDO
      WRITE(89,900)
      WRITE(89,953)
      DO NT=1,NTOX
        WRITE(89,950)NT,TOXBEG2TB(NT),TOXOUT2TB(NT),TOXBMO2TB(NT),
     &  TOXEND2TB(NT),TOXERR2TB(NT),RTOXERE2TB(NT),RTOXERO2TB(NT),
     &  TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
      ENDDO
      WRITE(89,900)
C
      WRITE(84,8888)TIME,RVERDE2T,RWVERDE2T,RBVERDE2T
C
      WRITE(89,*)' NEW TOTAL VOLUME (W_C+BED,WAT+SED) VOLUME BALANCE'
      WRITE(89,*)' _BEG_ = BEGINNING GLOBAL VOLUME'
      WRITE(89,*)' _OUT_ = NET VOLUME GOING OUT OF DOMAIN'
      WRITE(89,*)' _BMO_ = BEGINNING - OUT'
      WRITE(89,*)' _END_ = ENDING GLOBAL VOLUME'
      WRITE(89,*)' _ERR_ = END - BMO'
      WRITE(89,*)' _R_E_ = RELATIVE ERROR = ERR/END'
      WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
      WRITE(89,900)
      WRITE(89,954)
      WRITE(89,960)VOLBEG2T,VOLOUT,VOLBMO2T,
     &             VOLEND2T,VOLERR2T,RVERDE2T,RVERDO2T,SBLOUT2TT
      WRITE(89,900)
C      WRITE(89,*)' NEW WATER VOLUME (W_COL+BED WATER) VOLUME BALANCE'
      WRITE(89,*)' NEW WATER VOLUME (W_COL,WATER+SED) VOLUME BALANCE'
C        WRITE(89,*)' _BEG_ = BEGINNING GLOBAL WATER VOLUME'
      WRITE(89,*)' _BEG_ = BEGINNING GLOBAL W_COL VOLUME'
C       WRITE(89,*)' _OUT_ = NET WATER VOLUME GOING OUT OF DOMAIN'
      WRITE(89,*)' _OUT_ = NET W_COL VOLUME GOING OUT OF DOMAIN'
      WRITE(89,*)' _BMO_ = BEGINNING - OUT'
C         WRITE(89,*)' _END_ = ENDING GLOBAL WATER VOLUME'
      WRITE(89,*)' _END_ = ENDING GLOBAL W_COL VOLUME'
      WRITE(89,*)' _ERR_ = END - BMO'
      WRITE(89,*)' _R_E_ = RELATIVE ERROR = ERR/END'
      WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
      WRITE(89,900)
      WRITE(89,955)
      WRITE(89,960)WVOLBEG2T,WVOLOUT,WVOLBMO2T,
     &             WVOLEND2T,WVOLERR2T,RWVERDE2T,RWVERDO2T,VOLMORPH2T
C
      WRITE(89,900)
      WRITE(89,*)' SEE NOTES IN BAL2T5 FOR INTERPRETATION OF FOLLOWING'
      WRITE(89,900)
C
      DO NS=1,NSED
      WRITE(89,900)
      WRITE(89,8899)'SEDBEG2T(NS)    = ',SEDBEG2T(NS)
      WRITE(89,8899)'SEDBEG2TB(NS)   = ',SEDBEG2TB(NS)
      WRITE(89,8899)'SEDBEG2TW(NS)   = ',SEDBEG2TW(NS)
      WRITE(89,8899)'SEDEND2T(NS)    = ',SEDEND2T(NS)
      WRITE(89,8899)'SEDEND2TB(NS)   = ',SEDEND2TB(NS)
      WRITE(89,8899)'SEDEND2TW(NS)   = ',SEDEND2TW(NS) 
      WRITE(89,8899)'SEDOUT2T(NS)    = ',SEDOUT2T(NS)
      WRITE(89,8899)'SEDFLUX2T(NS)   = ',SEDFLUX2T(NS)
      ENDDO
C
      DO NS=1,NSND
C  MODIFY SNDOUT2T BACK TO ORGINAL DEFINITION
      SNDOUT2T(NS)=SNDOUT2T(NS)-SBLOUT2T(NS)
      WRITE(89,900)
      WRITE(89,8899)'SNDBEG2T(NS)    = ',SNDBEG2T(NS)
      WRITE(89,8899)'SNDBEG2TB(NS)   = ',SNDBEG2TB(NS)
      WRITE(89,8899)'SNDBEG2TW(NS)   = ',SNDBEG2TW(NS)
      WRITE(89,8899)'SNDEND2T(NS)    = ',SNDEND2T(NS)
      WRITE(89,8899)'SNDEND2TB(NS)   = ',SNDEND2TB(NS)
      WRITE(89,8899)'SNDEND2TW(NS)   = ',SNDEND2TW(NS) 
      WRITE(89,8899)'SNDOUT2T(NS)    = ',SNDOUT2T(NS)
      WRITE(89,8899)'SNDFLUX2T(NS)   = ',SNDFLUX2T(NS)
      WRITE(89,8899)'SNDFBL2T(NS)    = ',SNDFBL2T(NS) 
      WRITE(89,8899)'SBLOUT2T(NS)    = ',SBLOUT2T(NS)
      ENDDO
C
      DO NT=1,NTOX
C MODIFY TOXOUT2T BACK TO ORIGINAL DEFINITION
       TOXOUT2T(NT)=TOXOUT2T(NT)-TOXBLB2T(NT)
      WRITE(89,900)
      WRITE(89,8899)'TOXBEG2T(NT)    = ',TOXBEG2T(NT)
      WRITE(89,8899)'TOXBEG2TB(NT)   = ',TOXBEG2TB(NT)
      WRITE(89,8899)'TOXBEG2TW(NT)   = ',TOXBEG2TW(NT)
      WRITE(89,8899)'TOXEND2T(NT)    = ',TOXEND2T(NT)
      WRITE(89,8899)'TOXEND2TB(NT)   = ',TOXEND2TB(NT)
      WRITE(89,8899)'TOXEND2TW(NT)   = ',TOXEND2TW(NT)
      WRITE(89,8899)'TOXOUT2T(NT)    = ',TOXOUT2T(NT)
      WRITE(89,8899)'TOXFLUXW2T(NT)  = ',TOXFLUXW2T(NT)
      WRITE(89,8899)'TOXFLUXB2T(NT)  = ',TOXFLUXB2T(NT)
      WRITE(89,8899)'TOXFBL2T(NT)    = ',TOXFBL2T(NT)
      WRITE(89,8899)'TOXBLB2T(NT)    = ',TOXBLB2T(NT)
      WRITE(89,8899)'TADFLUX2T(NT)   = ',TADFLUX2T(NT)
      ENDDO
C
      CLOSE(89)
      CLOSE(81)
      CLOSE(82)
      CLOSE(83)
      CLOSE(84)
C
 8899 FORMAT(A18,E15.7)
  950 FORMAT(I5,12E17.9)
  960 FORMAT(5X,12E17.9)
  949 FORMAT('  NDYE DYEBEG2T         DYEOUT2T         DYEBMO2T     ',
     &'    DYEEND2T         DYEERR2T         RDYEERE2T',
     &'        RDYEERO2T    ',/)
  951 FORMAT('  NSED SEDBEG2T(NS)     SEDOUT2T(NS)      SEDBMO2T(NS)  ',
     &  '  SEDEND2T(NS)     SEDERR2T(NS)     RSEDERE2T(NS)',
     &  '    RSEDERO2T(NS)',
     &       '    SEDFLUX2T(NS)',/)
  952 FORMAT('  NSND SNDBEG2T(NS)     SNDOUT2T(NS)     SNDBMO2T(NS)  ',
     &  '   SNDEND2T(NS)     SNDERR2T(NS)     RSNDERE2T(NS)',
     &  '    RSNDERO2T(NS)',
     &       '    SNDFLUX2T(NS)    SBLOUT2T(NS)',/)
  953 FORMAT('   NT  TOXBEG2T(NT)     TOXOUT2T(NT)     TOXBMO2T(NT)  ',
     &  '   TOXEND2T(NT)     TOXERR2T(NT)     RTOXERE2T(NT)',
     &  '    RTOXERO2T(NT)',
     &       '    TOXFLUXW2T(NT)    TOXFLUXB2T(NT)   TADFLUX2T(NT)',/)
  954 FORMAT('       VOLBEG2T         VOLOUT           VOLBMO2T ',
     &     '       ',
     &     ' VOLEND2T         VOLERR2T         RVERDE2T',
     &     '         RVERDO2T       VOLBLOUT',/)
  955 FORMAT('       WVOLBEG2T        WVOLOUT          WVOLBMO2T ',
     &     '       ',
     &    ' WVOLEND2T       WVOLERR2T        RWVERDE2T',
     &    '        RWVERDO2T        VOLMORPH2T',/)
  890 FORMAT (' VOLUME, MASS, AND ENERGY BALANCE OVER',I12,' TIME STEPS'
     &,' ENDING AT TIME ',F12.4,//)
  891 FORMAT (' INITIAL VOLUME    INITIAL SALT    INITIAL DYE     '
     &,'INITIAL ENER    INITIAL UMO     INITIAL VMO     '
     &,'INITIAL AMO',/)
  892 FORMAT (1X,7(E14.6,2X))
  893 FORMAT (' VOLUME OUT        SALT OUT        DYE OUT         '
     &,'ENERGY OUT      UMO OUT         VMO OUT',/)
  894 FORMAT (' INITIAL-OUT VOL   INIT-OUT SALT   INIT-OUT DYE    '
     &,'INIT-OUT ENER   INIT-OUT UMO    INIT-OUT VMO',/)
  895 FORMAT (' FINAL VOLUME      FINAL SALT      FINAL DYE       '
     &,'FINAL ENERGY    FINAL UMO       FINAL VMO       '
     &,'FINAL AMO',/)
  896 FORMAT (' VOLUME ERR        SALT ERR        DYE ERR         '
     &,'ENERGY ERR      UMO ERR         VMO ERR',/)
  897 FORMAT (' R VOL/END ER      R SAL/END ER    R DYE/END ER    '
     &,'R ENE/END ER    R UMO/END ER    R VMO/END ER',/)
  898 FORMAT (' R VOL/OUT ER      R SAL/OUT ER    R DYE/OUT ER    '
     &,'R ENE/OUT ER    R UMO/OUT ER    R VMO/OUT ER',/)
  899 FORMAT (////)
  900 FORMAT (//)
  901 FORMAT(' UUEBEG =  ',E14.6)
  902 FORMAT(' UUEOUT =  ',E14.6)
  903 FORMAT(' UUEBMO =  ',E14.6)
  904 FORMAT(' UUEEND =  ',E14.6)
  905 FORMAT(' VVEBEG =  ',E14.6)
  906 FORMAT(' VVEOUT =  ',E14.6)
  907 FORMAT(' VVEBMO =  ',E14.6)
  908 FORMAT(' VVEEND =  ',E14.6)
  909 FORMAT(' PPEBEG =  ',E14.6)
  910 FORMAT(' PPEOUT =  ',E14.6)
  911 FORMAT(' PPEBMO =  ',E14.6)
  912 FORMAT(' PPEEND =  ',E14.6)
  913 FORMAT(' BBEBEG =  ',E14.6)
  914 FORMAT(' BBEOUT =  ',E14.6)
  915 FORMAT(' BBEBMO =  ',E14.6)
  916 FORMAT(' BBEEND =  ',E14.6)
 8888 FORMAT(F10.2,10E14.6)
C
C**********************************************************************C
C
C     RESET COUNTER
C
      NBAL=0
C
      ENDIF 
C
C**********************************************************************C
C
      RETURN
      END
