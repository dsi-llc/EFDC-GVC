C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION CSEDTAUS(DENBULK,TAUCO,VDRO,VDR,VDRC,IOPT,L)
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
      INCLUDE 'EFDC.PAR'
C
c#######################################################################
C  HQI Added, 11/18/2003, HAMRICK COMMENTED OUT SINCE NOT NEEDED FOR
C  FOR 11/24 VERSION OF IOPT= 99 THAT IS ACTIVE AS OF 01/08/2004
C
CJH      INCLUDE 'EFDC.CMN'
CJH      REAL rSEDPHI(4),rSIGPHI,D50SIG,D90SIG,rSNDBT,cohcon,rfracoh
CJH      REAL rTVAR3W,rTVAR3E,Z90
CJH      real rfrac1,rfrac2,rfrac3,ubnd,lbnd,usiz,lsiz
C
c#######################################################################
C
C **  CALCULATES CRITIAL STRESS FOR SURFACE EROSION OF COHESIVE 
C **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
C
C **  IOPT=1  BASED ON 
C **
C **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
C **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
C **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
C
C **  IOPT=2 & 3  BASED ON J. M. HAMRICK'S MODIFICATION OF
C **
C **  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATION
C **  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
C
C **  IOPT=4  BASED ON J. M. HAMRICK'S PARAMETERIZATION OF SEDFLUME
C     TEST DATA
C **
C
      IF(IOPT.EQ.1)THEN
        DENBULK=0.001*DENBULK
        IF(DENBULK.LE.1.065)THEN
          CSEDTAUS=0.0
        ELSE
          TMP=(DENBULK-1.065)**0.2
          CSEDTAUS=0.001*(0.883*TMP+0.05)
        ENDIF
      ENDIF
C
      IF(IOPT.EQ.2)THEN
        CSEDTAUS=TAUCO*(1.+VDRO)/(1.+VDR)
      ENDIF
C
      IF(IOPT.EQ.3)THEN
        CSEDTAUS=TAUCO*(1.+VDRO)/(1.+VDRC)
      ENDIF
C
c      IF(IOPT.EQ.3)THEN
c        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
c      ENDIF
C
      IF(IOPT.EQ.4)THEN
c        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
        CSEDTAUS=TAUCO
      ENDIF
C
      IF(IOPT.EQ.5)THEN
c        CSEDTAUS=TAUCO*EXP((1.+VDRO)/(1.+VDRC))
        CSEDTAUS=TAUCO
      ENDIF
C

c#######################################################################
C  HQI Change, 08/25/03, and 11/24/03  SO and RM
C  Change to implement critical shear stress option
C  CSEDTAUS is tau/rho with tau in dyne/cm^2
C  IWRSP(1) = 99 for shear stress as a function of bulk density from
C                sed-flume experiments
      IF(IOPT.GE.99)THEN
c     Bulk density in g/cc and Tau in N/m**2
         IF(L.LE.265) then
C     Woods Pond cells
c     CSEDTAUS=0.015*((DENBULK/1000.)**13.1)/1000
c     CSEDTAUS=1.2/1000.
            CSEDTAUS=0.2/1000.
         ELSE
C     North of Woods Pond cells
c     CSEDTAUS=8.0E-6*((DENBULK/1000.)**17.4)/1000.
c     CSEDTAUS=0.6/1000.
            CSEDTAUS=0.4/1000.
         ENDIF
      ENDIF
C
c#######################################################################
C
c#######################################################################
C  HQI Added, 11/18/2003
C  Compute the D90 of a cell bed from the fractions of each grain size
C  and the computed D50's (i.e. the Deff for the four classes.
c
c      IF(IOPT.EQ.9999)THEN
c      KTOP=KBT(L)
c      SEDDIA(1)=22.0*1.e-6 !convert micron to meter
c
c      rfracoh= SEDB(L,KTOP,1)/( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)
c     &                          +SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
c      rfrac1 = rfracoh + SNDB(L,KTOP,1)
c     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
c      rfrac2 = rfrac1 + SNDB(L,KTOP,2)
c     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
c      rfrac3 = rfrac2 + SNDB(L,KTOP,3)
c     &  /( SEDB(L,KTOP,1)+SNDB(L,KTOP,1)+SNDB(L,KTOP,2)+SNDB(L,KTOP,3) )
c
c      ubnd = 999.
c      lbnd = -999.
c
c 6400, 570, 160, 63 in microns
c
c      if ( rfrac3 .GT. 0.9 .AND. rfrac3 .LT. ubnd ) THEN
c        ubnd = rfrac3
c        usiz = SEDDIA(4)
c      END if
c      if ( rfrac2 .GT. 0.9 .AND. rfrac2 .LT. ubnd ) THEN
c        ubnd = rfrac2
c        usiz = SEDDIA(3)
c      END if
c      if ( rfrac1 .GT. 0.9 .AND. rfrac1 .LT. ubnd ) THEN
c        ubnd = rfrac1
c        usiz = SEDDIA(2)
c      END if
c      if ( rfracoh .GT. 0.9 .AND. rfracoh .LT. ubnd ) THEN
c        ubnd = rfracoh
c        usiz = SEDDIA(1)
c      END if
c
c      if ( rfracoh .LT. 0.9 .AND. rfracoh .GT. lbnd ) THEN
c        lbnd = rfracoh
c        lsiz = SEDDIA(1)
c      END if
c      if ( rfrac1 .LT. 0.9 .AND. rfrac1 .GT. lbnd ) THEN
c        lbnd = rfrac1
c        lsiz = SEDDIA(2)
c      END if
c      if ( rfrac2 .LT. 0.9 .AND. rfrac2 .GT. lbnd ) THEN
c        lbnd = rfrac2
c        lsiz = SEDDIA(3)
c      END if
c      if ( rfrac3 .LT. 0.9 .AND. rfrac3 .GT. lbnd ) THEN
c        lbnd = rfrac3
c        lsiz = SEDDIA(4)
c      END if
c
c      IF ( lbnd .LT. 0. ) THEN
c        lbnd = 0.
c        lsize = 0.
c      END IF
c      D90SIG = lsiz + (usiz-lsiz)*(0.9-lbnd)/(ubnd-lbnd)
c
c
C
C  SEDPHIC ** GP-CONVERT SEDIMENT DIAMETERS IN M TO MM AND SET PHI SIZE
C
c      SEDDIA(1)=22.0*1.e-6 !convert micron to meter
c      DO NS=1,NSED+NSND
c        rSEDPHI(NS)=-LOG(1000.*SEDDIA(NS))/LOG(2.)
c      END DO
C
C **  GP - SET MEAN PHI FOR TOP LAYER OF BED
C
c      rTVAR3W=0.
c      rTVAR3E=0.
c      rSIGPHI=0.
C
c      DO NX=1,NSND
c        NS=NX+NSED
c        rTVAR3W=rTVAR3W+rSEDPHI(NS)*SNDB(L,KTOP,NX)
c        rTVAR3E=rTVAR3E+SNDB(L,KTOP,NX)
c      END DO
c      rTVAR3W=rTVAR3W+rSEDPHI(1)*SEDB(L,KTOP,1)
c      rTVAR3E=rTVAR3E+SEDB(L,KTOP,1)
C
c      IF (rTVAR3E.LE.0.) THEN
c        rTVAR3E=1.
c      ELSE
c        rTVAR3W=rTVAR3W/rTVAR3E
c      END IF
C
c      DO NX=1,NSND
c        NS=NX+NSED
c        rSIGPHI = rSIGPHI + ((rSEDPHI(NS)-rTVAR3W)**2)
c     &                     *SNDB(L,KTOP,NX)
c      END DO
c      rSIGPHI = rSIGPHI + ((rSEDPHI(1)-rTVAR3W)**2)
c     &                     *SEDB(L,KTOP,1)
c      rSIGPHI = rSIGPHI/rTVAR3E
c     write(6,*) rSIGPHI, rTVAR3E
c     stop
C
c      rSIGPHI=SQRT(rSIGPHI)
c      rSIGPHI=2.**(rSIGPHI)
c
c
C **  SET MEAN D50
C
c      D50SIG=0.
c      rSNDBT=0.
C
c      DO NX=1,NSND
c       NS=NX+NSED
c       D50SIG=D50SIG+SNDB(L,KTOP,NX)*LOG(SEDDIA(NS))
c       !D50SIG=D50SIG+SNDB(L,KTOP,NX)*(SEDDIA(NS))
c       rSNDBT=rSNDBT+SNDB(L,KTOP,NX)
c      ENDDO
c      D50SIG=D50SIG+SEDB(L,KTOP,1)*LOG(SEDDIA(1))
c      !D50SIG=D50SIG+SEDB(L,KTOP,1)*(SEDDIA(1))
c      rSNDBT=rSNDBT+SEDB(L,KTOP,1)
C
c      D50SIG=EXP(D50SIG/rSNDBT)
c      !D50SIG=D50SIG/rSNDBT
c
c
C ** Compute the D90 from the standard deviation of grain size
C    distributionin bed
c
c
c      Z90 = 1.281551  !(Z-score for the 90th percentile)
C     D90SIG= (((D50SIG*1000.) + (Z90*rSIGPHI)))/1000.
c
C ** Cohesive concentration
c      cohcon= (SEDB(L,KTOP,1)*1e-6)/HBED(L,KTOP)
c
c      rfracoh = rfracoh*100.
c
c      CSEDTAUS=(0.36*((D90SIG/D50SIG)**0.948803))
c     &         *(0.613976*(rfracoh**(-0.109768)))
c     &         *(0.846585*(cohcon**0.217732))
c
c      CSEDTAUS = CSEDTAUS/1000.
C     Limit tau critical to 2 to 50 dynes/cm^2 based on Earl  Hayter
C     recommendation  09/12/03
c      if (CSEDTAUS.LT.(2./10000.)) then
c        CSEDTAUS = 2./10000.
c      end if
c
c     !write(6,'(f12.5)') D50SIG, D90SIG, rfracoh, COHCON, CSEDTAUS
c     write(6,'(f12.5)')CSEDTAUS
c     IF ( L .EQ. 4393 ) THEN
c       stop
c     END IF
c
c      ENDIF

c#######################################################################

      RETURN
      END
