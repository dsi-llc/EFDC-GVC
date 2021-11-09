C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
C      SUBROUTINE WWQTS(TINDAY)
      SUBROUTINE WWQTS
C
C**********************************************************************C
C
C **  SHEN'S MODIFICATION TO OUTPUT MACROALGAE
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
C WRITE TIME-SERIES OUTPUT: WQCHLX=1/WQCHLX
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      DIMENSION WQVOUT(NWQVM+6)
C
C      TINDAY=TINDAY
C
      OPEN(1,FILE='WQWCTS.OUT',STATUS='UNKNOWN',POSITION='APPEND')
C
      IF(ISDYNSTP.EQ.0)THEN
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN
        TIMTMP=TIMTMP/TCTMSR  
      ELSE
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF
C
      DO M=1,IWQTS
        DO K=1,KC
          LL=LWQTS(M)
          IZ=IWQZMAP(LL,K)
C          DO NW=1,NWQV
          DO NW=1,NWQV+1    ! Sen Bai
            WQVO(LL,K,NW) = WQVO(LL,K,NW)*0.5
          ENDDO
          NWQOUT=0
          IF(ISTRWQ(1).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,1)*WQCHLC
          ENDIF
          IF(ISTRWQ(2).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,2)*WQCHLD
          ENDIF
          IF(ISTRWQ(3).EQ.1)THEN
            NWQOUT=NWQOUT+1
            WQVOUT(NWQOUT)=WQVO(LL,K,3)*WQCHLG
          ENDIF
C          DO NW=4,NWQV
           DO NW=4,NWQV+1  ! Sen Bai 
            IF(ISTRWQ(NW).EQ.1)THEN
              NWQOUT=NWQOUT+1
              WQVOUT(NWQOUT)=WQVO(LL,K,NW)
            ENDIF
          ENDDO
C Sen Bai
	     NWQOUT=NWQOUT+1
           WQVOUT(NWQOUT)=WQVO(LL,K,IDNOTRVA)
C End of Sen Bai

C ** ADD TEMPERATURE AND SALINITY SATUARATION
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT)=TEM(LL,K)
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT)=SAL(LL,K)
C ** ADD DO SATUARATION
          NWQOUT=NWQOUT+1
          TVAL1=1./(TEM(LL,K)+273.15)
          TVAL2=TVAL1*TVAL1
          TVAL3=TVAL1*TVAL2
          TVAL4=TVAL2*TVAL2
          RLNSAT1=-139.3441+(1.575701E+5*TVAL1)-(6.642308E+7*TVAL2)
     &                   +(1.2438E+10*TVAL3)-(8.621949E+11*TVAL4)
          RLNSAT2=RLNSAT1-SAL(LL,K)*( 1.7674E-2-(1.0754E+1*TVAL1)
     &                           +(2.1407E+3*TVAL2) )
          WQVOUT(NWQOUT) = EXP(RLNSAT2)
C ** ADD FLOW INDUCED REAERATION
          NWQOUT=NWQOUT+1
          WQVOUT(NWQOUT) = 0.
          IF(IWQKA(IZ) .EQ. 2)THEN
            UMRM = 0.5*( U(LL,K) + U(LL+1   ,K) )
            VMRM = 0.5*( V(LL,K) + V(LNC(LL),K) )
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
            WQVOUT(NWQOUT) = WQKRO(IZ) * XMRM**0.5 / HP(LL)**0.5
          ENDIF
C ** ADD WIND INDUCED REAERATION
          NWQOUT=NWQOUT+1
          WINDREA = WINDST(LL)
          WQVOUT(NWQOUT)=0.728*SQRT(WINDREA)
     &                  +(0.0372*WINDREA-0.317)*WINDREA
          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP,
     &       (WQVOUT(NWOUT),NWOUT=1,NWQOUT)
        ENDDO
      ENDDO
C
   71 FORMAT(3I5,F11.5,1X, 25E11.3)
C
C     K=1
COLD      DO M=1,IWQTS
COLD         DO K=1,KC
COLD          LL=LWQTS(M)
COLD          DO NW=1,NWQV
COLD            WQVO(LL,K,NW) = WQVO(LL,K,NW)*0.5
COLD          ENDDO
COLD          CHLWQ = WQVO(LL,K,1)*WQCHLC + WQVO(LL,K,2)*WQCHLD
COLD     *      + WQVO(LL,K,3)*WQCHLG
COLD          CHLC = WQVO(LL,K,1)*WQCHLC
COLD          CHLD = WQVO(LL,K,2)*WQCHLD
COLDC          CHLCD = CHLC + CHLD
COLD          CHLG = WQVO(LL,K,3)*WQCHLG
COLD          TBWQ = WQVO(LL,K,1)+WQVO(LL,K,2)+WQVO(LL,K,3)
COLD          TOCWQ = WQVO(LL,K,4)+WQVO(LL,K,5)+WQVO(LL,K,6) + TBWQ
COLD          IF(IWQSRP.EQ.1)THEN
COLD            O2WQ = MAX(WQVO(LL,K,19), 0.0)
COLD            TAMDWQ = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQVO(LL,K,20) )
COLD            TAMPWQ = WQVO(LL,K,20) - TAMDWQ
COLD            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*TAMPWQ)
COLD            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*TAMPWQ)
COLD           ELSE IF(IWQSRP.EQ.2)THEN
COLD            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*SEDT(LL,K))
COLD            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*SEDT(LL,K))
COLD           ELSE
COLD            PO4DWQ = WQVO(LL,K,10)
COLD            SADWQ  = WQVO(LL,K,17)
COLD          ENDIF
COLD          XPO4DWQ = MAX(PO4DWQ,0.0)
COLD          APCWQ = 1.0 / (WQCP1PRM + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
COLD          TPWQ = WQVO(LL,K,7)+WQVO(LL,K,8)+WQVO(LL,K,9)+WQVO(LL,K,10)
COLD     *      + APCWQ*TBWQ
COLD          TNWQ = WQVO(LL,K,11)+WQVO(LL,K,12)+WQVO(LL,K,13)+WQVO(LL,K,14)
COLD     *      +WQVO(LL,K,15) + WQANCC*WQVO(LL,K,1)+WQANCD*WQVO(LL,K,2)
COLD     *      +WQANCG*WQVO(LL,K,3)
COLD          IF(IWQSI.EQ.1)THEN
COLD            TSIWQ = WQVO(LL,K,16)+WQVO(LL,K,17) + WQASCD*WQVO(LL,K,2)
COLD           ELSE
COLD            TSIWQ = 0.0
COLD          ENDIF
C
C M. MORTON 03/11/96 ADDED BOD5 TERM FOR OUTPUT:
COLD          IZ=IWQZMAP(LL,K)
COLD          BOD5 = 2.67 * ( WQVO(LL,K,5)*(1.0-EXP(-5.0*WQKLC))
COLD     +           + WQVO(LL,K,6)*(1.0-EXP(-5.0*WQKDC(IZ)))
COLD     +           + WQVO(LL,K,18)*(1.0-EXP(-5.0*WQKCD(IZ)))
COLD     +           + WQVO(LL,K,1)*(1.0-EXP(-5.0*WQBMRC(IZ)))
COLD     +           + WQVO(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(IZ)))
COLD     +           + WQVO(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(IZ))) )
COLD     +           + 4.57 * WQVO(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
COLD
CMRM          BOD5=0.0       ! JI, 10/2/97
C MRM MODIFIED OUTPUT FORMAT: REMOVED TAM AND TMP AND ADDED
C DIATOMS (CHLD) AND GREEN ALGAE (CHLG); REMOVED APCWQ AND
C ADDED CYANOBACTERIA (CHLC):
COLD          IF(IDNOTRVA.EQ.0)THEN
COLD          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP, CHLWQ,TOCWQ,WQVO(LL,K,6),
COLD     +      TPWQ, WQVO(LL,K,9), WQVO(LL,K,10), PO4DWQ, CHLC, TNWQ,
COLD     +      (WQVO(LL,K,NW),NW=13,15), TSIWQ, (WQVO(LL,K,NW),NW=16,17),
COLD     +      SADWQ, BOD5, WQVO(LL,K,19), CHLD, CHLG, WQVO(LL,K,NWQV)
COLD         ELSE
C HHTMP = WATER DEPTH (METERS)
C          HHTMP = GI*P(LL) - BELV(LL)
C CHLM  = MACROALGAE BIOMASS IN MICROGRAMS/SQUARE METER:
C          CHLM = WQVO(LL,K,NWQV+1) * WQCHLM * HHTMP
C CHLM IN UG/L AS FOLLOWS:
COLD          CHLM = WQVO(LL,K,NWQV+1) * WQCHLM
COLD          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP, CHLWQ,TOCWQ,WQVO(LL,K,6),
COLD     +      TPWQ, WQVO(LL,K,9), WQVO(LL,K,10), PO4DWQ, CHLC, TNWQ,
COLD     +      (WQVO(LL,K,NW),NW=13,15), TSIWQ, (WQVO(LL,K,NW),NW=16,17),
COLD     +      SADWQ, BOD5, WQVO(LL,K,19), CHLD, CHLG,WQVO(LL,K,NWQV),
COLD     +      CHLM
C MRM     +      WQVO(LL,K,NWQV+1)
COLD          ENDIF
C          WRITE(1,71) IL(LL),JL(LL),K, TIMTMP,CHLWQ,TOCWQ,WQVO(LL,K,6),
C     *      TPWQ,WQVO(LL,K,9),WQVO(LL,K,10),PO4DWQ,APCWQ,TNWQ,
C     *      (WQVO(LL,K,NW),NW=13,15),TSIWQ,(WQVO(LL,K,NW),NW=16,17),
C     *      SADWQ,(WQVO(LL,K,NW),NW=18,20),TAMPWQ,WQVO(LL,K,NWQV)
COLD         ENDDO
COLD      ENDDO
C
      CLOSE(1)
C
C   71 FORMAT(3I5,F11.5, 1P, 23E11.3)
C
      RETURN
      END
