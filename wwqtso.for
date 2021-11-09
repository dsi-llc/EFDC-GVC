C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQTSO
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
C WRITE TIME-SERIES OUTPUT: WQCHLX=1/WQCHLX
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      TINDAY=TINDAY
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
C     K=1
      DO M=1,IWQTS
         DO K=1,KC
          LL=LWQTS(M)
          DO NW=1,NWQV
            WQVO(LL,K,NW) = WQVO(LL,K,NW)*0.5
          ENDDO
          CHLWQ = WQVO(LL,K,1)*WQCHLC + WQVO(LL,K,2)*WQCHLD
     *      + WQVO(LL,K,3)*WQCHLG
          CHLC = WQVO(LL,K,1)*WQCHLC
          CHLD = WQVO(LL,K,2)*WQCHLD
C          CHLCD = CHLC + CHLD
          CHLG = WQVO(LL,K,3)*WQCHLG
          TBWQ = WQVO(LL,K,1)+WQVO(LL,K,2)+WQVO(LL,K,3)
          TOCWQ = WQVO(LL,K,4)+WQVO(LL,K,5)+WQVO(LL,K,6) + TBWQ
          IF(IWQSRP.EQ.1)THEN
            O2WQ = MAX(WQVO(LL,K,19), 0.0)
            TAMDWQ = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQVO(LL,K,20) )
            TAMPWQ = WQVO(LL,K,20) - TAMDWQ
            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*TAMPWQ)
            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*TAMPWQ)
           ELSE IF(IWQSRP.EQ.2)THEN
            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*SEDT(LL,K))
            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*SEDT(LL,K))
           ELSE
            PO4DWQ = WQVO(LL,K,10)
            SADWQ  = WQVO(LL,K,17)
          ENDIF
          XPO4DWQ = MAX(PO4DWQ,0.0)
          APCWQ = 1.0 / (WQCP1PRM + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
          TPWQ = WQVO(LL,K,7)+WQVO(LL,K,8)+WQVO(LL,K,9)+WQVO(LL,K,10)
     *      + APCWQ*TBWQ
          TNWQ = WQVO(LL,K,11)+WQVO(LL,K,12)+WQVO(LL,K,13)+WQVO(LL,K,14)
     *      +WQVO(LL,K,15) + WQANCC*WQVO(LL,K,1)+WQANCD*WQVO(LL,K,2)
     *      +WQANCG*WQVO(LL,K,3)
          IF(IWQSI.EQ.1)THEN
            TSIWQ = WQVO(LL,K,16)+WQVO(LL,K,17) + WQASCD*WQVO(LL,K,2)
           ELSE
            TSIWQ = 0.0
          ENDIF
C
C M. MORTON 03/11/96 ADDED BOD5 TERM FOR OUTPUT:
C
          IZ=IWQZMAP(LL,K)
          BOD5 = 2.67 * ( WQVO(LL,K,5)*(1.0-EXP(-5.0*WQKLC))
     +           + WQVO(LL,K,6)*(1.0-EXP(-5.0*WQKDC(IZ)))
     +           + WQVO(LL,K,18)*(1.0-EXP(-5.0*WQKCD(IZ)))
     +           + WQVO(LL,K,1)*(1.0-EXP(-5.0*WQBMRC(IZ)))
     +           + WQVO(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(IZ)))
     +           + WQVO(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(IZ))) )
     +           + 4.57 * WQVO(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
C
C MRM          BOD5=0.0       ! JI, 10/2/97
C MRM MODIFIED OUTPUT FORMAT: REMOVED TAM AND TMP AND ADDED
C DIATOMS (CHLD) AND GREEN ALGAE (CHLG); REMOVED APCWQ AND
C ADDED CYANOBACTERIA (CHLC):
          WRITE(1,71) IL(LL),JL(LL),K,TIMTMP, CHLWQ,TOCWQ,WQVO(LL,K,6),
     +      TPWQ, WQVO(LL,K,9), WQVO(LL,K,10), PO4DWQ, CHLC, TNWQ,
     +      (WQVO(LL,K,NW),NW=13,15), TSIWQ, (WQVO(LL,K,NW),NW=16,17),
     +      SADWQ, BOD5, WQVO(LL,K,19), CHLD, CHLG, WQVO(LL,K,NWQV)
C
C          WRITE(1,71) IL(LL),JL(LL),K, TIMTMP,CHLWQ,TOCWQ,WQVO(LL,K,6),
C     *      TPWQ,WQVO(LL,K,9),WQVO(LL,K,10),PO4DWQ,APCWQ,TNWQ,
C     *      (WQVO(LL,K,NW),NW=13,15),TSIWQ,(WQVO(LL,K,NW),NW=16,17),
C     *      SADWQ,(WQVO(LL,K,NW),NW=18,20),TAMPWQ,WQVO(LL,K,NWQV)
         ENDDO
      ENDDO
C
      CLOSE(1)
C
   71 FORMAT(3I5,F11.5, 1P, 21E11.3)
C
      RETURN
      END
