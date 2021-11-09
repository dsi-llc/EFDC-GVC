C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION FSBDLD(DIASED,GPDIASED,D50,DEP,PEXP,PHID,CSHIELDS,
     &                     SBDLDP,ISOPT)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C 12/03/2001        john hamrick       12/03/2001       john hamrick
C  modifed to include alternate formulas for bed load phi
C
C----------------------------------------------------------------------C
C
      INCLUDE 'EFDC.PAR'
C
C **  CALCULATES DIMNSIONLESS BED LOAD TRANSPORT COEFFICIENT
C
C     DIASED = SAND GRAIN DIAMETER 
C     GPDIASED    = SAND GRAIN SPECIFIC GRAVITY 
C     CSHIELDS     = SAND GRAIN SETTLING VELOCITY
C
C **  ISOPT=0  USE CONSTANT VALUE
C
      IF(ISOPT.EQ.0) FSBDLD=SBDLDP
C
C **  ISOPT=1  BASED ON 
C **
C **  VAN RIJN, L. C., 1984: SEDIMENT TRANSPORT, PART I: BED
C **  LOAD TRANSPORT, J. HYDRAULIC ENGINEERING, 110, 1431-1455.
C
      IF(ISOPT.EQ.1)THEN
        RD=DIASED*SQRT(GPDIASED)*1.E6
        RD=(1./RD)**0.2
        TMP=CSHIELDS**2.1
        FSBDLD=0.053*RD/TMP
      ENDIF
C
C **  ISOPT=2  BASED ON MODIFIED ENGULAND-HANSEN FORMULA
C **  REFERENCE TO BE ADDED
C
      IF(ISOPT.EQ.2)THEN
        TMP1=(DEP/D50)**0.33333
        TMP2=(PEXP/PHID)**1.125
        FSBDLD=2.0367*TMP1*TMP2
      ENDIF
C
C **  ISOPT=3  BASED ON WU, WANG AND JIA, J. HYDR. RES. V38, 3000
C
      IF(ISOPT.EQ.3)THEN
        TMP1=0.03*((PHID/PEXP)**0.6)
        TMP2=TMP1**2.2
        FSBDLD=0.0053/TMP2
      ENDIF
C

      RETURN
      END 
