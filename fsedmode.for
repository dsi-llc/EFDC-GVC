C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION FSEDMODE(WS,USTOT,USGRN,RSNDM,ISNDM1,ISNDM2,IMODE)
C
C **  FSEDMODE SETS BEDLOAD (IMODE=1) AND SUSPENDED LOAD (IMODE=2)
C **  TRANSPORT FRACTIONS
C
      INCLUDE 'EFDC.PAR'
C
C     CHOOSE BETWEEEN TOTAL STRESS SHEAR VELOCITY AND GRAIN STRESS 
C     SHEAR VELOCITY
C
      IF(WS.EQ.0)THEN  ! DSI
         FSEDMODE=1.0
         RETURN
      ENDIF
      IF(ISNDM2.EQ.0)THEN
        US=USTOT
      ELSE
        US=USGRN
      ENDIF
C
      USDWS=US/WS
C
C     CHOOSE BETWEEN MODE OPTIONS
C
C     ISNDM1=0 SET BOTH BEDLOAD AND SUSPENDED LOAD FRACTIONS TO 1.0
C     this option allows maximum bed load and suspened load transport
C     independent of each other
C     
      IF(ISNDM1.EQ.0) FSEDMODE=1.0
C
C     ISNDM1=1 SET BEDLOAD FRACTION TO 1. USE BINARY RELATIONSHIP FOR SUSPENDED
C     this option allows maximum bedload transport. suspended load is zero
C     for usdws.lt.rsndm or maximum for usdws.ge.rsndm
C     ie:  fsedmodeBL = 1 and fsedmodeSL = 0 or 1
C
      IF(ISNDM1.EQ.1)THEN
        FSEDMODE=0.
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.0
        ELSE
          IF(USDWS.GE.RSNDM)FSEDMODE=1.
        END IF
      END IF
C
C     ISNDM1=2 SET BEDLOAD FRACTION TO 1, USE LINEAR RELATIONSHIP FOR SUSPENDED
C     this option allows maximum bedload transport. suspended load ranges from zero
C     to maximum in a linear manner
C     ie:  fsedmodeBL = 1 and 0 < fsedmodeSL < 1
C
      IF(ISNDM1.EQ.2)THEN
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.0
        ELSE
          TMPVAL=((USDWS)-0.4)/9.6
          TMPVAL=MIN(TMPVAL,1.0)
          TMPVAL=MAX(TMPVAL,0.0)
          FSEDMODE=TMPVAL
        END IF
      END IF
C
C     ISNDM1=3 USE BINARY RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
C     this option does not allow simulataneous transport of bedload 
C     and suspended load.  for usdws.lt.rsndm maximum bedload transport
C     occurs.  otherwise maximum suspended load occurs
c     ie.  either: fsedmodeBL =1 fsedmodeSL = 0 
c              or: fsedmodeBL =0 fsedmodeSL = 1 
C
      IF(ISNDM1.EQ.3)THEN
        FSEDMODE=0.
        IF(IMODE.EQ.1)THEN
          IF(USDWS.LT.RSNDM)FSEDMODE=1.
        ELSE
          IF(USDWS.GE.RSNDM)FSEDMODE=1.
        END IF
      END IF
C
C     ISNDM1=4 USE LINEAR RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
C     this options partitions total load between suspended load and 
C     bedload using the linear suspended load relationships
c     ie fsedmodeBL + fsedmodeSL = 1
C
      IF(ISNDM1.EQ.4)THEN
        TMPVAL=((USDWS)-0.4)/9.6
        TMPVAL=MIN(TMPVAL,1.0)
        TMPVAL=MAX(TMPVAL,0.0)
        IF(IMODE.EQ.1)THEN
          FSEDMODE=1.-TMPVAL
        ELSE
          FSEDMODE=TMPVAL
        END IF
      END IF
C
      RETURN
      END
