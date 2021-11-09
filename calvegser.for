C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALVEGSER (ISTL)
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
C  NVEGSER = number of vegetation time series
C  NVEGSERV(NVEGTPM) = time series id for specific vegetation class 
C  MVEGTLAST(NVEGSERM) = place holder in interpolation table
C  TCVEGSER(NVEGSERM) = time conversion factor for time variable
C  TVEGSER(NDVEGSER,NVEGSERM) = time of data
C  VEGSERRT(NVEGSERM) = current value of RDLPSQ
C  VEGSERBT(NVEGSERM) = current value of BPVEG
C  VEGSERHT(NVEGSERM) = current value of HPVEG
C  VEGSERR(NDVEGSER,NVEGSERM) = time varying values of RDLPSQ
C  VEGSERB(NDVEGSER,NVEGSERM) = time varying values of BPVEG
C  VEGSERH(NDVEGSER,NVEGSERM) = time varying values of HPVEG
C
C**********************************************************************C
C
C ** SUBROUTINE CALVEGSER UPDATES TIME VARIABLE VEGETATION RESISTANCE 
C ** PARAMETERS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      IF(NVEGSER.GT.0)THEN
C
      DO NS=1,NVEGSER
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)/TCVEGSER(NS)+TBEGIN*(TCON/TCVEGSER(NS))
      ELSE
        TIME=TIMESEC/TCVEGSER(NS)
      ENDIF
      M1=MVEGTLAST(NS)
  100 CONTINUE
      M2=M1+1
      IF(TIME.GT.TVEGSER(M2,NS))THEN
       M1=M2
       GOTO 100
      ELSE
       MVEGTLAST(NS)=M1
      ENDIF      
C
      TDIFF=TVEGSER(M2,NS)-TVEGSER(M1,NS)
      WTM1=(TVEGSER(M2,NS)-TIME)/TDIFF
      WTM2=(TIME-TVEGSER(M1,NS))/TDIFF
      VEGSERRT(NS)=WTM1*VEGSERR(M1,NS)+WTM2*VEGSERR(M2,NS)
      VEGSERBT(NS)=WTM1*VEGSERB(M1,NS)+WTM2*VEGSERB(M2,NS)
      VEGSERHT(NS)=WTM1*VEGSERH(M1,NS)+WTM2*VEGSERH(M2,NS)
C
      ENDDO
C
      DO M=1,MVEGTYP
	  NSTMP=NVEGSERV(M)
	  IF(NSTMP.GT.0)THEN
          RDLPSQ(M)=VEGSERRT(NSTMP)
	    BPVEG(M)=VEGSERBT(NSTMP)
	    HPVEG(M)=VEGSERHT(NSTMP)
          BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
          PVEGX(M)=1.-BETVEG(M)*BDLTMP
          PVEGY(M)=1.-BETVEG(M)*BDLTMP
          PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
          BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
	  ENDIF
      ENDDO
C
      ENDIF
C
C 6000 FORMAT('N, VEGSERT = ',I6,4X,F12.4)
C
C**********************************************************************C
C
      RETURN
      END
