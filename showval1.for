C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SHOWVAL1
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
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER BLANK,ASTER,CSURF(32),CSALS(20),CSALB(20)
C
      DATA BLANK/' '/
      DATA ASTER/'*'/
C
C**********************************************************************C
C
      IF(NSHOWC.EQ.NSHOWR)THEN
        OPEN(1,FILE='SHOW.INP',STATUS='UNKNOWN')
         DO NSKIP=1,6
         READ(1,1)
         ENDDO
        READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW,ISHPRT
        READ(1,*)ZSSMIN,ZSSMAX,SSALMAX
        CLOSE(1)
        IF(NSHTYPE.EQ.1)THEN
          IZSMIN=NINT(ZSSMIN)
          IZSMAX=NINT(ZSSMAX)
          ISALMAX=NINT(SSALMAX)
          NSHOWC=0
          WRITE(6,2)   
          WRITE(6,8)   
          WRITE(6,9)   
          WRITE(6,10)IZSMIN,IZSMAX,ISALMAX,ISALMAX   
          WRITE(6,2)         
         ELSE
          IF(NSHTYPE.EQ.2)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,3)
            WRITE(6,4)
            WRITE(6,5)
            WRITE(6,6)
            WRITE(6,5)
            WRITE(6,2)
           ELSE
            IF(NSHTYPE.EQ.3)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,3)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,66)
            WRITE(6,5)
            WRITE(6,2)
            ENDIF
            IF(NSHTYPE.GE.5)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,33)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,67)
            WRITE(6,5)
            WRITE(6,2)
            ENDIF
            IF(NSHTYPE.EQ.4)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,34)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,68)
            WRITE(6,5)
            WRITE(6,2)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C     
      IMODTMP=MOD(N,ISHPRT)
      IF(IMODTMP.NE.0)THEN
        NSHOWC=NSHOWC+1
        RETURN
      ENDIF
C
      IF(NSHTYPE.EQ.1)THEN
        NSHOWC=NSHOWC+1
        DO M=1,32
        CSURF(M)=BLANK
        ENDDO
        DO M=1,20
        CSALS(M)=BLANK
        CSALB(M)=BLANK
        ENDDO
        L=LIJ(ICSHOW,JCSHOW)
	  KBP=KGVCP(L)
        ZSURF=(HP(L)+BELV(L))*100.
        IF(IS1DCHAN.GT.0) ZSURF=HP(L)*100.
        ZSTMP=(31.*(ZSURF-ZSSMIN)/(ZSSMAX-ZSSMIN))+1.
        IZSTMP=NINT(ZSTMP)
        IF(IZSTMP.GT.32)IZSTMP=32
        IF(IZSTMP.LT.1)IZSTMP=1
        CSURF(IZSTMP)=ASTER
        SSTMP=(19.*SAL(L,KC)/SSALMAX)+1.
        SBTMP=(19.*SAL(L,KBP)/SSALMAX)+1.
        ISSTMP=NINT(SSTMP)
        ISBTMP=NINT(SBTMP)
        IF(ISSTMP.GT.20)ISSTMP=20
        IF(ISSTMP.LT.1)ISSTMP=1
        IF(ISBTMP.GT.20)ISBTMP=20
        IF(ISBTMP.LT.1)ISBTMP=1
        CSALS(ISSTMP)=ASTER
        CSALB(ISBTMP)=ASTER
        WRITE(6,11)N,ICSHOW,JCSHOW,(CSURF(I),I=1,32),(CSALS(J),J=1,20),
     &                                               (CSALB(K),K=1,20)
       ELSE   
        NSHOWC=NSHOWC+1
        IF(ISDYNSTP.EQ.0)THEN
          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIME=TIMESEC/86400.
        ENDIF
	  IF(TIME.GE.9999.5) TIME=TIME-10000.0
	  IF(TIME.LT.0.0) TIME=ABS(TIME)
        L=LIJ(ICSHOW,JCSHOW)
        KBP=KGVCP(L)
        LN=LNC(L)
        ZSURF=(HP(L)+BELV(L))*100.
        IF(IS1DCHAN.GT.0) ZSURF=HP(L)*100.
        UTMP=0.5*STCUV(L)*(U(L+1,KC)+U(L,KC))*100.
        VTMP=0.5*STCUV(L)*(V(LN,KC)+V(L,KC))*100.
        VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
        UTMP=0.5*STCUV(L)*(U(L+1,KBP)+U(L,KBP))*100.
        VTMP=0.5*STCUV(L)*(V(LN,KBP)+V(L,KBP))*100.
        VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
        AVKS=AV(L,KS)*10000.*HP(L)
        AVKB=AV(L,KBP)*10000.*HP(L)
        ABKS=AB(L,KS)*10000.*HP(L)
        ABKB=AB(L,KBP)*10000.*HP(L)
        SALKC=SAL(L,KC)
        SALKB=SAL(L,KBP)
        IF(NSHTYPE.EQ.6)THEN
          SEDKC=SEDT(L,KC)/1000.
          SEDKB=SEDT(L,KBP)/1000.
          SNDKC=SNDT(L,KC)/1000.
          SNDKB=SNDT(L,KBP)/1000.
         ELSE
          SEDKC=SEDT(L,KC)
          SEDKB=SEDT(L,KBP)
          SNDKC=SNDT(L,KC)
          SNDKB=SNDT(L,KBP)
        ENDIF
        IZSURF=NINT(ZSURF)
        IVELEKC=NINT(VELEKC)
        IVELNKC=NINT(VELNKC)
        ISALKC=NINT(SALKC)
        ISEDKC=NINT(SEDKC)
        ISNDKC=NINT(SNDKC)
        ITEMKC=NINT(TEM(L,KC))
        IAVKS=NINT(AVKS)
        IABKS=NINT(ABKS)
        IVELEKB=NINT(VELEKB)
        IVELNKB=NINT(VELNKB)
        ISALKB=NINT(SALKB)
        ISEDKB=NINT(SEDKB)
        ISNDKB=NINT(SNDKB)
        ITEMKB=NINT(TEM(L,KBP))
        IAVKB=NINT(AVKB)
        IABKB=NINT(ABKB)
        IF(NSHTYPE.EQ.2)THEN
          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
     &                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     &                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
         ELSE
          IF(NSHTYPE.EQ.3)THEN
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     &                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     &                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
          ENDIF
          IF(NSHTYPE.GE.5)THEN
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
C     &                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,
C     &                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
C     &                            IVELEKC,IVELNKC,ISNDKC,IAVKS,IABKS,
C     &                            IVELEKB,IVELNKB,ISNDKB,IAVKB,IABKB
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     &                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,
     &                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB
          WRITE(6,79)ISNDKC,ISNDKB
          ENDIF
          IF(NSHTYPE.EQ.4)THEN
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     &                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,
     &                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C
    1 FORMAT (80X)
    2 FORMAT('----------------------------------------------------',
     &       '-------------------------------------------')
    3 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SAL |  AV ',
     &       1X,'|  AB  | VELE | VELN | SAL |  AV  |  AB  |')
   33 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SED |  AV ',
     &       1X,'|  AB  | VELE | VELN | SED |  AV  |  AB  |')
   34 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | TEM |  AV ',
     &       1X,'|  AB  | VELE | VELN | TEM |  AV  |  AB  |')
    4 FORMAT('| STEP |     |     |      | SURF | SURF | SUR | SURF',
     &       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
    5 FORMAT('|      |     |     |      |      |      |     |     ',
     &       1X,'|      |      |      |     |      |      |')
    6 FORMAT('|      |     |     |  CM  | CM/S | CM/S | PSU',
     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | PSU |CMSQ/S|CMSQ/S|')
    7 FORMAT('|',I6   ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
     &       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
    9 FORMAT('| STEP |     |     |              (CM)              |',
     &       1X,'       (PSU)       |       (PSU)        |')
   10 FORMAT('|      |     |     |',I3,'                          ',
     &        I3,'| 0               ',I3,'| 0               ',I3,'|')
   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
     &       20A1,'|')
   44 FORMAT('|      |     |     |      | SURF | SURF | SUR | SURF',
     &       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
   66 FORMAT('| DAYS |     |     |  CM  | CM/S | CM/S | PSU',
     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | PSU |CMSQ/S|CMSQ/S|')
   67 FORMAT('|      |     |     |  CM  | CM/S | CM/S |MG/L',
     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S |MG/L |CMSQ/S|CMSQ/S|')
   68 FORMAT('| DAYS |     |     |  CM  | CM/S | CM/S | D:C',
     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | D:C |CMSQ/S|CMSQ/S|')
   77 FORMAT('|',F6.1 ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
   79 FORMAT('|',6X ,'|',4X,1X,'|',4X,1X,'|',5X,1X,'|',5X,1X,
     &       '|',5X,1X,'|',I4,1X,'|',5X,1X,'|',5X,1X,'|',5X,1X,
     &       '|',5X,1X,'|',I4,1X,'|',5X,1X,'|',5X,1X,'|')
C
C**********************************************************************C
C
      RETURN
      END
