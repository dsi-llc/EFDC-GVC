C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SHOWVAL3
C
C  THIS VERSION OF SHOWVAL ARISES AS SHOWVAL1 DURING FSPLIT OF 
C  HOUSATONIC FILE SET UNDER SUBROUTINE ZZZ000.FOR
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0A
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
      CHARACTER BLANK,ASTER,CSURF(32),CSALS(20),CSALB(20)
C
      REAL*8 SBDLDT(5000)
      DATA BLANK/' '/
      DATA ASTER/'*'/
      DATA NSHOWTITLE/0/      
C
      IF(NSHOWC.EQ.NSHOWR)THEN
        OPEN(1,FILE='SHOW.INP',STATUS='UNKNOWN')
        DO NSKIP=1,6
          READ(1,1)
        ENDDO
        READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW,ISHPRT
        READ(1,*)ZSSMIN,ZSSMAX,SSALMAX
        CLOSE(1)
      ENDIF

      NSHOWTITLE=NSHOWTITLE+1
C      IF(NSHOWC.EQ.NSHOWR)THEN
      IF(MOD(NSHOWTITLE,NSHOWR).EQ.1)THEN
        IF(NSHTYPE.EQ.1)THEN
          ! *** SALINITY SPECIAL
          IZSMIN=NINT(ZSSMIN)
          IZSMAX=NINT(ZSSMAX)
          ISALMAX=NINT(SSALMAX)
          WRITE(6,2)
          WRITE(6,8)
          WRITE(6,9)
          WRITE(6,10)IZSMIN,IZSMAX,ISALMAX,ISALMAX
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.2)THEN
          ! *** N & SALINITY
          WRITE(6,2)
          WRITE(6,3)
          WRITE(6,4)
          WRITE(6,6)
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.3)THEN
          ! *** TIME & SALINITY
          WRITE(6,2)
          WRITE(6,3)
          WRITE(6,44)
          WRITE(6,66)
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.4)THEN
          ! *** TIME & TEMPERATURE
          WRITE(6,2)
          WRITE(6,34)
          WRITE(6,44)
          WRITE(6,68)
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.5)THEN
          ! *** TIME, SALINITY & SEDIMENT
          WRITE(6,2)
          WRITE(6,33)
          WRITE(6,44)
          WRITE(6,67)
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.6)THEN
          ! *** TIME, COHESIVE AND NONCOHESIVE SEDIMENT
          WRITE(6,2)
          WRITE(6,39)
          WRITE(6,49)
          WRITE(6,69)
          WRITE(6,2)
        ELSEIF(NSHTYPE.EQ.8)THEN
          ! *** TIME, COHESIVE AND NONCOHESIVE SEDIMENT FOR HOUSATONIC RIVER
          WRITE(6,2)
          WRITE(6,590)
          WRITE(6,591)
          WRITE(6,592)
          WRITE(6,2)
        ENDIF
      ENDIF
C
      IF(NSHTYPE.EQ.1)THEN
        DO M=1,32
          CSURF(M)=BLANK
        ENDDO
        DO M=1,20
          CSALS(M)=BLANK
          CSALB(M)=BLANK
        ENDDO
        L=LIJ(ICSHOW,JCSHOW)
        IF(IS1DCHAN.GT.0)THEN
          ZSURF=HP(L)
        ELSE
          ZSURF=HP(L)+BELV(L)
        ENDIF
        ZSTMP=(31.*(ZSURF-ZSSMIN)/(ZSSMAX-ZSSMIN))+1.
        IZSTMP=NINT(ZSTMP)
        IF(IZSTMP.GT.32)IZSTMP=32
        IF(IZSTMP.LT.1)IZSTMP=1
        CSURF(IZSTMP)=ASTER
        SSTMP=(19.*SAL(L,KC)/SSALMAX)+1.
        SBTMP=(19.*SAL(L,1)/SSALMAX)+1.
        ISSTMP=NINT(SSTMP)
        ISBTMP=NINT(SBTMP)
        IF(ISSTMP.GT.20)ISSTMP=20
        IF(ISSTMP.LT.1)ISSTMP=1
        IF(ISBTMP.GT.20)ISBTMP=20
        IF(ISBTMP.LT.1)ISBTMP=1
        CSALS(ISSTMP)=ASTER
        CSALB(ISBTMP)=ASTER
        WRITE(6,11)N,ICSHOW,JCSHOW,(CSURF(I),I=1,32),(CSALS(J),J=1,20),
     &      (CSALB(K),K=1,20)
      ELSE
       IF(MOD(NSHOWTITLE,ISHPRT).EQ.0)THEN   
        IF(ISDYNSTP.EQ.0)THEN
          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIME=TIMESEC/86400.
        ENDIF
        L=LIJ(ICSHOW,JCSHOW)
        LN=LNC(L)
        IF(IS1DCHAN.GT.0)THEN
          ZSURF=HP(L)
        ELSE
          ZSURF=HP(L)+BELV(L)
        ENDIF
        UTMP=0.5*STCUV(L)*(U(L+1,KC)+U(L,KC))*100.
        VTMP=0.5*STCUV(L)*(V(LN,KC)+V(L,KC))*100.
        VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
        UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))*100.
        VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))*100.
        VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
        VELLL=0.01*SQRT(VELEKB**2+VELNKB**2)
        FROUDE=VELLL/SQRT(9.81*HP(L))
        AVKB=AV(L,1)*10000.*HP(L)
        ABKB=AB(L,1)*10000.*HP(L)
        IF(KC.GT.1)THEN
          AVKS=AV(L,KS)*10000.*HP(L)
          ABKS=AB(L,KS)*10000.*HP(L)
        ELSE
          AVKS=AVKB
          ABKS=ABKB
        ENDIF
        SALKC=SAL(L,KC)
        SALKB=SAL(L,1)
        ! *** COMPUTE THE TSS
        SEDT(L,KC)=0.
        SEDT(L,1)=0.
        DO NS=1,NSED
          SEDT(L,KC)=SEDT(L,KC)+SED(L,KC,NS)
          SEDT(L,1)=SEDT(L,1)+SED(L,1,NS)
        ENDDO
        SNDT(L,KC)=0.
        SNDT(L,1)=0.
        IF(NSHTYPE.LT.8)THEN
          DO NS=1,NSND
            SNDT(L,KC)=SNDT(L,KC)+SND(L,KC,NS)
            SNDT(L,1)=SNDT(L,1)+SND(L,1,NS)
          ENDDO
        ENDIF
        ! *** COMPUTE THE TOTAL BEDLOAD
        SBDLDT(L)=0.
        IF(NSHTYPE.EQ.6)THEN
          DO NS=1,NSND
            SBDLDT(L)=SBDLDT(L)+QSBDLDX(L,NS)+QSBDLDY(L,NS)
          ENDDO
        ENDIF
        IF(NSHTYPE.EQ.5.OR.NSHTYPE.EQ.6)THEN
          SEDKC=SEDT(L,KC)
          SEDKB=SEDT(L,1)
        ENDIF
        IF(NSHTYPE.EQ.6)THEN
          SNDKC=SNDT(L,KC)
          SNDKB=SNDT(L,1)
          SBDLDTOT=SBDLDT(L)
        ENDIF
        IF(NSHTYPE.EQ.8)THEN
          SEDK1=SED(L,KC,1)
          SNDK1=SND(L,KC,1)
          SNDK2=SND(L,KC,2)
          SNDK3=SND(L,KC,3)
          SNDK4=SND(L,KC,4)
c         SNDK5=SND(L,KC,5)
          SBDLD1=QSBDLDX(L,1)+QSBDLDY(L,1)
          SBDLD2=QSBDLDX(L,2)+QSBDLDY(L,2)
          SBDLD3=QSBDLDX(L,3)+QSBDLDY(L,3)
          SBDLD4=QSBDLDX(L,4)+QSBDLDY(L,4)
c         SBDLD5=QSBDLDX(L,5)+QSBDLDY(L,5)
        ENDIF
        IVELEKC=NINT(VELEKC)
        IVELNKC=NINT(VELNKC)
        ISALKC=NINT(SALKC)
        ISEDKC=NINT(SEDKC)
        ISNDKC=NINT(SNDKC)
        ISBEDLD=NINT(SBDLDTOT)
        ITEMKC=NINT(TEM(L,KC))
        IAVKS=NINT(AVKS)
        IABKS=NINT(ABKS)
        IVELEKB=NINT(VELEKB)
        IVELNKB=NINT(VELNKB)
        ISALKB=NINT(SALKB)
        ISEDKB=NINT(SEDKB)
        ISNDKB=NINT(SNDKB)
        ITEMKB=NINT(TEM(L,1))
        IAVKB=NINT(AVKB)
        IABKB=NINT(ABKB)
        ISBDLD1=NINT(SBDLD1)
        ISBDLD2=NINT(SBDLD2)
        ISBDLD3=NINT(SBDLD3)
        ISBDLD4=NINT(SBDLD4)
        ISBDLD5=NINT(SBDLD5)
        ISEDK1=NINT(SEDK1)
        ISNDK1=NINT(SNDK1)
        ISNDK2=NINT(SNDK2)
        ISNDK3=NINT(SNDK3)
        ISNDK4=NINT(SNDK4)
        ISNDK5=NINT(SNDK5)
        TAUBB=TAUB(L)*1000.0
        TAUSE=TAUBSED(L)*1000.0
        TAUSN=TAUBSND(L)*1000.0
        IF(NSHTYPE.EQ.2)THEN
          WRITE(6,7)N,ZSURF,
     &        IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     &        IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
        ELSEIF(NSHTYPE.EQ.3)THEN
          WRITE(6,77)TIME,ZSURF,
     &        IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     &        IVELEKB,IVELNKB,ISALKB,IACTALL,NCORDRY
        ELSEIF(NSHTYPE.EQ.4)THEN
          WRITE(6,77)TIME,ZSURF,
     &        IVELEKC,IVELNKC,TEM(L,KC),IAVKS,IABKS,
     &        IVELEKB,IVELNKB,TEM(L,1),IAVKB,IABKB
        ELSEIF(NSHTYPE.EQ.5)THEN
          WRITE(6,77)TIME,ZSURF,
     &        IVELEKC,IVELNKC,ISALKC,ISEDKC,
     &        IVELEKB,IVELNKB,ISALKB,ISEDKB,NCORDRY
        ELSEIF(NSHTYPE.EQ.6)THEN
          WRITE(6,679)TIME,ZSURF,IVELEKC,IVELNKC,ISEDKC,ISNDKC,
     &               IVELEKB,IVELNKB,ISEDKB,ISNDKC,ISBEDLD,NCORDRY
        ELSEIF(NSHTYPE.EQ.8)THEN
          WRITE(6,593)TIME,ZSURF,IVELEKC,IVELNKC,ISEDK1,ISNDK1,ISNDK2,
     &ISNDK3,ISNDK4,ISNDK5,ISBDLD1,ISBDLD2,ISBDLD3,ISBDLD4,ISBDLD5,
     .TAUSE,TAUBB
        ENDIF
       ENDIF
      ENDIF
C
    1 FORMAT (80X)
    2 FORMAT('----------------------------------------------------'
     & ,'---------------------------')
    3 FORMAT(' TIME     ZSUR   VELE  VELN  SAL   AV  AB VELE VELN  SAL')
   33 FORMAT('    TIME    ZSUR   VELE  VELN   SAL   SED  VELE  VELN   SA
     &L   SED  # Dry')
   34 FORMAT('   TIME   ZSUR  VELE VELN   TEM   AV    AB  VELE VELN   TE
     .M   AV    AB')
   44 FORMAT('                SURF SURF  SURF  SURF  SURF BOTT BOTT  BOT
     .T  BOTT  BOTT')
   68 FORMAT('   DAYS     m   cm/s cm/s   D:C cm2/s cm2/s cm/s cm/s   D:
     .C cm2/s cm2/s')
   77 FORMAT(F7.3,F7.2,I6,I5,F6.1,2I6,2I5,F6.1,2I6)
   67 FORMAT('    days     cm    cm/s  cm/s   psu  mg/L  cm/s  cm/s   ps
     &u  mg/L')
C
C - COHESIVE AND NONCOHESIVE SEDIMENT OUTPUT - NSHTYPE=6
C
   39 FORMAT('    TIME    ZSUR VELE VELN COHSED NONSED VELE VELN COHSED  
     &NONSED  BDLD # Dry')
   49 FORMAT('                 SURF SURF  SURF   SURF  BOTT BOTT  BOTT   
     & BOTT        cells')
   69 FORMAT('    days     cm  cm/s cm/s  mg/L   mg/L  cm/s cm/s  mg/L   
     & mg/L')
  679 FORMAT(F9.4,F7.1,2I5,2I7,2I5,2I7,2I6)
C
C - COHESIVE AND NONCOHESIVE SEDIMENT OUTPUT for the HOUSATONIC & GREEN RIVERS - NSHTYPE=8
C
  590 FORMAT('    TIME    ZSUR  VELE  VELN COHSED  SED1  SND1  SND2  SND
     &3  SND4 BEDLD1 BEDLD2 BEDLD3 BEDLD4 BEDLD5  TAUBSED       TAUB')
  591 FORMAT('                                                          
     &           ')
  592 FORMAT('    days      m   cm/s  cm/s  mg/L   mg/L  mg/L  mg/L  mg/
     &L  mg/L   gm/s   gm/s   gm/s   gm/s   gm/s     Pa          Pa')
  593 FORMAT(F9.3,F7.1,I6,I6,I7,5I6,5I7,2E12.5,2F8.3)
             !1234567891234567 |234||234||234||234||234||234||234||234||234||234|
   35 FORMAT(' TIME     ZSUR   VELE VELN  TOX   AV   AB VELE VELN  TOX')
    4 FORMAT(' STEP            SURF SURF SURF BOTT BOTT BOTT BOTT BOTT',
     &       ' #ACT #DRY')
             !1234567891234567 |234||234||234||234||234||234||234||234||234||234|
    6 FORMAT('                  CM   CM/S  CM/S ',
     &     ' CM/S  CM/S  CM/S  PSU    PIPE  CELLS')
   66 FORMAT(' DAYS             CM   CM/S  CM/S ',
     &     ' CM/S  CM/S  CM/S  PSU    PIPE  CELLS')
             !1234567891234567 |234||234||234||234||234||234||234||234||234||234|
    7 FORMAT(I6,F7.2,11I5)

    8 FORMAT(' TIME    I   J       FREE SURFACE ELEVATION     ',1X
     &    ,' SURFACE SALINITY     BOTTOM SALINITY')
    9 FORMAT(' STEP                         (CM)              ',1X
     &    ,'       (PSU)              (PSU)')
   10 FORMAT('                ',I3,'                          ',I3
     &    ,' 0               ',I3,' 0               ',I3)
   11 FORMAT(I5,1X,I4,1X,I4,1X,32A1,20A1,20A1)
C
      NSHOWC=NSHOWC+1
      RETURN
      END
