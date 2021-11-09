C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SHOWVAL2
C
C **  LAST MODIFIED BY MIKE MORTON ON 8 AUGUST 2001
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
        READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW, ISHPRT                      
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
C            WRITE(6,5)                                                      
            WRITE(6,6)
C            WRITE(6,5)                                                      
            WRITE(6,2)
           ELSE
            IF(NSHTYPE.EQ.3)THEN
            NSHOWC=0        
            WRITE(6,2)
C            WRITE(6,3)                                                      
            WRITE(6,33)                                                      
            WRITE(6,44)
C            WRITE(6,5)                                                      
            WRITE(6,66)
C            WRITE(6,5)                                                      
            WRITE(6,2)
            ENDIF
            IF(NSHTYPE.EQ.4)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,33)
            WRITE(6,44)
C            WRITE(6,5)                                                      
            WRITE(6,67)
C            WRITE(6,5)                                                      
            WRITE(6,2)
            ENDIF
            IF(NSHTYPE.EQ.5)THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,34)
            WRITE(6,44)
C            WRITE(6,5)                                                      
            WRITE(6,68)
C            WRITE(6,5)                                                      
            WRITE(6,2)
            ENDIF
          ENDIF
        ENDIF
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
        ZSURF=(HP(L)+BELV(L))*100.
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
     &                                               (CSALB(K),K=1,20)
      ELSE
       IF(MOD(N,ISHPRT) .EQ. 0)THEN                                      
        NSHOWC=NSHOWC+1
        IF(ISDYNSTP.EQ.0)THEN
          TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        ELSE
          TIME=TIMESEC/86400.
        ENDIF
        L=LIJ(ICSHOW,JCSHOW)
        LN=LNC(L)
        ZSURF=(HP(L)+BELV(L))*100.
        UTMP=0.5*STCUV(L)*(U(L+1,KC)+U(L,KC))*100.
        VTMP=0.5*STCUV(L)*(V(LN,KC)+V(L,KC))*100.
        VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
        UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))*100.
        VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))*100.
        VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
        AVKS=AV(L,KS)*10000.*HP(L)
        AVKB=AV(L,1)*10000.*HP(L)
        ABKS=AB(L,KS)*10000.*HP(L)
        ABKB=AB(L,1)*10000.*HP(L)
        SALKC=SAL(L,KC)
        SALKB=SAL(L,1)
C       SEDKC=SEDT(L,KC)/1000.
C       SEDKB=SEDT(L,1)/1000.
C       SNDKC=SNDT(L,KC)/1000.
C       SNDKB=SNDT(L,1)/1000.
        SEDKC=SEDT(L,KC)
        SEDKB=SEDT(L,1)
        SNDKC=SNDT(L,KC)
        SNDKB=SNDT(L,1)
        TEMKC=TEM(L,KC)                                                      
        TEMKB=TEM(L,1)                                                       
C        IZSURF=NINT(ZSURF)
C        IVELEKC=NINT(VELEKC)
C        IVELNKC=NINT(VELNKC)
C        ISALKC=NINT(SALKC)
C        ISEDKC=NINT(SEDKC)
C        ISNDKC=NINT(SNDKC)
C        ITEMKC=NINT(TEM(L,KC))
C        IAVKS=NINT(AVKS)
C        IABKS=NINT(ABKS)
C        IVELEKB=NINT(VELEKB)
C        IVELNKB=NINT(VELNKB)
C        ISALKB=NINT(SALKB)
C        ISEDKB=NINT(SEDKB)
C        ISNDKB=NINT(SNDKB)
C        ITEMKB=NINT(TEM(L,1))
C        IAVKB=NINT(AVKB)
C        IABKB=NINT(ABKB)
        IF(NSHTYPE.EQ.2)THEN
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,                                 
C     &                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,       
C     &                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB        
C M. MORTON WRITES REAL VARIABLES TO PC SCREEN:                              
          WRITE(6,7)N,ICSHOW,JCSHOW, ZSURF,                                  
     +                             VELEKC, VELNKC, SALKC, AVKS, ABKS,        
     +                             VELEKB, VELNKB, SALKB, AVKB, ABKB         
         ELSE
          IF(NSHTYPE.EQ.3)THEN
C          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,                             
C     &                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,       
C     &                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB        
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,                        
     &                             VELEKC, VELNKC, SALKC,                    
     &                             VELEKB, VELNKB, SALKB                     
          ENDIF
          IF(NSHTYPE.EQ.4)THEN
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,                                 
C     &                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,       
C     &                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB        
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,                                 
C     &                            IVELEKC,IVELNKC,ISNDKC,IAVKS,IABKS,       
C     &                            IVELEKB,IVELNKB,ISNDKB,IAVKB,IABKB        
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,                        
     &                             VELEKC, VELNKC, SEDKC,                    
     &                             VELEKB, VELNKB, SEDKB                     
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,                        
     &                             VELEKC, VELNKC, SNDKC,                    
     &                             VELEKB, VELNKB, SNDKB                     
          ENDIF
          IF(NSHTYPE.EQ.5)THEN
C          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,                             
C     &                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,       
C     &                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB        
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,                        
     &                             VELEKC, VELNKC, TEMKC,                    
     &                             VELEKB, VELNKB, TEMKB                     
          ENDIF
        ENDIF
       ENDIF                                                                
      ENDIF
C
C**********************************************************************C
C
    1 FORMAT (80X)
C NOTE: M. MORTON CHANGED FORMAT STATEMENTS FOR CLEANER APPEARENCE ON PC:    
    2 FORMAT(' --------------------------------------------',                
     &       '-----------------------------------')                          
    3 FORMAT('      TIME  I  J  ZSUR  VELE  VELN  SAL    AV ',               
     &       '   AB  VELE  VELN  SAL    AV    AB')                           
    4 FORMAT('      STEP              SURF  SURF SURF  SURF ',               
     &       ' SURF  BOTT  BOTT BOTT  BOTT  BOTT')                           
    5 FORMAT(' ')                                                            
    6 FORMAT('                    CM  CM/S  CM/S  PPT CMSQS ',               
     &       'CMSQS  CM/S  CM/S  PPT CMSQS CMSQS')                           
    7 FORMAT(' ',I9,1X,I2,1X,I2,1X,F5.1,1X,F5.0,1X,                          
     &       F5.0,1X,F4.1,1X,F5.1,1X,F5.1,1X,F5.0,1X,                        
     &       F5.0,1X,F4.1,1X,F5.1,1X,F5.1)                                   

    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
     &       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
    9 FORMAT('| STEP |     |     |              (CM)              |',        
     &       1X,'       (PSU)       |       (PSU)        |')
   10 FORMAT('|      |     |     |',I3,'                          ',
     &        I3,'| 0               ',I3,'| 0               ',I3,'|')
   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
     &       20A1,'|')
   33 FORMAT('     MODEL   TIME   I   J    WSE   VELE   VELN    SAL',        
     &       '   VELE   VELN    SAL')                                        
   34 FORMAT('     MODEL   TIME   I   J   ZSUR   VELE   VELN    TEM',        
     +       '   VELE   VELN    TEM')                                        
   44 FORMAT('      TIME   STEP                  SURF   SURF   SURF',        
     &       '   BOTT   BOTT   BOTT')                                        
   66 FORMAT('      DAYS                    CM   CM/S   CM/S    PPT',        
     &       '   CM/S   CM/S    PPT')                                        
   67 FORMAT('      DAYS                    CM   CM/S   CM/S   MG/L',        
     +       '   CM/S   CM/S   MG/L')                                        
   68 FORMAT('      DAYS                    CM   CM/S   CM/S   DEGC',        
     +       '   CM/S   CM/S   DEGC')                                        
   77 FORMAT(' ',F9.5,1X,I6,1X,I3,1X,I3,1X,F6.1,1X,F6.1,1X,                  
     &       F6.1,1X,F6.1,1X,F6.1,1X,F6.1,1X,F6.1)                           

C NOTE: THE FOLLOWING FORMAT STATEMENTS WERE JOHN HAMRICK'S ORIGINAL 
C       CODE:   
C   2 FORMAT('----------------------------------------------------',
C     &       '-------------------------------------------')
C    3 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SAL |  AV ',
C     &       1X,'|  AB  | VELE | VELN | SAL |  AV  |  AB  |')
C   33 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SED |  AV ',
C     &       1X,'|  AB  | VELE | VELN | SED |  AV  |  AB  |')
C   34 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | TEM |  AV ',
C     &       1X,'|  AB  | VELE | VELN | TEM |  AV  |  AB  |')
C    4 FORMAT('| STEP |     |     |      | SURF | SURF | SUR | SURF',
C     &       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
C    5 FORMAT('|      |     |     |      |      |      |     |     ',
C     &       1X,'|      |      |      |     |      |      |')
C    6 FORMAT('|      |     |     |  CM  | CM/S | CM/S | PSU',
C     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | PSU |CMSQ/S|CMSQ/S|')
C    7 FORMAT('|',I6   ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
C     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
C     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
C    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
C     &       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
C    9 FORMAT('| STEP |     |     |              (CM)              |',
C     &       1X,'       (PSU)       |       (PSU)        |')
C   10 FORMAT('|      |     |     |',I3,'                          ',
C     &        I3,'| 0               ',I3,'| 0               ',I3,'|')
C   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
C     &       20A1,'|')
C   44 FORMAT('|      |     |     |      | SURF | SURF | SUR | SURF',
C     &       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
C   66 FORMAT('| DAYS |     |     |  CM  | CM/S | CM/S | PSU',
C     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | PSU |CMSQ/S|CMSQ/S|')
C   67 FORMAT('|      |     |     |  CM  | CM/S | CM/S |MG/L',
C     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S |MG/L |CMSQ/S|CMSQ/S|')
C   68 FORMAT('| DAYS |     |     |  CM  | CM/S | CM/S | D:C',
C     &       1X,'|CMSQ/S|CMSQ/S| CM/S | CM/S | D:C |CMSQ/S|CMSQ/S|')
C   77 FORMAT('|',F6.2 ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
C     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
C     &       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
C
C**********************************************************************C
C
      RETURN
      END
