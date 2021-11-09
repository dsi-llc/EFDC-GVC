      SUBROUTINE WASPHYDROLINKgvc
C------------------------------------------------------------------------------
C PURPOSE:
C
C   Subroutine WASPHYDROLINKgvc writes the hyd file for WINWASP water quality model
C                              when the generalized vertical grid (hybrid) is used.
C
C VARIABLE LIST:
C
C MODIFICATION HISTORY:
C
C   Date       Author         Comments
C   ---------- -------------- -------------------------------------------------
C   4/4/2006   Hugo Rodriguez Create the subroutine
C   4/11/2006  Hugo Rodriguez Add qfactor
c   5/23/2006  Hugo Rodriguez Add a limit maximum dispersion specific to certain cells
C------------------------------------------------------------------------------
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
      INCLUDE 'hydrolink_set.INT'
      parameter(nf=31000,nsg=12000)
C
      DIMENSION LAUX(ICM,JCM,KCM)
      INTEGER LU,LD,IU,ID,JU,JD,NAUX,nlay,SUR
      INTEGER istartyear,istartmonth,istartday,istarthour,istartminute
     *,istartsecond,Ihl_debug,Ihl_mode,inumsegconsts,j,
     *FLAGWASPBC(NQSERM,KCM)
      CHARACTER*6  sn(lcm),sn1
      CHARACTER*20 HYDFIL
      CHARACTER*23 segname
      CHARACTER*17 SN2
      character*80 DESCRIPTION,MODELERNAME
      character*3 Itext, Jtext, Ktext

	REAL*8  AUX,tmin,thour,tsec
        REAL  rinterval,crnu(nf),brintt(nf),flow(nf)
        REAL  SegVolume(nsg),SegDepth(nsg),SegVel(nsg)
        REAL  SegSalt(nsg),SegTemp(nsg),vol1,vol2,ad,adcoeff,abwmax                !hnr
        real  abwmx(nsg)                                                   !hnr
C
C
C**********************************************************************C
C
C **  READ CONTROL DATA FILE EFDC.WSP
C
C----------------------------------------------------------------------C
C
      SVPT=1.
      IF(NTSMMT.LT.NTSPTC)SVPT=0.
C
      IF(JSWASP.EQ.1) THEN

c      for jswasp=1 only first entry
        OPEN(1,FILE='EFDC.WSP',STATUS='UNKNOWN')
        WRITE(6,*)'EFDC.WSP opened'

C
C1**  READ CELL VOLUME PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
C
C2**  READ DIFFUSION PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) NRFLD,SCALR,CONVR,ISNKH,adcoeff,abwmax					!hnr

        DO LT=2,LALT                                      !hnr
          l=lij(illt(lt),jllt(lt))                        !hnr
          abwmx(l)=abwmax                                 !hnr
        END DO                                            !hnr
C
C3**  READ ADVECTION PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS
C
C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
C
        READ(1,1)
        READ(1,1)
        READ(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
C
        do i=1,5                                   !hnr
          read(1,*,err=11)                         !hnr
        end do                                     !hnr
11    continue                                   !hnr
        DO LT=2,LALT                               !hnr
          read(1,*,err=12)i,j,abwmax               !hnr
          l=lij(i,j)                               !hnr
          abwmx(l)=abwmax                          !hnr
        END DO                                     !hnr
12     continue                                   !hnr

        CLOSE(1)

        OPEN(1,FILE='ABmax.txt',STATUS='UNKNOWN')   !hnr
        write(1,*)'    I    J     ABmax'            !hnr
        DO LT=2,LALT                                      !hnr
          l=lij(illt(lt),jllt(lt))                        !hnr
          write(1,21)illt(lt),jllt(lt),abwmx(l)           !hnr
        END DO                                            !hnr
        close(1)
21      format(2I5,f10.6)
        WRITE(6,*)'EFDC.WSP read succesfully and ABmax.txt written'
C
    1   FORMAT (80X)
C
C     read qser file to check for flows only to some layers

        IF(NQSER.GE.1)THEN                                                !HNR
          OPEN(1,FILE='QSER.INP',STATUS='UNKNOWN')                        !HNR
C **  SKIP OVER TITLE AND AND HEADER LINES                              !HNR
C                                                                       !HNR
          DO IS=1,14                                                      !HNR
           READ(1,1)                                                      !HNR
          ENDDO                                                           !HNR
C                                                                       !HNR
          DO NS=1,NQSER                                                   !HNR
            READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),   !HNR
     &                   RMULADJ,ADDADJ,ICHGQS                          !HNR
            IF(ISO.GT.0) GOTO 860                                         !HNR
            IF(ISTYP.EQ.1)THEN                                            !HNR
              READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)                        !HNR
              do k=1,kc                                                   !HNR
                if(wkq(k).lt.1.e-6) flagwaspbC(ns,k)=1                    !HNR
              end do                                                      !HNR
              IF(ISO.GT.0) GOTO 860                                       !HNR
              DO M=1,MQSER(NS)                                            !HNR
                READ(1,*,IOSTAT=ISO)TQSER(M,NS),QSERTMP                   !HNR
                IF(ISO.GT.0) GOTO 860                                     !HNR
              END DO                                                      !HNR
            ELSE                                                          !HNR
              DO M=1,MQSER(NS)                                            !HNR
                READ(1,*,IOSTAT=ISO)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)    !HNR
                IF(ISO.GT.0) GOTO 860                                     !HNR
              END DO                                                      !HNR
            END IF                                                        !HNR
          END DO                                                          !HNR
        END IF                                                            !HNR
C                                                                       !HNR
        CLOSE(1)                                                          !HNR
        goto 862                                                          !HNR
860     WRITE(6,861)                                                      !HNR
861     FORMAT('  READ ERROR FOR FILE QSER.INP ')                         !HNR
862     continue                                                          !HNR
C
C**********************************************************************C
C
C **  DEFINE EFDC-WASP CORRESPONDENCE AND INITIALIZE FILES
C
C----------------------------------------------------------------------C
C
c=======================================================================
c      hlopen parameters
c=======================================================================
        Ihl_handle = 0
        Ihl_debug = 0
        Ihl_mode = 1	!Ihl_mode=0 to READ from dll hyd file, =1 to WRITE to dll hyd file
c=======================================================================
c     Set the debug flag 0=No debug 1=Debug (LOG.OUT)
c=======================================================================
        call hlsetdebug(Ihl_debug)
c=======================================================================
c     Open the file
c=======================================================================
        call hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
        if(ierror.gt.0) then
          call hlgetlasterror(errstring)
          write(6,6000)ierror,errstring,1
          stop
        end if
c=======================================================================
c     Set the language to FORTRAN
c=======================================================================
        call hlsetlanguage(Ihl_handle,1,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,2
          stop
        end if
c=======================================================================
c     Store a description string
c=======================================================================
        DESCRIPTION='   '//char(0)
        call hladddescription(Ihl_handle,0,DESCRIPTION,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,3
          stop
        end if
c=======================================================================
c     Store the modeler name
c=======================================================================
        MODELERNAME='Created by:  '//char(0)
        call hladddescription(Ihl_handle,1,MODELERNAME,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,4
          stop
        end if
c=======================================================================
c     Set the creator
c=======================================================================
        call hlsetcreator(Ihl_handle, 1, ierror)   !1 for fortran program
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,5
          stop
        end if
c=======================================================================
c     Set the seed moment (start date of the hyd file)
c=======================================================================
c     Calculate istartmonth, istartday,istartyear
        j=1461*(iyear+4800+int((imon-14)/12))/4+367*(imon-2-12*
     &    int((imon-14)/12))/12-(3*((iyear+4800+int((imon-14)/12)
     &    +100)/100))/4+iday-32075+tbegin+idays
	istartyear=100*(int(4*(j+68569)/146097)-49)+int(4000*
     &   (int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/
     &   1461001)+int(int(80*int(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*
     &   int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartmonth=int(80*int(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+
     &   68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)+2-12*int(int(80*
     &   int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartday=int(int(j+68569-(146097*int(4*int(j+68569)/146097)
     &   +3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)-2447*int(80*int(
     &   int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)+1)/1461001)/4+31)/2447)/80
	istarthour=0
	istartminute=0
	istartsecond=0

c     Adjust hour, minutes and seconds for the start day if necessary
	    IBEGIN=IDAYS*NTSPTC
	    AUX=FLOAT(IBEGIN)/FLOAT(NTSMMT)
	    NAUX=INT(AUX)
	      IF(AUX.GT.NAUX) then
	        AUX=(NAUX+1.-AUX)
                tsec=AUX*NTSMMT*TCON/NTSPTC
                istartsecond=INT(tsec)
                IF(tsec.GT.60) THEN
                  TMIN=tsec/60.0
                  istartminute=INT(TMIN)
                  istartsecond=INT((TMIN-istartminute)*60.)
                  IF(istartminute.GT.60) THEN
                    THOUR=FLOAT(istartminute)/60.0
                    istarthour=INT(THOUR)
                    istartminute=INT((THOUR-istarthour)*60.)
                  ENDIF
                ENDIF
              ENDIF
c=======================================================================

        call hlsetseedmoment(Ihl_handle,istartmonth,istartday,
     +     istartyear,istarthour,istartminute,istartsecond,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,6
          stop
        end if
c=======================================================================
c       instead of setting the number of layer, this sets the MAX number
C         of layers
	call hlsetnumlayers(Ihl_handle,kc,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,7
          stop
        end if

C        OPEN(1,FILE='hydcheck.txt',STATUS='UNKNOWN')                  !debug
C        write(1,*)'seed moment: m/d/y/h/min/s',istartmonth,istartday,  !debug 
C     +     istartyear,istarthour,istartminute,istartsecond            !debug
C        Write(1,*)'max number of layers',kc                              !debug
        
c=======================================================================
c     Set the number of segments
c=======================================================================
c        NJUN=KC*(LCLT-2)
        NJUN=0
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          nlay=kc-kgvcp(l)+1
          NJUN=NJUN+nlay
        END DO
        if(njun.gt.nsg) then
        write(6,500)NJUN,NSG
        STOP
        END IF
500     FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,
     +'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumsegments(Ihl_handle,njun,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,8
          stop
        end if
C        Write(1,*)'number of segments',njun                              !debug
c=======================================================================
c               Set WASP Segment Names (defaulft as I,J,K)
c=======================================================================
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          sn(l)='      '//char(0)
        END DO
        OPEN(94,FILE='segname.inp',iostat=ios,STATUS='old')
        IF(ios.eq.0) then
          do i=1,4
            read(94,*)
          end do
          do kk=1,la
            READ(94,*,err=111)I,J,SN1
            L=LIJ(I,J)
            sn(l)=sn1
          end do
        end if
111     CLOSE(94)
        I=0
        
        DO K=KC,1,-1
          DO LT=2,LALT
            l=lij(illt(lt),jllt(lt))
            if (LGVCP(L,K)) THEN
    	      I=I+1
	      LAUX(ILLT(LT),JLLT(LT),K)=I
              if(il(l).ge.100) then
                WRITE(itext,"(I3)")IL(L)
              elseif(il(l).ge.10) then
                WRITE(itext,"(I2)")IL(L)
              else
                WRITE(itext,"(I1)")IL(L)
              end if
              if(jl(l).ge.100) then
                WRITE(jtext,"(I3)")jL(L)
              elseif(jl(l).ge.10) then
                WRITE(jtext,"(I2)")jL(L)
              else
                WRITE(jtext,"(I1)")jL(L)
              end if
              if(k.ge.100) then
                WRITE(ktext,"(I3)")k
              elseif(k.ge.10) then
                WRITE(ktext,"(I2)")k
              else
                WRITE(ktext,"(I1)")k
              end if
	      sn2=' I='//itext//' J='//jtext//' K='//ktext
    	      segname=sn2//sn(l)//char(0)
	      call hlsetsegname(ihl_handle,i,segname,ierror)
C        Write(1,*)'segname',i,segname                              !debug
            END IF
          END DO
        END DO
c=======================================================================
c               Set WASP Segment Type and Segment Below
c=======================================================================
        
        OPEN(1,FILE='seginfo.wsp',STATUS='UNKNOWN')       
        DO K=KC,2,-1
          DO LT=2,LALT
            l=lij(illt(lt),jllt(lt))
            IF (LGVCP(L,K)) THEN
              LU=LAUX(ILLT(LT),JLLT(LT),K)
              IF (LGVCP(L,K-1)) THEN
	        LD=LAUX(ILLT(LT),JLLT(LT),K-1)
              ELSE
	        LD=0
              ENDIF
              SUR=2
              IF(K.EQ.KC) SUR=1
c	      call hlsetsegname(ihl_handle,i,segname,ierror)
C        Write(1,*)ILLT(LT),JLLT(LT),K,LU,LD,SUR                              !debug
              Write(1,*)LU,SUR,LD        
            END IF
          END DO
        END DO

        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          IF (LGVCP(L,1)) THEN
            LU=LAUX(ILLT(LT),JLLT(LT),1)
            LD=0
            SUR=2
c           call hlsetsegname(ihl_handle,i,segname,ierror)
C           Write(1,*)ILLT(LT),JLLT(LT),1,LU,LD,SUR                              !debug
            Write(1,*)LU,SUR,LD        
          END IF
        END DO
        close(1)
c=======================================================================
c     Set the number of flow paths
c=======================================================================
        NCHNH=0
        NCHNV=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCU(L,K)) THEN
              NCHNH=NCHNH+1
C              WRITE(1,*)I,J,K,'U'                         !debug
            END IF
            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (LGVCU(LIJ(I+1,J),K)) THEN
                NCHNH=NCHNH+1
C                WRITE(1,*)I,J,K,'U8'                         !debug
              END IF
            END IF
            IF (LGVCV(L,K)) THEN
              NCHNH=NCHNH+1
C             WRITE(1,*)I,J,K,'V'                          !debug
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              IF (LGVCV(LIJ(I,J+1),K)) THEN
                NCHNH=NCHNH+1
C                WRITE(1,*)I,J,K,'V8'                         !debug
              END IF
            END IF
          END DO
        END DO            
C        Write(1,*)'number of horizontal flows',nchnh                              !debug
C          NCHNH=NCHNH+INT(SUBO(L))
C          IF (IJCTLT(I+1,J).EQ.8) THEN
C            IF (SUBO(L+1).EQ.1.) NCHNH=NCHNH+1
C          END IF
C          NCHNH=NCHNH+INT(SVBO(L))
C          IF (IJCTLT(I,J+1).EQ.8) THEN
C            IF (SVBO(LNC(L)).EQ.1.) NCHNH=NCHNH+1
C          END IF

        DO LT=2,LALT
          DO k=1,KS
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCW(L,K)) THEN
              NCHNV=NCHNV+1
C              WRITE(1,*)I,J,K,'W'                            !debug
            END IF
          END DO
        END DO
C        Write(1,*)'number of vertical flows',nchnv                              !debug
        NCHN=NCHNH+NCHNV
C        Write(1,*)'number of hor+vert flows (w/o qser)',nchn                              !debug

C        NCHN=KC*NCHNH+(KC-1)*NCHNV

	DO LQ=1,NQSIJ
          I=IQS(LQ)                                     !HNR
          J=JQS(LQ)                                     !HNR
          L=LIJ(I,J)
	  IF(LIJLT(I,J).GT.0) THEN
            NS=NQSERQ(LQ)                               !HNR
	    DO K=KGVCP(L),KC
	      IF(flagwaspbC(ns,k).NE.1) THEN
	        NCHN=NCHN+1                            !HNR
	      END IF
	    END DO
	  END IF
	END DO
C        Write(1,*)'number of hor+vert flows (qser ADDED)',nchn                              !debug
c
	NQ=NQCTL
	DO LQ=1,NQCTL
	  I=IQCTLU(LQ)
	  J=JQCTLU(LQ)
          L=LIJ(I,J)
	  IF(LIJLT(I,J).GT.0) THEN
	    DO K=KGVCP(L),KC
              NCHN=NCHN+1                            !HNR
	    END DO
	  END IF
	END DO

        if(NCHN.gt.NF) then
        write(6,600)NCHN,NF
        STOP
        END IF
600     FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,
     +   'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
C        Write(1,*)'number of flows',nchn                              !debug
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,9
          stop
        end if
c=======================================================================
c     Set the flow path and direction
c=======================================================================
        LWASP=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCU(L,K)) THEN
              LWASP=LWASP+1
              LDTM=LAUX(I,J,K)
              LUTM=LDTM-1
              IF (IJCTLT(I-1,J).EQ.8) LUTM=0
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
C        Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM,1                              !debug
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,10
                  stop
               end if
            END IF
            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (LGVCU(L+1,K)) THEN
C              IF (SUBO(L+1).EQ.1.) THEN
	        LWASP=LWASP+1
                LDTM=0
                LUTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
C        Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM,1                              !debug
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,11
                  stop
               end if
              END IF
            END IF
          END DO

          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (LGVCV(L,K)) THEN
              LWASP=LWASP+1
              LSLT=LSCLT(LT)
              LDTM=LAUX(I,J,K)
              LUTM=LAUX(ILLT(LSLT),JLLT(LSLT),K)
              IF(IJCTLT(I,J-1).EQ.8) LUTM=0
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
C        Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM,2                              !debug
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,12
                  stop
               end if
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              LN=LNC(L)
              IF (LGVCV(LN,K)) THEN
c              IF (SVBO(LN).EQ.1.) THEN
                LWASP=LWASP+1
                LDTM=0
                LUTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
C        Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM,2                              !debug
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,13
                  stop
               end if
              END IF
            END IF
          END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQSIJ
            I=IQS(LT)
            J=JQS(LT)
            L=LIJ(I,J)
            IF(LIJLT(I,J).EQ.0) GOTO 100
            NS=NQSERQ(Lt)             
            IF(flagwaspbC(ns,k).EQ.1) GOTO 100                   !HNR
            IF (LGVCP(L,K)) THEN
              LWASP=LWASP+1
              LDTM=LAUX(I,J,K)
              LUTM=0
              call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
C              Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM                              !debug
              if(ierror .gt. 0)then
                call hlgetlasterror(errstring)
                write(6,6000) ierror, errstring,14
                stop
              end if
            END IF
100       END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQCTL
            I=IQCTLU(LT)
            J=JQCTLU(LT)
            LUTM=LAUX(I,J,K)
            I=IQCTLD(LT)
            J=JQCTLD(LT)
            LDTM=LAUX(I,J,K)
	      IF(LUTM.EQ.0.AND.LDTM.EQ.0) GOTO 200
              LWASP=LWASP+1
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,15
                  stop
               end if
200       END DO
        END DO

        IF(KC.GT.1)THEN
          DO K=KS,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCW(L,K)) THEN
C              IF (SWB(L).EQ.1.) THEN
                LWASP=LWASP+1
                LUTM=LAUX(I,J,K+1)
                LDTM=LAUX(I,J,K)
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,3,ierror)
C        Write(1,*)'FLOW No  FROM  TO',LWASP,LUTM,LDTM                              !debug
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring,16
                  stop
               end if
              END IF
            END DO
          END DO
        ENDIF
c=======================================================================
c     Set the number of segment constituents
c=======================================================================

	inumsegconsts=3  !volume,depth,velocity
	if(istran(2).GE.1) inumsegconsts=4   !temperature modeled in EFDC and transfered                       !hnr
	if(istran(1).GE.1) inumsegconsts=5   !salinity modeled in EFDC and transfered                          !hnr

        call hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,17
          stop
        end if
c=======================================================================
c     Set the number of flow path constituents
c=======================================================================
c when we add sed transport we need to add more
        call hlsetnumfpconsts(Ihl_handle, 3, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,18
          stop
        end if
c=======================================================================
c     Now we will set all the segment constituent types
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,19
          stop
        end if
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,20
          stop
        end if
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,21
          stop
        end if
c=======================================================================
c next only if temperature is transfered
        if(istran(2).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
        if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,22
          stop
        end if
        end if
c=======================================================================
c next only if salinity is transfered
        if(istran(1).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,23
          stop
        end if
        end if
c=======================================================================
c     Set all the flow constituent types
c=======================================================================
        call hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,24
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,25
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,26
          stop
        end if
c=======================================================================
	call hlsetvartimestep(Ihl_handle,0,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,27
          stop
        end if
c=======================================================================
	call hlsethydtimestep(Ihl_handle,dt,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,28
          stop
        end if
c=======================================================================
	rinterval=(dt*NTSMMT)/86400.
	call hlsetupdateint(Ihl_handle,rinterval,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,29
          stop
        end if
c=======================================================================
	call hlsethydtowaspratio(Ihl_handle,NTSMMT,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring,30
          stop
        end if
c=======================================================================
C
C     INITIAL CONDITIONS WHEN IDAYS=0
C
        IF(IDAYS.EQ.0) THEN
          IF(ISRESTI.EQ.0) THEN
C     INITIAL CONDITION FOR A COLD START
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                SegVel(LWASP)=0.0
                SegDepth(LWASP)=HP(L)*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
                IF(NTSMMT.LT.NTSPTC) THEN
                  SegDepth(LWASP)=HMP(L)*DZC(K)*GVCSCLP(L)
                  SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                  SegSalt(LWASP)=SALINIT(L,K)
                  SegTemp(LWASP)=TEMINIT(L,K)
                END IF
              END IF
            END DO
          END DO
          DO I=1,NCHN
	    FLOW(I)=0.0
	    CRNU(I)=0.0
	    BRINTT(I)=0.0
          END DO
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,31
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,32
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,33
            stop
          end if
c         next only if temperature is transfered
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,34
              stop
            end if
          end if
c         next only if salinity is transfered
          if(istran(1).GE.1) then
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,35
              stop
            end if
          end if
c=======================================================================
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,36
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,37
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,38
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,39
            stop
          end if
          ELSE
C     INITIAL CONDITIONS FROM A RESTART FILE
          LWASP=0
C     Segment Properties from the RESTART file
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                LN=LNC(L)
                VELX=0.5*(U(L,K)+U(L+1,K))
                VELY=0.5*(V(L,K)+V(LN,K))
                VELZ=0.5*(W(L,K-1)+W(L,K))
                VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
                SegVel(LWASP)=VELMAG
                SegDepth(LWASP)=HP(L)*DZC(K)*GVCSCLP(L)
                SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP)=SAL(L,K)
                SegTemp(LWASP)=TEM(L,K)
              END IF
            END DO
          END DO
C Advection and dispersion in the X-direction:
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLW=0.0
              IF (LGVCU(L,K)) THEN  
                vol1=DYP(L)*AH(L,K)*DZC(K)*
     +              GVCSCLP(L)*HP(L)/DXP(L)
                vol2=DYP(L-1)*AH(L-1,K)*DZC(K)*
     +               GVCSCLP(L-1)*HP(L-1)/DXP(L-1)
                ADDLW=0.5*(vol1+vol2)
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)         !hnr
                  vol2=DXYP(L-1)*HP(L-1)*DZC(K)*GVCSCLP(L-1)       !hnr
                  if (vol1.LT.vol2) vol2=vol1                    !hnr
                  ad=vol2/dt*adcoeff                             !hnr
                  if (addlw.GT.ad) then                          !hnr
                    addlw=ad                                     !hnr
                  end if                                         !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addlw=0.                                     !hnr
                  end if                                         !hnr
                LWASP=LWASP+1
                vol1=DYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                vol2=DYP(L-1)*HP(L-1)*DZC(K)*GVCSCLP(L-1)
                FLOW(LWASP)=U(L,K)*0.5*(vol1+vol2)
	        vol1=1.0/DXP(L)
	        vol2=1.0/DXP(L-1)
	        CRNU(LWASP)=U(L,K)*0.5*(vol1+vol2)
  	        BRINTT(LWASP)=ADDLW
              END IF
              IF (IJCTLT(I+1,J).EQ.8) THEN
                IF (LGVCU(L+1,K)) THEN  
C                IF (SUBO(L+1).EQ.1.) THEN
                  LWASP=LWASP+1
                  vol1=DYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DYP(L+1)*HP(L+1)*DZC(K)*GVCSCLP(L+1)
	          FLOW(LWASP)=U(L+1,K)*0.5*(vol1+vol2)
                  vol1=1.0/DXP(L)
	          vol2=1.0/DXP(L+1)
	          CRNU(LWASP)=U(L+1,K)*0.5*(vol1+vol2)
                  BRINTT(LWASP)=0.0                              !hnr
                END IF
              END IF
            END DO
C
C Advection and dispersion in the Y-direction:
C
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLS=0.0
              LS=LSC(L)
              IF (LGVCV(L,K)) THEN
c              IF (SVB(L).EQ.1.) THEN
                vol1=DXP(L)*AH(L,K)*DZC(K)*
     +              GVCSCLP(L)*HP(L)/DYP(L)
                vol2=DXP(LS)*AH(LS,K)*DZC(K)*
     +               GVCSCLP(LS)*HP(LS)/DYP(LS)
                ADDLS=0.5*(vol1+vol2)
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)         !hnr
                  vol2=DXYP(Ls)*HP(Ls)*DZC(K)*GVCSCLP(Ls)      !hnr
                  if (vol1.LT.vol2) vol2=vol1                    !hnr
                  ad=vol2/dt*adcoeff                             !hnr
                  if (addls.GT.ad) then                          !hnr
                    addls=ad                                     !hnr
                  end if                                         !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addls=0.                                     !hnr
                  end if                                         !hnr
	        LWASP=LWASP+1
                vol1=DXP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                vol2=DXP(LS)*HP(LS)*DZC(K)*GVCSCLP(LS)
                FLOW(LWASP)=V(L,K)*0.5*(vol1+vol2)
                  vol1=1.0/DYP(L)
	          vol2=1.0/DYP(LS)
	          CRNU(LWASP)=V(L,K)*0.5*(vol1+vol2)
                BRINTT(LWASP)=ADDLS
              END IF
              IF (IJCTLT(I,J+1).EQ.8) THEN
                LN=LNC(L)
                IF (LGVCV(LN,K)) THEN
c                IF (SVBO(LN).EQ.1.) THEN
  	          LWASP=LWASP+1
                  vol1=DXP(L)*HP(L)*DZC(K)*GVCSCLP(L)
                  vol2=DXP(LN)*HP(LN)*DZC(K)*GVCSCLP(LN)
                  FLOW(LWASP)=V(LN,K)*0.5*(vol1+vol2)
                  vol1=1.0/DYP(L)
	          vol2=1.0/DYP(LN)
	          CRNU(LWASP)=V(LN,K)*0.5*(vol1+vol2)
	          BRINTT(LWASP)=0.0                                !hnr
                END IF
              END IF
            END DO
          END DO

c Advection in input flows
          DO K=KC,1,-1
            DO LT=1,nqsij
              I=Iqs(LT)
              J=Jqs(LT)
              IF(LIJLT(I,J).EQ.0) GOTO 310
              NS=NQSERQ(Lt)
              L=LQS(Lt)
              IF(flagwaspbC(ns,k).EQ.1) GOTO 310                       !HNR
              IF (LGVCP(L,K)) THEN
                LWASP=LWASP+1
                FLOW(LWASP)=RQSMUL(Lt)*(QSS(K,Lt)+
     &                   qfactor(lt)*QSERT(K,NS))
                CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)
     $                     /(dzc(k)*HP(L)*GVCSCLP(L))
                BRINTT(LWASP)=0.0
              END IF
310         END DO
          END DO
c  Advection in structure flows
          DO K=KC,1,-1
            DO LT=1,nqctl
              Iu=Iqctlu(LT)
              Ju=Jqctlu(LT)
              Lu=Lij(iu,ju)
              Id=Iqctld(LT)
              Jd=Jqctld(LT)
              Ld=Lij(id,jd)
              IF(LU.EQ.0.AND.LD.EQ.0) GOTO 410
              flowx=RQCMUL(Lt)*QCTLT(K,Lt)
              UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HP(L)*GVCSCLP(L))
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
              BRINTT(LWASP)=0.0
410         END DO
          END DO
C
C Advection and dispersion in the Z-direction:
C
          IF (KC.GT.1) THEN
            DO K=KS,1,-1
              DO LT=2,LALT
                I=ILLT(LT)
                J=JLLT(LT)
                L=LIJ(I,J)
                ADDL=ab(l,k)*hp(l)                            !hnr  ab is Ev/H
                if (addl.gt.abwmx(l)) then
                  ab(l,k)=abwmx(l)/hp(l)                      !hnr  ab is Ev/H
                end if
                addl=0.0
                IF (LGVCW(L,K)) THEN
                  ADDL=DXYP(L)*AB(L,K)*DZIG(K)*GVCSCLPI(L)                    !hnr  ab is Ev/H
                  vol1=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)                           !hnr
                  vol2=DXYP(L)*HP(L)*DZC(K+1)*GVCSCLP(L)                         !hnr
                  if (vol1.LT.vol2) vol2=vol1                                      !hnr
                  ad=adcoeff*vol2/dt                                               !hnr
                  if (addl.GT.ad) then                                             !hnr
                    addl=ad                                                        !hnr
                  end if                                                           !hnr
                  LWASP=LWASP+1
	          FLOW(LWASP)=-DXYP(L)*W(L,K)
                  CRNU(LWASP)=W(L,K)*DZIG(K)*GVCSCLPI(L)/HP(L)
	          BRINTT(LWASP)=ADDL
                END IF
              END DO
            END DO
          END IF
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,40
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,41
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,42
            stop
          end if
c next only if temperature is transfered
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,43
              stop
            end if
          end if
c next only if salinity is transfered
          if(istran(1).GE.1) then
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring,44
              stop
            end if
          end if
c=======================================================================
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,45
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,46
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,47
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,48
            stop
          end if
        END IF
C      END INITIAL CONDITIONS WHEN IDAYS=0
      END IF
C     FINISH INITIALIZATION OF FILES (JSWASP=1)
      GOTO 3000
      END IF
      IF(N.LT.IBEGIN) GOTO 3000
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP DATA
C
C Advection and dispersion in the X-direction:
C
      LWASP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLW=0.0
          IF (LGVCU(L,K)) THEN
            LW=L-1
            vol1=DYP(L)*DZC(K)*HLPF(L)*GVCSCLP(L)
            vol2=DYP(LW)*DZC(K)*HLPF(LW)*GVCSCLP(LW)
            ADDLW=0.5*(vol1+vol2)*AHULPF(L,K)*DXIU(L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                !hnr
            vol2=DXYP(LW)*HLPF(LW)*DZC(K)*GVCSCLP(LW)             !hnr
            if (vol1.LT.vol2) vol2=vol1                           !hnr
            ad=adcoeff*vol2/dt                                    !hnr
            if (addlw.GT.ad) then                                 !hnr
              addlw=ad                                            !hnr
            end if                                                !hnr
            if(ishdmf.lt.2) then                                  !hnr
              addlw=0.                                            !hnr
            end if                                                !hnr
            vol1=DYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
            vol2=DYP(L-1)*HLPF(L-1)*DZC(K)*GVCSCLP(L-1)
c            FLOWX=0.5*(vol1+vol2)*ULPF(L,K)
            FLOWX=UHLPF(L,K)*DYU(L)*DZC(K)*GVCSCLU(L)
c            UDDXTMP=ULPF(L,K)*DXIU(L)
            UDDXTMP=UHLPF(L,K)*DXIU(L)/HU(L)
            LWASP=LWASP+1
c              write(6,*)lwasp,vol1,vol2,ulpf(l,k)
            FLOW(LWASP)=FLOWX
            CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=ADDLW
          END IF
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (LGVCU(L+1,K)) THEN
              vol1=DYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
              vol2=DYP(L+1)*HLPF(L+1)*DZC(K)*GVCSCLP(L+1)
c              FLOWX=0.5*(vol1+vol2)*ULPF(L+1,K)
              FLOWX=UHLPF(L+1,K)*DYU(L+1)*DZC(K)*GVCSCLU(L+1)
c              UDDXTMP=ULPF(L+1,K)*DXIU(L+1)
              UDDXTMP=UHLPF(L+1,K)*DXIU(L+1)/HU(L+1)
              vol1=DYP(L)*DZC(K)*HLPF(L)*GVCSCLP(L)
              vol2=DYP(L+1)*DZC(K)*HLPF(L+1)*GVCSCLP(L+1)
              ADDLW=0.5*(vol1+vol2)*AHULPF(L+1,K)*DXIU(L)
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                !hnr
              vol2=DXYP(L+1)*HLPF(L+1)*DZC(K)*GVCSCLP(L+1)             !hnr
              if (vol1.LT.vol2) vol2=vol1                           !hnr
              ad=adcoeff*vol2/dt                                    !hnr
              if (addlw.GT.ad) then                                 !hnr
                addlw=ad                                            !hnr
              end if                                                !hnr
              if(ishdmf.lt.2) then                                  !hnr
                addlw=0.                                            !hnr
              end if                                                !hnr
              LWASP=LWASP+1
c              write(6,*)lwasp,flowx
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
	      BRINTT(LWASP)=ADDLW
            END IF
          END IF
        END DO
C
C Advection and dispersion in the Y-direction:
C
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLS=0.0
          IF (LGVCV(L,K)) THEN
            LS=LSC(L)
            vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       !hnr
            vol2=DXP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)                     !hnr
            ADDLS=AHVLPF(L,K)*0.5*(vol1+vol2)*DYIV(L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       !hnr
            vol2=DXYP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)                    !hnr
            if (vol1.LT.vol2) vol2=vol1                                  !hnr
            ad=adcoeff*vol2/dt                                           !hnr
            if (addls.GT.ad) then                                        !hnr
              addls=ad                                                   !hnr
            end if                                                       !hnr
            if(ishdmf.lt.2) then                                         !hnr
              addls=0.                                                   !hnr
            end if                                                       !hnr
            vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
            vol2=DXP(LS)*HLPF(LS)*DZC(K)*GVCSCLP(LS)
c            FLOWY=0.5*(VOL1+VOL2)*VLPF(L,K)
            FLOWY=VHLPF(L,K)*DXV(L)*DZC(K)*GVCSCLV(L)
c            VDDYTMP=VLPF(L,K)*DYIV(L)
            VDDYTMP=VHLPF(L,K)*DYIV(L)/HV(L)
	    LWASP=LWASP+1
c              write(6,*)lwasp,flowy,vddytmp,addls
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
          END IF
          IF (IJCTLT(I,J+1).EQ.8) THEN
            LN=LNC(L)
            IF (LGVCV(LN,K)) THEN
              vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
              vol2=DXP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)
c              FLOWY=0.5*(VOL1+VOL2)*VLPF(LN,K)
              FLOWY=VHLPF(LN,K)*DXV(LN)*DZC(K)*GVCSCLV(LN)
c              VDDYTMP=VLPF(LN,K)*DYIV(LN)
              VDDYTMP=VHLPF(LN,K)*DYIV(LN)/HV(LN)
              vol1=DXP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       !hnr
              vol2=DXP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)                     !hnr
              ADDLS=AHVLPF(LN,K)*0.5*(vol1+vol2)*DYIV(LN)
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                       !hnr
              vol2=DXYP(LN)*HLPF(LN)*DZC(K)*GVCSCLP(LN)                    !hnr
              if (vol1.LT.vol2) vol2=vol1                                  !hnr
              ad=adcoeff*vol2/dt                                           !hnr
              if (addls.GT.ad) then                                        !hnr
                addls=ad                                                   !hnr
              end if                                                       !hnr
              if(ishdmf.lt.2) then                                         !hnr
                addls=0.                                                   !hnr
              end if                                                       !hnr
              LWASP=LWASP+1
c              write(6,*)lwasp,flowy
	      FLOW(LWASP)=FLOWY
              CRNU(LWASP)=VDDYTMP
	      BRINTT(LWASP)=ADDLS
            END IF
          END IF
        END DO
      END DO

c Advection in input flows
      DO K=KC,1,-1
        DO LT=1,nqsij
          I=Iqs(LT)
          J=Jqs(LT)
          IF(LIJLT(I,J).EQ.0) GOTO 300
          NS=NQSERQ(Lt)
          L=LQS(Lt)
          IF(flagwaspbC(ns,k).EQ.1) GOTO 300                       !HNR
          IF (LGVCP(L,K)) THEN
            flowx=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
            UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*GVCSCLP(L)*HLPF(L))
            LWASP=LWASP+1
c              write(6,*)lwasp,flowx
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=0.0
          END IF
300     END DO
      END DO
c  Advection and dispersion in structure flows
      DO K=KC,1,-1
        DO LT=1,nqctl
          Iu=Iqctlu(LT)
          Ju=Jqctlu(LT)
          Lu=Lij(iu,ju)
          Id=Iqctld(LT)
          Jd=Jqctld(LT)
          Ld=Lij(id,jd)
          IF(LU.EQ.0.AND.LD.EQ.0) GOTO 400
          flowx=RQCMUL(Lt)*QCTLT(K,Lt)
          UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*GVCSCLP(L)*HmP(L))
          if(ishdmf.lt.2) then                                            !hnr
            addls=0.                                                      !hnr
          end if                                                          !hnr
          LWASP=LWASP+1
	  FLOW(LWASP)=FLOWX
	  CRNU(LWASP)=UDDXTMP
	  BRINTT(LWASP)=0.0
400     END DO
      END DO
C
C Advection and dispersion in the Z-direction:
C
      IF (KC.GT.1) THEN
        DO K=KS,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDL=ablpf(l,k)*hlpf(l)                                                 !hnr  ab is Ev/H
            if(addl.gt.abwmx(l)) then
              ablpf(l,k)=abwmx(l)/hlpf(l)                                             !hnr
            end if
            ADDL=0.0
            IF (LGVCW(L,K)) THEN
              ADDL=DXYP(L)*ABLPF(L,K)*DZIG(K)*GVCSCLPI(L)             !hnr  ab is Ev/H
              vol1=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)                        !hnr
              vol2=DXYP(L)*HLPF(L)*DZC(K+1)*GVCSCLP(L)                      !hnr
              if (vol1.LT.vol2) vol2=vol1                                   !hnr
              ad=adcoeff*vol2/dt                                            !hnr
              if (addl.GT.ad) then                                          !hnr
                addl=ad                                                     !hnr
              end if                                                        !hnr
              FLOWZ=-DXYP(L)*wlpf(l,k)
              WDDZTMP=wlpf(l,k)*DZIG(K)*GVCSCLPI(L)/HLPF(L)
              KPTMP=K+1
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWZ
	      CRNU(LWASP)=WDDZTMP
	      BRINTT(LWASP)=ADDL
            END IF
          END DO
        END DO
      END IF
C
C Segment Properties:
C
      LCELTMP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          LN=LNC(L)
          IF (LGVCP(L,K)) THEN
            LCELTMP=LCELTMP+1
            VOLUM=DXYP(L)*HLPF(L)*DZC(K)*GVCSCLP(L)
            IF(NTSMMT.LT.NTSPTC) THEN
              VOLUM=DXYP(L)*HP(L)*DZC(K)*GVCSCLP(L)
            END IF
            DEPTH=HLPF(L)*DZC(K)*GVCSCLP(L)
            VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) +UHLPF(L+1,K)+SVPT*UVPT
     +      (L+1,K))/HLPF(L)
            VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) +VHLPF(LN,K)+SVPT*VVPT
     +      (LN,K) )/HLPF(L)
            VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) +WLPF(L,K)+SVPT*WVPT
     +      (L,K) )
            VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
            SegSalt(LCELTMP)=SALLPF(L,K)
            SegTemp(LCELTMP)=TEMLPF(L,K)
	    SegVolume(LCELTMP)=VOLUM
	    SegDepth(LCELTMP)=DEPTH
	    SegVel(LCELTMP)=VELMAG
          END IF
        END DO
      END DO
C
	   call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,49
            stop
         end if
	   call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,50
            stop
         end if
	   call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,51
            stop
         end if
c next only if temperature is transfered
       if(istran(2).GE.1) then
	   call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,52
            stop
         end if
       end if
c next only if salinity is transfered
       if(istran(1).GE.1) then
	   call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,53
            stop
         end if
       end if
c=======================================================================
c          do ii=1,nchn
c           write(6,*)ii,flow(ii),crnu(ii),brintt(ii)
c           end do
           
	   call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,54
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,55
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,56
            stop
         end if
	   call hlmomentcomplete(Ihl_Handle,ierror)
         if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring,57
            stop
         end if

C
3000  JSWASP=0
6000  format('Error ',I10, ' : ', A30,I3)
C----------------------------------------------------------------------c
C ** end SUBROUTINE WASPHYDROLINKGVC
C----------------------------------------------------------------------c
      RETURN
      END
