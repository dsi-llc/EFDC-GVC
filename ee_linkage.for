      SUBROUTINE EE_LINKAGE(JSEXPLORER)

      !----------------------------------------------------------------

      ! **  SUBROUTINE EE_LINKAGE (OLD EE_LINKAGE.FOR) WRITES BINARY OUTPUT FILES:
      ! **    EE_WC     - WATER COLUMN AND TOP LAYER OF SEDIMENTS
      ! **    EE_BED    - SEDIMENT BED LAYER INFORMATION
      ! **    EE_WQ     - WATER QUALITY INFORMATION FOR THE WATER COLUMN
      ! **    EE_SD     - SEDIMENT DIAGENSIS INFORMATION
      ! **    EE_ARRAYS - GENERAL/USER DEFINED ARRAY DUMP. LINKED TO  
      ! **                EFDC_EXPLORER FOR DISPLAY
      !----------------------------------------------------------------

      ! *** Notes:

      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
 
      ! *** EFDC_EXPLORER LINKAGE IS DESIGNED FOR SINGLE PRECISION VARIABLES
      INTEGER(4) VER
      INTEGER(4) IWQ(40), NACTIVE, NWQVAR, LC2
      INTEGER(4) JSEXPLORER,NS,NW,MW,NSEDSTEPS,NBEDSTEPS
      INTEGER(4) L,K,ISYS,NT,NX,N1,NSXD,HSIZE,NVAR,ITYPE,ITIMEVAR

      REAL(4)    TMPVAL,WQ
      REAl(4)    ZERO, SHEAR
      REAL(8)    EETIME
  
      LOGICAL(4) LWRITE
  
      CHARACTER*8 ARRAYNAME
  
      SAVE IWQ
      SAVE NSEDSTEPS, NBEDSTEPS, NWQVAR

      IF(ISDYNSTP == 0)THEN  
        DELT=DT  
      ELSE  
        DELT=DTDYN  
      ENDIF  
      NACTIVE=LA-1

      ! ***************************************************************************************
      ! **  INITIAL CALL TO WRITE HEADERS AND INITIAL CONDITIONS AFTER FULL INITIALIZATION
      ! ***************************************************************************************
      IF(JSEXPLORER == 1)THEN
        ! *** WATER DEPTHS (OLD SURFPLT)
        OPEN(95,FILE='EE_WS.OUT',STATUS='UNKNOWN')
        CLOSE(95,STATUS='DELETE')
        OPEN(95,FILE='EE_WS.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &       FORM='BINARY')
        VER=102
        WRITE(95)VER,INT(IC,4),INT(JC,4),NACTIVE
        CLOSE(95,STATUS='KEEP')

        ! *** VELOCITIES (OLD VELPLTH)
        OPEN(95,FILE='EE_VEL.OUT',STATUS='UNKNOWN')
        CLOSE(95,STATUS='DELETE')
        OPEN(95,FILE='EE_VEL.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &      FORM='BINARY')
        VER=104
        LC2 = LC+1
        WRITE(95)VER,INT(IC,4),INT(JC,4),INT(KC,4),NACTIVE
        WRITE(95)REAL(RSSBCE(1:LC2),4),REAL(RSSBCW(1:LC2),4),
     &           REAL(RSSBCS(1:LC2),4),REAL(RSSBCN(1:LC2),4)
        CLOSE(95,STATUS='KEEP')
        
        IF(ISSPH(8) >= 1)THEN
          ! *** WATER COLUMN AND TOP LAYER OF SEDIMENT
          OPEN(95,FILE='EE_WC.OUT',STATUS='UNKNOWN')
          CLOSE(95,STATUS='DELETE')
          OPEN(95,FILE='EE_WC.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &         FORM='BINARY')
          !------ Header contains Support Information for Output Extraction ------
          VER=109
          HSIZE = 80 + (NSED+NSND)*4    !2*4+3*4+7*4+3*4+5*4 + (NSED+NSND)*4
          WRITE(95) VER,HSIZE
          WRITE(95) INT(KC,4),INT(KB,4),NACTIVE
          WRITE(95) (INT(ISTRAN(I),4),I=1,7)
          WRITE(95) INT(NSED,4),INT(NSND,4),INT(NTOX,4)
          LWRITE = .FALSE.
          WRITE(95) INT(ISWAVE,4),INT(ISBEDSTR,4),LWRITE,INT(ISBDLDBC,4)
     &             ,REAL(TBEDIT(2),4)

          NSXD=NSED+NSND
          DO NS=1,NSXD
            WRITE(95)REAL(SEDDIA(NS),4)
          ENDDO
          CLOSE(95,STATUS='KEEP')
        ENDIF

        ! *** SEDIMENT BED LAYERS
        IF(ISBEXP >= 1 .AND. KB > 1)THEN
          IF(ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1)THEN
            OPEN(10,FILE='EE_BED.OUT',STATUS='UNKNOWN')
            CLOSE(10,STATUS='DELETE')
            OPEN(10,FILE='EE_BED.OUT',STATUS='UNKNOWN',
     &           ACCESS='SEQUENTIAL',FORM='BINARY')
            !------ Header contains Support Information for Output Extraction ------
            VER=108
            HSIZE = 60 + (NSED+NSND)*4    !2*4+3*4+7*4+3*4 + (NSED+NSND)*4
            WRITE(10) VER,HSIZE
            WRITE(10) INT(KC,4),INT(KB,4),NACTIVE
            WRITE(10) (INT(ISTRAN(I),4),I=1,7)
            WRITE(10) INT(NSED,4),INT(NSND,4),INT(NTOX,4)

            DO NS=1,NSXD
              WRITE(10)REAL(SEDDIA(NS),4)
            ENDDO
            CLOSE(10,STATUS='KEEP')
            NBEDSTEPS=ISBEXP-1
          ENDIF
        ENDIF

        ! *** WATER QUALITY MODEL (HEM3D) RESULTS
        IF(ISTRAN(8) > 0)THEN
          OPEN(95,FILE='EE_WQ.OUT',STATUS='UNKNOWN')
          CLOSE(95,STATUS='DELETE')
          OPEN(95,FILE='EE_WQ.OUT',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &      FORM='BINARY')
          NWQVAR=NWQV
          IF(IDNOTRVA > 0)NWQVAR = NWQVAR+1   ! *** Macroalgae always needs to be the last WQ variable since it is not advected
          !------ Header contains Support Information for Output Extraction ------
          VER=102
          HSIZE = 24 + NWQVAR*4    !2*4+2*4+1*4 + NWQVAR*4
          WRITE(95) VER,HSIZE
          WRITE(95) INT(KC,4),INT(KB,4),NACTIVE
          WRITE(95) INT(NWQVAR,4)

          IWQ=0
          DO MW=1,NWQVAR
            IWQ(MW)=ISTRWQ(MW)
          ENDDO
          IF(IDNOTRVA > 0)IWQ(NWQVAR)=1

          WRITE(95)(INT(IWQ(NW),4),NW=1,NWQVAR)
          CLOSE(95,STATUS='KEEP')

          ! *** SAVE SEDIMENT DIAGENESIS RESULTS
          IF(ISSDBIN < 0)THEN
            OPEN(95,FILE='EE_SD.OUT',STATUS='UNKNOWN')
            CLOSE(95,STATUS='DELETE')
            OPEN(95,FILE='EE_SD.OUT',STATUS='UNKNOWN',
     &          ACCESS='SEQUENTIAL',FORM='BINARY')
            !------ Header contains Support Information for Output Extraction ------
            VER=101
            HSIZE = 20              !2*4+3*4+1*4
            WRITE(95) VER,HSIZE
            WRITE(95) INT(KC,4),INT(KB,4),NACTIVE
            CLOSE(95,STATUS='KEEP')
            NSEDSTEPS=-1
          ENDIF
      
          IF(ISRPEM > 0)THEN
            OPEN(95,FILE='EE_RPEM.OUT',STATUS='UNKNOWN')
            CLOSE(95,STATUS='DELETE')
            OPEN(95,FILE='EE_RPEM.OUT',STATUS='UNKNOWN',
     &            ACCESS='SEQUENTIAL',FORM='BINARY')
            !------ Header contains Support Information for Output Extraction ------
            VER=101
            HSIZE = 28 + NACTIVE*4             !2*4+3*4+2*4 + (LA-1)*4
            WRITE(95) VER,HSIZE
            WRITE(95) INT(KC,4),INT(KB,4),NACTIVE
            WRITE(95) INT(NRPEM,4),INT(INITRPEM,4)
            DO L=2,LA
              LWRITE =.false. ! LMASKRPEM(L)
              WRITE(95) LWRITE
            ENDDO
            CLOSE(95,STATUS='KEEP')
          ENDIF
      
        ENDIF

        IF(ISINWV == 2)THEN
          OPEN(95,FILE='EE_ARRAYS.OUT',STATUS='UNKNOWN')
          CLOSE(95,STATUS='DELETE')
          OPEN(95,FILE='EE_ARRAYS.OUT',STATUS='UNKNOWN',
     &          ACCESS='SEQUENTIAL',FORM='BINARY')
          VER=100
          WRITE(95)VER
          NVAR = 3
          WRITE(95)NVAR    ! # OF TIME VARYING ARRAYS

          ! FLAGS: ARRAY TYPE, TIME VARIABLE
          ! ARRAY TYPE:    0 = L            DIM'D
          !                1 = L,KC         DIM'D
          !                2 = L,0:KC       DIM'D
          !                3 = L,KB         DIM'D
          !                4 = L,KC,NCLASS  DIM'D
          ! TIME VARIABLE: 0 = NOT CHANGING
          !                1 = TIME VARYING

          ! *****************************************************************
          ! *** THE FOLLOWING ARRAYS ARE TIME INVARIANT SO ONLY WRITTEN ONCE
          ! *****************************************************************

          !WRITE(95)0,0
          !ARRAYNAME='SUB'
          !WRITE(95)ARRAYNAME
          !DO L=2,LA
          !  WRITE(95)SUB3D(L,K)
          !ENDDO

          !WRITE(95)0,0
          !ARRAYNAME='SVB'
          !WRITE(95)ARRAYNAME
          !DO L=2,LA
          !  WRITE(95)SVB(L)
          !ENDDO
          !ITIMEVAR = 0
          !ITYPE = 1
          !WRITE(95)ITYPE,ITIMEVAR
          !ARRAYNAME='FXWAVE'
          !WRITE(95)ARRAYNAME
          !DO K=1,KC
          !  DO L=2,LA
          !    WRITE(95)REAL(FXWAVE(L,K),4)
          !  ENDDO
          !ENDDO
          !
          !ITYPE = 1
          !WRITE(95)ITYPE,ITIMEVAR
          !ARRAYNAME='FYWAVE'
          !WRITE(95)ARRAYNAME
          !DO K=1,KC
          !  DO L=2,LA
          !    WRITE(95)REAL(FYWAVE(L,K),4)
          !  ENDDO
          !ENDDO

          CLOSE(95,STATUS='KEEP')

        ENDIF

        !------ END OF HEADER SECTION ------
    
      ELSEIF(JSEXPLORER == -1)THEN
        ! *** FORCE ALL OUTPUT
        NSEDSTEPS=32000
        NBEDSTEPS=32000
     
      ENDIF

      ! ***************************************************************************************
      ! *** EFDC_EXPLORER LINKAGE OUTPUT SECTION
      ! ***************************************************************************************

      ! *** SET TIMING FOR INITIAL CONDITION
      IF(ISDYNSTP == 0)THEN
        EETIME=DT*FLOAT(N)+TCON*TBEGIN
      ELSE
        EETIME=TIMESEC
      ENDIF
      EETIME=EETIME/86400.
      WRITE(6, '( A,F10.4 )' )'EFDC_EXPLORER LINKAGE SNAPSHOT: ',EETIME
      
      ! *** WATER DEPTHS (OLD SURFPLT)
      OPEN(95,FILE='EE_WS.OUT',POSITION='APPEND',STATUS='OLD',
     &     FORM='BINARY')
      WRITE (95)INT(N,4),REAL(EETIME,4),REAL(DELT,4)
      DO L=2,LA
        WRITE(95)REAL(HP(L),4)
      ENDDO
      CALL FLUSH(95)
      CLOSE(95,STATUS='KEEP')

      ! *** VELOCITIES (OLD VELPLTH)
      OPEN(95,FILE='EE_VEL.OUT',POSITION='APPEND',STATUS='OLD',
     &        FORM='BINARY')
      WRITE (95)INT(N,4),REAL(EETIME,4),REAL(DELT,4)

      ! *** Write the UVW Instantaneous Velocity Field (Unrotated)
      DO L=2,LA
        WRITE(95)(REAL(U(L,K),4),REAL(V(L,K),4),REAL(W(L,K),4),K=1,KC)
      ENDDO
      CALL FLUSH(95)
      CLOSE(95,STATUS='KEEP')

      IF(ISSPH(8) >= 1)THEN
        ! *** WATER COLUMN AND TOP LAYER OF SEDIMENT
        OPEN(95,FILE='EE_WC.OUT',STATUS='OLD',POSITION='APPEND',
     &          FORM='BINARY')
        WRITE(95)REAL(EETIME,4),NACTIVE

        ! *** WRITE THE TOP LAYER INDEX
        IF(ISTRAN(6) == 1 .OR. ISTRAN(7) >= 1)THEN
          DO L=2,LA
            WRITE(95)INT(KBT(L),4)
          ENDDO
        ENDIF
    
        ! *** WRITE THE WATER COLUMN AND TOP LAYER OF SEDIMENT DATA, IF NEEDED
        DO L=2,LA
          N1=KBT(L)
          IF(ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)THEN
            IF(ISBEDSTR >= 1)THEN
              WRITE(95)REAL(TAUBSED(L),4)
              IF(ISBEDSTR == 1)THEN
                WRITE(95)REAL(TAUBSND(L),4)
              ENDIF
            ELSE
              ! *** TOTAL BED SHEAR STRESS
              WRITE(95)REAL(TAUB(L),4)
            ENDIF
          ELSE
            ! *** TOTAL BED SHEAR STRESS
            SHEAR=MAX(QQ(L,0),QQMIN)/CTURB2
            WRITE(95)SHEAR
          ENDIF
          IF(ISWAVE >= 1)THEN
            ! *** Shear due to Current Only
            SHEAR = (RSSBCE(L)*TBX(L+1)    + RSSBCW(L)*TBX(L))**2  + 
     &              (RSSBCN(L)*TBY(LNC(L)) + RSSBCS(L)*TBY(L))**2  
            SHEAR = 0.5*SQRT(SHEAR)  
            WRITE(95)REAL(QQWV3(L),4)  ! *** Bed Shear due to Waves Only
            WRITE(95)SHEAR             ! *** Bed Shear due to Current Only
            IF(ISWAVE >= 3)THEN
              WRITE(95)REAL(WVWHA(L),4),REAL(WVFRQL(L),4),
     &                 REAL(WACCWE(L),4)
              IF(ISWAVE == 4)THEN
                ! ***         DISSIPATION        SXX               SYY               SXY   (M3/S2)
                WRITE(95)REAL(WVDISP(L,KC),4),REAL(WVHUU(L,KC),4),
     &                   REAL(WVHVV(L,KC),4), REAL(WVHUV(L,KC),4)
              ENDIF
            ENDIF
          ENDIF
          IF(ISTRAN(1) == 1)WRITE(95)(REAL(SAL(L,K),4),K=1,KC)
          IF(ISTRAN(2) == 1)THEN
            WRITE(95)(REAL(TEM(L,K),4),K=1,KC)
            IF(TBEDIT(2) > 0.)WRITE(95)REAL(TEMB(L,KBHM),4)
          ENDIF
          IF(ISTRAN(3) == 1)WRITE(95,ERR=999,IOSTAT=ISYS)
     &                            (REAL(DYE(L,K),4),K=1,KC)
          IF(ISTRAN(4) == 1)WRITE(95)(REAL(SFL(L,K),4),K=1,KC)
          IF(ISTRAN(5) == 1)THEN
            WRITE(95)(REAL(TOXB(L,N1,NT),4),NT=1,NTOX)
            WRITE(95)((REAL(TOX(L,K,NT),4),K=1,KC),NT=1,NTOX)
          ENDIF
          IF(ISTRAN(6) == 1 .OR. ISTRAN(7) >= 1)THEN
            WRITE(95)REAL(BELV(L),4),REAL(HBED(L,N1),4),
     &               REAL(BDENBED(L,N1),4),REAL(PORBED(L,N1),4)
            IF(ISTRAN(6) == 1)THEN
              WRITE(95) (REAL(SEDB(L,N1,NS),4),
     &                   REAL(VFRBED(L,N1,NS),4),NS=1,NSED)
              WRITE(95)((REAL(SED(L,K,NS),4),K=1,KC),NS=1,NSED)
            ENDIF
            IF(ISTRAN(7) == 1)THEN
              WRITE(95)(REAL(SNDB(L,N1,NX),4), 
     &                  REAL(VFRBED(L,N1,NX+NSED),4),NX=1,NSND)
              WRITE(95)((REAL(SND(L,K,NX),4),K=1,KC),NX=1,NSND)
              IF(ISBDLDBC > 0)THEN
                WRITE(95)(REAL(CQBEDLOADX(L,NX),4),
     &                    REAL(CQBEDLOADY(L,NX),4),NX=1,NSND)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        CALL FLUSH(95)
        CLOSE(95,STATUS='KEEP')
      ENDIF

      ! *** NOW OUTPUT ALL THE BEDINFO TO A SINGLE FILE
      IF(ISBEXP >= 1 .AND. KB > 1)THEN
        NBEDSTEPS=NBEDSTEPS+1
        IF((((ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1) .AND. KB > 1) .AND.
     &                  NBEDSTEPS >= ISBEXP) .OR. JSEXPLORER == 1)THEN
          OPEN(10,FILE='EE_BED.OUT',STATUS='UNKNOWN',POSITION='APPEND',
     &            FORM='BINARY')
          WRITE(10)REAL(EETIME,4),NACTIVE
          DO L=2,LA
            WRITE(10)INT(KBT(L),4)
          ENDDO
          DO L=2,LA
            DO K=1,KB   
              WRITE(10)REAL(HBED(L,K),4),REAL(BDENBED(L,K),4),
     &                 REAL(PORBED(L,K),4)
              IF(ISTRAN(6) >= 1)THEN
                DO NS=1,NSED
                  WRITE(10)REAL(SEDB(L,K,NS),4)
                ENDDO
              ENDIF
              IF(ISTRAN(7) >= 1)THEN
                DO NX=1,NSND
                  NS=NSED+NX
                  WRITE(10)REAL(SNDB(L,K,NX),4)
                ENDDO
              ENDIF
              IF(ISTRAN(5) >= 1)THEN
                DO NT=1,NTOX
                  WRITE(10)REAL(TOXB(L,K,NT),4)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          CALL FLUSH(10)
          CLOSE(10,STATUS='KEEP')
          NBEDSTEPS=0
        ENDIF
      ENDIF

      ! *** WATER QUALITY
      IF(ISTRAN(8) > 0)THEN
        !  1) CHC - cyanobacteria 
        !  2) CHD - diatom algae 
        !  3) CHG - green algae 
        !  4) ROC - refractory particulate organic carbon 
        !  5) LOC - labile particulate organic carbon 
        !  6) DOC - dissolved organic carbon 
        !  7) ROP - refractory particulate organic phosphorus 
        !  8) LOP - labile particulate organic phosphorus 
        !  9) DOP - dissolved organic phosphorus 
        ! 10) P4D - total phosphate
        ! 11) RON - refractory particulate organic nitrogen 22) macroalgae
        ! 12) LON - labile particulate organic nitrogen
        ! 13) DON - dissolved organic nitrogen
        ! 14) NHX - ammonia nitrogen
        ! 15) NOX - nitrate nitrogen
        ! 16) SUU - particulate biogenic silica
        ! 17) SAA - dissolved available silica
        ! 18) COD - chemical oxygen demand
        ! 19) DOX - dissolved oxygen
        ! 20) TAM - total active metal
        ! 21) FCB - fecal coliform bacteria
        OPEN(95,FILE='EE_WQ.OUT',STATUS='UNKNOWN',POSITION='APPEND',
     &          FORM='BINARY')
        WRITE(95)REAL(EETIME,4)
        DO L=2,LA
          DO K=1,KC
            DO NW=1,NWQVAR
              IF(IWQ(NW) > 0)THEN
                WQ=WQV(L,K,NW)
                WRITE(95)WQ
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CALL FLUSH(95)
        CLOSE(95,STATUS='KEEP')

        ! *** SAVE SEDIMENT DIAGENESIS RESULTS
        IF(IWQBEN > 0 .AND. ISSDBIN < 0)THEN
          ! *** IF JSEXPLORER=1 THEN WRITE THE ARRAYS (I.E. IC'S)
          NSEDSTEPS=NSEDSTEPS+1
          IF(NSEDSTEPS >= ABS(ISSDBIN) .OR. JSEXPLORER == 1)THEN
            OPEN(95,FILE='EE_SD.OUT',STATUS='UNKNOWN',POSITION='APPEND',
     &              FORM='BINARY')
            WRITE(95)REAL(EETIME,4)
            DO L=2,LA

              !   SMPON = Conc. Particulate Org. Nitrogen   in G-class 1, 2 & 3  (g/m3) dim(LA,NSMGM)
              !   SMPOP = Conc. Particulate Org. Phosphorus in G-class 1, 2 & 3  (g/m3) dim(LA,NSMGM)
              !   SMPOC = Conc. Particulate Org. Carbon     in G-class 1, 2 & 3  (g/m3) dim(LA,NSMGM)

              ! *** DEPOSITION FLUXES
              ! SMDFN(LL,?) = Sediment Flux To The Sediment Bed From PON Into G1, G2, & G3
              ! SMDFP(LL,?) = Sediment Flux To The Sediment Bed From POP Into G1, G2, & G3
              ! SMDFC(LL,?) = Sediment Flux To The Sediment Bed From POC Into G1, G2, & G3

              !  SM1NH4 = Conc. NH4-N in layer 1 (g/m3)  dim(LA)
              !  SM2NH4 = Conc. NH4-N in layer 2 (g/m3)
              !  SM1NO3 = Conc. NO3-N in layer 1 (g/m3)
              !  SM2NO3 = Conc. NO3-N in layer 2 (g/m3)
              !  SM1PO4 = Conc. PO4-P in layer 1 (g/m3)
              !  SM2PO4 = Conc. PO4-P in layer 2 (g/m3)
              !  SM1H2S = Conc. Sulfide (H2S) in layer 1 (g/m3)
              !  SM2H2S = Conc. Sulfide (H2S) in layer 2 (g/m3)
              !   SMPSI = Conc. Particulate biogenic silica in layer 2 (g/m3)
              !   SM1SI = Conc. Dissolved available silica in layer 1 (g/m3)
              !   SM2SI = Conc. Dissolved available silica in layer 2 (g/m3)
              !   SMBST = Accumulated benthic stress (days)
              !     SMT = Sediment temperature (degC)

              ! *** SEDIMENT OXYGEN DEMANDS
              !  SMCSOD = CARBONACEOUS SOD
              !  SMNSOD = NITROGENOUS SOD

              ! *** BENTHIC FLUXES
              !  WQBFNH4 = AMMONIUM FLUX
              !  WQBFNO3 = NITRATE FLUX
              !   WQBFO2 = O2 SEDIMENT FLUX (SOD)
              !  WQBFCOD = COD FLUX
              ! WQBFPO4D = PO4 FLUX
              !  WQBFSAD = SILICA FLUX

              WRITE(95)(REAL(SMPON(L,K),4),K=1,3)
              WRITE(95)(REAL(SMPOP(L,K),4),K=1,3)
              WRITE(95)(REAL(SMPOC(L,K),4),K=1,3)
              WRITE(95)(REAL(SMDFN(L,K),4),K=1,3)
              WRITE(95)(REAL(SMDFP(L,K),4),K=1,3)
              WRITE(95)(REAL(SMDFC(L,K),4),K=1,3)
              WRITE(95)REAL(SM1NH4(L),4),REAL(SM2NH4(L),4)
              WRITE(95)REAL(SM1NO3(L),4),REAL(SM2NO3(L),4)
              WRITE(95)REAL(SM1PO4(L),4),REAL(SM2PO4(L),4)
              WRITE(95)REAL(SM1H2S(L),4),REAL(SM2H2S(L),4)
              WRITE(95)REAL(SM1SI(L),4), REAL(SM2SI(L),4)
              WRITE(95)REAL(SMPSI(L),4)
              WRITE(95)REAL(SMBST(L),4),REAL(SMT(L),4)
              WRITE(95)REAL(SMCSOD(L),4),REAL(SMNSOD(L),4)
              WRITE(95)REAL(WQBFNH4(L),4),REAL(WQBFNO3(L),4),
     &                 REAL(WQBFO2(L),4)
              WRITE(95)REAL(WQBFCOD(L),4),REAL(WQBFPO4D(L),4),
     &                 REAL(WQBFSAD(L),4)
            ENDDO

            CALL FLUSH(95)
            CLOSE(95,STATUS='KEEP')
            NSEDSTEPS=0
          ENDIF
        ENDIF
    
      ENDIF

      ! *** INTERNAL ARRAYS
      IF(ISINWV == 2)THEN
        ZERO=0.0
        OPEN(95,FILE='EE_ARRAYS.OUT',STATUS='UNKNOWN', POSITION='APPEND'
     &         ,FORM='BINARY')

        ! *** TIME VARIABLE USER SPECIFIED ARRAYS IN THIS SECTION
        ITIMEVAR = 1
        ITYPE = 1
        WRITE(95)ITYPE,ITIMEVAR
        ARRAYNAME='AH'
        WRITE(95)ARRAYNAME
        DO K=1,KC
          DO L=2,LA
            WRITE(95)REAL(AH(L,K),4)
          ENDDO
        ENDDO

        ITYPE = 1
        WRITE(95)ITYPE,ITIMEVAR
        ARRAYNAME='AV'
        WRITE(95)ARRAYNAME
        DO K=1,KC
          DO L=2,LA
            TMPVAL=AV(L,K)*HP(L)
            WRITE(95)TMPVAL
          ENDDO
        ENDDO

        ITYPE = 2
        WRITE(95)ITYPE,ITIMEVAR
        ARRAYNAME='QQ'
        WRITE(95)ARRAYNAME
        DO K=0,KC
          DO L=2,LA
            WRITE(95)REAL(QQ(L,K),4)
          ENDDO
        ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FMDUX' 
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(FMDUX(L,k),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FMDUY' 
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(FMDUY(L,k),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FMDVX' 
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(FMDVX(L,k),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FMDVY' 
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(FMDVY(L,k),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FXVEG'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(SFL(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FYVEG'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(SFL2(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FXVEGE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  WRITE(95)REAL(FXVEGE(L),4)
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FYVEGE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  WRITE(95)REAL(FYVEGE(L),4)
        !ENDDO

        ! *** THE FOLLOWING ARE PROVIDED FOR REFERENCE IF THE USER WANTS DIRECT ACCESS TO THESE ARRAYS

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='UHDYE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  WRITE(95)REAL(UHDYE(L),4)
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='VHDXE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  WRITE(95)REAL(VHDXE(L),4)
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FXE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  TMPVAL=FXE(L)*DELT*DXIU(L)
        !  WRITE(95)TMPVAL
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FYE'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  TMPVAL=FYE(L)*DELT*DYIV(L)
        !  WRITE(95)TMPVAL
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FUHX'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  TMPVAL=AHC(L,1)*DELT*DXIU(L)
        !  WRITE(95)TMPVAL
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FVHY'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  TMPVAL=AHC(L,2)*DELT*DYIV(L)
        !  WRITE(95)TMPVAL
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FUHU'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(AHU(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FVHU'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(AMCU(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FVHV'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(AMCV(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 1
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='FUHV'
        !WRITE(95)ARRAYNAME
        !DO K=1,KC
        !  DO L=2,LA
        !    WRITE(95)REAL(AMSU(L,K),4)
        !  ENDDO
        !ENDDO

        !ITYPE = 0
        !WRITE(95)ITYPE,ITIMEVAR
        !ARRAYNAME='TATMT'
        !WRITE(95)ARRAYNAME
        !DO L=2,LA
        !  WRITE(95)TATMT(L)
        !ENDDO

        CALL FLUSH(95)
        CLOSE(95,STATUS='KEEP')

      ENDIF

      RETURN

  999 STOP ' Error writing SNAPSHOT file'

      END
