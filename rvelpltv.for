C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RVELPLTV
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
C **  SUBROUTINE VELPLTV WRITES A FILE FOR VERTICAL PLANE CONTOURING
C **  OF VELOCITY NORMAL TO AN ARBITARY SEQUENCE OF (I,J) POINTS AND
C **  AND VERTICAL PLANE TANGENTIAL-VERTICAL VELOCITY VECTORS
C
C**********************************************************************C
C
      INCLUDE 'EFDC.PAR'
      INCLUDE 'EFDC.CMN'
C
C**********************************************************************C
C
      REAL RVELN (KCM,100)
      REAL RVELT (KCM,100)
      REAL RW (KCM,100)
      REAL PVELN (KCM,100)
      REAL PVELT (KCM,100)
      REAL PWX(KCM,100)
      REAL RLVELN (KCM,100)
      REAL RLVELT (KCM,100)
      REAL RLW (KCM,100)
      CHARACTER*80 TITLE10,TITLE20,TITLE30
      CHARACTER*80 TITLE40,TITLE50,TITLE60
      CHARACTER*80 TITLE70,TITLE80,TITLE90
C
C**********************************************************************C
C
      IF(JSRVPV.NE.1) GOTO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      TITLE10='NORMAL ERT VELOCITY CONTOURS'
      TITLE20='NORMAL VPT VELOCITY CONTOURS'
      TITLE30='NORMAL MMT VELOCITY CONTOURS'
      TITLE40='TANGENTIAL ERT VELOCITY CONTOURS'
      TITLE50='TANGENTIAL VPT VELOCITY CONTOURS'
      TITLE60='TANGENTIAL MMT VELOCITY CONTOURS'
      TITLE70='TANGENTIAL ERT VELOCITY VECTORS'
      TITLE80='TANGENTIAL VPT VELOCITY VECTORS'
      TITLE90='TANGENTIAL MMT VELOCITY VECTORS'
C
      LEVELS=KC
C
      IF(ISECVPV.GE.1)THEN
        OPEN(11,FILE='RVLCNV1.OUT')
        OPEN(21,FILE='PVLCNV1.OUT')
        OPEN(31,FILE='MVLCNV1.OUT')
        OPEN(41,FILE='RVLCVT1.OUT')
        OPEN(51,FILE='PVLCVT1.OUT')
        OPEN(61,FILE='MVLCVT1.OUT')
        OPEN(71,FILE='RVLVCV1.OUT')
        OPEN(81,FILE='PVLVCV1.OUT')
        OPEN(91,FILE='MVLVCV1.OUT')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        CLOSE(31,STATUS='DELETE')
        CLOSE(41,STATUS='DELETE')
        CLOSE(51,STATUS='DELETE')
        CLOSE(61,STATUS='DELETE')
        CLOSE(71,STATUS='DELETE')
        CLOSE(81,STATUS='DELETE')
        CLOSE(91,STATUS='DELETE')
        OPEN(11,FILE='RVLCNV1.OUT')
        OPEN(21,FILE='PVLCNV1.OUT')
        OPEN(31,FILE='MVLCNV1.OUT')
        OPEN(41,FILE='RVLCVT1.OUT')
        OPEN(51,FILE='PVLCVT1.OUT')
        OPEN(61,FILE='MVLCVT1.OUT')
        OPEN(71,FILE='RVLVCV1.OUT')
        OPEN(81,FILE='PVLVCV1.OUT')
        OPEN(91,FILE='MVLVCV1.OUT')
      ENDIF 
      IF(ISECVPV.GE.2)THEN
        OPEN(12,FILE='RVLCNV2.OUT')
        OPEN(22,FILE='PVLCNV2.OUT')
        OPEN(32,FILE='MVLCNV2.OUT')
        OPEN(42,FILE='RVLCVT2.OUT')
        OPEN(52,FILE='PVLCVT2.OUT')
        OPEN(62,FILE='MVLCVT2.OUT')
        OPEN(72,FILE='RVLVCV2.OUT')
        OPEN(82,FILE='PVLVCV2.OUT')
        OPEN(92,FILE='MVLVCV2.OUT')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        CLOSE(32,STATUS='DELETE')
        CLOSE(42,STATUS='DELETE')
        CLOSE(52,STATUS='DELETE')
        CLOSE(62,STATUS='DELETE')
        CLOSE(72,STATUS='DELETE')
        CLOSE(82,STATUS='DELETE')
        CLOSE(92,STATUS='DELETE')
        OPEN(12,FILE='RVLCNV2.OUT')
        OPEN(22,FILE='PVLCNV2.OUT')
        OPEN(32,FILE='MVLCNV2.OUT')
        OPEN(42,FILE='RVLCVT2.OUT')
        OPEN(52,FILE='PVLCVT2.OUT')
        OPEN(62,FILE='MVLCVT2.OUT')
        OPEN(72,FILE='RVLVCV2.OUT')
        OPEN(82,FILE='PVLVCV2.OUT')
        OPEN(92,FILE='MVLVCV2.OUT')
      ENDIF
      IF(ISECVPV.GE.3)THEN
        OPEN(13,FILE='RVLCNV3.OUT')
        OPEN(23,FILE='PVLCNV3.OUT')
        OPEN(33,FILE='MVLCNV3.OUT')
        OPEN(43,FILE='RVLCVT3.OUT')
        OPEN(53,FILE='PVLCVT3.OUT')
        OPEN(63,FILE='MVLCVT3.OUT')
        OPEN(73,FILE='RVLVCV3.OUT')
        OPEN(83,FILE='PVLVCV3.OUT')
        OPEN(93,FILE='MVLVCV3.OUT')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        CLOSE(33,STATUS='DELETE')
        CLOSE(43,STATUS='DELETE')
        CLOSE(53,STATUS='DELETE')
        CLOSE(63,STATUS='DELETE')
        CLOSE(73,STATUS='DELETE')
        CLOSE(83,STATUS='DELETE')
        CLOSE(93,STATUS='DELETE')
        OPEN(13,FILE='RVLCNV3.OUT')
        OPEN(23,FILE='PVLCNV3.OUT')
        OPEN(33,FILE='MVLCNV3.OUT')
        OPEN(43,FILE='RVLCVT3.OUT')
        OPEN(53,FILE='PVLCVT3.OUT')
        OPEN(63,FILE='MVLCVT3.OUT')
        OPEN(73,FILE='RVLVCV3.OUT')
        OPEN(83,FILE='PVLVCV3.OUT')
        OPEN(93,FILE='MVLVCV3.OUT')
      ENDIF
      IF(ISECVPV.GE.4)THEN
        OPEN(14,FILE='RVLCNV4.OUT')
        OPEN(24,FILE='PVLCNV4.OUT')
        OPEN(34,FILE='MVLCNV4.OUT')
        OPEN(44,FILE='RVLCVT4.OUT')
        OPEN(54,FILE='PVLCVT4.OUT')
        OPEN(64,FILE='MVLCVT4.OUT')
        OPEN(74,FILE='RVLVCV4.OUT')
        OPEN(84,FILE='PVLVCV4.OUT')
        OPEN(94,FILE='MVLVCV4.OUT')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        CLOSE(34,STATUS='DELETE')
        CLOSE(44,STATUS='DELETE')
        CLOSE(54,STATUS='DELETE')
        CLOSE(64,STATUS='DELETE')
        CLOSE(74,STATUS='DELETE')
        CLOSE(84,STATUS='DELETE')
        CLOSE(94,STATUS='DELETE')
        OPEN(14,FILE='RVLCNV4.OUT')
        OPEN(24,FILE='PVLCNV4.OUT')
        OPEN(34,FILE='MVLCNV4.OUT')
        OPEN(44,FILE='RVLCVT4.OUT')
        OPEN(54,FILE='PVLCVT4.OUT')
        OPEN(64,FILE='MVLCVT4.OUT')
        OPEN(74,FILE='RVLVCV4.OUT')
        OPEN(84,FILE='PVLVCV4.OUT')
        OPEN(94,FILE='MVLVCV4.OUT')
      ENDIF 
      IF(ISECVPV.GE.5)THEN
        OPEN(15,FILE='RVLCNV5.OUT')
        OPEN(25,FILE='PVLCNV5.OUT')
        OPEN(35,FILE='MVLCNV5.OUT')
        OPEN(45,FILE='RVLCVT5.OUT')
        OPEN(55,FILE='PVLCVT5.OUT')
        OPEN(65,FILE='MVLCVT5.OUT')
        OPEN(75,FILE='RVLVCV5.OUT')
        OPEN(85,FILE='PVLVCV5.OUT')
        OPEN(95,FILE='MVLVCV5.OUT')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        CLOSE(35,STATUS='DELETE')
        CLOSE(45,STATUS='DELETE')
        CLOSE(55,STATUS='DELETE')
        CLOSE(65,STATUS='DELETE')
        CLOSE(75,STATUS='DELETE')
        CLOSE(85,STATUS='DELETE')
        CLOSE(95,STATUS='DELETE')
        OPEN(15,FILE='RVLCNV5.OUT')
        OPEN(25,FILE='PVLCNV5.OUT')
        OPEN(35,FILE='MVLCNV5.OUT')
        OPEN(45,FILE='RVLCVT5.OUT')
        OPEN(55,FILE='PVLCVT5.OUT')
        OPEN(65,FILE='MVLCVT5.OUT')
        OPEN(75,FILE='RVLVCV5.OUT')
        OPEN(85,FILE='PVLVCV5.OUT')
        OPEN(95,FILE='MVLVCV5.OUT')
      ENDIF
      IF(ISECVPV.GE.6)THEN
        OPEN(16,FILE='RVLCNV6.OUT')
        OPEN(26,FILE='PVLCNV6.OUT')
        OPEN(36,FILE='MVLCNV6.OUT')
        OPEN(46,FILE='RVLCVT6.OUT')
        OPEN(56,FILE='PVLCVT6.OUT')
        OPEN(66,FILE='MVLCVT6.OUT')
        OPEN(76,FILE='RVLVCV6.OUT')
        OPEN(86,FILE='PVLVCV6.OUT')
        OPEN(96,FILE='MVLVCV6.OUT')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        CLOSE(36,STATUS='DELETE')
        CLOSE(46,STATUS='DELETE')
        CLOSE(56,STATUS='DELETE')
        CLOSE(66,STATUS='DELETE')
        CLOSE(76,STATUS='DELETE')
        CLOSE(86,STATUS='DELETE')
        CLOSE(96,STATUS='DELETE')
        OPEN(16,FILE='RVLCNV6.OUT')
        OPEN(26,FILE='PVLCNV6.OUT')
        OPEN(36,FILE='MVLCNV6.OUT')
        OPEN(46,FILE='RVLCVT6.OUT')
        OPEN(56,FILE='PVLCVT6.OUT')
        OPEN(66,FILE='MVLCVT6.OUT')
        OPEN(76,FILE='RVLVCV6.OUT')
        OPEN(86,FILE='PVLVCV6.OUT')
        OPEN(96,FILE='MVLVCV6.OUT')
      ENDIF
      IF(ISECVPV.GE.7)THEN
        OPEN(17,FILE='RVLCNV7.OUT')
        OPEN(27,FILE='PVLCNV7.OUT')
        OPEN(37,FILE='MVLCNV7.OUT')
        OPEN(47,FILE='RVLCVT7.OUT')
        OPEN(57,FILE='PVLCVT7.OUT')
        OPEN(67,FILE='MVLCVT7.OUT')
        OPEN(77,FILE='RVLVCV7.OUT')
        OPEN(87,FILE='PVLVCV7.OUT')
        OPEN(97,FILE='MVLVCV7.OUT')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        CLOSE(37,STATUS='DELETE')
        CLOSE(47,STATUS='DELETE')
        CLOSE(57,STATUS='DELETE')
        CLOSE(67,STATUS='DELETE')
        CLOSE(77,STATUS='DELETE')
        CLOSE(87,STATUS='DELETE')
        CLOSE(97,STATUS='DELETE')
        OPEN(17,FILE='RVLCNV7.OUT')
        OPEN(27,FILE='PVLCNV7.OUT')
        OPEN(37,FILE='MVLCNV7.OUT')
        OPEN(47,FILE='RVLCVT7.OUT')
        OPEN(57,FILE='PVLCVT7.OUT')
        OPEN(67,FILE='MVLCVT7.OUT')
        OPEN(77,FILE='RVLVCV7.OUT')
        OPEN(87,FILE='PVLVCV7.OUT')
        OPEN(97,FILE='MVLVCV7.OUT')
      ENDIF 
      IF(ISECVPV.GE.8)THEN
        OPEN(18,FILE='RVLCNV8.OUT')
        OPEN(28,FILE='PVLCNV8.OUT')
        OPEN(38,FILE='MVLCNV8.OUT')
        OPEN(48,FILE='RVLCVT8.OUT')
        OPEN(58,FILE='PVLCVT8.OUT')
        OPEN(68,FILE='MVLCVT8.OUT')
        OPEN(78,FILE='RVLVCV8.OUT')
        OPEN(88,FILE='PVLVCV8.OUT')
        OPEN(98,FILE='MVLVCV8.OUT')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        CLOSE(38,STATUS='DELETE')
        CLOSE(48,STATUS='DELETE')
        CLOSE(58,STATUS='DELETE')
        CLOSE(68,STATUS='DELETE')
        CLOSE(78,STATUS='DELETE')
        CLOSE(88,STATUS='DELETE')
        CLOSE(98,STATUS='DELETE')
        OPEN(18,FILE='RVLCNV8.OUT')
        OPEN(28,FILE='PVLCNV8.OUT')
        OPEN(38,FILE='MVLCNV8.OUT')
        OPEN(48,FILE='RVLCVT8.OUT')
        OPEN(58,FILE='PVLCVT8.OUT')
        OPEN(68,FILE='MVLCVT8.OUT')
        OPEN(78,FILE='RVLVCV8.OUT')
        OPEN(88,FILE='PVLVCV8.OUT')
        OPEN(98,FILE='MVLVCV8.OUT')
      ENDIF
      IF(ISECVPV.GE.9)THEN
        OPEN(19,FILE='RVLCNV9.OUT')
        OPEN(29,FILE='PVLCNV9.OUT')
        OPEN(39,FILE='MVLCNV9.OUT')
        OPEN(49,FILE='RVLCVT9.OUT')
        OPEN(59,FILE='PVLCVT9.OUT')
        OPEN(69,FILE='MVLCVT9.OUT')
        OPEN(79,FILE='RVLVCV9.OUT')
        OPEN(89,FILE='PVLVCV9.OUT')
        OPEN(99,FILE='MVLVCV9.OUT')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        CLOSE(39,STATUS='DELETE')
        CLOSE(49,STATUS='DELETE')
        CLOSE(59,STATUS='DELETE')
        CLOSE(69,STATUS='DELETE')
        CLOSE(79,STATUS='DELETE')
        CLOSE(89,STATUS='DELETE')
        CLOSE(99,STATUS='DELETE')
        OPEN(19,FILE='RVLCNV9.OUT')
        OPEN(29,FILE='PVLCNV9.OUT')
        OPEN(39,FILE='MVLCNV9.OUT')
        OPEN(49,FILE='RVLCVT9.OUT')
        OPEN(59,FILE='PVLCVT9.OUT')
        OPEN(69,FILE='MVLCVT9.OUT')
        OPEN(79,FILE='RVLVCV9.OUT')
        OPEN(89,FILE='PVLVCV9.OUT')
        OPEN(99,FILE='MVLVCV9.OUT')
      ENDIF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LUN4=40+IS
      LUN5=50+IS
      LUN6=60+IS
      LUN7=70+IS
      LUN8=80+IS
      LUN9=90+IS
      LINES=NIJVPV(IS)
      WRITE (LUN1,99) TITLE10,CVTITLE(LUN1)
      WRITE (LUN2,99) TITLE20,CVTITLE(LUN2)
      WRITE (LUN3,99) TITLE30,CVTITLE(LUN3)
      WRITE (LUN4,99) TITLE40,CVTITLE(LUN4)
      WRITE (LUN5,99) TITLE50,CVTITLE(LUN5)
      WRITE (LUN6,99) TITLE60,CVTITLE(LUN6)
      WRITE (LUN7,99) TITLE70,CVTITLE(LUN7)
      WRITE (LUN8,99) TITLE80,CVTITLE(LUN8)
      WRITE (LUN9,99) TITLE90,CVTITLE(LUN9)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN3,101)LINES,LEVELS
      WRITE (LUN4,101)LINES,LEVELS
      WRITE (LUN5,101)LINES,LEVELS
      WRITE (LUN6,101)LINES,LEVELS
      WRITE (LUN7,101)LINES,LEVELS
      WRITE (LUN8,101)LINES,LEVELS
      WRITE (LUN9,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      WRITE (LUN3,250)(ZZ(K),K=1,KC)
      WRITE (LUN4,250)(ZZ(K),K=1,KC)
      WRITE (LUN5,250)(ZZ(K),K=1,KC)
      WRITE (LUN6,250)(ZZ(K),K=1,KC)
      WRITE (LUN7,250)(ZZ(K),K=1,KC)
      WRITE (LUN8,250)(ZZ(K),K=1,KC)
      WRITE (LUN9,250)(ZZ(K),K=1,KC)
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      CLOSE(LUN4)
      CLOSE(LUN5)
      CLOSE(LUN6)
      CLOSE(LUN7)
      CLOSE(LUN8)
      CLOSE(LUN9)
      ENDDO
C
      JSRVPV=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
      ELSE
        TIME=TIMESEC/TCON
      ENDIF
C
      IF(ISECVPV.GE.1)THEN
        OPEN(11,FILE='RVLCNV1.OUT',POSITION='APPEND')
        OPEN(21,FILE='PVLCNV1.OUT',POSITION='APPEND')
        OPEN(31,FILE='MVLCNV1.OUT',POSITION='APPEND')
        OPEN(41,FILE='RVLCVT1.OUT',POSITION='APPEND')
        OPEN(51,FILE='PVLCVT1.OUT',POSITION='APPEND')
        OPEN(61,FILE='MVLCVT1.OUT',POSITION='APPEND')
        OPEN(71,FILE='RVLVCV1.OUT',POSITION='APPEND')
        OPEN(81,FILE='PVLVCV1.OUT',POSITION='APPEND')
        OPEN(91,FILE='MVLVCV1.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.2)THEN
        OPEN(12,FILE='RVLCNV2.OUT',POSITION='APPEND')
        OPEN(22,FILE='PVLCNV2.OUT',POSITION='APPEND')
        OPEN(32,FILE='MVLCNV2.OUT',POSITION='APPEND')
        OPEN(42,FILE='RVLCVT2.OUT',POSITION='APPEND')
        OPEN(52,FILE='PVLCVT2.OUT',POSITION='APPEND')
        OPEN(62,FILE='MVLCVT2.OUT',POSITION='APPEND')
        OPEN(72,FILE='RVLVCV2.OUT',POSITION='APPEND')
        OPEN(82,FILE='PVLVCV2.OUT',POSITION='APPEND')
        OPEN(92,FILE='MVLVCV2.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.3)THEN
        OPEN(13,FILE='RVLCNV3.OUT',POSITION='APPEND')
        OPEN(23,FILE='PVLCNV3.OUT',POSITION='APPEND')
        OPEN(33,FILE='MVLCNV3.OUT',POSITION='APPEND')
        OPEN(43,FILE='RVLCVT3.OUT',POSITION='APPEND')
        OPEN(53,FILE='PVLCVT3.OUT',POSITION='APPEND')
        OPEN(63,FILE='MVLCVT3.OUT',POSITION='APPEND')
        OPEN(73,FILE='RVLVCV3.OUT',POSITION='APPEND')
        OPEN(83,FILE='PVLVCV3.OUT',POSITION='APPEND')
        OPEN(93,FILE='MVLVCV3.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.4)THEN
        OPEN(14,FILE='RVLCNV4.OUT',POSITION='APPEND')
        OPEN(24,FILE='PVLCNV4.OUT',POSITION='APPEND')
        OPEN(34,FILE='MVLCNV4.OUT',POSITION='APPEND')
        OPEN(44,FILE='RVLCVT4.OUT',POSITION='APPEND')
        OPEN(54,FILE='PVLCVT4.OUT',POSITION='APPEND')
        OPEN(64,FILE='MVLCVT4.OUT',POSITION='APPEND')
        OPEN(74,FILE='RVLVCV4.OUT',POSITION='APPEND')
        OPEN(84,FILE='PVLVCV4.OUT',POSITION='APPEND')
        OPEN(94,FILE='MVLVCV4.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.5)THEN
        OPEN(15,FILE='RVLCNV5.OUT',POSITION='APPEND')
        OPEN(25,FILE='PVLCNV5.OUT',POSITION='APPEND')
        OPEN(35,FILE='MVLCNV5.OUT',POSITION='APPEND')
        OPEN(45,FILE='RVLCVT5.OUT',POSITION='APPEND')
        OPEN(55,FILE='PVLCVT5.OUT',POSITION='APPEND')
        OPEN(65,FILE='MVLCVT5.OUT',POSITION='APPEND')
        OPEN(75,FILE='RVLVCV5.OUT',POSITION='APPEND')
        OPEN(85,FILE='PVLVCV5.OUT',POSITION='APPEND')
        OPEN(95,FILE='MVLVCV5.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.6)THEN
        OPEN(16,FILE='RVLCNV6.OUT',POSITION='APPEND')
        OPEN(26,FILE='PVLCNV6.OUT',POSITION='APPEND')
        OPEN(36,FILE='MVLCNV6.OUT',POSITION='APPEND')
        OPEN(46,FILE='RVLCVT6.OUT',POSITION='APPEND')
        OPEN(56,FILE='PVLCVT6.OUT',POSITION='APPEND')
        OPEN(66,FILE='MVLCVT6.OUT',POSITION='APPEND')
        OPEN(76,FILE='RVLVCV6.OUT',POSITION='APPEND')
        OPEN(86,FILE='PVLVCV6.OUT',POSITION='APPEND')
        OPEN(96,FILE='MVLVCV6.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.7)THEN
        OPEN(17,FILE='RVLCNV7.OUT',POSITION='APPEND')
        OPEN(27,FILE='PVLCNV7.OUT',POSITION='APPEND')
        OPEN(37,FILE='MVLCNV7.OUT',POSITION='APPEND')
        OPEN(47,FILE='RVLCVT7.OUT',POSITION='APPEND')
        OPEN(57,FILE='PVLCVT7.OUT',POSITION='APPEND')
        OPEN(67,FILE='MVLCVT7.OUT',POSITION='APPEND')
        OPEN(77,FILE='RVLVCV7.OUT',POSITION='APPEND')
        OPEN(87,FILE='PVLVCV7.OUT',POSITION='APPEND')
        OPEN(97,FILE='MVLVCV7.OUT',POSITION='APPEND')
      ENDIF 
      IF(ISECVPV.GE.8)THEN
        OPEN(18,FILE='RVLCNV8.OUT',POSITION='APPEND')
        OPEN(28,FILE='PVLCNV8.OUT',POSITION='APPEND')
        OPEN(38,FILE='MVLCNV8.OUT',POSITION='APPEND')
        OPEN(48,FILE='RVLCVT8.OUT',POSITION='APPEND')
        OPEN(58,FILE='PVLCVT8.OUT',POSITION='APPEND')
        OPEN(68,FILE='MVLCVT8.OUT',POSITION='APPEND')
        OPEN(78,FILE='RVLVCV8.OUT',POSITION='APPEND')
        OPEN(88,FILE='PVLVCV8.OUT',POSITION='APPEND')
        OPEN(98,FILE='MVLVCV8.OUT',POSITION='APPEND')
      ENDIF
      IF(ISECVPV.GE.9)THEN
        OPEN(19,FILE='RVLCNV9.OUT',POSITION='APPEND')
        OPEN(29,FILE='PVLCNV9.OUT',POSITION='APPEND')
        OPEN(39,FILE='MVLCNV9.OUT',POSITION='APPEND')
        OPEN(49,FILE='RVLCVT9.OUT',POSITION='APPEND')
        OPEN(59,FILE='PVLCVT9.OUT',POSITION='APPEND')
        OPEN(69,FILE='MVLCVT9.OUT',POSITION='APPEND')
        OPEN(79,FILE='RVLVCV9.OUT',POSITION='APPEND')
        OPEN(89,FILE='PVLVCV9.OUT',POSITION='APPEND')
        OPEN(99,FILE='MVLVCV9.OUT',POSITION='APPEND')
      ENDIF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LUN4=40+IS
      LUN5=50+IS
      LUN6=60+IS
      LUN7=70+IS
      LUN8=80+IS
      LUN9=90+IS
      WRITE (LUN1,100)N,TIME
      WRITE (LUN2,100)N,TIME
      WRITE (LUN3,100)N,TIME
      WRITE (LUN4,100)N,TIME
      WRITE (LUN5,100)N,TIME
      WRITE (LUN6,100)N,TIME
      WRITE (LUN7,100)N,TIME
      WRITE (LUN8,100)N,TIME
      WRITE (LUN9,100)N,TIME
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       LN=LNC(L)
       LS=LSC(L)
        DO K=1,KC
        RVELN(K,NN)=50.*((UHLPF(L+1,K)+UHLPF(L,K))*COSC
     &            +(VHLPF(LN,K)+VHLPF(L,K))*SINC)/HLPF(L)
        RVELT(K,NN)=-50.*((UHLPF(L+1,K)+UHLPF(L,K))*SINC
     &            -(VHLPF(LN,K)+VHLPF(L,K))*COSC)/HLPF(L)
        RW(K,NN)=50.*(WLPF(L,K)+WLPF(L,K-1))
        PVELN(K,NN)=50.*((UVPT(L+1,K)+UVPT(L,K))*COSC
     &            +(VVPT(LN,K)+VVPT(L,K))*SINC)/HLPF(L)
        PVELT(K,NN)=-50.*((UVPT(L+1,K)+UVPT(L,K))*SINC
     &            -(VVPT(LN,K)+VVPT(L,K))*COSC)/HLPF(L)
        PWX(K,NN)=50.*(WVPT(L,K)+WVPT(L,K-1))
        RLVELN(K,NN)=RVELN(K,NN)+PVELN(K,NN)
        RLVELT(K,NN)=RVELT(K,NN)+PVELT(K,NN)
        RLW(K,NN)=RW(K,NN)+PWX(K,NN)
        ENDDO
       ENDDO
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN6,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN9,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN6,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN9,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(RVELN(K,NN),K=1,KC)
       WRITE(LUN2,250)(PVELN(K,NN),K=1,KC)
       WRITE(LUN3,250)(RLVELN(K,NN),K=1,KC)
       WRITE(LUN4,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN5,250)(PVELT(K,NN),K=1,KC)
       WRITE(LUN6,250)(RLVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN8,250)(PVELT(K,NN),K=1,KC)
       WRITE(LUN9,250)(RLVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RW(K,NN),K=1,KC)
       WRITE(LUN8,250)(PWX(K,NN),K=1,KC)
       WRITE(LUN9,250)(RLW(K,NN),K=1,KC)
       ENDDO
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      CLOSE(LUN4)
      CLOSE(LUN5)
      CLOSE(LUN6)
      CLOSE(LUN7)
      CLOSE(LUN8)
      CLOSE(LUN9)
      ENDDO
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
CMRM  200 FORMAT(2I5,1X,1P,6E13.5) 
CMRM  250 FORMAT(1P,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
