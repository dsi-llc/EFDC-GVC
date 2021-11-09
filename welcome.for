C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE WELCOME
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
C**********************************************************************C
C
      WRITE(6,100)
      WRITE(6,101)
      WRITE(6,102)
      WRITE(6,103)
      WRITE(6,104)
      WRITE(6,105)
      WRITE(6,106)
      WRITE(6,107)
      WRITE(6,108)
      WRITE(6,109)
      WRITE(6,110)
      WRITE(6,111)
      WRITE(6,112)
      WRITE(6,113)
      WRITE(6,114)
      WRITE(6,115)
      WRITE(6,116)
      WRITE(6,117)
      WRITE(6,118)
      WRITE(6,119)
      WRITE(6,120)
      WRITE(6,121)
      WRITE(6,122)
      WRITE(6,123)
      WRITE(6,124)
      WRITE(6,125)
      WRITE(6,126)
      WRITE(6,127)
      WRITE(6,128)
C
  100 FORMAT(/)
  101 FORMAT(5X,'***********************************',
     +'************************************')
  102 FORMAT(5X,'*                                  ',
     +'                                   *')
  103 FORMAT(5X,'*                                  ',
     +'                                   *')
  104 FORMAT(5X,'*            EEEEEEEEE    FFFFFFFFF',
     +'    DDDDDDDD       CCCCCCCC        *')
  105 FORMAT(5X,'*           EEE          FFF       ',
     +'   DDD     DD    CCC      CC       *')
  106 FORMAT(5X,'*          EEE          FFF        ',
     +'  DDD     DD    CCC                *')
  107 FORMAT(5X,'*         EEEEEEEE     FFFFFFFF    ',
     +' DDD     DD    CCC                 *')
  108 FORMAT(5X,'*        EEE          FFF          ',
     +'DDD     DD    CCC                  *')
  109 FORMAT(5X,'*       EEE          FFF          D',
     +'DD     DD    CCC      CC           *')
  110 FORMAT(5X,'*      EEEEEEEEE    FFF          DD',
     +'DDDDDDDD      CCCCCCCCC            *')
  111 FORMAT(5X,'*                                  ',
     +'                                   *')
  112 FORMAT(5X,'*                                  ',
     +'                                   *')
  113 FORMAT(5X,'*                                  ',
     +'                                   *')
  114 FORMAT(5X,'*   TTTTTTTTTTTT         ENVIRONMEN',
     +'TAL FLUID DYNAMICS CODE            *')
  115 FORMAT(5X,'*       TTT   TT                   ',
     +'                                   *')
  116 FORMAT(5X,'*       TTT   TT         DEVELOPED ',
     +'BY JOHN M. HAMRICK                 *')
  117 FORMAT(5X,'*       TTTTTTTTTTTT               ',
     +'                                   *')
  118 FORMAT(5X,'*       TTT   TT         MAINTAINED',
     +' AND SUPPORTED BY TETRA TECH,INC.  *')
  119 FORMAT(5X,'*       TTT   TT                   ',
     +'                                   *')
  120 FORMAT(5X,'*       TTT   TT  TT     RELEASE DA',
     +'TE: 14 JUNE 2007                   *')
  121 FORMAT(5X,'*       TTT    TTTT                ',
     +'                                   *')
  122 FORMAT(5X,'*                                  ',
     +'                                   *')
  123 FORMAT(5X,'*                                  ',
     +'                                   *')
  124 FORMAT(5X,'*     EFDC AND "ENVIRONMENTAL FLUID',
     +' DYNAMICS CODE" ARE TRADEMARKS     *')
  125 FORMAT(5X,'*     OF JOHN M. HAMRICK           ',
     +'                                   *')
  126 FORMAT(5X,'*                                  ',
     +'                                   *')
  127 FORMAT(5X,'***********************************',
     +'************************************')
  128 FORMAT(//)
C
      RETURN
      END
