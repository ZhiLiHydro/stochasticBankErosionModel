!                    *****************
                     SUBROUTINE CONDIN
!                    *****************
!
!
!***********************************************************************
! TELEMAC2D   V6P3                                   21/08/2010
!***********************************************************************
!
!brief    INITIALISES THE PHYSICAL PARAMETERS H, U, V ETC.
!
!history  J-M HERVOUET (LNHE)
!+        30/08/2007
!+        V6P0
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  M.S.TURNBULL (HRW), N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        06/12/2011
!+        V6P2
!+   Addition of the Tsunami displacement (based on Okada's model)
!+   by calling CONDI_OKADA and of the TPXO tidal model by calling
!+   CONDI_TPXO (the TPXO model being coded in module TPXO)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE TPXO
      USE OKADA
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER ITRAC
      DOUBLE PRECISION BID
      LOGICAL OK
!
!-----------------------------------------------------------------------
!
!   INITIALISES THE TIME
!
      AT = 0.D0
      OK = .TRUE.
!
!-----------------------------------------------------------------------
!
!   INITIALISES THE VELOCITIES: ZERO VELOCITIES
!
      CALL OS('X=0     ',X=U)
      CALL OS('X=0     ',X=V)
!=======================================================================
!    VARIABLE 23
!=======================================================================
      CALL OS( 'X=C     ',X=PRIVE%ADR(1)%P,C=0.D0 )
      CALL FIND_IN_SEL(PRIVE%ADR(1)%P,'VARIABLE 23     ',
     &     T2D_FILES(T2DGEO)%LU,T2D_FILES(T2DGEO)%FMT,W,OK,TIME=BID)
!=======================================================================
!=======================================================================
!    VARIABLE 24
!=======================================================================
      CALL OS( 'X=C     ',X=PRIVE%ADR(2)%P,C=0.D0 )
      CALL FIND_IN_SEL(PRIVE%ADR(2)%P,'VARIABLE 24     ',
     &     T2D_FILES(T2DGEO)%LU,T2D_FILES(T2DGEO)%FMT,W,OK,TIME=BID)
!=======================================================================
!
!-----------------------------------------------------------------------
!
!   INITIALISES THE WATER DEPTH H
!
      IF(CDTINI(1:10).EQ.'COTE NULLE'.OR.
     &   CDTINI(1:14).EQ.'ZERO ELEVATION') THEN
        CALL OS( 'X=0     ' , X=H )
        CALL OS( 'X=X-Y   ' , X=H , Y=ZF )
      ELSEIF(CDTINI(1:14).EQ.'COTE CONSTANTE'.OR.
     &       CDTINI(1:18).EQ.'CONSTANT ELEVATION') THEN
        CALL OS( 'X=C     ' , H , H  , H , COTINI )
        CALL OS( 'X=X-Y   ' , H , ZF , H , 0.D0   )
      ELSEIF(CDTINI(1:13).EQ.'HAUTEUR NULLE'.OR.
     &       CDTINI(1:10).EQ.'ZERO DEPTH') THEN
        CALL OS( 'X=C     ' , H , H  , H , 0.D0  )
      ELSEIF(CDTINI(1:17).EQ.'HAUTEUR CONSTANTE'.OR.
     &       CDTINI(1:14).EQ.'CONSTANT DEPTH') THEN
        CALL OS( 'X=C     ' , H , H  , H , HAUTIN )
      ELSEIF(CDTINI(1:25).EQ.'ALTIMETRIE SATELLITE TPXO'.OR.
     &       CDTINI(1:24).EQ.'TPXO SATELLITE ALTIMETRY') THEN
        CALL OS('X=-Y    ',X=H,Y=ZF)
        CALL CONDI_TPXO(NPOIN,MESH%NPTFR,MESH%NBOR%I,
     &                  X,Y,H%R,U%R,V%R,
     &                  LIHBOR%I,LIUBOR%I,KENT,KENTU,
     &                  GEOSYST,NUMZONE,LAMBD0,PHI0,
     &                  T2D_FILES,T2DBB1,T2DBB2,
     &                  MARDAT,MARTIM,INTMICON,MSL)
      ELSEIF(CDTINI(1:13).EQ.'PARTICULIERES'.OR.
     &       CDTINI(1:10).EQ.'PARTICULAR'.OR.
     &       CDTINI(1:07).EQ.'SPECIAL') THEN
!
!
!

        IF(LNG.EQ.1) WRITE(LU,10)
        IF(LNG.EQ.2) WRITE(LU,11)
10      FORMAT(1X,'CONDIN : AVEC DES CONDITIONS INITIALES PARTICULIERES'
     &         ,/,'         VOUS DEVEZ MODIFIER CONDIN')
11      FORMAT(1X,'CONDIN : WITH SPECIAL INITIAL CONDITIONS'
     &         ,/,'         YOU HAVE TO MODIFY CONDIN')
        CALL PLANTE(1)
        STOP
!
!  END OF CODE TO BE MODIFIED BY USER
!
      ELSE
        IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'CONDIN : CONDITION INITIALE NON PREVUE : ',CDTINI
        ENDIF
        IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'CONDIN: INITIAL CONDITION UNKNOWN: ',CDTINI
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
!   INITIALISES TSUNAMI DISPLACEMENT
!
      IF(OPTTSUNAMI.EQ.1) THEN
        CALL CONDI_OKADA(NPOIN,X,Y,H%R,COETSUNAMI,LAMBD0,PHI0)
      ENDIF
!
!-----------------------------------------------------------------------
!
!   INITIALISES THE TRACERS
!
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC
          CALL OS('X=C     ',X=T%ADR(ITRAC)%P,C=TRAC0(ITRAC))
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
! INITIALISES THE VISCOSITY
!
      CALL OS('X=C     ',X=VISC,C=PROPNU)
!
!-----------------------------------------------------------------------
!
      RETURN
      END





!                    *****************
                     SUBROUTINE NOEROD
!                    *****************
!
     & (H , ZF , ZR , Z , X , Y , NPOIN , CHOIX , NLISS )
!
!***********************************************************************
! SISYPHE   V6P1                                   21/07/2011
!***********************************************************************
!
!brief    FIXES THE NON-ERODABLE BED ELEVATION ZR.
!
!note     METHODS OF TREATMENT OF NON-ERODABLE BEDS CAN LEAD TO ZF.
!note  CHOOSE TO SMOOTH THE SOLUTION WITH NLISS > 0.
!
!history  C. LENORMANT
!+
!+        V5P1
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CHOIX          |-->| SELECTED METHOD FOR THE TREATMENT OF RIGID BEDS
!| H              |-->| WATER DEPTH
!| NLISS          |<->| NUMBER OF SMOOTHINGS
!| NPOIN          |-->| NUMBER OF 2D POINTS
!| X,Y            |-->| 2D COORDINATES
!| Z              |-->| FREE SURFACE
!| ZF             |-->| BED LEVEL
!| ZR             |<--| RIGID BED LEVEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D, ONLY:T2D_FILES, T2DGEO, NPRIV,W
      USE DECLARATIONS_SISYPHE , ONLY : IELMT, MESH,PRIVE
  
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN):: NPOIN , CHOIX
      INTEGER, INTENT(INOUT):: NLISS
!
      DOUBLE PRECISION, INTENT(IN)::  Z(NPOIN) , ZF(NPOIN)
      DOUBLE PRECISION , INTENT(IN)::  X(NPOIN) , Y(NPOIN), H(NPOIN)
      DOUBLE PRECISION , INTENT(INOUT)::  ZR(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
!      INTEGER ERR
      DOUBLE PRECISION BID
!      REAL, ALLOCATABLE :: WSEB(:)
      LOGICAL OK
      OK = .TRUE.
!=======================================================================
!    VARIABLE 23
!=======================================================================
      CALL OS( 'X=C     ',X=PRIVE%ADR(1)%P,C=0.D0 )
      CALL FIND_IN_SEL(PRIVE%ADR(1)%P,'VARIABLE 23     ',
     &     T2D_FILES(T2DGEO)%LU,T2D_FILES(T2DGEO)%FMT,W,OK,TIME=BID)
!=======================================================================
!=======================================================================
!    VARIABLE 24
!=======================================================================
      CALL OS( 'X=C     ',X=PRIVE%ADR(2)%P,C=0.D0 )
      CALL FIND_IN_SEL(PRIVE%ADR(2)%P,'VARIABLE 24     ',
     &     T2D_FILES(T2DGEO)%LU,T2D_FILES(T2DGEO)%FMT,W,OK,TIME=BID)
!=======================================================================
!
! LOAD DATA FROM SELAFIN FILE
!
!      ALLOCATE(WSEB(NPOIN),STAT=ERR)
!      IF(ERR.NE.0) STOP
!
!      IF( NPRIV.GE.1 ) THEN
!       CALL OS( 'X=C     ',X=PRIVE%ADR(1)%P,C=0.D0 )
!       CALL OS( 'X=C     ',PRIVE%ADR(1)%P, ZF,ZF,(0.D0) )
!       CALL FIND_IN_SEL(PRIVE%ADR(1)%P,'VARIABLE 23     ',
!       CALL FIND_IN_SEL(PRIVE%ADR(1)%P,'BOTTOM          ',
!     &     T2D_FILES(T2DGEO)%LU,WSEB,OK,TIME=BID)
!      ENDIF
!      DEALLOCATE(WSEB)
!      IF( OK ) THEN
!         WRITE(LU,'(//,1X,A)') 'SUCCESS'
!      ELSE
!         WRITE(LU,'(//,1X,A)') 'FAIL'
!      ENDIF
!
!--------------------
! RIGID BEDS POSITION
!---------------------
!
!     DEFAULT VALUE:       ZR=ZF-100 EVERYWHERE
!
      IF(NPRIV.GE.1) THEN
       DO I = 1, NPOIN
        IF (PRIVE%ADR(1)%P%R(I).LT.1.D0) THEN !RIGID IF GT.0
         ZR(I)=-100.D0
        ELSE
         ZR(I)=ZF(I)
        ENDIF
       ENDDO
      ENDIF
!      CALL OV( 'X=C       ',ZR, ZF, ZF, -100.D0, NPOIN )
!
!------------------
! SMOOTHING OPTION
!------------------
!
!     NLISS : NUMBER OF SMOOTHING IF  (ZF - ZR ) NEGATIVE
!             DEFAULT VALUE : NLISS = 0 (NO SMOOTHING)
!
      NLISS = 0
!
!-----------------------------------------------------------------------
!
      RETURN
      END



!                    *****************
                     SUBROUTINE CORSTR
!                    *****************
!
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    CORRECTS THE FRICTION COEFFICIENT ON THE BOTTOM
!+                WHEN IT IS VARIABLE IN TIME.
!
!warning  USER SUBROUTINE; MUST BE CODED BY THE USER; THIS IS MERELY AN EXAMPLE
!code
!+2D   DO I = 1 , NPOIN
!+2D     IF(AT.GT.1200.D0) THEN
!+2D       CHESTR%R(I) = 40.D0
!+2D     ELSE
!+2D       CHESTR%R(I) = 60.D0
!+2D     ENDIF
!+2D   ENDDO
!
!history  J-M HERVOUET (LNHE)
!+        17/08/1994
!+        V5P6
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
!     C2D: EXAMPLE FOR TELEMAC-2D
!     C3D: EXAMPLE FOR TELEMAC-3D
!
      USE DECLARATIONS_TELEMAC2D
!3D   USE DECLARATIONS_TELEMAC3D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      INTEGER I
      DOUBLE PRECISION ZL,NLOW,NLRG,XMIN,XMAX
      LOGICAL OK
      DOUBLE PRECISION BID
!      ZL=106.9D0
      NLOW=0.03D0
      NLRG=0.2D0
!      XMIN=-2.2D0
!      XMAX=4.D0

!3D   INTEGER I
!
!-----------------------------------------------------------------------
!
!      DO I = 1 , NPOIN
!        ZL = 0.9D0*ZL0+(X(I)-XMIN)*ZL0/(XMAX-XMIN)
!        ZL = ZL0
!        IF(ZF%R(I).GT.ZL) THEN
!          CHESTR%R(I) = NLRG
!        ELSE
!          CHESTR%R(I) = NLOW
!        ENDIF
!      ENDDO



!=======================================================================
!    VARIABLE 23
!      OK = .TRUE.
!=======================================================================
!      CALL OS( 'X=C     ',X=PRIVE%ADR(1)%P,C=0.D0 )
!      CALL FIND_IN_SEL(PRIVE%ADR(1)%P,'VARIABLE 23     ',
!     &     T2D_FILES(T2DGEO)%LU,T2D_FILES(T2DGEO)%FMT,W,OK,TIME=BID)
!=======================================================================

      DO I=1,NPOIN
        IF (PRIVE%ADR(1)%P%R(I).LT.1.D0) THEN
          CHESTR%R(I) = NLOW
        ELSE
          CHESTR%R(I) = NLRG
        ENDIF
      ENDDO



!      IF(NCSIZE.GT.1) THEN
!        CALL PARCOM(PRIVE%ADR(1)%P,4,MESH)
!      ENDIF

!
!-----------------------------------------------------------------------
!
!3D   DO I = 1 , NPOIN2
!3D     IF(AT.GT.1200.D0) THEN
!3D       RUGOF%R(I) = 40.D0
!3D     ELSE
!3D       RUGOF%R(I) = 60.D0
!3D     ENDIF
!3D   ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END








!                    *******************
                     SUBROUTINE MAXSLOPE
!                    *******************
!
     &(SLOPE,ZF,ZR,XEL,YEL,NELEM,NELMAX,NPOIN,IKLE,EVOL,UNSV2D,MESH,
     & ZFCL_MS,AVAIL,NOMBLAY,NSICLA)
!
!***********************************************************************
! SISYPHE   V6P3                                   21/07/2011
!***********************************************************************
!
!brief    COLLAPSE OF SAND WITH A SLOPE GREATER THAN A
!+                STABILITY CRITERION.
!+        For more explanation see release notes 5.8
!
!history  J-M HERVOUET (LNH)
!+        16/11/2007
!+        V5P8
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        08/03/2013
!+        V6P3
!+   Now possible with several classes of sediment.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| EVOL           |<->| WORK ARRAY, THEN EVOLUTION DUE TO SLIDE
!| IKLE           |-->| CONNECTIVITY TABLE
!| MESH           |-->| MESH STRUCTURE
!| NELEM          |-->| NUMBER OF ELEMENTS
!| NELMAX      -   |-->| MAXIMUM NUMBER OF ELEMENTS
!| NPOIN          |-->| NUMBER OF POINTS IN THE MESH
!| SLOPE          |-->| MAXIMUM SLOPE IN DEGREES
!| UNSV2D         |-->| INVERSE OF INTEGRAL OF BASES
!| XEL,YEL        |-->| MESH COORDINATES PER ELEMENT
!| ZF             |<->| BOTTOM THAT WILL BE MODIFIED
!| ZR             |-->| NON ERODABLE BED
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
      USE INTERFACE_SISYPHE, EX_MAXSLOPE => MAXSLOPE
      USE DECLARATIONS_TELEMAC2D, ONLY : PRIVE,AT,DT,LISPRD,LT,
     &                                   T2D_FILES,T2DGEO,W
      USE DECLARATIONS_SISYPHE,   ONLY : MOFAC,X,Y,TOB
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NELEM,NELMAX,NPOIN,NOMBLAY,NSICLA
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3)
!
      DOUBLE PRECISION, INTENT(IN   ) :: SLOPE
      DOUBLE PRECISION, INTENT(INOUT) :: ZF(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: ZR(NPOIN)

      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,3),YEL(NELMAX,3)
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,NOMBLAY,NSICLA)
!
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: EVOL,ZFCL_MS
      TYPE(BIEF_OBJ), INTENT(IN)      :: UNSV2D
      TYPE(BIEF_MESH) :: MESH
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IELEM,I1,I2,I3,I,IG1,IG2,IR1,IR2,J
      DOUBLE PRECISION X2,X3,Y2,Y3,Z2,Z3,A,B,L,ZC,DEUXSURF,TANSL
      DOUBLE PRECISION Q(3),QG1,QG2,QR1,QR2
!
      LOGICAL CASE2
!
      INTRINSIC SQRT,SIN,COS,TAN,ATAN,ERF,ERFC,MIN,MAX,MINVAL,MAXVAL

      INTEGER N,CLOCK,BER,SWAPNUM,SWAPTMP,SHUFFLE(NELEM),INTEREST(NELEM)
      INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
      DOUBLE PRECISION DICE,RISK,ZL,LOCALTOB,MAXTOB,CRITTOB,D,AVETOB
      DOUBLE PRECISION MAXLOCALTOB
      INTEGER LINEAR_OR_GAUSS,SS,K,COUPLINGPERIOD,PRINTOUTPERIOD
      DOUBLE PRECISION BID
      LOGICAL OK
!
!-----------------------------------------------------------------------
!
!      DO IELEM=1,NELEM
!        SHUFFLE(IELEM)=IELEM
!      ENDDO

      ZL=108.D0
      CRITTOB=1.D-2
      LINEAR_OR_GAUSS=2
      COUPLINGPERIOD=1
      PRINTOUTPERIOD=3600
      MAXTOB=MAXVAL(TOB%R)
      
!      DO I=1,NPOIN
!        IF (PRIVE%ADR(1)%P%R(I).LT.1.D0) THEN
!          ZR(I)=-100.D0
!        ELSE
!          ZR(I)=ZF(I)
!        ENDIF
!      ENDDO

      TANSL=TAN(4.D0*ATAN(1.D0)*SLOPE/180.D0)
!
!     INITIALISES THE RIGHT-HAND SIDE EVOL TO ZERO
!
      CALL CPSTVC(UNSV2D,EVOL)
      CALL OS('X=0     ',X=EVOL)

      IF(MOD(LT,COUPLINGPERIOD).EQ.0) THEN

      DO IELEM=1,NELEM
        SHUFFLE(IELEM)=IELEM
      ENDDO

      DO I=1,NPOIN
        IF (PRIVE%ADR(1)%P%R(I).LT.1.D0) THEN
          ZR(I)=ZF(I)-100.D0
        ELSE
          ZR(I)=ZF(I)
        ENDIF
      ENDDO

!
!     ONE CLASS VERSION
!
      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL SYSTEM_CLOCK(COUNT=CLOCK)
      SEED=CLOCK+37*(/(I-1,I=1,N)/)
      CALL RANDOM_SEED(PUT=SEED)

      DO IELEM=1,NELEM
        CALL RANDOM_NUMBER(DICE)
        SWAPNUM=INT(DICE*NELEM)
        IF(SWAPNUM.EQ.0) SWAPNUM=1
        IF(SWAPNUM.EQ.(NELEM+1)) SWAPNUM=NELEM
        SWAPTMP=SHUFFLE(IELEM)
        SHUFFLE(IELEM)=SHUFFLE(SWAPNUM)
        SHUFFLE(SWAPNUM)=SWAPTMP
      ENDDO

      J=1
      DO IELEM=1,NELEM
        I=SHUFFLE(IELEM)
        I1=IKLE(I,1)
        I2=IKLE(I,2)
        I3=IKLE(I,3)
        IF((ZF(I1).GT.ZL.OR.ZF(I2).GT.ZL.OR.ZF(I3).GT.ZL).AND.
     &     (ZF(I1).LT.ZL.OR.ZF(I2).LT.ZL.OR.ZF(I3).LT.ZL).AND.
     &     (PRIVE%ADR(1)%P%R(I1).GT.0.D0.OR.
     &      PRIVE%ADR(1)%P%R(I2).GT.0.D0.OR.
     &      PRIVE%ADR(1)%P%R(I3).GT.0.D0).AND.
     &     (PRIVE%ADR(1)%P%R(I1).LT.1.D0.OR.
     &      PRIVE%ADR(1)%P%R(I2).LT.1.D0.OR.
     &      PRIVE%ADR(1)%P%R(I3).LT.1.D0).AND.
     &    PRIVE%ADR(2)%P%R(I1).LT.1.D0.AND.
     &    PRIVE%ADR(2)%P%R(I2).LT.1.D0.AND.
     &    PRIVE%ADR(2)%P%R(I3).LT.1.D0) THEN
            INTEREST(J)=I
            J=J+1
        ENDIF
      ENDDO
      J=J-1

      AVETOB=0.D0
      MAXLOCALTOB=0.D0
      SS=0
      IF(J.NE.0) THEN
        DO IELEM=1,J
          CALL RANDOM_NUMBER(DICE)
          I=INTEREST(IELEM)
          I1=IKLE(I,1)
          I2=IKLE(I,2)
          I3=IKLE(I,3)
          LOCALTOB=MAX(TOB%R(I1),TOB%R(I2),TOB%R(I3))
          IF(LOCALTOB.GT.MAXLOCALTOB) MAXLOCALTOB=LOCALTOB
          AVETOB=AVETOB+LOCALTOB
          IF(LINEAR_OR_GAUSS.EQ.1) THEN
            IF(LOCALTOB.GT.MAXTOB) THEN
              RISK=1.D0
            ELSE
              RISK=LOCALTOB/MAXTOB
            ENDIF
          ELSEIF(LINEAR_OR_GAUSS.EQ.2) THEN
            D=4*(LOCALTOB-CRITTOB)/SQRT(LOCALTOB**2.D0+CRITTOB**2.D0)
            RISK=0.5D0*(1.D0+ERF(D/SQRT(2.D0)))
          ENDIF
        
          BER=0
          IF(DICE.LT.RISK) THEN
            BER=1
            SS=SS+1
          ELSE
            BER=0
          ENDIF
        
!        IF(DICE.LT.0.333333D0) THEN
!          K=1
!        ELSEIF(DICE.LT.0.666667D0) THEN
!          K=2
!        ELSE
!          K=3
!        ENDIF

          IF(PRIVE%ADR(1)%P%R(I1).EQ.1.D0.AND.BER.EQ.1) THEN
            PRIVE%ADR(1)%P%R(I1)=0.D0
!            PRINT *,'BER1_YES!'
          ENDIF
          IF(PRIVE%ADR(1)%P%R(I2).EQ.1.D0.AND.BER.EQ.1) THEN
            PRIVE%ADR(1)%P%R(I2)=0.D0
!            PRINT *,'BER2_YES!'
          ENDIF
          IF(PRIVE%ADR(1)%P%R(I3).EQ.1.D0.AND.BER.EQ.1) THEN
            PRIVE%ADR(1)%P%R(I3)=0.D0
!            PRINT *,'BER3_YES!'
          ENDIF
        ENDDO

        AVETOB=AVETOB/J
      ENDIF

      IF(MOD(LT,PRINTOUTPERIOD).EQ.0) THEN

      PRINT *,''
      PRINT *,'============= SHEAR STRESS INFORMATION =============='
      PRINT *,'==             ',LT/COUPLINGPERIOD,'         
     &                       =='
      PRINT *,'== ----------------------------------------------- =='
      PRINT *,'== CRITICAL TAU = ',CRITTOB,'        =='
      PRINT *,'== ----------------------------------------------- =='
      PRINT *,'== AVERAGE  TAU = ',AVETOB,'        =='
      PRINT *,'== RATIO        = ',AVETOB/CRITTOB,'        =='
      D=4*(AVETOB-CRITTOB)/SQRT(AVETOB**2.D0+CRITTOB**2.D0)
      RISK=0.5D0*(1.D0+ERF(D/SQRT(2.D0)))
      PRINT *,'== AVERAGE RISK = ',RISK,'        =='
      PRINT *,'== ----------------------------------------------- =='
      PRINT *,'== MAXIMUM  TAU = ',MAXLOCALTOB,'        =='
      PRINT *,'== RATIO        = ',MAXLOCALTOB/CRITTOB,'        =='
      D=4*(MAXLOCALTOB-CRITTOB)/SQRT(MAXLOCALTOB**2.D0+CRITTOB**2.D0)
      RISK=0.5D0*(1.D0+ERF(D/SQRT(2.D0)))
      PRINT *,'== MAXIMUM RISK = ',RISK,'        =='
      PRINT *,'== ----------------------------------------------- =='
      PRINT *,'== BER SUCC #   = ',SS,'                    =='
      PRINT *,'====================================================='
      PRINT *,''

      ENDIF

      DEALLOCATE(SEED)

      ENDIF

      IF(NSICLA.EQ.1) THEN
!
!       LOOP ON ELEMENTS
!
        DO IELEM=1,NELEM
!
          I1=IKLE(IELEM,1)
          I2=IKLE(IELEM,2)
          I3=IKLE(IELEM,3)
!
          X2=XEL(IELEM,2)
          X3=XEL(IELEM,3)
          Y2=YEL(IELEM,2)
          Y3=YEL(IELEM,3)
          Z2=ZF(I2)-ZF(I1)
          Z3=ZF(I3)-ZF(I1)
!
!         TWICE THE TRIANGLE AREA
!
          DEUXSURF=X2*Y3-X3*Y2
!
!         AVERAGE BOTTOM IN THE ELEMENT
!
          ZC=(ZF(I1)+ZF(I2)+ZF(I3))/3.D0
!
!         COMPONENTS OF BOTTOM GRADIENT
!
          A=(Z2*Y3-Z3*Y2)/DEUXSURF
          B=(Z3*X2-Z2*X3)/DEUXSURF
!
!         CORRECTING FACTOR ON SLOPE
!
          L=MIN(1.D0,TANSL/MAX(SQRT(A**2+B**2),1.D-8))
!
!         L LIMITED DUE TO NON-ERODABLE BEDS : ZF MUST NOT GO BELOW ZR
!
          IF(ZF(I1).GT.ZC) L=MAX(L,(ZR(I1)-ZC)/MAX(ZF(I1)-ZC,1.D-8))
          IF(ZF(I2).GT.ZC) L=MAX(L,(ZR(I2)-ZC)/MAX(ZF(I2)-ZC,1.D-8))
          IF(ZF(I3).GT.ZC) L=MAX(L,(ZR(I3)-ZC)/MAX(ZF(I3)-ZC,1.D-8))
!
!         BUILDS THE RIGHT-HAND SIDE
!
!         HERE THE EVOLUTIONS ARE MULTIPLIED BY SURFAC/3
!         BECAUSE THE REAL EVOLUTION TAKING INTO ACCOUNT OTHER ELEMENTS
!         WILL NEED A FACTOR (SURFAC/3)/(INTEGRAL OF BASIS)
!
          EVOL%R(I1)=EVOL%R(I1)+(1.D0-L)*(ZC-ZF(I1))*DEUXSURF/6.D0
          EVOL%R(I2)=EVOL%R(I2)+(1.D0-L)*(ZC-ZF(I2))*DEUXSURF/6.D0
          EVOL%R(I3)=EVOL%R(I3)+(1.D0-L)*(ZC-ZF(I3))*DEUXSURF/6.D0
!
        ENDDO
!
      ELSE
!
!       MULTI-CLASS VERSION
!
!       INITIALING TO 0. THE EVOLUTIONS DUE TO SLIDE FOR EACH CLASS
!
        DO I=1,NSICLA
          CALL OS('X=0     ',X=ZFCL_MS%ADR(I)%P)
        ENDDO
!
!       LOOP ON ELEMENTS
!
        DO IELEM=1,NELEM
!
          I1=IKLE(IELEM,1)
          I2=IKLE(IELEM,2)
          I3=IKLE(IELEM,3)
!
          X2=XEL(IELEM,2)
          X3=XEL(IELEM,3)
          Y2=YEL(IELEM,2)
          Y3=YEL(IELEM,3)
          Z2=ZF(I2)-ZF(I1)
          Z3=ZF(I3)-ZF(I1)
!
!         TWICE THE TRIANGLE AREA
!
          DEUXSURF=X2*Y3-X3*Y2
!
!         AVERAGE BOTTOM IN THE ELEMENT
!
          ZC=(ZF(I1)+ZF(I2)+ZF(I3))/3.D0
!
!         COMPONENTS OF BOTTOM GRADIENT
!
          A=(Z2*Y3-Z3*Y2)/DEUXSURF
          B=(Z3*X2-Z2*X3)/DEUXSURF
!
!         CORRECTING FACTOR ON SLOPE
!
          L=MIN(1.D0,TANSL/MAX(SQRT(A**2+B**2),1.D-8))
!
!         L LIMITED DUE TO NON-ERODABLE BEDS : ZF MUST NOT GO BELOW ZR
!
          IF(ZF(I1).GT.ZC) L=MAX(L,(ZR(I1)-ZC)/MAX(ZF(I1)-ZC,1.D-8))
          IF(ZF(I2).GT.ZC) L=MAX(L,(ZR(I2)-ZC)/MAX(ZF(I2)-ZC,1.D-8))
          IF(ZF(I3).GT.ZC) L=MAX(L,(ZR(I3)-ZC)/MAX(ZF(I3)-ZC,1.D-8))
!
!         BUILDS THE RIGHT-HAND SIDE
!
!         HERE THE EVOLUTIONS ARE MULTIPLIED BY SURFAC/3
!         BECAUSE THE REAL EVOLUTION TAKING INTO ACCOUNT OTHER ELEMENTS
!         WILL NEED A FACTOR (SURFAC/3)/(INTEGRAL OF BASIS)
!
!         FIRST IN TERMS OF QUANTITIES BROUGHT TO POINTS
!
          Q(1)=(1.D0-L)*(ZC-ZF(I1))*DEUXSURF/6.D0
          Q(2)=(1.D0-L)*(ZC-ZF(I2))*DEUXSURF/6.D0
          Q(3)=(1.D0-L)*(ZC-ZF(I3))*DEUXSURF/6.D0
!
          EVOL%R(I1)=EVOL%R(I1)+Q(1)
          EVOL%R(I2)=EVOL%R(I2)+Q(2)
          EVOL%R(I3)=EVOL%R(I3)+Q(3)
!
!         TAKING INTO ACCOUNT THE QUANTITIES TO UPDATE ZFCL_MS
!         IG1 AND IG2 : POINTS THAT GIVE
!         IR1 AND IR2 : POINTS THAT RECEIVE
!         CASE2: TWO POINTS GIVE TO THE THIRD ONE (THE OTHER CASE IS
!                ONE POINT GIVES TO THE TWO OTHERS)
          CASE2=.FALSE.
!
!         PARAMETERISING TO REDUCE THE 6 CASES TO 2
!
          IF(Q(1).GE.0.D0) THEN
            IF(Q(2).GE.0.D0) THEN
!             3 GIVES TO 1 AND 2
              IG1=I3
              QG1=Q(3)
              IR1=I1
              QR1=Q(1)
              IR2=I2
              QR2=Q(2)
            ELSE
              IF(Q(3).GE.0.D0) THEN
!               2 GIVES TO 1 AND 3
                IG1=I2
                QG1=Q(2)
                IR1=I1
                QR1=Q(1)
                IR2=I3
                QR2=Q(3)               
              ELSE
!               2 AND 3 GIVE TO 1
                IG1=I2
                QG1=Q(2)
                IG2=I3
                QG2=Q(3)
                IR1=I1
                QR1=Q(1)
                CASE2=.TRUE.
              ENDIF
            ENDIF
          ELSE
            IF(Q(2).GT.0.D0) THEN
              IF(Q(3).GT.0.D0) THEN
!               1 GIVES TO 2 AND 3
                IG1=I1
                QG1=Q(1)
                IR1=I2
                QR1=Q(2)
                IR2=I3
                QR2=Q(3)
              ELSE
!               1 AND 3 GIVE TO 2
                IG1=I1
                QG1=Q(1)
                IG2=I3
                QG2=Q(3)
                IR1=I2
                QR1=Q(2)
                CASE2=.TRUE.
              ENDIF
            ELSE
!             1 AND 2 GIVE TO 3
              IG1=I1
              QG1=Q(1)
              IG2=I2
              QG2=Q(2)
              IR1=I3
              QR1=Q(3)
              CASE2=.TRUE.
            ENDIF
          ENDIF
!
          IF(CASE2) THEN
!
!           THE TWO DONNORS CASE : IG1 AND IG2 GIVE TO IR1
!           ZFCL_MS IS HERE VOLUMES
!
            DO I=1,NSICLA
              ZFCL_MS%ADR(I)%P%R(IG1)=ZFCL_MS%ADR(I)%P%R(IG1)
     &                               +QG1*AVAIL(IG1,1,I)
              ZFCL_MS%ADR(I)%P%R(IG2)=ZFCL_MS%ADR(I)%P%R(IG2)
     &                               +QG2*AVAIL(IG2,1,I)
              ZFCL_MS%ADR(I)%P%R(IR1)=ZFCL_MS%ADR(I)%P%R(IR1)
     &                               -QG1*AVAIL(IG1,1,I)
     &                               -QG2*AVAIL(IG2,1,I)
            ENDDO            
!
          ELSE
!
!           THE ONE DONNOR CASE : IG1 GIVES TO IR1 AND IR2
!           ZFCL_MS IS HERE VOLUMES
!
            DO I=1,NSICLA
              ZFCL_MS%ADR(I)%P%R(IG1)=ZFCL_MS%ADR(I)%P%R(IG1)
     &                               +QG1*AVAIL(IG1,1,I)
              ZFCL_MS%ADR(I)%P%R(IR1)=ZFCL_MS%ADR(I)%P%R(IR1)
     &                               +QR1*AVAIL(IG1,1,I)
              ZFCL_MS%ADR(I)%P%R(IR2)=ZFCL_MS%ADR(I)%P%R(IR2)
     &                               +QR2*AVAIL(IG1,1,I)
            ENDDO                    
!
          ENDIF
!
        ENDDO       
!
!       ADDING VOLUMES IN PARALLEL
!
        IF(NCSIZE.GT.1) THEN
          DO I=1,NSICLA
            CALL PARCOM(ZFCL_MS%ADR(I)%P,2,MESH)
          ENDDO
        ENDIF
!
!       FINAL DIVISION BY THE INTEGRAL OF BASES: VOLUMES CHANGED INTO
!       BED VARIATIONS (LIKE EVOL BELOW)
!
        DO I=1,NSICLA
          DO J=1,NPOIN
            ZFCL_MS%ADR(I)%P%R(J)=ZFCL_MS%ADR(I)%P%R(J)
     &                               *UNSV2D%R(J)/MOFAC
          ENDDO
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!     FINAL RESOLUTION
!
!      PRINT *,MAXVAL(TOB%R),'---TOB BEFOREPARCOM---'
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM(EVOL,2,MESH)
        CALL PARCOM(PRIVE%ADR(1)%P,4,MESH)
      ENDIF
!      PRINT *,MAXVAL(TOB%R),'---TOB AFTERPARCOM---'

!      DO I=1,NPOIN
!        IF(PRIVE%ADR(1)%P%R(I).EQ.0.D0.OR.
!     &     PRIVE%ADR(1)%P%R(I).EQ.1.D0) THEN
!        ELSE
!          PRINT *,'FOUND V23 NOT 0 OR 1 ERROR!!!!!!!!!!!!!!!!!!!!!!!!!'
!        ENDIF
!      ENDDO
!
!     FINAL DIVISION BY THE INTEGRAL OF BASES: QUANTITIES CHANGED INTO
!     ELEVATIONS
!
      DO I=1,NPOIN
        EVOL%R(I)=EVOL%R(I)*UNSV2D%R(I)/MOFAC
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END





