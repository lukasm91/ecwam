! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TAUT_Z0(IUSFG,          &
&                  HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC, &
&                  USTAR, Z0, Z0B, CHRNCK)
!$loki routine seq

! ----------------------------------------------------------------------

!**** *TAUT_Z0* - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.


!**   INTERFACE.
!     ----------

!       *CALL* *TAUT_Z0(KIJS, KIJL, IUSFG, FL1, WAVNUM,
!                       UTOP, UDIR, TAUW, TAUWDIR, RNFAC,
!                       USTAR, Z0, Z0B, CHRNCK)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                           0 DO NOT USE THE FIELD USTAR
!          *FL1*     - 2D-SPECTRA
!          *WAVNUM*  - WAVE NUMBER
!          *HALP*    - 1/2 PHILLIPS PARAMETER
!          *UTOP*    - WIND SPEED AT REFERENCE LEVEL XNLEV
!          *UDIR*    - WIND SPEED DIRECTION AT REFERENCE LEVEL XNLEV
!          *TAUW*    - WAVE STRESS.
!          *TAUWDIR* - WAVE STRESS DIRECTION.
!          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *USTAR*   - FRICTION VELOCITY
!          *Z0*      - ROUGHNESS LENGTH
!          *Z0B*     - BACKGROUND ROUGHNESS LENGTH
!          *CHRNCK*  - CHARNOCK COEFFICIENT

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESS LENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!       IT WAS EXTENDED TO INCLUDE THE GRAVITY-CAPILLARY MODEL FOR THE CALCULATION
!       OF THE BACKGROUND ROUGHNESS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK, LLGCBZ0
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G, GM1, EPSUS, EPSMIN, ACD, BCD, CDMAX
      USE YOWPHYS  , ONLY : XKAPPA, XNLEV, RNU, RNUM, ALPHA, ALPHAMIN, ALPHAMAX, &
     &                      ANG_GC_A, ANG_GC_B, ANG_GC_C
      USE YOWTABL  , ONLY : EPS1 

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "chnkmin.intfb.h"
#include "stress_gc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IUSFG
      REAL(KIND=JWRB), INTENT(IN) :: HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC
      REAL(KIND=JWRB), INTENT(INOUT) :: USTAR
      REAL(KIND=JWRB), INTENT(OUT) :: Z0, Z0B, CHRNCK


      INTEGER(KIND=JWIM), PARAMETER :: NITER=17

      REAL(KIND=JWRB), PARAMETER :: TWOXMP1=3.0_JWRB

      INTEGER(KIND=JWIM) :: IDX, ITER
      INTEGER(KIND=JWIM) :: IFRPH

      ! Cd and Z0 from Hersbach 2010, ECMWF Tech Memo (without the viscous part)
!     CD = ACDLIN + BCDLIN*SQRT(PCHAR) * U10
      REAL(KIND=JWRB), PARAMETER :: ACDLIN=0.0008_JWRB
      REAL(KIND=JWRB), PARAMETER :: BCDLIN=0.00047_JWRB
      REAL(KIND=JWRB) :: ALPHAGM1

      REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
      REAL(KIND=JWRB) :: PCE_GC
      REAL(KIND=JWRB) :: Z0MINRST
      REAL(KIND=JWRB) :: CHARNOCK_MIN
      REAL(KIND=JWRB) :: COSDIFF 
      REAL(KIND=JWRB) :: ZCHAR
      REAL(KIND=JWRB) :: US2TOTAUW, USMAX
      REAL(KIND=JWRB) :: XLOGXL, XKUTOP, XOLOGZ0
      REAL(KIND=JWRB) :: USTOLD, USTNEW, TAUOLD, TAUNEW, X, F, DELF, CDFG
      REAL(KIND=JWRB) :: USTM1, Z0TOT, Z0CH, Z0VIS, HZ0VISO1MX, ZZ
      REAL(KIND=JWRB) :: CONST, TAUV, DEL
      REAL(KIND=JWRB) :: RNUEFF, RNUKAPPAM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ALPHAOG, XMIN
      REAL(KIND=JWRB) :: W1
      REAL(KIND=JWRB) :: TAUWACT, TAUWEFF 
      REAL(KIND=JWRB) :: ANG_GC, TAUUNR

      LOGICAL :: LLCOSDIFF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',0,ZHOOK_HANDLE)

      XLOGXL=LOG(XNLEV)
      US2TOTAUW=1.0_JWRB+EPS1

!     ONLY take the contribution of TAUW that is in the wind direction
        COSDIFF = COS(UDIR-TAUWDIR)
        TAUWACT = MAX(TAUW*COSDIFF, EPSMIN )
        LLCOSDIFF = (COSDIFF > 0.9_JWRB )

!  USING THE CG MODEL:
       IF (LLGCBZ0) THEN

        IF (LLCAPCHNK) THEN
          CHARNOCK_MIN = CHNKMIN(UTOP)
          ALPHAOG = CHARNOCK_MIN*GM1
        ELSE
        ALPHAOG= 0.0_JWRB
        ENDIF

        USMAX = MAX(-0.21339_JWRB + 0.093698_JWRB*UTOP -0.0020944_JWRB*UTOP**2 + 5.5091E-5_JWRB*UTOP**3, 0.03_JWRB)
        TAUWEFF = MIN(TAUWACT*US2TOTAUW, USMAX**2 )

        RNUEFF = 0.04_JWRB*RNU

        RNUKAPPAM1 = RNUEFF/XKAPPA

        PCE_GC = 0.001_JWRB * IUSFG + (1-IUSFG) * 0.005_JWRB

        IF (IUSFG == 0 ) THEN
          ALPHAGM1 = ALPHA*GM1
          IF ( UTOP < 1.0_JWRB ) THEN
            CDFG = 0.002_JWRB
          ELSEIF ( LLCOSDIFF ) THEN
            X = MIN(TAUWACT/MAX(USTAR,EPSUS)**2,0.99_JWRB)
            ZCHAR = MIN( ALPHAGM1 * USTAR**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(UTOP-35._JWRB)) )
            ZCHAR = MIN(ZCHAR,ALPHAMAX)
            CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR) * UTOP
          ELSE
            CDFG = CDM(UTOP)
          ENDIF
          USTAR = UTOP*SQRT(CDFG)
        ENDIF

        W1 = 0.85_JWRB - 0.05_JWRB*( TANH(10.0_JWRB*(UTOP-5.0_JWRB)) + 1.0_JWRB )

        XKUTOP = XKAPPA * UTOP

        USTOLD = USTAR
        TAUOLD = USTOLD**2

        DO ITER=1,NITER
!         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
          Z0 = MAX(XNLEV/(EXP(MIN(XKUTOP/USTOLD, 50.0_JWRB))-1.0_JWRB), Z0MIN)
          ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
          TAUV = RNUKAPPAM1*USTOLD/Z0

          ANG_GC = ANG_GC_A + ANG_GC_B * TANH(ANG_GC_C * TAUOLD)

          TAUUNR = STRESS_GC(ANG_GC, USTAR, Z0, Z0MIN, HALP, RNFAC)

!         TOTAL kinematic STRESS:
          TAUNEW = TAUWEFF + TAUV + TAUUNR
          USTNEW = SQRT(TAUNEW)
          USTAR = W1*USTOLD+(1.0_JWRB-W1)*USTNEW

!         CONVERGENCE ?
          DEL = USTAR-USTOLD
          IF (ABS(DEL) < PCE_GC*USTAR) EXIT 
          TAUOLD = USTAR**2
          USTOLD = USTAR
        ENDDO
        ! protection just in case there is no convergence
        IF (ITER > NITER ) THEN
          CDFG = CDM(UTOP)
          USTAR = UTOP*SQRT(CDFG)
          Z0MINRST = USTAR**2 * ALPHA*GM1
          Z0 = MAX(XNLEV/(EXP(XKUTOP/USTAR)-1.0_JWRB), Z0MINRST)
          Z0B = Z0MINRST
        ELSE
          Z0 = MAX(XNLEV/(EXP(XKUTOP/USTAR)-1.0_JWRB), Z0MIN)
          Z0B = Z0*SQRT(TAUUNR/TAUOLD)
        ENDIF

!       Refine solution
        X = TAUWEFF/TAUOLD

        IF (X < 0.99_JWRB) THEN
          USTOLD = USTAR
          TAUOLD = MAX(USTOLD**2,TAUWEFF)

          DO ITER=1,NITER
            X = MIN(TAUWEFF/TAUOLD, 0.99_JWRB)
            USTM1 = 1.0_JWRB/MAX(USTOLD,EPSUS)
            !!!! Limit how small z0 could become
            !!!! This is a bit of a compromise to limit very low Charnock for intermediate high winds (15 -25 m/s)
            !!!! It is not ideal !!!
            Z0 = MAX(XNLEV/(EXP(MIN(XKUTOP/USTOLD, 50.0_JWRB))-1.0_JWRB), Z0MIN)

            TAUUNR = STRESS_GC(ANG_GC, USTOLD, Z0, Z0MIN, HALP, RNFAC)

            Z0B = MAX( Z0*SQRT(TAUUNR/TAUOLD), ALPHAOG*TAUOLD)
            Z0VIS = RNUM*USTM1
            HZ0VISO1MX = 0.5_JWRB*Z0VIS/(1.0_JWRB-X)
            Z0 = HZ0VISO1MX+SQRT(HZ0VISO1MX**2+Z0B**2/(1.0_JWRB-X))

            XOLOGZ0= 1.0_JWRB/(XLOGXL-LOG(Z0))
            F = USTOLD-XKUTOP*XOLOGZ0
            ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B**2+0.5_JWRB*Z0VIS*Z0-Z0**2) &
&                / (2.0_JWRB*Z0**2*(1.0_JWRB-X)-Z0VIS*Z0)

            DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ
            IF (DELF /= 0.0_JWRB) USTAR = USTOLD-F/DELF

!           CONVERGENCE ?
            DEL = USTAR-USTOLD

            IF (ABS(DEL) < PCE_GC*USTAR) EXIT 
            USTOLD = USTAR
            TAUOLD = MAX(USTOLD**2,TAUWEFF)
          ENDDO
          ! protection just in case there is no convergence
          IF (ITER > NITER ) THEN
            CDFG = CDM(UTOP)
            USTAR = UTOP*SQRT(CDFG)
            Z0MINRST = USTAR**2 * ALPHA*GM1
            Z0 = MAX(XNLEV/(EXP(XKUTOP/USTAR)-1.0_JWRB), Z0MINRST)
            Z0B = Z0MINRST
            CHRNCK = MAX(G*Z0/USTAR**2, ALPHAMIN)
          ELSE
            CHRNCK = MAX( G*(Z0B/SQRT(1.0_JWRB-X))/MAX(USTAR,EPSUS)**2, ALPHAMIN)
          ENDIF

        ELSE
          USTM1 = 1.0_JWRB/MAX(USTAR, EPSUS)
          Z0VIS = RNUM*USTM1
          CHRNCK = MAX(G*(Z0-Z0VIS) * USTM1**2, ALPHAMIN)
        ENDIF


       ELSE

        TAUWEFF = TAUWACT*US2TOTAUW

        IF (LLCAPCHNK) THEN
          CHARNOCK_MIN = CHNKMIN(UTOP)
          XMIN = 0.15_JWRB*(ALPHA-CHARNOCK_MIN)
          ALPHAOG = CHARNOCK_MIN*GM1
        ELSE
          XMIN= 0.0_JWRB
          ALPHAOG= ALPHA*GM1
        ENDIF

        XKUTOP = XKAPPA * UTOP

        USTOLD = (1-IUSFG)*UTOP*SQRT(MIN(ACD+BCD*UTOP,CDMAX)) + IUSFG*USTAR
        TAUOLD = MAX(USTOLD**2,TAUWEFF)
        USTAR = SQRT(TAUOLD)
        USTM1 = 1.0_JWRB/MAX(USTAR,EPSUS) 

        DO ITER=1,NITER
          X = MAX(TAUWACT/TAUOLD,XMIN)
          Z0CH = ALPHAOG*TAUOLD/SQRT(1.0_JWRB-X)
          Z0VIS = RNUM*USTM1
          Z0TOT = Z0CH+Z0VIS

          XOLOGZ0= 1.0_JWRB/(XLOGXL-LOG(Z0TOT))
          F = USTAR-XKUTOP*XOLOGZ0
          ZZ = USTM1*(Z0CH*(2.0_JWRB-TWOXMP1*X)/(1.0_JWRB-X)-Z0VIS)/Z0TOT
          DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ

          IF (DELF /= 0.0_JWRB) USTAR = USTAR-F/DELF
          TAUNEW = MAX(USTAR**2,TAUWEFF)
          USTAR = SQRT(TAUNEW)
          IF (TAUNEW == TAUOLD) EXIT
          USTM1 = 1.0_JWRB/MAX(USTAR,EPSUS)
          TAUOLD = TAUNEW
        ENDDO

        Z0 = Z0CH
        Z0B = ALPHAOG*TAUOLD
        CHRNCK = MAX(G*Z0*USTM1**2, ALPHAMIN)


       ENDIF

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',1,ZHOOK_HANDLE)

CONTAINS

!  INLINE FUNCTION.
!  ----------------

!  Simple empirical fit to model drag coefficient
   FUNCTION CDM(U10)
      !$loki routine seq
      REAL(KIND=JWRB), INTENT(IN) :: U10 
      REAL(KIND=JWRB) :: CDM

      CDM = MAX(MIN(0.0006_JWRB+0.00008_JWRB*U10, 0.001_JWRB+0.0018_JWRB*EXP(-0.05_JWRB*(U10-33._JWRB))),0.001_JWRB)
   END FUNCTION CDM

END SUBROUTINE TAUT_Z0
