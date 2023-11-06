! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE Z0WAVE (US, TAUW, UTOP, Z0, Z0B, CHRNCK)
!$loki routine seq

! ----------------------------------------------------------------------

!**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.

!*    PURPOSE.
!     --------

!       COMPUTE ROUGHNESS LENGTH. 

!**   INTERFACE.
!     ----------

!       *CALL* *Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
!          *KIJS* - INDEX OF FIRST GRIDPOINT.
!          *KIJL* - INDEX OF LAST GRIDPOINT.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *UTOP* - WIND SPEED.
!          *Z0*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *Z0B*  - BACKGROUND ROUGHNESS LENGTH.
!          *CHRNCK- CHARNOCK COEFFICIENT

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK
      USE YOWPCONS , ONLY : G, GM1
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWTABL  , ONLY : EPS1

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "chnkmin.intfb.h"

      REAL(KIND=JWRB),INTENT(IN)  ::  US, TAUW, UTOP
      REAL(KIND=JWRB),INTENT(OUT) ::  Z0, Z0B, CHRNCK


      REAL(KIND=JWRB) :: UST2, UST3, ARG
      REAL(KIND=JWRB) :: ALPHAOG

! ----------------------------------------------------------------------

      IF (LLCAPCHNK) THEN
        ALPHAOG= CHNKMIN(UTOP)*GM1
      ELSE
        ALPHAOG= ALPHA*GM1
      ENDIF

      UST2 = US**2
      UST3 = US**3
      ARG = MAX(UST2-TAUW,EPS1)
      Z0 = ALPHAOG*UST3/SQRT(ARG)
      Z0B = ALPHAOG*UST2
      CHRNCK = G*Z0/UST2

      END SUBROUTINE Z0WAVE
