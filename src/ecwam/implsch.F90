! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE IMPLSCH_FIRST_PART (KIJS, KIJL, FL1,                         &
 &                  WAVNUM, CGROUP, CIWA, &
 &                  EMAXDPT, &
 &                  WDWAVE, CICOVER, &
 &                  EMEAN, FMEAN, &
 &                  F1MEAN, AKMEAN, XKMEAN, &
 &                  FLM, &
 &                  COSWDIF, SINWDIF2, &
 &                  CIREDUC)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.


!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (KIJS, KIJL, FL1,
!    &                   WVPRPT,
!    &                   WVENVI, FF_NOW,
!    &                   INTFLDS, WAM2NEMO,
!    &                   MIJ,  XLLWS)
!      *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!      *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!      *FL1*     - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *WVPRPT*  - WAVE PROPERTIES FIELDS
!      *WVENVI*  - WAVE ENVIRONMENT  
!      *FF_NOW*    FORCING FIELDS
!      *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!      *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!      *MIJ*       LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *XLLWS*     TOTAL WINDSEA MASK FROM INPUT SOURCE TERM


!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL1 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.


!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,   &
 &                             INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU, LWNEMOCOUSTRN 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR  ,LICERUN   ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWPCONS , ONLY : WSEMEAN_MIN
      USE YOWSTAT  , ONLY : IDELT    ,LBIWBK
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
#include "ciwabr.intfb.h"
#include "fkmean.intfb.h"
#include "sdepthlim.intfb.h"
! ----------------------------------------------------------------------


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CIWA

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: EMAXDPT

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE, CICOVER


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: DELT, DELTM, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: F1MEAN, AKMEAN, XKMEAN

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: CIREDUC

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
          SINWDIF2(IJ,K) = SIN(TH(K)-WDWAVE(IJ))**2
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
      IF (LBIWBK) THEN
         CALL SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
      ENDIF

!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                    &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          FLM(IJ,K)=FLMIN*MAX(0.0_JWRB, COSWDIF(IJ,K))**2
        ENDDO
      ENDDO

!     COMPUTE DAMPING COEFFICIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
!!! testing sea ice attenuation (might need to restrict usage when needed)
      IF (LCIWABR) THEN
        CALL CIWABR(KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIREDUC)
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)*CIREDUC(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)


END SUBROUTINE IMPLSCH_FIRST_PART
SUBROUTINE IMPLSCH_SDISSIP (KIJS, KIJL, FL1,                         &
 &                  WAVNUM, CINV, XK2CG, &
 &                  INDEP, &
 &                  AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR, &
 &                  UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, &
 &                  MIJ, XLLWS, &
 &                  FMEANWS, EMEAN, FMEAN, &
 &                  F1MEAN, XKMEAN, &
 &                  PHIWA, &
 &                  FLM, &
 &                  COSWDIF, SINWDIF2, &
 &                  RHOWGDFTH, &
 &                  FLD, SL, SPOS, &
 &                  SSOURCE)


      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,   &
 &                             INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU, LWNEMOCOUSTRN 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR  ,LICERUN   ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWPCONS , ONLY : WSEMEAN_MIN
      USE YOWSTAT  , ONLY : IDELT    ,LBIWBK
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
#include "sdissip.intfb.h"
! ----------------------------------------------------------------------


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: XK2CG

      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: INDEP

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE, CICOVER, AIRD, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, WSWAVE
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: DELT, DELTM, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL) :: HALP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: FMEANWS, EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: PHIWA

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(INOUT) :: RHOWGDFTH
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: SSOURCE

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      DELTM = 1.0_JWRB/DELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------

      CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,                 &
     &              INDEP, WAVNUM, XK2CG,&
     &              EMEAN, F1MEAN, XKMEAN,                    &
     &              UFRIC, COSWDIF, AIRD)

!     Save source term contributions relevant for the calculation of ocean fluxes
      IF (LCFLX .AND. .NOT.LWVFLX_SNL) THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              SSOURCE(IJ,K,M) = SL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)


END SUBROUTINE IMPLSCH_SDISSIP
SUBROUTINE IMPLSCH_SNONLIN (KIJS, KIJL, FL1, FLD, SL,                &
 &                  WAVNUM, DEPTH, AKMEAN)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
#include "snonlin.intfb.h"
! ----------------------------------------------------------------------


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: AKMEAN
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)
      CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)
      IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

END SUBROUTINE IMPLSCH_SNONLIN
SUBROUTINE IMPLSCH_AFTER_SNONLIN (KIJS, KIJL, FL1,                         &
 &                  WAVNUM, CINV, XK2CG, STOKFAC, &
 &                  EMAXDPT, DEPTH, IOBND, IODP,      &
 &                  AIRD, WDWAVE, CICOVER, WSWAVE, &
 &                  UFRIC, CITHICK, &
 &                  NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, &
 &                  NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, &
 &                  NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
 &                  WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS, &
 &                  TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD, &
 &                  PHIEPS, PHIAW, &
 &                  MIJ, XLLWS, &
 &                  FMEANWS, EMEAN, FMEAN, &
 &                  F1MEAN, AKMEAN, XKMEAN, &
 &                  PHIWA, &
 &                  FLM, &
 &                  COSWDIF, &
 &                  RHOWGDFTH, &
 &                  FLD, SL, SPOS, &
 &                  CIREDUC, &
 &                  SSOURCE)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.


!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (KIJS, KIJL, FL1,
!    &                   WVPRPT,
!    &                   WVENVI, FF_NOW,
!    &                   INTFLDS, WAM2NEMO,
!    &                   MIJ,  XLLWS)
!      *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!      *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!      *FL1*     - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *WVPRPT*  - WAVE PROPERTIES FIELDS
!      *WVENVI*  - WAVE ENVIRONMENT  
!      *FF_NOW*    FORCING FIELDS
!      *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!      *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!      *MIJ*       LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *XLLWS*     TOTAL WINDSEA MASK FROM INPUT SOURCE TERM


!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL1 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.


!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,   &
 &                             INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU, LWNEMOCOUSTRN 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR  ,LICERUN   ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWPCONS , ONLY : WSEMEAN_MIN
      USE YOWSTAT  , ONLY : IDELT    ,LBIWBK
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sbottom.intfb.h"
#include "imphftail.intfb.h"
#include "sdiwbk.intfb.h"
#include "setice.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"
! ----------------------------------------------------------------------


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: STOKFAC

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IODP
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IOBND

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE, CICOVER, AIRD, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC, WSWAVE
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: DELT, DELTM, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB, FL1_LOC, SSOURCE_LOC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: EMEANWS, USFM
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: FMEANWS, EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: PHIWA

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(INOUT) :: COSWDIF
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(INOUT) :: RHOWGDFTH
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: CIREDUC
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: SSOURCE

      LOGICAL :: LCFLX

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      DELTM = 1.0_JWRB/DELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*    2.3.1 ITERATIVELY UPDATE STRESS AND COMPUTE WIND INPUT TERMS. 
!           -------------------------------------------------------

      DO IJ=KIJS,KIJL
        IF (LCFLX .AND. LWVFLX_SNL) THEN
!     Save source term contributions relevant for the calculation of ocean fluxes
!!!!!!  SL must only contain contributions contributed to fluxes into the oceans
!       MODULATE SL BY IMPLICIT FACTOR
          DO M=1,NFRE
            DO K=1,NANG
              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              SSOURCE(IJ,K,M) = SL(IJ,K,M)/GTEMP1
            ENDDO
          ENDDO
        ENDIF


        IF (LBIWBK) THEN
          CALL SDIWBK_PW(IJ, KIJL, FL1 ,FLD, SL, DEPTH(IJ), EMAXDPT(IJ), EMEAN(IJ), F1MEAN(IJ))
        ENDIF

        CALL SBOTTOM_PW(IJ, KIJL, FL1, FLD, SL, WAVNUM, DEPTH(IJ))

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!     INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!     FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

        USFM = UFRIC(IJ)*MAX(FMEANWS(IJ), FMEAN(IJ))

        IF (LLUNSTR) THEN
          DO K=1,NANG
            DO M=1,NFRE
              FL1_LOC = FL1(IJ,K,M)
              SSOURCE_LOC = SSOURCE(IJ,K,M)

              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
              FLHAB = ABS(GTEMP2)
              FLHAB = MIN(FLHAB,USFM*COFRM4(M)*DELT)
              FL1_LOC = FL1_LOC + IOBND(IJ)*SIGN(FLHAB,GTEMP2)
              FL1_LOC = MAX(IODP(IJ)*CIREDUC(IJ,K,M)*FL1_LOC,FLM(IJ,K))
              SSOURCE_LOC = SSOURCE_LOC + DELTM * MIN(FLMAX(M)-FL1_LOC,0.0_JWRB)
              FL1_LOC = MIN(FL1_LOC,FLMAX(M))

              FL1(IJ,K,M) = FL1_LOC
              SSOURCE(IJ,K,M) = SSOURCE_LOC
            ENDDO
          ENDDO
        ELSE
          DO K=1,NANG
            DO M=1,NFRE
              FL1_LOC = FL1(IJ,K,M)
              SSOURCE_LOC = SSOURCE(IJ,K,M)

              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
              FLHAB = ABS(GTEMP2)
              FLHAB = MIN(FLHAB,USFM*COFRM4(M)*DELT)
              FL1_LOC = FL1_LOC + SIGN(FLHAB,GTEMP2)
              FL1_LOC = MAX(CIREDUC(IJ,K,M)*FL1_LOC,FLM(IJ,K))
              SSOURCE_LOC = SSOURCE_LOC + DELTM * MIN(FLMAX(M)-FL1_LOC,0.0_JWRB)
              FL1_LOC = MIN(FL1_LOC,FLMAX(M))

              FL1(IJ,K,M) = FL1_LOC
              SSOURCE(IJ,K,M) = SSOURCE_LOC
            ENDDO
          ENDDO
        ENDIF

        IF (LCFLX) THEN
          CALL WNFLUXES_PW (IJ, KIJL,                              &
     &                 MIJ(IJ), RHOWGDFTH,                          &
     &                 CINV,                             &
     &                 SSOURCE, CICOVER(IJ),                 &
     &                 PHIWA(IJ),                                   &
     &                 EMEAN(IJ), F1MEAN(IJ), WSWAVE(IJ),            &
     &                 WDWAVE(IJ), UFRIC(IJ), AIRD(IJ),&
     &                 NPHIEPS(IJ), NTAUOC(IJ), NSWH(IJ), NMWP(IJ), NEMOTAUX(IJ), &
     &                 NEMOTAUY(IJ), NEMOWSWAVE(IJ), NEMOPHIF(IJ), &
     &                 TAUXD(IJ), TAUYD(IJ), TAUOCXD(IJ), TAUOCYD(IJ), TAUOC(IJ), &
     &                 PHIOCD(IJ), PHIEPS(IJ), PHIAW(IJ), &
     &                 .TRUE.)
        ENDIF
! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------

        CALL FKMEAN_PW(IJ, KIJL, FL1, WAVNUM,                      &
     &              EMEAN(IJ), FMEAN(IJ), F1MEAN(IJ), AKMEAN(IJ), XKMEAN(IJ))

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
        CALL FEMEANWS_PW(IJ, KIJL, FL1, XLLWS, FMEANWS(IJ), EMEANWS)

        CALL IMPHFTAIL_PW(IJ, KIJL, MIJ(IJ), FLM, WAVNUM, XK2CG, FL1)


!     UPDATE WINDSEA VARIANCE AND MEAN FREQUENCY IF PASSED TO ATMOSPHERE
!     ------------------------------------------------------------------
        IF (LWFLUX) THEN
          IF (EMEANWS < WSEMEAN_MIN) THEN
            WSEMEAN(IJ) = WSEMEAN_MIN 
            WSFMEAN(IJ) = 2._JWRB*FR(NFRE)
          ELSE
            WSEMEAN(IJ) = EMEANWS
            WSFMEAN(IJ) = FMEANWS(IJ) 
          ENDIF
        ENDIF


!*    2.6 SET FL1 ON ICE POINTS TO ZERO
!         -----------------------------

        IF (LICERUN .AND. LMASKICE) THEN
          CALL SETICE_PW(IJ, KIJL, FL1, CICOVER(IJ), COSWDIF)
        ENDIF


!*    2.7 SURFACE STOKES DRIFT AND STRAIN IN SEA ICE
!         ------------------------------------------

        CALL STOKESTRN_PW(IJ, KIJL, FL1, WAVNUM, STOKFAC, DEPTH(IJ), WSWAVE(IJ), WDWAVE(IJ), CICOVER(IJ), CITHICK(IJ), &
 &                     USTOKES(IJ), VSTOKES(IJ), STRNMS(IJ), NEMOUSTOKES(IJ), NEMOVSTOKES(IJ), NEMOSTRN(IJ))
      ENDDO

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)


END SUBROUTINE IMPLSCH_AFTER_SNONLIN

SUBROUTINE IMPLSCH (KIJS, KIJL, FL1,                         &
 &                  WAVNUM, CGROUP, CIWA, CINV, XK2CG, STOKFAC, &
 &                  EMAXDPT, INDEP, DEPTH, IOBND, IODP,      &
 &                  AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR, &
 &                  UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, CITHICK, &
 &                  NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, &
 &                  NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, &
 &                  NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
 &                  WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS, &
 &                  TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD, &
 &                  PHIEPS, PHIAW, &
 &                  MIJ, XLLWS)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.


!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (KIJS, KIJL, FL1,
!    &                   WVPRPT,
!    &                   WVENVI, FF_NOW,
!    &                   INTFLDS, WAM2NEMO,
!    &                   MIJ,  XLLWS)
!      *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!      *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!      *FL1*     - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *WVPRPT*  - WAVE PROPERTIES FIELDS
!      *WVENVI*  - WAVE ENVIRONMENT  
!      *FF_NOW*    FORCING FIELDS
!      *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!      *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!      *MIJ*       LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *XLLWS*     TOTAL WINDSEA MASK FROM INPUT SOURCE TERM


!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL1 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.


!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,   &
 &                             INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU, LWNEMOCOUSTRN 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR  ,LICERUN   ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWPCONS , ONLY : WSEMEAN_MIN
      USE YOWSTAT  , ONLY : IDELT    ,LBIWBK
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
#include "ciwabr.intfb.h"
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sbottom.intfb.h"
#include "imphftail.intfb.h"
#include "sdepthlim.intfb.h"
#include "sdissip.intfb.h"
#include "sdiwbk.intfb.h"
#include "setice.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"
! ----------------------------------------------------------------------


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CIWA
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: STOKFAC

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: INDEP
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IODP
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IOBND

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE, CICOVER, AIRD, WSTAR, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, WSWAVE
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: DELT, DELTM, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL) :: EMEAN, FMEAN, HALP
      REAL(KIND=JWRB), DIMENSION(KIJL) :: EMEANWS, FMEANWS, USFM
      REAL(KIND=JWRB), DIMENSION(KIJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL(KIND=JWRB), DIMENSION(KIJL) :: PHIWA

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: RHOWGDFTH
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: CIREDUC 
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: SSOURCE 

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      DELTM = 1.0_JWRB/DELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
          SINWDIF2(IJ,K) = SIN(TH(K)-WDWAVE(IJ))**2
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
      IF (LBIWBK) THEN
         CALL SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
      ENDIF

!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                    &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          FLM(IJ,K)=FLMIN*MAX(0.0_JWRB, COSWDIF(IJ,K))**2
        ENDDO
      ENDDO

!     COMPUTE DAMPING COEFFICIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
!!! testing sea ice attenuation (might need to restrict usage when needed)
      IF (LCIWABR) THEN
        CALL CIWABR(KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIREDUC)
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)*CIREDUC(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*    2.3.1 ITERATIVELY UPDATE STRESS AND COMPUTE WIND INPUT TERMS. 
!           -------------------------------------------------------

      LUPDTUS = .TRUE.
      NCALL = 2
      DO ICALL = 1, NCALL 
        CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                      &
     &               LUPDTUS,                                       &
     &               FL1,                                           &
     &               WAVNUM, CINV, XK2CG,      &
     &               WSWAVE, WDWAVE, AIRD,     &
     &               WSTAR, CICOVER,           &
     &               COSWDIF, SINWDIF2,                             &
     &               FMEAN, HALP, FMEANWS,                          &
     &               FLM,                                           &
     &               UFRIC, TAUW, TAUWDIR,     &
     &               Z0M, Z0B, CHRNCK, PHIWA,  &
     &               FLD, SL, SPOS,                                 &
     &               MIJ, RHOWGDFTH, XLLWS)

      ENDDO

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------

      CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,                 &
     &              INDEP, WAVNUM, XK2CG,&
     &              EMEAN, F1MEAN, XKMEAN,                    &
     &              UFRIC, COSWDIF, AIRD)

!     Save source term contributions relevant for the calculation of ocean fluxes
      IF (LCFLX .AND. .NOT.LWVFLX_SNL) THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              SSOURCE(IJ,K,M) = SL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

      IF (LCFLX .AND. LWVFLX_SNL) THEN
!     Save source term contributions relevant for the calculation of ocean fluxes
!!!!!!  SL must only contain contributions contributed to fluxes into the oceans
!       MODULATE SL BY IMPLICIT FACTOR
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              SSOURCE(IJ,K,M) = SL(IJ,K,M)/GTEMP1
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      CALL SDIWBK(KIJS, KIJL, FL1 ,FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)

      CALL SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!     INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!     FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DO IJ=KIJS,KIJL
        USFM(IJ) = UFRIC(IJ)*MAX(FMEANWS(IJ), FMEAN(IJ))
      ENDDO

      IF (LLUNSTR) THEN
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=KIJS,KIJL
              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
              FLHAB = ABS(GTEMP2)
              FLHAB = MIN(FLHAB,USFM(IJ)*COFRM4(M)*DELT)
              FL1(IJ,K,M) = FL1(IJ,K,M) + IOBND(IJ)*SIGN(FLHAB,GTEMP2)
              FL1(IJ,K,M) = MAX(IODP(IJ)*CIREDUC(IJ,K,M)*FL1(IJ,K,M),FLM(IJ,K))
              SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M) + DELTM * MIN(FLMAX(M)-FL1(IJ,K,M),0.0_JWRB)
              FL1(IJ,K,M) = MIN(FL1(IJ,K,M),FLMAX(M))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=KIJS,KIJL
              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
              FLHAB = ABS(GTEMP2)
              FLHAB = MIN(FLHAB,USFM(IJ)*COFRM4(M)*DELT)
              FL1(IJ,K,M) = FL1(IJ,K,M) + SIGN(FLHAB,GTEMP2)
              FL1(IJ,K,M) = MAX(CIREDUC(IJ,K,M)*FL1(IJ,K,M),FLM(IJ,K))
              SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M) + DELTM * MIN(FLMAX(M)-FL1(IJ,K,M),0.0_JWRB)
              FL1(IJ,K,M) = MIN(FL1(IJ,K,M),FLMAX(M))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (LCFLX) THEN
        CALL WNFLUXES (KIJS, KIJL,                              &
     &                 MIJ, RHOWGDFTH,                          &
     &                 CINV,                             &
     &                 SSOURCE, CICOVER,                 &
     &                 PHIWA,                                   &
     &                 EMEAN, F1MEAN, WSWAVE,            &
     &                 WDWAVE, UFRIC, AIRD,&
     &                 NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, &
     &                 NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
     &                 TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, &
     &                 PHIOCD, PHIEPS, PHIAW, &
     &                 .TRUE.)
      ENDIF
! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                      &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, FMEANWS, EMEANWS)

      CALL IMPHFTAIL(KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)


!     UPDATE WINDSEA VARIANCE AND MEAN FREQUENCY IF PASSED TO ATMOSPHERE
!     ------------------------------------------------------------------
      IF (LWFLUX) THEN
        DO IJ=KIJS,KIJL
          IF (EMEANWS(IJ) < WSEMEAN_MIN) THEN
            WSEMEAN(IJ) = WSEMEAN_MIN 
            WSFMEAN(IJ) = 2._JWRB*FR(NFRE)
          ELSE
            WSEMEAN(IJ) = EMEANWS(IJ)
            WSFMEAN(IJ) = FMEANWS(IJ) 
          ENDIF
        ENDDO
      ENDIF


!*    2.6 SET FL1 ON ICE POINTS TO ZERO
!         -----------------------------

      IF (LICERUN .AND. LMASKICE) THEN
        CALL SETICE(KIJS, KIJL, FL1, CICOVER, COSWDIF)
      ENDIF


!*    2.7 SURFACE STOKES DRIFT AND STRAIN IN SEA ICE
!         ------------------------------------------

      CALL STOKESTRN(KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, &
 &                   USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)


END SUBROUTINE IMPLSCH
