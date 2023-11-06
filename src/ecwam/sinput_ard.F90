! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
  SUBROUTINE SINPUT_ARD (IDX, NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE, UFRIC, Z0M, COSWDIF,  &
  & SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS)
  !$loki routine seq
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.
    
    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
    !        *Z0M* - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
    !       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLCAPCHNK, LLNORMAGAM
    USE YOWFRED, ONLY: FR, TH, DFIM, COSTH, SINTH, ZPIFR, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE, NANG_PARAM
    USE YOWPCONS, ONLY: G, GM1, EPSMIN, EPSUS, ZPI
    USE YOWPHYS, ONLY: ZALP, TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2, RN1_RN, RNU, RNUM, SWELLF, SWELLF2, SWELLF3, SWELLF4,  &
    & SWELLF5, SWELLF6, SWELLF7, SWELLF7M1, Z0RAT, Z0TUBMAX, ABMIN, ABMAX
    USE YOWTEST, ONLY: IU06
    USE YOWTABL, ONLY: IAB, SWELLFT
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    #include "wsigstar.intfb.h"
    INTEGER(KIND=JWIM), INTENT(IN) :: NGST
    LOGICAL, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, IDX
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN) :: WDWAVE, WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN) :: RAORW, WSTAR, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: XLLWS
    
    
    INTEGER(KIND=JWIM) :: K, M, IND, IGST
    
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
    REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
    REAL(KIND=JWRB) :: FU, FUD, NU_AIR, SMOOTH, HFTSWELLF6, Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
    REAL(KIND=JWRB) :: ARG, DELABM1
    REAL(KIND=JWRB) :: TAUPX, TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    REAL(KIND=JWRB) :: CONST, SIG, SIG2, COEF, COEF5, DFIM_SIG2
    REAL(KIND=JWRB) :: CONSTF, CONST11, CONST22
    REAL(KIND=JWRB) :: Z0VIS, Z0NOZ, FWW
    REAL(KIND=JWRB) :: PVISC, PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N
    REAL(KIND=JWRB) :: UORBT
    REAL(KIND=JWRB) :: AORB
    REAL(KIND=JWRB) :: TEMP
    REAL(KIND=JWRB) :: RE
    REAL(KIND=JWRB) :: RE_C
    REAL(KIND=JWRB) :: ZORB
    REAL(KIND=JWRB) :: CNSN, SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG, SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR, AIRD_PVISC
    REAL(KIND=JWRB), DIMENSION(2) :: XSTRESS, YSTRESS, FLP, SLP
    REAL(KIND=JWRB), DIMENSION(2) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
    REAL(KIND=JWRB), DIMENSION(2) :: UCNZALPD
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: DSTAB1, TEMP1, TEMP2
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM, 2) :: GAM0, DSTAB
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM) :: COSLP
    
    LOGICAL :: LTAUWSHELTER
    ! ----------------------------------------------------------------------
    
    
    AVG_GST = 1.0_JWRB / NGST
    CONST1 = BETAMAXOXKAPPA2
    CONSTN = DELTH / (XKAPPA*ZPI)
    
    ABS_TAUWSHELTER = ABS(TAUWSHELTER)
    IF (ABS_TAUWSHELTER == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
    ELSE
      LTAUWSHELTER = .true.
    END IF
    
    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    IF (NGST > 1) CALL WSIGSTAR(WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
  
    IF (LLNORMAGAM) THEN
      CSTRNFAC = CONSTN*RNFAC / RAORW
    END IF
    ! ----------------------------------------------------------------------
    
    IF (LLSNEG) THEN
      !!!!  only for the negative sinput
      NU_AIR = RNU
      FACM1_NU_AIR = 4.0_JWRB / NU_AIR
      
      FAC_NU_AIR = RNUM
      
      FU = ABS(SWELLF3)
      FUD = SWELLF2
      DELABM1 = REAL(IAB) / (ABMAX - ABMIN)
      
      
      !       computation of Uorb and Aorb
      UORBT = EPSMIN
      AORB = EPSMIN
      
      DO M=1,NFRE
        SIG = ZPIFR(M)
        SIG2 = SIG**2
        DFIM_SIG2 = DFIM(M)*SIG2
        
        K = 1
        TEMP = FL1(IDX, K, M)
        DO K=2,NANG
          TEMP = TEMP + FL1(IDX, K, M)
        END DO
        
        UORBT = UORBT + DFIM_SIG2*TEMP
        AORB = AORB + DFIM(M)*TEMP
      END DO
      
      UORBT = 2.0_JWRB*SQRT(UORBT)          ! this is the significant orbital amplitude
      AORB = 2.0_JWRB*SQRT(AORB)          ! this 1/2 Hs
      RE = FACM1_NU_AIR*UORBT*AORB          ! this is the Reynolds number
      Z0VIS = FAC_NU_AIR / MAX(UFRIC, 0.0001_JWRB)
      Z0TUB = Z0RAT*MIN(Z0TUBMAX, Z0M)
      Z0NOZ = MAX(Z0VIS, Z0TUB)
      ZORB = AORB / Z0NOZ
      
      !         compute fww
      XI = (LOG10(MAX(ZORB, 3.0_JWRB)) - ABMIN)*DELABM1
      IND = MIN(IAB - 1, INT(XI))
      DELI1 = MIN(1.0_JWRB, XI - REAL(IND, kind=JWRB))
      DELI2 = 1.0_JWRB - DELI1
      FWW = SWELLFT(IND)*DELI2 + SWELLFT(IND + 1)*DELI1
      TEMP2 = FWW*UORBT
      
      !       Define the critical Reynolds number
      IF (SWELLF6 == 1.0_JWRB) THEN
        RE_C = SWELLF4
      ELSE
        HFTSWELLF6 = 1.0_JWRB - SWELLF6
        RE_C = SWELLF4*(2.0_JWRB / AORB)**HFTSWELLF6
      END IF
      
      !       Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7 > 0.0_JWRB) THEN
        SMOOTH = 0.5_JWRB*TANH((RE - RE_C)*SWELLF7M1)
        PTURB = 0.5_JWRB + SMOOTH
        PVISC = 0.5_JWRB - SMOOTH
      ELSE
        IF (RE <= RE_C) THEN
          PTURB = 0.0_JWRB
          PVISC = 0.5_JWRB
        ELSE
          PTURB = 0.5_JWRB
          PVISC = 0.0_JWRB
        END IF
      END IF
      
      AIRD_PVISC = PVISC*RAORW
      
    END IF
    
    
    
    ! Initialisation
    
    IF (NGST == 1) THEN
      USTP(1) = UFRIC
    ELSE IF (NGST == 2) THEN
      USTP(1) = UFRIC*(1.0_JWRB + SIG_N)
      USTP(2) = UFRIC*(1.0_JWRB - SIG_N)
    END IF
    
    DO IGST=1,NGST
      USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
    END DO
    
    IF (LTAUWSHELTER) THEN
      DO IGST=1,NGST
        XSTRESS(IGST) = 0.0_JWRB
        YSTRESS(IGST) = 0.0_JWRB
        USG2(IGST) = USTP(IGST)**2
        TAUX(IGST) = USG2(IGST)*SIN(WDWAVE)
        TAUY(IGST) = USG2(IGST)*COS(WDWAVE)
      END DO
      
      ROGOROAIR = G / RAORW
      
    ELSE
      DO K=1,NANG
        COSLP(K) = COSWDIF(IDX, K)
      END DO
    END IF
    
    
    !*    2. MAIN LOOP OVER FREQUENCIES.
    !        ---------------------------
    
    IF (.not.LLNORMAGAM) THEN
      GAMNORMA(:) = 1.0_JWRB
    END IF
    
    IF (.not.LLSNEG) THEN
      DSTAB(:, :) = 0.0_JWRB
    END IF
    
    DO M=1,NFRE
      SIG = ZPIFR(M)
      SIG2 = SIG**2
      CONST = SIG*CONST1
      
      COEF = -SWELLF*16._JWRB*SIG2 / G
      COEF5 = -SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG)
      
      IF (LTAUWSHELTER) THEN
        DO IGST=1,NGST
          TAUPX = TAUX(IGST) - ABS_TAUWSHELTER*XSTRESS(IGST)
          TAUPY = TAUY(IGST) - ABS_TAUWSHELTER*YSTRESS(IGST)
          USDIRP(IGST) = ATAN2(TAUPX, TAUPY)
          USTP(IGST) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
          USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
        END DO
        
        CONSTF = ROGOROAIR*CINV(IDX, M)*DFIM(M)
      END IF
      
      
      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------
      
      DO IGST=1,NGST
        UCN(IGST) = USTP(IGST)*CINV(IDX, M)
        UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP)
      END DO
      ZCN = LOG(WAVNUM(IDX, M)*Z0M)
      CNSN = CONST*RAORW
      
      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------
      
      DO K=1,NANG
        XLLWS(IDX, K, M) = 0.0_JWRB
      END DO
      
      IF (LLNORMAGAM) THEN
        XNGAMCONST = CSTRNFAC*XK2CG(IDX, M)
      END IF
      
      IF (LLSNEG) THEN
        !       SWELL DAMPING:
        DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IDX, M)
        TEMP1 = COEF*RAORW
      END IF
      
      DO IGST=1,NGST
        DO K=1,NANG
          IF (LTAUWSHELTER) THEN
            COSLP(K) = COS(TH(K) - USDIRP(IGST))
          END IF
          GAM0(K, IGST) = 0.0_JWRB
          IF (COSLP(K) > 0.01_JWRB) THEN
            X = COSLP(K)*UCN(IGST)
            ZLOG = ZCN + UCNZALPD(IGST) / COSLP(K)
            IF (ZLOG < 0.0_JWRB) THEN
              ZLOG2X = ZLOG*ZLOG*X
              GAM0(K, IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
              XLLWS(IDX, K, M) = 1.0_JWRB
            END IF
          END IF
        END DO
        
        IF (LLNORMAGAM) THEN
          
          SUMF = 0.0_JWRB
          SUMFSIN2 = 0.0_JWRB
          DO K=1,NANG
            SUMF = SUMF + GAM0(K, IGST)*FL1(IDX, K, M)
            SUMFSIN2 = SUMFSIN2 + GAM0(K, IGST)*FL1(IDX, K, M)*SINWDIF2(IDX, K)
          END DO
          
          ZNZ = XNGAMCONST*USTPM1(IGST)
          GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
          
        END IF
        
        IF (LLSNEG) THEN
          DO K=1,NANG
            DSTAB2 = TEMP1*(TEMP2 + (FU + FUD*COSLP(K))*USTP(IGST))
            DSTAB(K, IGST) = DSTAB1 + PTURB*DSTAB2
          END DO
        END IF
      END DO
      
      
      !*    2.2 UPDATE THE SHELTERING STRESS (in any),
      !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ---------------------------------------------------------
      
      DO K=1,NANG
        
        DO IGST=1,NGST
          ! SLP: only the positive contributions
          SLP(IGST) = GAM0(K, IGST)*GAMNORMA(IGST)
          FLP(IGST) = SLP(IGST) + DSTAB(K, IGST)
        END DO
        
        DO IGST=1,NGST
          SLP(IGST) = SLP(IGST)*FL1(IDX, K, M)
        END DO
        
        IF (LTAUWSHELTER) THEN
          CONST11 = CONSTF*SINTH(K)
          CONST22 = CONSTF*COSTH(K)
          DO IGST=1,NGST
            XSTRESS(IGST) = XSTRESS(IGST) + SLP(IGST)*CONST11
            YSTRESS(IGST) = YSTRESS(IGST) + SLP(IGST)*CONST22
          END DO
        END IF
        
        IGST = 1
        SLP_AVG = SLP(IGST)
        FLP_AVG = FLP(IGST)
        DO IGST=2,NGST
          SLP_AVG = SLP_AVG + SLP(IGST)
          FLP_AVG = FLP_AVG + FLP(IGST)
        END DO
        
        SPOS(IDX, K, M) = AVG_GST*SLP_AVG
        
        FLD(IDX, K, M) = AVG_GST*FLP_AVG
        SL(IDX, K, M) = FLD(IDX, K, M)*FL1(IDX, K, M)
        
      END DO
      
    END DO
    
    ! END LOOP OVER FREQUENCIES
    
    
    
  END SUBROUTINE SINPUT_ARD
