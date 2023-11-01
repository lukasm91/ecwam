#include <cstdio>
#define CUDA_CHECK(e)                                                          \
  {                                                                            \
    cudaError_t err = (e);                                                     \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, cudaGetErrorString(err));                                    \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

#include "device_ptr.h"
#include "field.h"

#include <cassert>
#include <iostream>

using real_t = float;
__host__ __device__ constexpr real_t operator"" _r(long double ld) {
  return ld;
}
__host__ __device__ constexpr real_t operator"" _r(unsigned long long int li) {
  return li;
}

__global__ void
snonlin(int const kijs, int const kijl, int const isnonlin, int mfrstlw,
        int const mlsthg, int const kfrh, int const nfre, int const nang, //
        real_t const dal1, real_t const dal2,
        // Fields
        DevicePtr<real_t, 4> fld, DevicePtr<real_t, 4> sl,
        DevicePtr<real_t const, 4> fl1, DevicePtr<real_t const, 3> wavnum,
        DevicePtr<real_t const, 2> depth, DevicePtr<real_t const, 2> akmean,
        DevicePtr<real_t const, 2> fr, DevicePtr<real_t const, 1> zpifr,
        DevicePtr<int const, 1> ikp, DevicePtr<int const, 1> ikp1,
        DevicePtr<int const, 1> ikm, DevicePtr<int const, 1> ikm1,
        DevicePtr<int const, 2> k1w, DevicePtr<int const, 2> k2w,
        DevicePtr<int const, 2> k11w, DevicePtr<int const, 2> k21w,
        DevicePtr<real_t const, 1> af11, DevicePtr<real_t const, 1> fklap,
        DevicePtr<real_t const, 1> fklap1, DevicePtr<real_t const, 1> fklam,
        DevicePtr<real_t const, 1> fklam1, DevicePtr<int const, 2> inlcoef,
        DevicePtr<real_t const, 2> rnlcoef) {

  int const ichnk = blockIdx.x + 1;
  int const ij = threadIdx.x + kijs;

  real_t enhfr;
  if (isnonlin == 0) {
    enhfr = max(0.75_r * depth(ij, ichnk) * akmean(ij, ichnk), 0.5_r);
    enhfr = 1.0_r +
            (5.5_r / enhfr) * (1.0_r - .833_r * enhfr) * exp(-1.25_r * enhfr);
  } else {
    // TODO
    assert(false);
  }

  ////
  // 2. FREQUENCY LOOP.
  int const mfr1stfr = -mfrstlw + 1;
  int const mfrlstfr = nfre - kfrh + mfr1stfr;

  for (int mc = 1; mc <= mlsthg; ++mc) {
    real_t enh;
    if (isnonlin == 0) {
      enh = enhfr;
    } else {
      // TODO
      assert(false);
    }

    int const mp = ikp(mc - mfrstlw + 1);
    int const mp1 = ikp1(mc - mfrstlw + 1);
    int const mm = ikm(mc - mfrstlw + 1);
    int const mm1 = ikm1(mc - mfrstlw + 1);
    int const ic = inlcoef(1, mc);
    int const ip = inlcoef(2, mc);
    int const ip1 = inlcoef(3, mc);
    int const im = inlcoef(4, mc);
    int const im1 = inlcoef(5, mc);

    real_t const ftail = rnlcoef(1, mc);

    real_t const fklamp = fklap(mc - mfrstlw + 1);
    // real_t const fklamp1 = fklap1(mc - mfrstlw + 1);
    real_t const gw1 = rnlcoef(2, mc);
    real_t const gw2 = rnlcoef(3, mc);
    real_t const gw3 = rnlcoef(4, mc);
    real_t const gw4 = rnlcoef(5, mc);
    real_t const fklampa = rnlcoef(6, mc);
    real_t const fklampb = rnlcoef(7, mc);
    real_t const fklamp2 = rnlcoef(8, mc);
    real_t const fklamp1 = rnlcoef(9, mc);
    real_t const fklapa2 = rnlcoef(10, mc);
    real_t const fklapb2 = rnlcoef(11, mc);
    real_t const fklap12 = rnlcoef(12, mc);
    real_t const fklap22 = rnlcoef(13, mc);

    real_t const fklamm = fklam(mc - mfrstlw + 1);
    // real_t const fklamm1 = fklam1(mc - mfrstlw + 1);
    real_t const gw5 = rnlcoef(14, mc);
    real_t const gw6 = rnlcoef(15, mc);
    real_t const gw7 = rnlcoef(16, mc);
    real_t const gw8 = rnlcoef(17, mc);
    real_t const fklamma = rnlcoef(18, mc);
    real_t const fklammb = rnlcoef(19, mc);
    real_t const fklamm2 = rnlcoef(20, mc);
    real_t const fklamm1 = rnlcoef(21, mc);
    real_t const fklama2 = rnlcoef(22, mc);
    real_t const fklamb2 = rnlcoef(23, mc);
    real_t const fklam12 = rnlcoef(24, mc);
    real_t const fklam22 = rnlcoef(25, mc);

    real_t const ftemp = af11(mc - mfrstlw + 1) * enh;

    if (mc > mfr1stfr && mc < mfrlstfr) {
      //       the interactions for MC are all within the fully resolved
      //       spectral domain

      for (int k = 1; k <= nang; ++k) {
        for (int kh = 1; kh <= 2; ++kh) {
          int const k1 = k1w(k, kh);
          int const k2 = k2w(k, kh);
          int const k11 = k11w(k, kh);
          int const k21 = k21w(k, kh);
          assert(kh == 1 ? k1 == ((k11 + 1 - 1) % nang) + 1 : ((k1 + 1 - 1) % nang) + 1 == k11);
          assert(kh == 2 ? k2 == ((k21 + 1 - 1) % nang) + 1 : ((k2 + 1 - 1) % nang) + 1 == k21);

          ///
          // 2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
          //         DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
          real_t const sap =
              gw1 * fl1(ij, k1, ip, ichnk) + gw2 * fl1(ij, k11, ip, ichnk) +
              gw3 * fl1(ij, k1, ip1, ichnk) + gw4 * fl1(ij, k11, ip1, ichnk);
          real_t const sam =
              gw5 * fl1(ij, k2, im, ichnk) + gw6 * fl1(ij, k21, im, ichnk) +
              gw7 * fl1(ij, k2, im1, ichnk) + gw8 * fl1(ij, k21, im1, ichnk);
          //// not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC
          ///)*FTAIL
          real_t const fij = fl1(ij, k, ic, ichnk);
          real_t fad1 = fij * (sap + sam);
          real_t const fad2 = fad1 - 2_r * sap * sam;
          fad1 = fad1 + fad2;
          real_t const fcen = ftemp * fij;
          real_t const ad = fad2 * fcen;
          real_t const delad = fad1 * ftemp;
          real_t const delap = (fij - 2_r * sam) * dal1 * fcen;
          real_t const delam = (fij - 2_r * sap) * dal2 * fcen;

          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 9>
              sl_assign = {{
                  {true, k, mc, -2_r * ad},
                  {true, k2, mm, +ad * fklamm1},
                  {true, k21, mm, +ad * fklamm2},
                  {true, k2, mm1, +ad * fklamma},
                  {true, k21, mm1, +ad * fklammb},
                  {true, k1, mp, +ad * fklamp1},
                  {true, k11, mp, +ad * fklamp2},
                  {true, k1, mp1, +ad * fklampa},
                  {true, k11, mp1, +ad * fklampb},
              }};
          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 9>
              fld_assign = {{
                  {true, k, mc, -2_r * delad},
                  {true, k2, mm, +delam * fklam12},
                  {true, k21, mm, +delam * fklam22},
                  {true, k2, mm1, +delam * fklama2},
                  {true, k21, mm1, +delam * fklamb2},
                  {true, k1, mp, +delap * fklap12},
                  {true, k11, mp, +delap * fklap22},
                  {true, k1, mp1, +delap * fklapa2},
                  {true, k11, mp1, +delap * fklapb2},
              }};

          int i;
          real_t sl_values[sl_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl_values[i] = sl(ij, i1, i2, ichnk) + v;
            ++i;
          }
          real_t fld_values[fld_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld_values[i] = fld(ij, i1, i2, ichnk) + v;
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl(ij, i1, i2, ichnk) = sl_values[i];
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld(ij, i1, i2, ichnk) = fld_values[i];
            ++i;
          }
        }
      }
    } else if (mc >= mfrlstfr) {
      for (int k = 1; k <= nang; ++k) {
        for (int kh = 1; kh <= 2; ++kh) {
          int const k1 = k1w(k, kh);
          int const k2 = k2w(k, kh);
          int const k11 = k11w(k, kh);
          int const k21 = k21w(k, kh);

          real_t const sap =
              gw1 * fl1(ij, k1, ip, ichnk) + gw2 * fl1(ij, k11, ip, ichnk) +
              gw3 * fl1(ij, k1, ip1, ichnk) + gw4 * fl1(ij, k11, ip1, ichnk);
          real_t const sam =
              gw5 * fl1(ij, k2, im, ichnk) + gw6 * fl1(ij, k21, im, ichnk) +
              gw7 * fl1(ij, k2, im1, ichnk) + gw8 * fl1(ij, k21, im1, ichnk);
          real_t const fij = fl1(ij, k, ic, ichnk) * ftail;
          real_t fad1 = fij * (sap + sam);
          real_t const fad2 = fad1 - 2_r * sap * sam;
          fad1 = fad1 + fad2;
          real_t const fcen = ftemp * fij;
          real_t const ad = fad2 * fcen;
          real_t const delad = fad1 * ftemp;
          real_t const delap = (fij - 2_r * sam) * dal1 * fcen;
          real_t const delam = (fij - 2_r * sap) * dal2 * fcen;

          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 9>
              sl_assign = {{
                  {true, k2, mm, +ad * fklamm1},
                  {true, k21, mm, +ad * fklamm2},
                  {mm1 <= nfre, k2, mm1, +ad * fklamma},
                  {mm1 <= nfre, k21, mm1, +ad * fklammb},
                  {mc <= nfre, k, mc, -2_r * ad},
                  {mp <= nfre, k1, mp, +ad * fklamp1},
                  {mp <= nfre, k11, mp, +ad * fklamp2},
                  {mp1 <= nfre, k1, mp1, +ad * fklampa},
                  {mp1 <= nfre, k11, mp1, +ad * fklampb},
              }};
          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 9>
              fld_assign = {{
                  {true, k2, mm, +delam * fklam12},
                  {true, k21, mm, +delam * fklam22},
                  {mm1 <= nfre, k2, mm1, +delam * fklama2},
                  {mm1 <= nfre, k21, mm1, +delam * fklamb2},
                  {mc <= nfre, k, mc, -2_r * delad},
                  {mc <= nfre && mp <= nfre, k1, mp, +delap * fklap12},
                  {mp <= nfre, k11, mp, +delap * fklap22},
                  {mp1 <= nfre, k1, mp1, +delap * fklapa2},
                  {mp1 <= nfre, k11, mp1, +delap * fklapb2},
              }};

          int i;
          real_t sl_values[sl_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl_values[i] = sl(ij, i1, i2, ichnk) + v;
            ++i;
          }
          real_t fld_values[fld_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld_values[i] = fld(ij, i1, i2, ichnk) + v;
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl(ij, i1, i2, ichnk) = sl_values[i];
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld(ij, i1, i2, ichnk) = fld_values[i];
            ++i;
          }
        }
      }

    } else {
      for (int k = 1; k <= nang; ++k) {
        for (int kh = 1; kh <= 2; ++kh) {
          int const k1 = k1w(k, kh);
          int const k2 = k2w(k, kh);
          int const k11 = k11w(k, kh);
          int const k21 = k21w(k, kh);

          real_t const sap =
              gw1 * fl1(ij, k1, ip, ichnk) + gw2 * fl1(ij, k11, ip, ichnk) +
              gw3 * fl1(ij, k1, ip1, ichnk) + gw4 * fl1(ij, k11, ip1, ichnk);
          real_t const sam =
              gw5 * fl1(ij, k2, im, ichnk) + gw6 * fl1(ij, k21, im, ichnk) +
              gw7 * fl1(ij, k2, im1, ichnk) + gw8 * fl1(ij, k21, im1, ichnk);
          real_t const fij = fl1(ij, k, ic, ichnk) * ftail;
          real_t fad1 = fij * (sap + sam);
          real_t const fad2 = fad1 - 2_r * sap * sam;
          fad1 = fad1 + fad2;
          real_t const fcen = ftemp * fij;
          real_t const ad = fad2 * fcen;
          real_t const delad = fad1 * ftemp;
          real_t const delap = (fij - 2_r * sam) * dal1 * fcen;
          real_t const delam = (fij - 2_r * sap) * dal2 * fcen;

          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 7>
              sl_assign = {{
                  {mm >= 1, k2, mm1, +ad * fklamma},
                  {mm >= 1, k21, mm1, +ad * fklammb},
                  {true, k, mc, -2_r * ad},
                  {true, k1, mp, +ad * fklamp1},
                  {true, k11, mp, +ad * fklamp2},
                  {true, k1, mp1, +ad * fklampa},
                  {true, k11, mp1, +ad * fklampb},
              }};
          cuda::std::array<cuda::std::tuple<bool, int, int, real_t>, 7>
              fld_assign = {{
                  {mm >= 1, k2, mm1, +delam * fklama2},
                  {mm >= 1, k21, mm1, +delam * fklamb2},
                  {true, k, mc, -2_r * delad},
                  {true, k1, mp, +delap * fklap12},
                  {true, k11, mp, +delap * fklap22},
                  {true, k1, mp1, +delap * fklapa2},
                  {true, k11, mp1, +delap * fklapb2},
              }};
          int i;
          real_t sl_values[sl_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl_values[i] = sl(ij, i1, i2, ichnk) + v;
            ++i;
          }
          real_t fld_values[fld_assign.size()] = {};
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld_values[i] = fld(ij, i1, i2, ichnk) + v;
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : sl_assign) {
            if (b)
              sl(ij, i1, i2, ichnk) = sl_values[i];
            ++i;
          }
          i = 0;
#pragma unroll
          for (auto &&[b, i1, i2, v] : fld_assign) {
            if (b)
              fld(ij, i1, i2, ichnk) = fld_values[i];
            ++i;
          }
        }
      }
    }
  }
}

extern "C" {
void snonlin_cuda_ext(
    int const kijs, int const kijl, int const nchnk, int const nang,
    int const nfre, real_t *fl1_ptr, real_t *fld_ptr, real_t *sl_ptr,
    real_t *wavnum_ptr, real_t *depth_ptr, real_t *akmean_ptr,
    real_t const *fr_ptr, real_t const *zpifr_ptr, int const *ikp_ptr,
    int const *ikp1_ptr, int const *ikm_ptr, int const *ikm1_ptr,
    int const *k1w_ptr, int const *k2w_ptr, int const *k11w_ptr,
    int const *k21w_ptr, real_t const *af11_ptr, real_t const *fklap_ptr,
    real_t const *fklap1_ptr, real_t const *fklam_ptr, real_t const *fklam1_ptr,
    real_t const dal1, real_t const dal2, int const mfrstlw, int const mlsthg,
    int const kfrh, int const *inlcoef_ptr, real_t const *rnlcoef_ptr,
    int const isnonlin) {

  int const nproma = kijl - kijs + 1;

  auto fl1 = Field(fl1_ptr, nproma, nang, nfre, nchnk);
  auto fld = Field(fld_ptr, nproma, nang, nfre, nchnk);
  auto sl = Field(sl_ptr, nproma, nang, nfre, nchnk);
  auto wavnum = Field(wavnum_ptr, nproma, nfre, nchnk);
  auto depth = Field(depth_ptr, nproma, nchnk);
  auto akmean = Field(akmean_ptr, nproma, nchnk);
  auto fr = Field(static_cast<real_t *>(nullptr), nproma, nchnk);
  auto zpifr = Field(static_cast<real_t *>(nullptr), nfre);
  auto ikp = Field(ikp_ptr, mlsthg - mfrstlw + 1);
  auto ikp1 = Field(ikp1_ptr, mlsthg - mfrstlw + 1);
  auto ikm = Field(ikm_ptr, mlsthg - mfrstlw + 1);
  auto ikm1 = Field(ikm1_ptr, mlsthg - mfrstlw + 1);
  auto k1w = Field(k1w_ptr, nang, 2);
  auto k2w = Field(k2w_ptr, nang, 2);
  auto k11w = Field(k11w_ptr, nang, 2);
  auto k21w = Field(k21w_ptr, nang, 2);
  auto af11 = Field(af11_ptr, mlsthg - mfrstlw + 1);
  auto fklap = Field(fklap_ptr, mlsthg - mfrstlw + 1);
  auto fklap1 = Field(fklap1_ptr, mlsthg - mfrstlw + 1);
  auto fklam = Field(fklam_ptr, mlsthg - mfrstlw + 1);
  auto fklam1 = Field(fklam1_ptr, mlsthg - mfrstlw + 1);
  auto inlcoef = Field(inlcoef_ptr, ninl, mlsthg);
  auto rnlcoef = Field(rnlcoef_ptr, nrnl, mlsthg);

  assert(isnonlin == 0);
  snonlin<<<nchnk, nproma>>>(
      kijs, kijl, isnonlin, mfrstlw, mlsthg, kfrh, nfre, nang, //
      dal1, dal2,                                              //
      fld, sl, fl1, wavnum, depth, akmean, fr, zpifr, ikp, ikp1, ikm, ikm1, k1w,
      k2w, k11w, k21w, af11, fklap, fklap1, fklam, fklam1, inlcoef, rnlcoef);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
}
}
