#pragma once

#include "field.h"

#include <cuda/std/array>
#include <cuda/std/tuple>

namespace device_ptr_impl {

constexpr int FORTRAN_TO_C_OFFSET = -1;
constexpr int C_TO_FORTRAN_OFFSET = 1;

struct placeholder_t {};
__device__ inline void assert_size_helper2(int const &arg, placeholder_t) {}
__device__ inline void assert_size_helper2(int const &arg, int new_size) {
  assert(arg == new_size);
  __builtin_assume(arg == new_size);
}

template <typename T>
__device__ T const &at(T const *__restrict__ arr, size_t i) {
  return arr[i];
}
template <typename T> __device__ T &at(T *__restrict__ arr, size_t i) {
  return arr[i];
}

template <typename T, int Dim> class DevicePtr {
public:
  template <typename... Args>
  __device__ auto operator()(Args &&...args) const -> decltype(auto) {
    static_assert(sizeof...(args) == Dim,
                  "Access DevicePtr with incorrect number of arguments");
    cuda::std::array<int, Dim> indices{std::forward<decltype(args)>(args)...};
    for (int i = 0; i < Dim; ++i) {
      assert(indices[i] > 0);
      assert(indices[i] <= sizes_[i]);
      assert(indices[i] <= sizes_[i] && indices[i] > 0);
    }

    int ret = 0;
#pragma unroll
    for (int i = 0; i < Dim; ++i) {
      int this_i = indices[i] + FORTRAN_TO_C_OFFSET;
#pragma unroll
      for (int ii = 0; ii < i; ++ii) {
        this_i *= sizes_[ii];
      }
      ret += this_i;
    }

    int total_size = 1;
    for (int i = 0; i < Dim; ++i)
      total_size *= sizes_[i];
    assert(ret < total_size);
    return at(pointer_, ret);
  };

  template <typename... Args>
  __device__ void assert_size(Args &&...args) const {
    static_assert(sizeof...(args) == Dim);
    assert_size_helper(cuda::std::forward_as_tuple(std::forward<Args>(args)...),
                       std::make_index_sequence<sizeof...(Args)>());
  }

  // Create const DevicePtr from non-const
  __device__ DevicePtr(DevicePtr<std::remove_const_t<T>, Dim> const &other)
      : pointer_(other.pointer_), sizes_(other.sizes_) {}
  // Create DevicePtr from Field of same tzpe
  __host__ DevicePtr(Field<T, Dim> const &f)
      : pointer_(f.get(on_device)), sizes_(f.sizes()) {}
  // Create DevicePtr<const T> from Field<T>
  template <typename TT, typename = std::enable_if_t<!std::is_const_v<TT> &&
                                                     std::is_const_v<T>>>
  __host__ DevicePtr(Field<TT, Dim> const &f)
      : pointer_(f.get(on_device)), sizes_(f.sizes()) {}

  friend class DevicePtr<std::add_const_t<T>, Dim>;

private:
  template <typename Args, std::size_t... Is>
  __device__ void assert_size_helper(Args &&args,
                                     std::index_sequence<Is...>) const {
    (assert_size_helper2(sizes_[Is], cuda::std::get<Is>(args)), ...);
  }

  T *pointer_;
  cuda::std::array<int, Dim> sizes_;
};
} // namespace device_ptr_impl

using device_ptr_impl::DevicePtr;
static constexpr __device__ device_ptr_impl::placeholder_t _placeholder{};
