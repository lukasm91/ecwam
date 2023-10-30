#pragma once

#include <cassert>
#include <cuda/std/array>
#include <memory>

#include "openacc.h"

constexpr int nrnl = 25;
constexpr int ninl = 5;

namespace field_impl {
template <typename> struct Temporary {};
struct OnDevice {};
struct OnHost {};

template <typename T, int Dim> class Field {
public:
  // not owning pointer
  template <typename... Sizes>
  explicit Field(T *pointer, Sizes &&...sizes)
      : pointer_(pointer, [](T *ptr) {}),
        sizes_({int(std::forward<Sizes>(sizes))...}) {
    // Although it is valid to have a Fortran array allocated with size = 0,
    // calling OpenACC `copyin` will not allocate it on the device, and thus the
    // array will not be "present" on the device (e.g. `pg_edgeidx, pg_edgeblk,
    // pg_vertidx, pg_exdist` in some cases for certain ranks).
    if (nelems() > 0) {
      assert(
          acc_is_present(const_cast<std::remove_const_t<T> *>(pointer_.get()),
                         nelems() * sizeof(T)));
    }
  }
  // owning pointer
  template <typename... Sizes>
  explicit Field(std::shared_ptr<T[]> &&pointer, Sizes &&...sizes)
      : pointer_(std::move(pointer)),
        sizes_({int(std::forward<Sizes>(sizes))...}) {
    // Although it is valid to have a Fortran array allocated with size = 0,
    // calling OpenACC `copyin` will not allocate it on the device, and thus the
    // array will not be "present" on the device (e.g. `pg_edgeidx, pg_edgeblk,
    // pg_vertidx, pg_exdist` in some cases for certain ranks).
    if (nelems() > 0) {
      assert(
          acc_is_present(const_cast<std::remove_const_t<T> *>(pointer_.get()),
                         nelems() * sizeof(T)));
    }
  }
  // temporary pointer
  template <typename... Sizes>
  explicit Field(Temporary<T>, Sizes &&...sizes)
      : sizes_({std::forward<Sizes>(sizes)...}) {
    T *ptr = new T[nelems()];
    acc_create(ptr, nelems() * sizeof(T));
    pointer_ =
        std::shared_ptr<T[]>{ptr, [nbytes = nelems() * sizeof(T)](T *ptr) {
                               acc_delete(ptr, nbytes);
                               delete[] ptr;
                             }};
    // Although it is valid to have a Fortran array allocated with size = 0,
    // calling OpenACC `copyin` will not allocate it on the device, and thus the
    // array will not be "present" on the device (e.g. `pg_edgeidx, pg_edgeblk,
    // pg_vertidx, pg_exdist` in some cases for certain ranks).
    if (nelems() > 0) {
      assert(
          acc_is_present(const_cast<std::remove_const_t<T> *>(pointer_.get()),
                         nelems() * sizeof(T)));
    }
  }
  // creating const field from non-const
  Field(Field<std::remove_const_t<T>, Dim> const &other)
      : pointer_(other.pointer_), sizes_(other.sizes_) {}
  Field() = default;

  friend class Field<std::add_const_t<T>, Dim>;

  T *get(OnDevice) const {
    return static_cast<T *>(
        acc_deviceptr(const_cast<std::remove_const_t<T> *>(get(OnHost()))));
  }
  T *get(OnHost) const { return pointer_.get(); }

  cuda::std::array<int, Dim> sizes() const { return sizes_; }
  int nelems() const {
    int ret = 1;
    for (const auto &s : sizes_)
      ret *= s;
    return ret;
  }

  using underlying_type = T;
  static constexpr int dimension = Dim;

private:
  std::shared_ptr<T[]> pointer_;
  cuda::std::array<int, Dim> sizes_;
};
template <typename T, typename... Sizes>
Field(T *pointer, Sizes &&...) -> Field<T, sizeof...(Sizes)>;
template <typename T, typename... Sizes>
Field(Temporary<T>, Sizes &&...) -> Field<T, sizeof...(Sizes)>;
template <typename T, typename... Sizes>
Field(std::shared_ptr<T[]> &&pointer, Sizes &&...sizes)
    -> Field<T, sizeof...(Sizes)>;
} // namespace field_impl

using field_impl::Field;

template <typename T>
static constexpr auto temporary = field_impl::Temporary<T>{};
static constexpr auto on_device = field_impl::OnDevice{};
static constexpr auto on_host = field_impl::OnHost{};
