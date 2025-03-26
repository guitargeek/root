#ifndef MyROOT_RVEC
#define MyROOT_RVEC

#define R__RVEC_NODISCARD [[nodiscard]]

#include <Rtypes.h> // R__CLING_PTRCHECK
#include <TError.h> // R__ASSERT

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits> // for numeric_limits
#include <memory> // uninitialized_value_construct
#include <new>
#include <numeric> // for inner_product
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace MyROOT {

namespace VecOps {
template<typename T>
class RVec;
}

namespace Internal {
namespace VecOps {

template<typename T>
using RVec = MyROOT::VecOps::RVec<T>;

// clang-format off
template <typename>
struct IsRVec : std::false_type {};

template <typename T>
struct IsRVec<MyROOT::VecOps::RVec<T>> : std::true_type {};
// clang-format on

constexpr bool All(const bool *vals, std::size_t size)
{
   for (auto i = 0u; i < size; ++i)
      if (!vals[i])
         return false;
   return true;
}

template <typename... T>
std::size_t GetVectorsSize(const std::string &id, const RVec<T> &... vs)
{
   constexpr const auto nArgs = sizeof...(T);
   const std::size_t sizes[] = {vs.size()...};
   if (nArgs > 1) {
      for (auto i = 1UL; i < nArgs; i++) {
         if (sizes[0] == sizes[i])
            continue;
         std::string msg(id);
         msg += ": input RVec instances have different lengths!";
         throw std::runtime_error(msg);
      }
   }
   return sizes[0];
}

template <typename F, typename... RVecs>
auto MapImpl(F &&f, RVecs &&... vs) -> RVec<decltype(f(vs[0]...))>
{
   const auto size = GetVectorsSize("Map", vs...);
   RVec<decltype(f(vs[0]...))> ret(size);

   for (auto i = 0UL; i < size; i++)
      ret[i] = f(vs[i]...);

   return ret;
}

template <typename Tuple_t, std::size_t... Is>
auto MapFromTuple(Tuple_t &&t, std::index_sequence<Is...>)
   -> decltype(MapImpl(std::get<std::tuple_size<Tuple_t>::value - 1>(t), std::get<Is>(t)...))
{
   constexpr const auto tupleSizeM1 = std::tuple_size<Tuple_t>::value - 1;
   return MapImpl(std::get<tupleSizeM1>(t), std::get<Is>(t)...);
}

/// Return the next power of two (in 64-bits) that is strictly greater than A.
/// Return zero on overflow.
inline uint64_t NextPowerOf2(uint64_t A)
{
   A |= (A >> 1);
   A |= (A >> 2);
   A |= (A >> 4);
   A |= (A >> 8);
   A |= (A >> 16);
   A |= (A >> 32);
   return A + 1;
}

/// This is all the stuff common to all SmallVectors.
class R__CLING_PTRCHECK(off) SmallVectorBase {
public:
   // This limits the maximum size of an RVec<char> to ~4GB but we don't expect this to ever be a problem,
   // and we prefer the smaller Size_T to reduce the size of each RVec object.
   using Size_T = int32_t;

protected:
   void *fBeginX;
   /// Always >= 0.
   // Type is signed only for consistency with fCapacity.
   Size_T fSize = 0;
   /// Always >= -1. fCapacity == -1 indicates the RVec is in "memory adoption" mode.
   Size_T fCapacity;

   /// The maximum value of the Size_T used.
   static constexpr size_t SizeTypeMax() { return std::numeric_limits<Size_T>::max(); }

   SmallVectorBase() = delete;
   SmallVectorBase(void *FirstEl, size_t TotalCapacity) : fBeginX(FirstEl), fCapacity(TotalCapacity) {}

   /// This is an implementation of the grow() method which only works
   /// on POD-like data types and is out of line to reduce code duplication.
   /// This function will report a fatal error if it cannot increase capacity.
   void grow_pod(void *FirstEl, size_t MinSize, size_t TSize);

   /// Report that MinSize doesn't fit into this vector's size type. Throws
   /// std::length_error or calls report_fatal_error.
   static void report_size_overflow(size_t MinSize);
   /// Report that this vector is already at maximum capacity. Throws
   /// std::length_error or calls report_fatal_error.
   static void report_at_maximum_capacity();

   /// If false, the RVec is in "memory adoption" mode, i.e. it is acting as a view on a memory buffer it does not own.
   bool Owns() const { return fCapacity != -1; }

public:
   size_t size() const { return fSize; }
   size_t capacity() const noexcept { return Owns() ? fCapacity : fSize; }

   R__RVEC_NODISCARD bool empty() const { return !fSize; }

   /// Set the array size to \p N, which the current array must have enough
   /// capacity for.
   ///
   /// This does not construct or destroy any elements in the vector.
   ///
   /// Clients can use this in conjunction with capacity() to write past the end
   /// of the buffer when they know that more elements are available, and only
   /// update the size later. This avoids the cost of value initializing elements
   /// which will only be overwritten.
   void set_size(size_t N)
   {
      if (N > capacity()) {
         throw std::runtime_error("Setting size to a value greater than capacity.");
      }
      fSize = N;
   }
};

/// Used to figure out the offset of the first element of an RVec
template <class T>
struct SmallVectorAlignmentAndSize {
   alignas(SmallVectorBase) char Base[sizeof(SmallVectorBase)];
   alignas(T) char FirstEl[sizeof(T)];
};

/// This is the part of SmallVectorTemplateBase which does not depend on whether the type T is a POD.
template <typename T>
class R__CLING_PTRCHECK(off) SmallVectorTemplateCommon : public SmallVectorBase {
   using Base = SmallVectorBase;

   /// Find the address of the first element.  For this pointer math to be valid
   /// with small-size of 0 for T with lots of alignment, it's important that
   /// SmallVectorStorage is properly-aligned even for small-size of 0.
   void *getFirstEl() const
   {
      return const_cast<void *>(reinterpret_cast<const void *>(reinterpret_cast<const char *>(this) +
                                                               offsetof(SmallVectorAlignmentAndSize<T>, FirstEl)));
   }
   // Space after 'FirstEl' is clobbered, do not add any instance vars after it.

protected:
   SmallVectorTemplateCommon(size_t Size) : Base(getFirstEl(), Size) {}

   void grow_pod(size_t MinSize, size_t TSize) { Base::grow_pod(getFirstEl(), MinSize, TSize); }

   /// Return true if this is a smallvector which has not had dynamic
   /// memory allocated for it.
   bool isSmall() const { return this->fBeginX == getFirstEl(); }

   /// Put this vector in a state of being small.
   void resetToSmall()
   {
      this->fBeginX = getFirstEl();
      // from the original LLVM implementation:
      // FIXME: Setting fCapacity to 0 is suspect.
      this->fSize = this->fCapacity = 0;
   }

public:
   // note that fSize is a _signed_ integer, but we expose it as an unsigned integer for consistency with STL containers
   // as well as backward-compatibility
   using size_type = size_t;
   using difference_type = ptrdiff_t;
   using value_type = T;
   using iterator = T *;
   using const_iterator = const T *;

   using const_reverse_iterator = std::reverse_iterator<const_iterator>;
   using reverse_iterator = std::reverse_iterator<iterator>;

   using reference = T &;
   using const_reference = const T &;
   using pointer = T *;
   using const_pointer = const T *;

   using Base::capacity;
   using Base::empty;
   using Base::size;

   // forward iterator creation methods.
   iterator begin() noexcept { return (iterator)this->fBeginX; }
   const_iterator begin() const noexcept { return (const_iterator)this->fBeginX; }
   const_iterator cbegin() const noexcept { return (const_iterator)this->fBeginX; }
   iterator end() noexcept { return begin() + size(); }
   const_iterator end() const noexcept { return begin() + size(); }
   const_iterator cend() const noexcept { return begin() + size(); }

   // reverse iterator creation methods.
   reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
   const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
   const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
   reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
   const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }
   const_reverse_iterator crend() const noexcept { return const_reverse_iterator(begin()); }

   size_type size_in_bytes() const { return size() * sizeof(T); }
   size_type max_size() const noexcept { return std::min(this->SizeTypeMax(), size_type(-1) / sizeof(T)); }

   size_t capacity_in_bytes() const { return capacity() * sizeof(T); }

   /// Return a pointer to the vector's buffer, even if empty().
   pointer data() noexcept { return pointer(begin()); }
   /// Return a pointer to the vector's buffer, even if empty().
   const_pointer data() const noexcept { return const_pointer(begin()); }

   reference front()
   {
      if (empty()) {
         throw std::runtime_error("`front` called on an empty RVec");
      }
      return begin()[0];
   }

   const_reference front() const
   {
      if (empty()) {
         throw std::runtime_error("`front` called on an empty RVec");
      }
      return begin()[0];
   }

   reference back()
   {
      if (empty()) {
         throw std::runtime_error("`back` called on an empty RVec");
      }
      return end()[-1];
   }

   const_reference back() const
   {
      if (empty()) {
         throw std::runtime_error("`back` called on an empty RVec");
      }
      return end()[-1];
   }
};

/// SmallVectorTemplateBase<TriviallyCopyable = false> - This is where we put
/// method implementations that are designed to work with non-trivial T's.
///
/// We approximate is_trivially_copyable with trivial move/copy construction and
/// trivial destruction. While the standard doesn't specify that you're allowed
/// copy these types with memcpy, there is no way for the type to observe this.
/// This catches the important case of std::pair<POD, POD>, which is not
/// trivially assignable.
template <typename T, bool = (std::is_trivially_copy_constructible<T>::value) &&
                             (std::is_trivially_move_constructible<T>::value) &&
                             std::is_trivially_destructible<T>::value>
class R__CLING_PTRCHECK(off) SmallVectorTemplateBase : public SmallVectorTemplateCommon<T> {
protected:
   SmallVectorTemplateBase(size_t Size) : SmallVectorTemplateCommon<T>(Size) {}

   static void destroy_range(T *S, T *E)
   {
      while (S != E) {
         --E;
         E->~T();
      }
   }

   /// Move the range [I, E) into the uninitialized memory starting with "Dest",
   /// constructing elements as needed.
   template <typename It1, typename It2>
   static void uninitialized_move(It1 I, It1 E, It2 Dest)
   {
      std::uninitialized_copy(std::make_move_iterator(I), std::make_move_iterator(E), Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory starting with "Dest",
   /// constructing elements as needed.
   template <typename It1, typename It2>
   static void uninitialized_copy(It1 I, It1 E, It2 Dest)
   {
      std::uninitialized_copy(I, E, Dest);
   }

   /// Grow the allocated memory (without initializing new elements), doubling
   /// the size of the allocated memory. Guarantees space for at least one more
   /// element, or MinSize more elements if specified.
   void grow(size_t MinSize = 0);

public:
   void push_back(const T &Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(Elt);
      this->set_size(this->size() + 1);
   }

   void push_back(T &&Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(::std::move(Elt));
      this->set_size(this->size() + 1);
   }

   void pop_back()
   {
      this->set_size(this->size() - 1);
      this->end()->~T();
   }
};

// Define this out-of-line to dissuade the C++ compiler from inlining it.
template <typename T, bool TriviallyCopyable>
void R__CLING_PTRCHECK(off) SmallVectorTemplateBase<T, TriviallyCopyable>::grow(size_t MinSize)
{
   // Ensure we can fit the new capacity.
   // This is only going to be applicable when the capacity is 32 bit.
   if (MinSize > this->SizeTypeMax())
      this->report_size_overflow(MinSize);

   // Ensure we can meet the guarantee of space for at least one more element.
   // The above check alone will not catch the case where grow is called with a
   // default MinSize of 0, but the current capacity cannot be increased.
   // This is only going to be applicable when the capacity is 32 bit.
   if (this->capacity() == this->SizeTypeMax())
      this->report_at_maximum_capacity();

   // Always grow, even from zero.
   size_t NewCapacity = size_t(NextPowerOf2(this->capacity() + 2));
   NewCapacity = std::min(std::max(NewCapacity, MinSize), this->SizeTypeMax());
   T *NewElts = static_cast<T *>(malloc(NewCapacity * sizeof(T)));
   R__ASSERT(NewElts != nullptr);

   // Move the elements over.
   this->uninitialized_move(this->begin(), this->end(), NewElts);

   if (this->Owns()) {
      // Destroy the original elements.
      destroy_range(this->begin(), this->end());

      // If this wasn't grown from the inline copy, deallocate the old space.
      if (!this->isSmall())
         free(this->begin());
   }

   this->fBeginX = NewElts;
   this->fCapacity = NewCapacity;
}

/// SmallVectorTemplateBase<TriviallyCopyable = true> - This is where we put
/// method implementations that are designed to work with trivially copyable
/// T's. This allows using memcpy in place of copy/move construction and
/// skipping destruction.
template <typename T>
class R__CLING_PTRCHECK(off) SmallVectorTemplateBase<T, true> : public SmallVectorTemplateCommon<T> {
   using SuperClass = SmallVectorTemplateCommon<T>;

protected:
   SmallVectorTemplateBase(size_t Size) : SmallVectorTemplateCommon<T>(Size) {}

   // No need to do a destroy loop for POD's.
   static void destroy_range(T *, T *) {}

   /// Move the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename It1, typename It2>
   static void uninitialized_move(It1 I, It1 E, It2 Dest)
   {
      // Just do a copy.
      uninitialized_copy(I, E, Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename It1, typename It2>
   static void uninitialized_copy(It1 I, It1 E, It2 Dest)
   {
      // Arbitrary iterator types; just use the basic implementation.
      std::uninitialized_copy(I, E, Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename T1, typename T2>
   static void uninitialized_copy(
      T1 *I, T1 *E, T2 *Dest,
      typename std::enable_if<std::is_same<typename std::remove_const<T1>::type, T2>::value>::type * = nullptr)
   {
      // Use memcpy for PODs iterated by pointers (which includes SmallVector
      // iterators): std::uninitialized_copy optimizes to memmove, but we can
      // use memcpy here. Note that I and E are iterators and thus might be
      // invalid for memcpy if they are equal.
      if (I != E)
         memcpy(reinterpret_cast<void *>(Dest), I, (E - I) * sizeof(T));
   }

   /// Double the size of the allocated memory, guaranteeing space for at
   /// least one more element or MinSize if specified.
   void grow(size_t MinSize = 0)
   {
      this->grow_pod(MinSize, sizeof(T));
   }

public:
   using iterator = typename SuperClass::iterator;
   using const_iterator = typename SuperClass::const_iterator;
   using reference = typename SuperClass::reference;
   using size_type = typename SuperClass::size_type;

   void push_back(const T &Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      memcpy(reinterpret_cast<void *>(this->end()), &Elt, sizeof(T));
      this->set_size(this->size() + 1);
   }

   void pop_back() { this->set_size(this->size() - 1); }
};

/// Storage for the SmallVector elements.  This is specialized for the N=0 case
/// to avoid allocating unnecessary storage.
template <typename T, unsigned N>
struct SmallVectorStorage {
   alignas(T) char InlineElts[N * sizeof(T)]{};
};

/// We need the storage to be properly aligned even for small-size of 0 so that
/// the pointer math in \a SmallVectorTemplateCommon::getFirstEl() is
/// well-defined.
template <typename T>
struct alignas(T) SmallVectorStorage<T, 0> {
};

/// The size of the inline storage of an RVec.
/// Our policy is to allocate at least 8 elements (or more if they all fit into one cacheline)
/// unless the size of the buffer with 8 elements would be over a certain maximum size.
template <typename T>
struct RVecInlineStorageSize {
private:
#ifdef R__HAS_HARDWARE_INTERFERENCE_SIZE
   static constexpr std::size_t cacheLineSize = std::hardware_destructive_interference_size;
#else
   // safe bet: assume the typical 64 bytes
   static constexpr std::size_t cacheLineSize = 64;
#endif
   static constexpr unsigned elementsPerCacheLine = (cacheLineSize - sizeof(SmallVectorBase)) / sizeof(T);
   static constexpr unsigned maxInlineByteSize = 1024;

public:
   static constexpr unsigned value =
      elementsPerCacheLine >= 8 ? elementsPerCacheLine : (sizeof(T) * 8 > maxInlineByteSize ? 0 : 8);
};

// A C++14-compatible implementation of std::uninitialized_value_construct
template <typename ForwardIt>
void UninitializedValueConstruct(ForwardIt first, ForwardIt last)
{
#if __cplusplus < 201703L
   for (; first != last; ++first)
      new (static_cast<void *>(std::addressof(*first))) typename std::iterator_traits<ForwardIt>::value_type();
#else
   std::uninitialized_value_construct(first, last);
#endif
}

/// An unsafe function to reset the buffer for which this RVec is acting as a view.
///
/// \note This is a low-level method that _must_ be called on RVecs that are already non-owning:
/// - it does not put the RVec in "non-owning mode" (fCapacity == -1)
/// - it does not free any owned buffer
template <typename T>
void ResetView(RVec<T> &v, T* addr, std::size_t sz)
{
   v.fBeginX = addr;
   v.fSize = sz;
}

} // namespace VecOps
} // namespace Internal

namespace Detail {
namespace VecOps {

/// This class consists of common code factored out of the SmallVector class to
/// reduce code duplication based on the SmallVector 'N' template parameter.
template <typename T>
class R__CLING_PTRCHECK(off) RVecImpl : public Internal::VecOps::SmallVectorTemplateBase<T> {
   using SuperClass = Internal::VecOps::SmallVectorTemplateBase<T>;

public:
   using iterator = typename SuperClass::iterator;
   using const_iterator = typename SuperClass::const_iterator;
   using reference = typename SuperClass::reference;
   using size_type = typename SuperClass::size_type;

protected:
   // Default ctor - Initialize to empty.
   explicit RVecImpl(unsigned N) : MyROOT::Internal::VecOps::SmallVectorTemplateBase<T>(N) {}

public:
   RVecImpl(const RVecImpl &) = delete;

   ~RVecImpl()
   {
      // Subclass has already destructed this vector's elements.
      // If this wasn't grown from the inline copy, deallocate the old space.
      if (!this->isSmall() && this->Owns())
         free(this->begin());
   }

   // also give up adopted memory if applicable
   void clear()
   {
      if (this->Owns()) {
         this->destroy_range(this->begin(), this->end());
         this->fSize = 0;
      } else {
         this->resetToSmall();
      }
   }

   void resize(size_type N)
   {
      if (N < this->size()) {
         if (this->Owns())
            this->destroy_range(this->begin() + N, this->end());
         this->set_size(N);
      } else if (N > this->size()) {
         if (this->capacity() < N)
            this->grow(N);
         for (auto I = this->end(), E = this->begin() + N; I != E; ++I)
            new (&*I) T();
         this->set_size(N);
      }
   }

   void resize(size_type N, const T &NV)
   {
      if (N < this->size()) {
         if (this->Owns())
            this->destroy_range(this->begin() + N, this->end());
         this->set_size(N);
      } else if (N > this->size()) {
         if (this->capacity() < N)
            this->grow(N);
         std::uninitialized_fill(this->end(), this->begin() + N, NV);
         this->set_size(N);
      }
   }

   void reserve(size_type N)
   {
      if (this->capacity() < N)
         this->grow(N);
   }

   void pop_back_n(size_type NumItems)
   {
      if (this->size() < NumItems) {
         throw std::runtime_error("Popping back more elements than those available.");
      }
      if (this->Owns())
         this->destroy_range(this->end() - NumItems, this->end());
      this->set_size(this->size() - NumItems);
   }

   R__RVEC_NODISCARD T pop_back_val()
   {
      T Result = ::std::move(this->back());
      this->pop_back();
      return Result;
   }

   void swap(RVecImpl &RHS);

   /// Add the specified range to the end of the SmallVector.
   template <typename in_iter,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<in_iter>::iterator_category, std::input_iterator_tag>::value>::type>
   void append(in_iter in_start, in_iter in_end)
   {
      size_type NumInputs = std::distance(in_start, in_end);
      if (NumInputs > this->capacity() - this->size())
         this->grow(this->size() + NumInputs);

      this->uninitialized_copy(in_start, in_end, this->end());
      this->set_size(this->size() + NumInputs);
   }

   /// Append \p NumInputs copies of \p Elt to the end.
   void append(size_type NumInputs, const T &Elt)
   {
      if (NumInputs > this->capacity() - this->size())
         this->grow(this->size() + NumInputs);

      std::uninitialized_fill_n(this->end(), NumInputs, Elt);
      this->set_size(this->size() + NumInputs);
   }

   void append(std::initializer_list<T> IL) { append(IL.begin(), IL.end()); }

   // from the original LLVM implementation:
   // FIXME: Consider assigning over existing elements, rather than clearing &
   // re-initializing them - for all assign(...) variants.

   void assign(size_type NumElts, const T &Elt)
   {
      clear();
      if (this->capacity() < NumElts)
         this->grow(NumElts);
      this->set_size(NumElts);
      std::uninitialized_fill(this->begin(), this->end(), Elt);
   }

   template <typename in_iter,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<in_iter>::iterator_category, std::input_iterator_tag>::value>::type>
   void assign(in_iter in_start, in_iter in_end)
   {
      clear();
      append(in_start, in_end);
   }

   void assign(std::initializer_list<T> IL)
   {
      clear();
      append(IL);
   }

   iterator erase(const_iterator CI)
   {
      // Just cast away constness because this is a non-const member function.
      iterator I = const_cast<iterator>(CI);

      if (I < this->begin() || I >= this->end()) {
         throw std::runtime_error("The iterator passed to `erase` is out of bounds.");
      }

      iterator N = I;
      // Shift all elts down one.
      std::move(I + 1, this->end(), I);
      // Drop the last elt.
      this->pop_back();
      return (N);
   }

   iterator erase(const_iterator CS, const_iterator CE)
   {
      // Just cast away constness because this is a non-const member function.
      iterator S = const_cast<iterator>(CS);
      iterator E = const_cast<iterator>(CE);

      if (S < this->begin() || E > this->end() || S > E) {
         throw std::runtime_error("Invalid start/end pair passed to `erase` (out of bounds or start > end).");
      }

      iterator N = S;
      // Shift all elts down.
      iterator I = std::move(E, this->end(), S);
      // Drop the last elts.
      if (this->Owns())
         this->destroy_range(I, this->end());
      this->set_size(I - this->begin());
      return (N);
   }

   iterator insert(iterator I, T &&Elt)
   {
      if (I == this->end()) { // Important special case for empty vector.
         this->push_back(::std::move(Elt));
         return this->end() - 1;
      }

      if (I < this->begin() || I > this->end()) {
         throw std::runtime_error("The iterator passed to `insert` is out of bounds.");
      }

      if (this->size() >= this->capacity()) {
         size_t EltNo = I - this->begin();
         this->grow();
         I = this->begin() + EltNo;
      }

      ::new ((void *)this->end()) T(::std::move(this->back()));
      // Push everything else over.
      std::move_backward(I, this->end() - 1, this->end());
      this->set_size(this->size() + 1);

      // If we just moved the element we're inserting, be sure to update
      // the reference.
      T *EltPtr = &Elt;
      if (I <= EltPtr && EltPtr < this->end())
         ++EltPtr;

      *I = ::std::move(*EltPtr);
      return I;
   }

   iterator insert(iterator I, const T &Elt)
   {
      if (I == this->end()) { // Important special case for empty vector.
         this->push_back(Elt);
         return this->end() - 1;
      }

      if (I < this->begin() || I > this->end()) {
         throw std::runtime_error("The iterator passed to `insert` is out of bounds.");
      }

      if (this->size() >= this->capacity()) {
         size_t EltNo = I - this->begin();
         this->grow();
         I = this->begin() + EltNo;
      }
      ::new ((void *)this->end()) T(std::move(this->back()));
      // Push everything else over.
      std::move_backward(I, this->end() - 1, this->end());
      this->set_size(this->size() + 1);

      // If we just moved the element we're inserting, be sure to update
      // the reference.
      const T *EltPtr = &Elt;
      if (I <= EltPtr && EltPtr < this->end())
         ++EltPtr;

      *I = *EltPtr;
      return I;
   }

   iterator insert(iterator I, size_type NumToInsert, const T &Elt)
   {
      // Convert iterator to elt# to avoid invalidating iterator when we reserve()
      size_t InsertElt = I - this->begin();

      if (I == this->end()) { // Important special case for empty vector.
         append(NumToInsert, Elt);
         return this->begin() + InsertElt;
      }

      if (I < this->begin() || I > this->end()) {
         throw std::runtime_error("The iterator passed to `insert` is out of bounds.");
      }

      // Ensure there is enough space.
      reserve(this->size() + NumToInsert);

      // Uninvalidate the iterator.
      I = this->begin() + InsertElt;

      // If there are more elements between the insertion point and the end of the
      // range than there are being inserted, we can use a simple approach to
      // insertion.  Since we already reserved space, we know that this won't
      // reallocate the vector.
      if (size_t(this->end() - I) >= NumToInsert) {
         T *OldEnd = this->end();
         append(std::move_iterator<iterator>(this->end() - NumToInsert), std::move_iterator<iterator>(this->end()));

         // Copy the existing elements that get replaced.
         std::move_backward(I, OldEnd - NumToInsert, OldEnd);

         std::fill_n(I, NumToInsert, Elt);
         return I;
      }

      // Otherwise, we're inserting more elements than exist already, and we're
      // not inserting at the end.

      // Move over the elements that we're about to overwrite.
      T *OldEnd = this->end();
      this->set_size(this->size() + NumToInsert);
      size_t NumOverwritten = OldEnd - I;
      this->uninitialized_move(I, OldEnd, this->end() - NumOverwritten);

      // Replace the overwritten part.
      std::fill_n(I, NumOverwritten, Elt);

      // Insert the non-overwritten middle part.
      std::uninitialized_fill_n(OldEnd, NumToInsert - NumOverwritten, Elt);
      return I;
   }

   template <typename ItTy,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<ItTy>::iterator_category, std::input_iterator_tag>::value>::type>
   iterator insert(iterator I, ItTy From, ItTy To)
   {
      // Convert iterator to elt# to avoid invalidating iterator when we reserve()
      size_t InsertElt = I - this->begin();

      if (I == this->end()) { // Important special case for empty vector.
         append(From, To);
         return this->begin() + InsertElt;
      }

      if (I < this->begin() || I > this->end()) {
         throw std::runtime_error("The iterator passed to `insert` is out of bounds.");
      }

      size_t NumToInsert = std::distance(From, To);

      // Ensure there is enough space.
      reserve(this->size() + NumToInsert);

      // Uninvalidate the iterator.
      I = this->begin() + InsertElt;

      // If there are more elements between the insertion point and the end of the
      // range than there are being inserted, we can use a simple approach to
      // insertion.  Since we already reserved space, we know that this won't
      // reallocate the vector.
      if (size_t(this->end() - I) >= NumToInsert) {
         T *OldEnd = this->end();
         append(std::move_iterator<iterator>(this->end() - NumToInsert), std::move_iterator<iterator>(this->end()));

         // Copy the existing elements that get replaced.
         std::move_backward(I, OldEnd - NumToInsert, OldEnd);

         std::copy(From, To, I);
         return I;
      }

      // Otherwise, we're inserting more elements than exist already, and we're
      // not inserting at the end.

      // Move over the elements that we're about to overwrite.
      T *OldEnd = this->end();
      this->set_size(this->size() + NumToInsert);
      size_t NumOverwritten = OldEnd - I;
      this->uninitialized_move(I, OldEnd, this->end() - NumOverwritten);

      // Replace the overwritten part.
      for (T *J = I; NumOverwritten > 0; --NumOverwritten) {
         *J = *From;
         ++J;
         ++From;
      }

      // Insert the non-overwritten middle part.
      this->uninitialized_copy(From, To, OldEnd);
      return I;
   }

   void insert(iterator I, std::initializer_list<T> IL) { insert(I, IL.begin(), IL.end()); }

   template <typename... ArgTypes>
   reference emplace_back(ArgTypes &&...Args)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(std::forward<ArgTypes>(Args)...);
      this->set_size(this->size() + 1);
      return this->back();
   }

   RVecImpl &operator=(const RVecImpl &RHS);

   RVecImpl &operator=(RVecImpl &&RHS);
};

template <typename T>
void RVecImpl<T>::swap(RVecImpl<T> &RHS)
{
   if (this == &RHS)
      return;

   // We can only avoid copying elements if neither vector is small.
   if (!this->isSmall() && !RHS.isSmall()) {
      std::swap(this->fBeginX, RHS.fBeginX);
      std::swap(this->fSize, RHS.fSize);
      std::swap(this->fCapacity, RHS.fCapacity);
      return;
   }

   // This block handles the swap of a small and a non-owning vector
   // It is more efficient to first move the non-owning vector, hence the 2 cases
   if (this->isSmall() && !RHS.Owns()) { // the right vector is non-owning
      RVecImpl<T> temp(0);
      temp = std::move(RHS);
      RHS = std::move(*this);
      *this = std::move(temp);
      return;
   } else if (RHS.isSmall() && !this->Owns()) { // the left vector is non-owning
      RVecImpl<T> temp(0);
      temp = std::move(*this);
      *this = std::move(RHS);
      RHS = std::move(temp);
      return;
   }

   if (RHS.size() > this->capacity())
      this->grow(RHS.size());
   if (this->size() > RHS.capacity())
      RHS.grow(this->size());

   // Swap the shared elements.
   size_t NumShared = this->size();
   if (NumShared > RHS.size())
      NumShared = RHS.size();
   for (size_type i = 0; i != NumShared; ++i)
      std::iter_swap(this->begin() + i, RHS.begin() + i);

   // Copy over the extra elts.
   if (this->size() > RHS.size()) {
      size_t EltDiff = this->size() - RHS.size();
      this->uninitialized_copy(this->begin() + NumShared, this->end(), RHS.end());
      RHS.set_size(RHS.size() + EltDiff);
      if (this->Owns())
         this->destroy_range(this->begin() + NumShared, this->end());
      this->set_size(NumShared);
   } else if (RHS.size() > this->size()) {
      size_t EltDiff = RHS.size() - this->size();
      this->uninitialized_copy(RHS.begin() + NumShared, RHS.end(), this->end());
      this->set_size(this->size() + EltDiff);
      if (RHS.Owns())
         this->destroy_range(RHS.begin() + NumShared, RHS.end());
      RHS.set_size(NumShared);
   }
}

template <typename T>
RVecImpl<T> &RVecImpl<T>::operator=(const RVecImpl<T> &RHS)
{
   // Avoid self-assignment.
   if (this == &RHS)
      return *this;

   // If we already have sufficient space, assign the common elements, then
   // destroy any excess.
   size_t RHSSize = RHS.size();
   size_t CurSize = this->size();
   if (CurSize >= RHSSize) {
      // Assign common elements.
      iterator NewEnd;
      if (RHSSize)
         NewEnd = std::copy(RHS.begin(), RHS.begin() + RHSSize, this->begin());
      else
         NewEnd = this->begin();

      // Destroy excess elements.
      if (this->Owns())
         this->destroy_range(NewEnd, this->end());

      // Trim.
      this->set_size(RHSSize);
      return *this;
   }

   // If we have to grow to have enough elements, destroy the current elements.
   // This allows us to avoid copying them during the grow.
   // From the original LLVM implementation:
   // FIXME: don't do this if they're efficiently moveable.
   if (this->capacity() < RHSSize) {
      if (this->Owns()) {
         // Destroy current elements.
         this->destroy_range(this->begin(), this->end());
      }
      this->set_size(0);
      CurSize = 0;
      this->grow(RHSSize);
   } else if (CurSize) {
      // Otherwise, use assignment for the already-constructed elements.
      std::copy(RHS.begin(), RHS.begin() + CurSize, this->begin());
   }

   // Copy construct the new elements in place.
   this->uninitialized_copy(RHS.begin() + CurSize, RHS.end(), this->begin() + CurSize);

   // Set end.
   this->set_size(RHSSize);
   return *this;
}

template <typename T>
RVecImpl<T> &RVecImpl<T>::operator=(RVecImpl<T> &&RHS)
{
   // Avoid self-assignment.
   if (this == &RHS)
      return *this;

   // If the RHS isn't small, clear this vector and then steal its buffer.
   if (!RHS.isSmall()) {
      if (this->Owns()) {
         this->destroy_range(this->begin(), this->end());
         if (!this->isSmall())
            free(this->begin());
      }
      this->fBeginX = RHS.fBeginX;
      this->fSize = RHS.fSize;
      this->fCapacity = RHS.fCapacity;
      RHS.resetToSmall();
      return *this;
   }

   // If we already have sufficient space, assign the common elements, then
   // destroy any excess.
   size_t RHSSize = RHS.size();
   size_t CurSize = this->size();
   if (CurSize >= RHSSize) {
      // Assign common elements.
      iterator NewEnd = this->begin();
      if (RHSSize)
         NewEnd = std::move(RHS.begin(), RHS.end(), NewEnd);

      // Destroy excess elements and trim the bounds.
      if (this->Owns())
         this->destroy_range(NewEnd, this->end());
      this->set_size(RHSSize);

      // Clear the RHS.
      RHS.clear();

      return *this;
   }

   // If we have to grow to have enough elements, destroy the current elements.
   // This allows us to avoid copying them during the grow.
   // From the original LLVM implementation:
   // FIXME: this may not actually make any sense if we can efficiently move
   // elements.
   if (this->capacity() < RHSSize) {
      if (this->Owns()) {
         // Destroy current elements.
         this->destroy_range(this->begin(), this->end());
      }
      this->set_size(0);
      CurSize = 0;
      this->grow(RHSSize);
   } else if (CurSize) {
      // Otherwise, use assignment for the already-constructed elements.
      std::move(RHS.begin(), RHS.begin() + CurSize, this->begin());
   }

   // Move-construct the new elements in place.
   this->uninitialized_move(RHS.begin() + CurSize, RHS.end(), this->begin() + CurSize);

   // Set end.
   this->set_size(RHSSize);

   RHS.clear();
   return *this;
}

template <typename T>
bool IsSmall(const MyROOT::VecOps::RVec<T> &v)
{
   return v.isSmall();
}

template <typename T>
bool IsAdopting(const MyROOT::VecOps::RVec<T> &v)
{
   return !v.Owns();
}

} // namespace VecOps
} // namespace Detail

namespace VecOps {
// Note that we open here with @{ the Doxygen group vecops and it is
// closed again at the end of the C++ namespace VecOps
/**
  * \defgroup vecops VecOps
  * A "std::vector"-like collection of values implementing handy operation to analyse them
  * @{
*/

// From the original SmallVector code:
// This is a 'vector' (really, a variable-sized array), optimized
// for the case when the array is small.  It contains some number of elements
// in-place, which allows it to avoid heap allocation when the actual number of
// elements is below that threshold.  This allows normal "small" cases to be
// fast without losing generality for large inputs.
//
// Note that this does not attempt to be exception safe.

template <typename T, unsigned int N>
class R__CLING_PTRCHECK(off) RVecN : public Detail::VecOps::RVecImpl<T>, Internal::VecOps::SmallVectorStorage<T, N> {
public:
   RVecN() : Detail::VecOps::RVecImpl<T>(N) {}

   ~RVecN()
   {
      if (this->Owns()) {
         // Destroy the constructed elements in the vector.
         this->destroy_range(this->begin(), this->end());
      }
   }

   explicit RVecN(size_t Size, const T &Value) : Detail::VecOps::RVecImpl<T>(N) { this->assign(Size, Value); }

   explicit RVecN(size_t Size) : Detail::VecOps::RVecImpl<T>(N)
   {
      if (Size > N)
         this->grow(Size);
      this->fSize = Size;
      MyROOT::Internal::VecOps::UninitializedValueConstruct(this->begin(), this->end());
   }

   template <typename ItTy,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<ItTy>::iterator_category, std::input_iterator_tag>::value>::type>
   RVecN(ItTy S, ItTy E) : Detail::VecOps::RVecImpl<T>(N)
   {
      this->append(S, E);
   }

   RVecN(std::initializer_list<T> IL) : Detail::VecOps::RVecImpl<T>(N) { this->assign(IL); }

   RVecN(const RVecN &RHS) : Detail::VecOps::RVecImpl<T>(N)
   {
      if (!RHS.empty())
         Detail::VecOps::RVecImpl<T>::operator=(RHS);
   }

   RVecN &operator=(const RVecN &RHS)
   {
      Detail::VecOps::RVecImpl<T>::operator=(RHS);
      return *this;
   }

   RVecN(RVecN &&RHS) : Detail::VecOps::RVecImpl<T>(N)
   {
      if (!RHS.empty())
         Detail::VecOps::RVecImpl<T>::operator=(::std::move(RHS));
   }

   RVecN(Detail::VecOps::RVecImpl<T> &&RHS) : Detail::VecOps::RVecImpl<T>(N)
   {
      if (!RHS.empty())
         Detail::VecOps::RVecImpl<T>::operator=(::std::move(RHS));
   }

   RVecN(const std::vector<T> &RHS) : RVecN(RHS.begin(), RHS.end()) {}

   RVecN &operator=(RVecN &&RHS)
   {
      Detail::VecOps::RVecImpl<T>::operator=(::std::move(RHS));
      return *this;
   }

   RVecN(T* p, size_t n) : Detail::VecOps::RVecImpl<T>(N)
   {
      this->fBeginX = p;
      this->fSize = n;
      this->fCapacity = -1;
   }

   RVecN &operator=(Detail::VecOps::RVecImpl<T> &&RHS)
   {
      Detail::VecOps::RVecImpl<T>::operator=(::std::move(RHS));
      return *this;
   }

   RVecN &operator=(std::initializer_list<T> IL)
   {
      this->assign(IL);
      return *this;
   }

   using reference = typename Internal::VecOps::SmallVectorTemplateCommon<T>::reference;
   using const_reference = typename Internal::VecOps::SmallVectorTemplateCommon<T>::const_reference;
   using size_type = typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type;
   using value_type = typename Internal::VecOps::SmallVectorTemplateCommon<T>::value_type;
   using Internal::VecOps::SmallVectorTemplateCommon<T>::begin;
   using Internal::VecOps::SmallVectorTemplateCommon<T>::size;

   reference operator[](size_type idx)
   {
      return begin()[idx];
   }

   const_reference operator[](size_type idx) const
   {
      return begin()[idx];
   }

   template <typename V, unsigned M, typename = std::enable_if<std::is_convertible<V, bool>::value>>
   RVecN operator[](const RVecN<V, M> &conds) const
   {
      const size_type n = conds.size();

      if (n != this->size()) {
         std::string msg = "Cannot index RVecN of size " + std::to_string(this->size()) +
                           " with condition vector of different size (" + std::to_string(n) + ").";
         throw std::runtime_error(msg);
      }

      size_type n_true = 0ull;
      for (auto c : conds)
         n_true += c; // relies on bool -> int conversion, faster than branching

      RVecN ret;
      ret.reserve(n_true);
      size_type j = 0u;
      for (size_type i = 0u; i < n; ++i) {
         if (conds[i]) {
            ret.push_back(this->operator[](i));
            ++j;
         }
      }
      return ret;
   }

   // conversion
   template <typename U, unsigned M, typename = std::enable_if<std::is_convertible<T, U>::value>>
   operator RVecN<U, M>() const
   {
      return RVecN<U, M>(this->begin(), this->end());
   }

   reference at(size_type pos)
   {
      if (pos >= size_type(this->fSize)) {
         std::string msg = "RVecN::at: size is " + std::to_string(this->fSize) + " but out-of-bounds index " +
                           std::to_string(pos) + " was requested.";
         throw std::out_of_range(msg);
      }
      return this->operator[](pos);
   }

   const_reference at(size_type pos) const
   {
      if (pos >= size_type(this->fSize)) {
         std::string msg = "RVecN::at: size is " + std::to_string(this->fSize) + " but out-of-bounds index " +
                           std::to_string(pos) + " was requested.";
         throw std::out_of_range(msg);
      }
      return this->operator[](pos);
   }

   /// No exception thrown. The user specifies the desired value in case the RVecN is shorter than `pos`.
   value_type at(size_type pos, value_type fallback)
   {
      if (pos >= size_type(this->fSize))
         return fallback;
      return this->operator[](pos);
   }

   /// No exception thrown. The user specifies the desired value in case the RVecN is shorter than `pos`.
   value_type at(size_type pos, value_type fallback) const
   {
      if (pos >= size_type(this->fSize))
         return fallback;
      return this->operator[](pos);
   }
};

// clang-format off
/**
\class MyROOT::VecOps::RVec
\brief A "std::vector"-like collection of values implementing handy operation to analyse them
\tparam T The type of the contained objects

A RVec is a container designed to make analysis of values' collections fast and easy.
Its storage is contiguous in memory and its interface is designed such to resemble to the one
of the stl vector. In addition the interface features methods and
[external functions](https://root.cern/doc/master/namespaceMyROOT_1_1VecOps.html) to ease the manipulation and analysis
of the data in the RVec.

\note MyROOT::VecOps::RVec can also be spelled simply MyROOT::RVec. Shorthand aliases such as MyROOT::RVecI or MyROOT::RVecD
are also available as template instantiations of RVec of fundamental types. The full list of available aliases:
- RVecB (`bool`)
- RVecC (`char`)
- RVecD (`double`)
- RVecF (`float`)
- RVecI (`int`)
- RVecL (`long`)
- RVecLL (`long long`)
- RVecU (`unsigned`)
- RVecUL (`unsigned long`)
- RVecULL (`unsigned long long`)

\note RVec does not attempt to be exception safe. Exceptions thrown by element constructors during insertions, swaps or
other operations will be propagated potentially leaving the RVec object in an invalid state.

\note RVec methods (e.g. `at` or `size`) follow the STL naming convention instead of the MyROOT naming convention in order
to make RVec a drop-in replacement for `std::vector`.

\htmlonly
<a href="https://doi.org/10.5281/zenodo.1253756"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1253756.svg" alt="DOI"></a>
\endhtmlonly

## Table of Contents
- [Example](\ref example)
- [Arithmetic operations, logical operations and mathematical functions](\ref operationsandfunctions)
- [Owning and adopting memory](\ref owningandadoptingmemory)
- [Sorting and manipulation of indices](\ref sorting)
- [Usage in combination with RDataFrame](\ref usagetdataframe)
- [Reference for the RVec class](\ref RVecdoxyref)
- [Reference for RVec helper functions](https://root.cern/doc/master/namespaceMyROOT_1_1VecOps.html)

\anchor example
## Example
Suppose to have an event featuring a collection of muons with a certain pseudorapidity,
momentum and charge, e.g.:
~~~{.cpp}
std::vector<short> mu_charge {1, 1, -1, -1, -1, 1, 1, -1};
std::vector<float> mu_pt {56, 45, 32, 24, 12, 8, 7, 6.2};
std::vector<float> mu_eta {3.1, -.2, -1.1, 1, 4.1, 1.6, 2.4, -.5};
~~~
Suppose you want to extract the transverse momenta of the muons satisfying certain
criteria, for example consider only negatively charged muons with a pseudorapidity
smaller or equal to 2 and with a transverse momentum greater than 10 GeV.
Such a selection would require, among the other things, the management of an explicit
loop, for example:
~~~{.cpp}
std::vector<float> goodMuons_pt;
const auto size = mu_charge.size();
for (size_t i=0; i < size; ++i) {
   if (mu_pt[i] > 10 && abs(mu_eta[i]) <= 2. &&  mu_charge[i] == -1) {
      goodMuons_pt.emplace_back(mu_pt[i]);
   }
}
~~~
These operations become straightforward with RVec - we just need to *write what
we mean*:
~~~{.cpp}
auto goodMuons_pt = mu_pt[ (mu_pt > 10.f && abs(mu_eta) <= 2.f && mu_charge == -1) ]
~~~
Now the clean collection of transverse momenta can be used within the rest of the data analysis, for
example to fill a histogram.

\anchor operationsandfunctions
## Arithmetic operations, logical operations and mathematical functions
Arithmetic operations on RVec instances can be performed: for example, they can be added, subtracted, multiplied.
~~~{.cpp}
RVec<double> v1 {1.,2.,3.,4.};
RVec<float> v2 {5.f,6.f,7.f,8.f};
auto v3 = v1+v2;
auto v4 = 3 * v1;
~~~
The supported operators are 
 - +, -, *, /
 - +=, -=, *=, /=
 - <, >, ==, !=, <=, >=, &&, ||
 - ~, !
 - &, |, ^
 - &=, |=, ^=
 - <<=, >>=

The most common mathematical functions are supported. It is possible to invoke them passing 
RVecs as arguments.
 - abs, fdim, fmod, remainder
 - floor, ceil, trunc, round, lround, llround
 - exp, exp2, expm1
 - log, log10, log2, log1p
 - pow
 - sqrt, cbrt
 - sin, cos, tan, asin, acos, atan, atan2, hypot
 - sinh, cosh, tanh, asinh, acosh
 - erf, erfc
 - lgamma, tgamma

\anchor owningandadoptingmemory
## Owning and adopting memory
RVec has contiguous memory associated to it. It can own it or simply adopt it. In the latter case,
it can be constructed with the address of the memory associated to it and its length. For example:
~~~{.cpp}
std::vector<int> myStlVec {1,2,3};
RVec<int> myRVec(myStlVec.data(), myStlVec.size());
~~~
In this case, the memory associated to myStlVec and myRVec is the same, myRVec simply "adopted it".
If any method which implies a re-allocation is called, e.g. *emplace_back* or *resize*, the adopted
memory is released and new one is allocated. The previous content is copied in the new memory and
preserved.

\anchor sorting
## Sorting and manipulation of indices

### Sorting
RVec complies to the STL interfaces when it comes to iterations. As a result, standard algorithms
can be used, for example sorting:
~~~{.cpp}
RVec<double> v{6., 4., 5.};
std::sort(v.begin(), v.end());
~~~

For convenience, helpers are provided too:
~~~{.cpp}
auto sorted_v = Sort(v);
auto reversed_v = Reverse(v);
~~~

### Manipulation of indices

It is also possible to manipulated the RVecs acting on their indices. For example,
the following syntax
~~~{.cpp}
RVecD v0 {9., 7., 8.};
auto v1 = Take(v0, {1, 2, 0});
~~~
will yield a new RVec<double> the content of which is the first, second and zeroth element of
v0, i.e. `{7., 8., 9.}`.

The `Argsort` and `StableArgsort` helper extracts the indices which order the content of a `RVec`.
For example, this snippet accomplishes in a more expressive way what we just achieved:
~~~{.cpp}
auto v1_indices = Argsort(v0); // The content of v1_indices is {1, 2, 0}.
v1 = Take(v0, v1_indices);
~~~

The `Take` utility allows to extract portions of the `RVec`. The content to be *taken*
can be specified with an `RVec` of indices or an integer. If the integer is negative,
elements will be picked starting from the end of the container:
~~~{.cpp}
RVecF vf {1.f, 2.f, 3.f, 4.f};
auto vf_1 = Take(vf, {1, 3}); // The content is {2.f, 4.f}
auto vf_2 = Take(vf, 2); // The content is {1.f, 2.f}
auto vf_3 = Take(vf, -3); // The content is {2.f, 3.f, 4.f}
~~~

\anchor usagetdataframe
## Usage in combination with RDataFrame
RDataFrame leverages internally RVecs. Suppose to have a dataset stored in a
TTree which holds these columns (here we choose C arrays to represent the
collections, they could be as well std::vector instances):
~~~{.bash}
  nPart            "nPart/I"            An integer representing the number of particles
  px               "px[nPart]/D"        The C array of the particles' x component of the momentum
  py               "py[nPart]/D"        The C array of the particles' y component of the momentum
  E                "E[nPart]/D"         The C array of the particles' Energy
~~~
Suppose you'd like to plot in a histogram the transverse momenta of all particles
for which the energy is greater than 200 MeV.
The code required would just be:
~~~{.cpp}
RDataFrame d("mytree", "myfile.root");
auto cutPt = [](RVecD &pxs, RVecD &pys, RVecD &Es) {
   auto all_pts = sqrt(pxs * pxs + pys * pys);
   auto good_pts = all_pts[Es > 200.];
   return good_pts;
   };

auto hpt = d.Define("pt", cutPt, {"px", "py", "E"})
            .Histo1D("pt");
hpt->Draw();
~~~
And if you'd like to express your selection as a string:
~~~{.cpp}
RDataFrame d("mytree", "myfile.root");
auto hpt = d.Define("pt", "sqrt(pxs * pxs + pys * pys)[E>200]")
            .Histo1D("pt");
hpt->Draw();
~~~
\anchor RVecdoxyref
**/
// clang-format on

template <typename T>
class R__CLING_PTRCHECK(off) RVec : public RVecN<T, Internal::VecOps::RVecInlineStorageSize<T>::value> {
   using SuperClass = RVecN<T, Internal::VecOps::RVecInlineStorageSize<T>::value>;

   friend void Internal::VecOps::ResetView<>(RVec<T> &v, T *addr, std::size_t sz);

public:
   using reference = typename SuperClass::reference;
   using const_reference = typename SuperClass::const_reference;
   using size_type = typename SuperClass::size_type;
   using value_type = typename SuperClass::value_type;
   using SuperClass::begin;
   using SuperClass::size;

   RVec() {}

   explicit RVec(size_t Size, const T &Value) : SuperClass(Size, Value) {}

   explicit RVec(size_t Size) : SuperClass(Size) {}

   template <typename ItTy,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<ItTy>::iterator_category, std::input_iterator_tag>::value>::type>
   RVec(ItTy S, ItTy E) : SuperClass(S, E)
   {
   }

   RVec(std::initializer_list<T> IL) : SuperClass(IL) {}

   RVec(const RVec &RHS) : SuperClass(RHS) {}

   RVec &operator=(const RVec &RHS)
   {
      SuperClass::operator=(RHS);
      return *this;
   }

   RVec(RVec &&RHS) : SuperClass(std::move(RHS)) {}

   RVec &operator=(RVec &&RHS)
   {
      SuperClass::operator=(std::move(RHS));
      return *this;
   }

   RVec(Detail::VecOps::RVecImpl<T> &&RHS) : SuperClass(std::move(RHS)) {}

   template <unsigned N>
   RVec(RVecN<T, N> &&RHS) : SuperClass(std::move(RHS)) {}

   template <unsigned N>
   RVec(const RVecN<T, N> &RHS) : SuperClass(RHS) {}

   RVec(const std::vector<T> &RHS) : SuperClass(RHS) {}

   RVec(T* p, size_t n) : SuperClass(p, n) {}

   // conversion
   template <typename U, typename = std::enable_if<std::is_convertible<T, U>::value>>
   operator RVec<U>() const
   {
      return RVec<U>(this->begin(), this->end());
   }

   using SuperClass::operator[];

   template <typename V, typename = std::enable_if<std::is_convertible<V, bool>::value>>
   RVec operator[](const RVec<V> &conds) const
   {
      return RVec(SuperClass::operator[](conds));
   }

   using SuperClass::at;

   friend bool MyROOT::Detail::VecOps::IsSmall<T>(const RVec<T> &v);

   friend bool MyROOT::Detail::VecOps::IsAdopting<T>(const RVec<T> &v);
};

template <typename T, unsigned N>
inline size_t CapacityInBytes(const RVecN<T, N> &X)
{
   return X.capacity_in_bytes();
}

template <typename T> struct PromoteTypeImpl;

template <> struct PromoteTypeImpl<float>       { using Type = float;       };
template <> struct PromoteTypeImpl<double>      { using Type = double;      };
template <> struct PromoteTypeImpl<long double> { using Type = long double; };

template <typename T> struct PromoteTypeImpl { using Type = double; };

template <typename T>
using PromoteType = typename PromoteTypeImpl<T>::Type;

template <typename U, typename V>
using PromoteTypes = decltype(PromoteType<U>() + PromoteType<V>());


template <typename T>
RVec<T> Take(const RVec<T> &v, const RVec<typename RVec<T>::size_type> &i)
{
   using size_type = typename RVec<T>::size_type;
   const size_type isize = i.size();
   RVec<T> r(isize);
   for (size_type k = 0; k < isize; k++)
      r[k] = v[i[k]];
   return r;
}

/// Take version that defaults to (user-specified) output value if some index is out of range
template <typename T>
RVec<T> Take(const RVec<T> &v, const RVec<typename RVec<T>::size_type> &i, const T default_val)
{
   using size_type = typename RVec<T>::size_type;
   const size_type isize = i.size();
   RVec<T> r(isize);
   for (size_type k = 0; k < isize; k++)
   {
      if (i[k] < v.size() && i[k]>=0){
         r[k] = v[i[k]];
      }
      else {
         r[k] = default_val;
      }
   }
   return r;
}

/// Return first `n` elements of an RVec if `n > 0` and last `n` elements if `n < 0`.
///
/// Example code, at the MyROOT prompt:
/// ~~~{.cpp}
/// using namespace MyROOT::VecOps;
/// RVecD v {2., 3., 1.};
/// auto firstTwo = Take(v, 2);
/// firstTwo
/// // (MyROOT::VecOps::RVec<double>) { 2.0000000, 3.0000000 }
/// auto lastOne = Take(v, -1);
/// lastOne
/// // (MyROOT::VecOps::RVec<double>) { 1.0000000 }
/// ~~~
template <typename T>
RVec<T> Take(const RVec<T> &v, const int n)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v.size();
   const size_type absn = std::abs(n);
   if (absn > size) {
      const auto msg = std::to_string(absn) + " elements requested from Take but input contains only " +
                       std::to_string(size) + " elements.";
      throw std::runtime_error(msg);
   }
   RVec<T> r(absn);
   if (n < 0) {
      for (size_type k = 0; k < absn; k++)
         r[k] = v[size - absn + k];
   } else {
      for (size_type k = 0; k < absn; k++)
         r[k] = v[k];
   }
   return r;
}

/// Return first `n` elements of an RVec if `n > 0` and last `n` elements if `n < 0`.
///
/// This Take version defaults to a user-specified value
/// `default_val` if the absolute value of `n` is
/// greater than the size of the RVec `v`
///
/// Example code, at the MyROOT prompt:
/// ~~~{.cpp}
/// using MyROOT::VecOps::RVec;
/// RVec<int> x{1,2,3,4};
/// Take(x,-5,1)
/// // (MyROOT::VecOps::RVec<int>) { 1, 1, 2, 3, 4 }
/// Take(x,5,20)
/// // (MyROOT::VecOps::RVec<int>) { 1, 2, 3, 4, 20 }
/// Take(x,-1,1)
/// // (MyROOT::VecOps::RVec<int>) { 4 }
/// Take(x,4,1)
/// // (MyROOT::VecOps::RVec<int>) { 1, 2, 3, 4 }
/// ~~~
template <typename T>
RVec<T> Take(const RVec<T> &v, const int n, const T default_val)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v.size();
   const size_type absn = std::abs(n);
   // Base case, can be handled by another overload of Take
   if (absn <= size) {
      return Take(v, n);
   }
   RVec<T> temp = v;
   // Case when n is positive and n > v.size()
   if (n > 0) {
      temp.resize(n, default_val);
      return temp;
   }
   // Case when n is negative and abs(n) > v.size()
   const auto num_to_fill = absn - size;
   MyROOT::VecOps::RVec<T> fill_front(num_to_fill, default_val);
   return Concatenate(fill_front, temp);
}

////////////////////////////////////////////////////////////////////////////////
/// Print a RVec at the prompt:
template <class T>
std::ostream &operator<<(std::ostream &os, const RVec<T> &v)
{
   // In order to print properly, convert to 64 bit int if this is a char
   constexpr bool mustConvert = std::is_same<char, T>::value || std::is_same<signed char, T>::value ||
                                std::is_same<unsigned char, T>::value || std::is_same<wchar_t, T>::value ||
                                std::is_same<char16_t, T>::value || std::is_same<char32_t, T>::value;
   using Print_t = typename std::conditional<mustConvert, long long int, T>::type;
   os << "{ ";
   auto size = v.size();
   if (size) {
      for (std::size_t i = 0; i < size - 1; ++i) {
         os << (Print_t)v[i] << ", ";
      }
      os << (Print_t)v[size - 1];
   }
   os << " }";
   return os;
}


/** @} */ // end of Doxygen group vecops

} // End of VecOps NS

} // End of MyROOT NS

inline void MyROOT::Internal::VecOps::SmallVectorBase::report_size_overflow(size_t MinSize)
{
   std::string Reason = "RVec unable to grow. Requested capacity (" + std::to_string(MinSize) +
                        ") is larger than maximum value for size type (" + std::to_string(SizeTypeMax()) + ")";
   throw std::length_error(Reason);
}

inline void MyROOT::Internal::VecOps::SmallVectorBase::report_at_maximum_capacity()
{
   std::string Reason = "RVec capacity unable to grow. Already at maximum size " + std::to_string(SizeTypeMax());
   throw std::length_error(Reason);
}

// Note: Moving this function into the header may cause performance regression.
inline void MyROOT::Internal::VecOps::SmallVectorBase::grow_pod(void *FirstEl, size_t MinSize, size_t TSize)
{
   // Ensure we can fit the new capacity.
   // This is only going to be applicable when the capacity is 32 bit.
   if (MinSize > SizeTypeMax())
      report_size_overflow(MinSize);

   // Ensure we can meet the guarantee of space for at least one more element.
   // The above check alone will not catch the case where grow is called with a
   // default MinSize of 0, but the current capacity cannot be increased.
   // This is only going to be applicable when the capacity is 32 bit.
   if (capacity() == SizeTypeMax())
      report_at_maximum_capacity();

   // In theory 2*capacity can overflow if the capacity is 64 bit, but the
   // original capacity would never be large enough for this to be a problem.
   size_t NewCapacity = 2 * capacity() + 1; // Always grow.
   NewCapacity = std::min(std::max(NewCapacity, MinSize), SizeTypeMax());

   void *NewElts;
   if (fBeginX == FirstEl || !this->Owns()) {
      NewElts = malloc(NewCapacity * TSize);
      R__ASSERT(NewElts != nullptr);

      // Copy the elements over.  No need to run dtors on PODs.
      memcpy(NewElts, this->fBeginX, size() * TSize);
   } else {
      // If this wasn't grown from the inline copy, grow the allocated space.
      NewElts = realloc(this->fBeginX, NewCapacity * TSize);
      R__ASSERT(NewElts != nullptr);
   }

   this->fBeginX = NewElts;
   this->fCapacity = NewCapacity;
}

#endif // MyROOT_RVEC
