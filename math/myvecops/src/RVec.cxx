#include "MyROOT/RVec.hxx"
using namespace MyROOT::VecOps;

// Check that no bytes are wasted and everything is well-aligned.
namespace {
struct Struct16B {
   alignas(16) void *X;
};
struct Struct32B {
   alignas(32) void *X;
};
}
static_assert(alignof(RVec<Struct16B>) >= alignof(Struct16B), "wrong alignment for 16-byte aligned T");
static_assert(alignof(RVec<Struct32B>) >= alignof(Struct32B), "wrong alignment for 32-byte aligned T");
static_assert(sizeof(RVec<Struct16B>) >= alignof(Struct16B), "missing padding for 16-byte aligned T");
static_assert(sizeof(RVec<Struct32B>) >= alignof(Struct32B), "missing padding for 32-byte aligned T");
static_assert(sizeof(RVecN<void *, 0>) == sizeof(int32_t) * 2 + sizeof(void *), "wasted space in RVec");
static_assert(sizeof(RVecN<void *, 1>) == sizeof(int32_t) * 2 + sizeof(void *) * 2,
              "wasted space in SmallVector size 1");
static_assert(sizeof(RVec<void *>) ==
                 sizeof(int32_t) * 2 +
                    (1 + MyROOT::Internal::VecOps::RVecInlineStorageSize<void *>::value) * sizeof(void *),
              "wasted space in RVec");

void MyROOT::Internal::VecOps::SmallVectorBase::report_size_overflow(size_t MinSize)
{
   std::string Reason = "RVec unable to grow. Requested capacity (" + std::to_string(MinSize) +
                        ") is larger than maximum value for size type (" + std::to_string(SizeTypeMax()) + ")";
   throw std::length_error(Reason);
}

void MyROOT::Internal::VecOps::SmallVectorBase::report_at_maximum_capacity()
{
   std::string Reason = "RVec capacity unable to grow. Already at maximum size " + std::to_string(SizeTypeMax());
   throw std::length_error(Reason);
}

// Note: Moving this function into the header may cause performance regression.
void MyROOT::Internal::VecOps::SmallVectorBase::grow_pod(void *FirstEl, size_t MinSize, size_t TSize)
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

#if (_VECOPS_USE_EXTERN_TEMPLATES)

namespace MyROOT {
namespace VecOps {

#define RVEC_DECLARE_FLOAT_TEMPLATE(T)  \
   template class RVec<T>;

#define RVEC_DECLARE_INTEGER_TEMPLATE(T) \
   template class RVec<T>;

RVEC_DECLARE_INTEGER_TEMPLATE(char)
RVEC_DECLARE_INTEGER_TEMPLATE(short)
RVEC_DECLARE_INTEGER_TEMPLATE(int)
RVEC_DECLARE_INTEGER_TEMPLATE(long)
RVEC_DECLARE_INTEGER_TEMPLATE(long long)

RVEC_DECLARE_INTEGER_TEMPLATE(unsigned char)
RVEC_DECLARE_INTEGER_TEMPLATE(unsigned short)
RVEC_DECLARE_INTEGER_TEMPLATE(unsigned int)
RVEC_DECLARE_INTEGER_TEMPLATE(unsigned long)
RVEC_DECLARE_INTEGER_TEMPLATE(unsigned long long)

RVEC_DECLARE_FLOAT_TEMPLATE(float)
RVEC_DECLARE_FLOAT_TEMPLATE(double)


} // namespace VecOps
} // namespace MyROOT

#endif // _VECOPS_USE_EXTERN_TEMPLATES
