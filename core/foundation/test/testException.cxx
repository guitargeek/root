#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "ROOT/TestSupport.hxx"

#include <ROOT/RError.hxx>

#include <memory>
#include <stdexcept>

namespace {

/// Used to verify that wrapped return values are not unnecessarily copied
struct ComplexReturnType {
   static int gNCopies;
   ComplexReturnType() { gNCopies++; }
   ComplexReturnType(const ComplexReturnType &) { gNCopies++; }
   ComplexReturnType(ComplexReturnType &&other) = default;
   ComplexReturnType &operator=(const ComplexReturnType &) { return *this; }
   ComplexReturnType &operator=(ComplexReturnType &&other) = default;
};
int ComplexReturnType::gNCopies = 0;

static ROOT::RResult<void> TestFailure()
{
   return R__FAIL("test failure");
}

static ROOT::RResult<void> TestSuccess()
{
   return ROOT::RResult<void>::Success();
}

static ROOT::RResult<void> TestSuccessOrFailure(bool succeed)
{
   if (succeed)
      return R__FORWARD_RESULT(TestSuccess());
   return R__FORWARD_RESULT(TestFailure());
}

static ROOT::RResult<int> TestSyscall(bool succeed)
{
   if (succeed)
      return 42;
   return R__FAIL("failure");
}

static ROOT::RResult<int> TestChain(bool succeed)
{
   auto rv = TestSyscall(succeed);
   return R__FORWARD_RESULT(rv);
}

static ROOT::RResult<int> TestChainMultiTypes(bool succeed)
{
   auto rv = TestSuccessOrFailure(succeed);
   if (!rv)
      return R__FORWARD_ERROR(rv);
   return 0;
}

static ROOT::RResult<ComplexReturnType> TestComplex()
{
   return ComplexReturnType();
}

class ExceptionX : public std::runtime_error {
public:
   explicit ExceptionX(const std::string &what) : std::runtime_error(what) {}
};

} // anonymous namespace

TEST(Exception, Report)
{
   try {
      TestChain(false);
      EXPECT_TRUE(false) << "Above line should have thrown!";
   } catch (const ROOT::RException &e) {
      ASSERT_EQ(2U, e.GetError().GetStackTrace().size());
      EXPECT_THAT(e.GetError().GetStackTrace().at(0).fFunction, ::testing::HasSubstr("TestSyscall(bool)"));
      EXPECT_THAT(e.GetError().GetStackTrace().at(1).fFunction, ::testing::HasSubstr("TestChain(bool)"));
   }
}

TEST(Exception, ForwardResult)
{
   auto res = TestChain(true);
   ASSERT_TRUE(static_cast<bool>(res));
   EXPECT_EQ(42, res.Inspect());
}

TEST(Exception, ForwardError)
{
   EXPECT_THROW(TestSuccessOrFailure(false), ROOT::RException);
   EXPECT_NO_THROW(TestSuccessOrFailure(true));

   auto res = TestChainMultiTypes(true);
   ASSERT_TRUE(static_cast<bool>(res));
   EXPECT_EQ(0, res.Inspect());

   EXPECT_THROW(TestChainMultiTypes(false), ROOT::RException);
}

TEST(Exception, DiscardReturnValue)
{
   EXPECT_THROW(TestFailure(), ROOT::RException);
   EXPECT_NO_THROW(TestSuccess());
}

TEST(Exception, CheckReturnValue)
{
   auto rv = TestFailure();
   EXPECT_FALSE(rv);
   // No exception / crash when the scope closes
}

TEST(Exception, DoubleThrow)
{
   // We have to suppress a warning because of the double throw.
   ROOT::TestSupport::FilterDiagsRAII filterDiags{
      [](int /*severity*/, bool, const char * /*location*/, const char *msg) {
         EXPECT_STREQ(msg, "unhandled RResult exception during stack unwinding");
      }};

   try {
      auto rv = TestFailure();
      // Throwing ExceptionX will destruct rv along the way. Since rv carries an error state, it would normally
      // throw an exception itself. In this test, we verify that rv surpresses throwing an exception if another
      // exception is currently active.
      throw ExceptionX("something else went wrong");
   } catch (const ExceptionX &) {
      // This will only catch ExceptionX but not RException. In case rv mistakenly throws an exception,
      // we would notice the test failure by a crash of the unit test.
   }
}

TEST(Exception, VoidThrowOnError)
{
   // no throw on success
   TestSuccess().ThrowOnError();
   // throw on failure
   EXPECT_THROW(TestFailure().ThrowOnError(), ROOT::RException);
}

TEST(Exception, Syscall)
{
   auto fd = TestSyscall(true);
   if (!fd) {
      // In production code, we would expect error handling code other than throw
      EXPECT_THROW(fd.Throw(), ROOT::RException);
   }
   EXPECT_EQ(42, fd.Inspect());

   EXPECT_THROW(TestSyscall(false).Inspect(), ROOT::RException);
}

TEST(Exception, ComplexReturnType)
{
   auto res = TestComplex();
   EXPECT_EQ(1, ComplexReturnType::gNCopies);
}

TEST(Exception, MoveOnlyReturnType)
{
   auto TestMoveOnly = []() -> ROOT::RResult<std::unique_ptr<int>> { return std::make_unique<int>(1); };
   auto res = TestMoveOnly();

   // Using Inspect to make a copy won't compile
   // auto copy_inner = res.Inspect();

   // This will compile, but we only have read-only access
   const auto &copy_inner = res.Inspect();
   EXPECT_EQ(1, *copy_inner);

   // Instead, Unwrap is required to get ownership of the move-only type
   auto move_inner = res.Unwrap();
   EXPECT_EQ(1, *move_inner);
   move_inner.reset();
   move_inner = std::make_unique<int>(2);
   EXPECT_EQ(2, *move_inner);
}
