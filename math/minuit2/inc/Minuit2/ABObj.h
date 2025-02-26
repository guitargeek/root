// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ABObj
#define ROOT_Minuit2_ABObj

namespace ROOT {

namespace Minuit2 {

class gen {};
class sym {};
class vec {};

template <class A, class B>
class AlgebraicSumType {
public:
   typedef gen Type;
};

template <class T>
class AlgebraicSumType<T, T> {
public:
   typedef T Type;
};

template <>
class AlgebraicSumType<vec, gen> {
   typedef gen Type;
};

template <>
class AlgebraicSumType<gen, vec> {
   typedef gen Type;
};

template <>
class AlgebraicSumType<vec, sym> {
   typedef gen Type;
};

template <>
class AlgebraicSumType<sym, vec> {
   typedef gen Type;
};

//

template <class A, class B>
class AlgebraicProdType {
   typedef gen Type;
};

template <class T>
class AlgebraicProdType<T, T> {
   typedef T Type;
};

template <>
class AlgebraicProdType<gen, gen> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, sym> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, gen> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<gen, sym> {
public:
   typedef gen Type;
};

template <>
class AlgebraicProdType<vec, gen> {
   typedef gen Type;
};

template <>
class AlgebraicProdType<gen, vec> {
public:
   typedef vec Type;
};

template <>
class AlgebraicProdType<vec, sym> {
   typedef gen Type;
};

template <>
class AlgebraicProdType<sym, vec> {
public:
   typedef vec Type;
};

// Helper base class to delete assignment.
//
// Note that base classes without any data members or virtual functions don't
// cause any runtime overhead.
//
// Assignment is often historically deleted (that was done probably to avoid
// mistakes from accidental re-assignment). Also define destructor and copy
// constructor in this case, according to the rule of five.
class DeleteAssignment {
public:
   DeleteAssignment() = default;
   ~DeleteAssignment() = default;
   DeleteAssignment(const DeleteAssignment &) = default;
   DeleteAssignment(DeleteAssignment &&) = default;
   DeleteAssignment &operator=(const DeleteAssignment &) = delete;
   DeleteAssignment &operator=(DeleteAssignment &&) = delete;
};

template <class Type, class M, class T=double>
class ABObj : public DeleteAssignment {

public:
   ABObj(const M &obj, T factor = 1.) : fObject(obj), fFactor(factor) {}

   const M &Obj() const { return fObject; }

   T f() const { return fFactor; }

private:
   M fObject;
   T fFactor;
};

// templated scaling operator *
template <class mt, class M, class T>
ABObj<mt, M, T> operator*(T f, const M &obj)
{
   return {obj, f};
}

// templated operator /
template <class mt, class M, class T>
ABObj<mt, M, T> operator/(const M &obj, T f)
{
   return {obj, T(1.) / f};
}

// templated unary operator -
template <class mt, class M, class T>
ABObj<mt, M, T> operator-(const M &obj)
{
   return {obj, -1.};
}

// factor * ABObj
template <class mt, class M, class T>
ABObj<mt, M, T> operator*(T f, const ABObj<mt, M, T> &obj)
{
   return {obj.Obj(), obj.f() * f};
}

// ABObj / factor
template <class mt, class M, class T>
ABObj<mt, M, T> operator/(const ABObj<mt, M, T> &obj, T f)
{
   return {obj.Obj(), obj.f() / f};
}

// -ABObj
template <class mt, class M, class T>
ABObj<mt, M, T> operator-(const ABObj<mt, M, T> &obj)
{
   return {obj.Obj(), T(-1.) * obj.f()};
}

template <class M1, class M2>
class ABSum : public DeleteAssignment {
public:
   ABSum(const M1 &a, const M2 &b) : fA(a), fB(b) {}

   const M1 &A() const { return fA; }
   const M2 &B() const { return fB; }

private:
   M1 fA;
   M2 fB;
};

// ABObj + ABObj
template <class atype, class A, class btype, class B, class T>
ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>
operator+(const ABObj<atype, A, T> &a, const ABObj<btype, B, T> &b)
{

   return {ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>(a, b)};
}

// ABObj - ABObj
template <class atype, class A, class btype, class B, class T>
ABObj<typename AlgebraicSumType<atype, btype>::Type, ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>
operator-(const ABObj<atype, A, T> &a, const ABObj<btype, B, T> &b)
{

   return {ABSum<ABObj<atype, A, T>, ABObj<btype, B, T>>(a, ABObj<btype, B, T>(b.Obj(), T(-1.) * b.f()))};
}

template <class M1, class M2>
class ABProd : public DeleteAssignment {
public:
   ABProd(const M1 &a, const M2 &b) : fA(a), fB(b) {}

   const M1 &A() const { return fA; }
   const M2 &B() const { return fB; }

private:
   M1 fA;
   M2 fB;
};

// ABObj * ABObj
template <class atype, class A, class btype, class B, class T>
ABObj<typename AlgebraicProdType<atype, btype>::Type, ABProd<ABObj<atype, A, T>, ABObj<btype, B, T>>, T>
operator*(const ABObj<atype, A, T> &a, const ABObj<btype, B, T> &b)
{

   return {ABProd<ABObj<atype, A, T>, ABObj<btype, B, T>>(a, b)};
}

template <class M, class T>
class VectorOuterProduct {

public:
   VectorOuterProduct(const M &obj) : fObject(obj) {}

   const M &Obj() const { return fObject; }

private:
   M fObject;
};

template <class M, class T>
ABObj<sym, VectorOuterProduct<ABObj<vec, M, T>, T>, T> Outer_product(const ABObj<vec, M, T> &obj)
{
   return {VectorOuterProduct<ABObj<vec, M, T>, T>(obj)};
}

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ABObj
