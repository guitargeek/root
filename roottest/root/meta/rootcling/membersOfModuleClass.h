#ifndef ROOTTEST_MEMBERSOFMODULECLASS_H
#define ROOTTEST_MEMBERSOFMODULECLASS_H

#include "Rtypes.h"

// The classes in this header go into the C++ module of
// membersOfModuleClassDict. The dictionaries membersOfModuleClassPatternDict,
// membersOfModuleClassSpecDict and membersOfModuleClassRedeclDict select some
// of their members, which they don't use otherwise, so these are not loaded
// from the module when rootcling scans for the selected classes.
namespace MembersOfModuleClass {

struct Outer {
   struct Inner1 {
      int fX = 0;
   };
   struct Inner2 {
      int fY = 0;
   };
};

struct OuterWithSpecialization {
   template <class T>
   struct Tmpl {
      T fT{};
   };
#ifdef __clang__
   // GCC doesn't accept explicit specializations in class scope.
   template <>
   struct Tmpl<char> {
      char fC = 0;
   };
#endif
};

// Member class templates that are declared more than once, and a template
// that is first declared as a friend. The instantiations of a template are
// reached from its first declaration, which is in the class.
struct OuterWithRedecls {
   template <class T>
   struct DefinedOutside;

   template <class T>
   struct Redeclared;
   template <class T>
   struct Redeclared {
      T fT{};
   };

   template <class T>
   friend struct FriendTmplFirstDeclaredAsFriend;
};

template <class T>
struct OuterWithRedecls::DefinedOutside {
   T fT{};
};

template <class T>
struct FriendTmplFirstDeclaredAsFriend {
   T fT{};
};

// Nested classes that are declared before they are defined. The loaded members
// of the class might contain the forward declaration but not the definition.
template <class T>
struct TmplWithForwardDecl {
   struct Member;
   struct Member {
      T fM{};
   };
};

struct UsesTmplWithForwardDecl {
   TmplWithForwardDecl<int> fT;
};

struct BaseWithForwardDecl {
   struct Inner;
   struct Inner {
      int fX = 0;
   };
};

struct DerivedWithUsing : BaseWithForwardDecl {
   using BaseWithForwardDecl::Inner;
};

// A member class template that is declared before a typedef to one of its
// instantiations and defined after it. If the instantiation is selected both by
// its name and by the typedef, the rule that is matched first wins. A full
// traversal visits the instantiation first, at the first declaration of the
// template.
struct OuterWithTypedef {
   template <class T>
   struct Tmpl;
   typedef Tmpl<int> TmplInt;
   template <class T>
   struct Tmpl {
      T fT{};
      ClassDefNV(Tmpl, 1);
   };
};

struct UsesOuterWithTypedef {
   OuterWithTypedef::Tmpl<int> fT;
};

// A class that is defined in the declaration of a typedef to it. If it is
// selected both by its name and by the typedef, the rule that is matched first
// wins. A full traversal visits the class first.
struct OuterWithTypedefStruct {
   typedef struct Inner {
      int fX = 0;
      ClassDefNV(Inner, 1);
   } Alias;
   int fA = 0;
};

// A class with static data members that are excluded by a pattern. The rule
// must match them, otherwise rootcling warns that it is unused.
struct OuterWithStatic {
   static int sCount;
   int fA = 0;
};

} // namespace MembersOfModuleClass

#endif
