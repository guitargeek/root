#pragma once

namespace N {
struct S {
   int fI;
};
} // namespace N

struct Dep1 {
   N::S fS;
};
