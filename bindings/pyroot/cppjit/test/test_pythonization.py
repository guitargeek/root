import sys

import py
from pytest import mark, raises
from support import (
    IS_CLANG_REPL,
    IS_LINUX_ARM,
    IS_MAC,
    IS_VALGRIND,
    setup_make,
)

currpath = py.path.local(__file__).dirpath()
test_dct = "pythonizables_cxx"


def setup_module(mod):
    setup_make("pythonizables")


class TestClassPYTHONIZATION:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppjit

        cls.pyzables = cppjit.load_reflection_info(cls.test_dct)

    @mark.xfail(condition=IS_MAC, reason="Fails on OSX")
    def test00_api(self):
        """Test basic semantics of the pythonization API"""

        import cppjit

        raises(TypeError, cppjit.py.add_pythonization, 1)

        def pythonizor1(klass, name):
            pass

        def pythonizor2(klass, name):
            pass

        pythonizor3 = pythonizor1

        cppjit.py.add_pythonization(pythonizor1)
        assert cppjit.py.remove_pythonization(pythonizor2) == False
        assert cppjit.py.remove_pythonization(pythonizor3) == True

        def pythonizor(klass, name):
            if name == "pyzables::SomeDummy1":
                klass.test = 1

        cppjit.py.add_pythonization(pythonizor)
        assert cppjit.gbl.pyzables.SomeDummy1.test == 1

        def pythonizor(klass, name):
            if name == "SomeDummy2":
                klass.test = 2

        cppjit.py.add_pythonization(pythonizor, "pyzables")

        # global pythonizors are still run even if namespaced ones available
        def pythonizor(klass, name):
            if name == "pyzables::SomeDummy2":
                klass.test = 3

        cppjit.py.add_pythonization(pythonizor)

        assert cppjit.gbl.pyzables.SomeDummy2.test == 3

        cppjit.cppdef("""
            namespace pyzables {
                class TObjString {
                public:
                    std::string s;
                    TObjString(std::string ss) : s(ss) {}
                    size_t Sizeof() { return s.size() + 1; }
                };
            }
        """)

        def root_pythonizor(klass, name):
            if name == "pyzables::TObjString":
                klass.__len__ = klass.Sizeof

        cppjit.py.add_pythonization(root_pythonizor)

        assert len(cppjit.gbl.pyzables.TObjString("aap")) == 4  # include '\0'

    def test01_size_mapping(self):
        """Use composites to map GetSize() onto buffer returns"""

        import cppjit

        def set_size(self, buf):
            buf.reshape((self.GetN(),))
            return buf

        cppjit.py.add_pythonization(
            cppjit.py.compose_method("NakedBuffers$", "Get[XY]$", set_size), "pyzables"
        )

        bsize, xval, yval = 3, 2, 5
        m = cppjit.gbl.pyzables.NakedBuffers(bsize, xval, yval)

        x = m.GetX()
        assert len(x) == bsize
        assert list(x) == list(map(lambda x: x * xval, range(bsize)))

        y = m.GetY()
        assert len(y) == bsize
        assert list(y) == list(map(lambda x: x * yval, range(bsize)))

    def test02_size_mapping_of_templated_method(self):
        """Use composites to map GetSize() onto buffer returns"""

        import cppjit

        def set_size(self, buf):
            buf.reshape((self.GetN(),))
            return buf

        cppjit.py.add_pythonization(
            cppjit.py.compose_method("NakedBuffers2.*Vector.*", "Get[XY]$", set_size),
            "pyzables",
        )

        bsize, xval, yval = 3, 2, 5
        m = cppjit.gbl.pyzables.NakedBuffers2[cppjit.gbl.pyzables.Vector](
            bsize, xval, yval
        )

        x = m.GetX()
        assert len(x) == bsize
        assert list(x) == list(map(lambda x: x * xval, range(bsize)))

        y = m.GetY()
        assert len(y) == bsize
        assert list(y) == list(map(lambda x: x * yval, range(bsize)))

    def test03_type_pinning(self):
        """Verify pinnability of returns"""

        import cppjit

        cppjit.gbl.pyzables.GimeDerived.__creates__ = True

        result = cppjit.gbl.pyzables.GimeDerived()
        assert type(result) == cppjit.gbl.pyzables.MyDerived

        cppjit.py.pin_type(cppjit.gbl.pyzables.MyBase)
        assert type(result) == cppjit.gbl.pyzables.MyDerived

    def test04_transparency(self):
        """Transparent use of smart pointers"""

        import cppjit

        Countable = cppjit.gbl.pyzables.Countable
        mine = cppjit.gbl.pyzables.mine

        assert type(mine) == Countable
        assert mine.m_check == 0xCDCDCDCD
        assert type(mine.__smartptr__()) == cppjit.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xCDCDCDCD
        assert mine.say_hi() == "Hi!"

    @mark.xfail(
        condition=IS_VALGRIND and IS_LINUX_ARM and IS_CLANG_REPL,
        run=False,
        reason="Crashes on Valgind Clang-Repl-ARM",
    )
    def test05_converters(self):
        """Smart pointer argument passing"""

        import cppjit

        pz = cppjit.gbl.pyzables
        mine = pz.mine

        assert 0xCDCDCDCD == pz.pass_mine_rp_ptr(mine)
        assert 0xCDCDCDCD == pz.pass_mine_rp_ref(mine)
        assert 0xCDCDCDCD == pz.pass_mine_rp(mine)

        assert 0xCDCDCDCD == pz.pass_mine_sp_ptr(mine)
        assert 0xCDCDCDCD == pz.pass_mine_sp_ref(mine)

        assert 0xCDCDCDCD == pz.pass_mine_sp_ptr(mine.__smartptr__())
        assert 0xCDCDCDCD == pz.pass_mine_sp_ref(mine.__smartptr__())

        assert 0xCDCDCDCD == pz.pass_mine_sp(mine)
        assert 0xCDCDCDCD == pz.pass_mine_sp(mine.__smartptr__())

        # TODO:
        # cppjit.gbl.mine = mine
        pz.renew_mine()

    @mark.xfail(
        condition=IS_VALGRIND and IS_LINUX_ARM and IS_CLANG_REPL,
        run=False,
        reason="Fails with Valgrind with Clang-Repl ARM",
    )
    def test06_executors(self):
        """Smart pointer return types"""

        import cppjit

        pz = cppjit.gbl.pyzables
        Countable = pz.Countable

        mine = pz.gime_mine_ptr()
        assert type(mine) == Countable
        assert mine.m_check == 0xCDCDCDCD
        assert type(mine.__smartptr__()) == cppjit.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xCDCDCDCD
        assert mine.say_hi() == "Hi!"

        mine = pz.gime_mine_ref()
        assert type(mine) == Countable
        assert mine.m_check == 0xCDCDCDCD
        assert type(mine.__smartptr__()) == cppjit.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xCDCDCDCD
        assert mine.say_hi() == "Hi!"

        mine = pz.gime_mine()
        assert type(mine) == Countable
        assert mine.m_check == 0xCDCDCDCD
        assert type(mine.__smartptr__()) == cppjit.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xCDCDCDCD
        assert mine.say_hi() == "Hi!"

    def test07_creates_flag(self):
        """Effect of creates flag on return type"""

        import gc

        import cppjit

        pz = cppjit.gbl.pyzables
        Countable = pz.Countable

        gc.collect()
        oldcount = Countable.sInstances  # there's eg. one global variable

        pz.gime_naked_countable.__creates__ = True
        for i in range(10):
            cnt = pz.gime_naked_countable()
            gc.collect()
            assert Countable.sInstances == oldcount + 1
            del cnt
        gc.collect()

        assert Countable.sInstances == oldcount

    def test08_base_class_pythonization(self):
        """Derived class should not re-pythonize base class pythonization"""

        import cppjit

        d = cppjit.gbl.pyzables.IndexableDerived()

        assert d[0] == 42
        assert d[-1] == 42
        # skip the IndexErorr test: pythonization for __getitem__[index] < size()
        # can not be applied strict enough (instead of an index, this could be an
        # associative  container, with 'index' a key, not a counter
        # raises(IndexError, d.__getitem__, 1)

    @mark.xfail(condition=IS_MAC, reason="Fails on OS X")
    def test09_cpp_side_pythonization(self):
        """Use of C++ side pythonizations"""

        import cppjit

        # explicit pythonization
        for kls in [
            cppjit.gbl.pyzables.WithCallback1,
            cppjit.gbl.pyzables.WithCallback2,
            # backwards compatibility with the old cppyy names
            cppjit.gbl.pyzables.WithCppyyCallback1,
            cppjit.gbl.pyzables.WithCppyyCallback2,
        ]:
            w = kls(42)
            assert hasattr(w, "GetInt")
            assert not hasattr(w, "get_int")
            assert w.GetInt() == 42

            assert hasattr(w, "SetInt")
            assert not hasattr(w, "set_int")
            w.SetInt(17)
            assert w.GetInt() == 17

            assert kls.klass_name == kls.__cpp_name__

        # up-the-hierarchy pythonization
        w = cppjit.gbl.pyzables.WithCallback3(42)
        assert hasattr(w, "GetInt")
        assert not hasattr(w, "get_int")
        assert w.GetInt() == 2 * 42

        assert hasattr(w, "SetInt")
        assert not hasattr(w, "set_int")
        w.SetInt(17)
        assert w.GetInt() == 4 * 17

        assert cppjit.gbl.pyzables.WithCallback2.klass_name == "pyzables::WithCallback3"

    @mark.xfail(condition=IS_MAC, reason="Fails on OS X")
    def test10_shared_ptr_reset(self):
        """Checks that smart pointer types are Pythonized with the special
        __smartptr__ member that can also be used to reset the underlying smart
        pointer."""

        import cppjit

        optr = cppjit.gbl.std.make_shared["std::string"]("hello smart pointer")
        o2 = cppjit.gbl.std.string()
        cppjit._backend.SetOwnership(
            o2, False
        )  # This object will be owned by the smart pointer
        optr.__smartptr__().reset(o2)
        assert optr == o2

    def test11_kwargs_struct(self):
        """A parameter named kwargs that takes an aggregate collects the
        keyword arguments that match no other parameter"""

        import cppjit

        cppjit.cppdef("""\
        namespace kwargs_struct {
        struct Range {
          std::string desc{"unset"};
          Range() = default;
          Range(const char* name, bool adjust = true) : desc(std::string{name} + (adjust ? "+" : "")) {}
          Range(double lo, double hi) : desc(std::to_string(int(lo)) + ":" + std::to_string(int(hi))) {}
          Range(std::nullptr_t) : desc("none") {}
        };
        struct Inner { bool Save = false; int Level = 7; };
        struct Nested {
          std::string desc{"unset"};
          Nested() = default;
          Nested(Inner const& kwargs) : desc((kwargs.Save ? "save," : "") + std::to_string(kwargs.Level)) {}
        };
        struct Options { Range Range; int PrintLevel = 1; Nested Fit; Inner Agg; };

        std::string str(Options const& o) {
          return o.Range.desc + "|" + std::to_string(o.PrintLevel) + "|" + o.Fit.desc + "|"
                 + std::to_string(o.Agg.Save) + "," + std::to_string(o.Agg.Level);
        }
        std::string func(int a, Options const& kwargs, const char* extra = "") {
          return std::to_string(a) + "|" + str(kwargs) + "|" + extra;
        }
        std::string func(int a, double b) { return "other"; }

        std::string pick(int a) { return "plain"; }
        std::string pick(int a, Options const& kwargs) { return "kwargs"; }

        template <class... Args>
        std::string tmpl(int a, Args const&... args) { return "template"; }
        std::string tmpl(int a, Options const& kwargs, int b = 0) { return "kwargs|" + str(kwargs); }

        struct Klass {
          std::string desc;
          Klass() = default;
          Klass(int a, Options kwargs) : desc(func(a, kwargs)) {}
          std::string meth(Options const& kwargs, int b = 0) { return std::to_string(b) + "|" + str(kwargs); }
        };

        struct Name { std::string s; Name(const char* n) : s(n) {} };
        std::string named(int a) { return "int"; }
        std::string named(Name const& name, Options const& kwargs) {
          if (kwargs.PrintLevel < 0) throw std::invalid_argument("negative PrintLevel");
          return name.s + "|" + str(kwargs);
        }

        struct Undefined;
        std::string incomplete(int a, Undefined const& kwargs) { return "incomplete"; }
        }""")

        ns = cppjit.gbl.kwargs_struct

        # no keywords: default-constructed struct, with default member initializers
        assert ns.func(1) == "1|unset|1|unset|0,7|"

        # members of class type follow C++ initializer syntax
        assert ns.func(1, Range="r1", PrintLevel=-1) == "1|r1+|-1|unset|0,7|"
        assert ns.func(1, Range=("r1", False)) == "1|r1|1|unset|0,7|"
        assert ns.func(1, Range=(0, 10)) == "1|0:10|1|unset|0,7|"
        assert ns.func(1, Range=None) == "1|none|1|unset|0,7|"
        assert ns.func(1, Fit=dict(Save=True, Level=3)) == "1|unset|1|save,3|0,7|"
        assert ns.func(1, Agg={"Save": True}) == "1|unset|1|unset|1,7|"

        # the kwargs parameter is keyword-only, other parameters take keywords as usual
        assert ns.func(1, "x", PrintLevel=2) == "1|unset|2|unset|0,7|x"
        assert ns.func(1, extra="x", PrintLevel=2) == "1|unset|2|unset|0,7|x"
        assert ns.func(1, kwargs=ns.Options()) == "1|unset|1|unset|0,7|"
        assert ns.func(1, 2.0) == "other"

        # methods and constructors
        assert ns.Klass().meth(5, PrintLevel=3) == "5|unset|3|unset|0,7"
        assert ns.Klass().meth(b=5) == "5|unset|1|unset|0,7"
        assert ns.Klass(1, Range="r").desc == "1|r+|1|unset|0,7|"

        # calls without keywords prefer the other overloads
        assert ns.pick(1) == "plain"
        assert ns.pick(1, PrintLevel=2) == "kwargs"

        # overload sets with templates
        assert ns.tmpl(1, 2.0) == "template"
        assert ns.tmpl(1, PrintLevel=3) == "kwargs|unset|3|unset|0,7"

        # errors of the overload that took the keywords have priority
        with raises(AttributeError, match="no data member"):
            ns.func(1, Bogus=True)
        with raises(AttributeError, match="no data member"):
            ns.tmpl(1, Bogus=True)
        with raises(AttributeError, match="no data member"):
            ns.func(1, Fit=dict(Bogus=True))
        with raises(TypeError, match="keyword argument 'Range'"):
            ns.func(1, Range=object())
        with raises(TypeError):
            ns.func(1, "x", "y")

        # also if the overload takes the arguments only with implicit conversions
        assert ns.named("n", PrintLevel=2) == "n|unset|2|unset|0,7"
        with raises(AttributeError, match="no data member"):
            ns.named("n", Bogus=True)
        with raises(cppjit.gbl.std.invalid_argument, match="negative PrintLevel"):
            ns.named("n", PrintLevel=-1)

        # a kwargs parameter of incomplete type is an ordinary parameter
        with raises(TypeError):
            ns.incomplete(1, PrintLevel=2)


## actual test run
if __name__ == "__main__":
    result = run_pytest(__file__)  # noqa: F821
    sys.exit(result)
