#ifndef CPYCPPYY_INTERPRETERLOCK_H
#define CPYCPPYY_INTERPRETERLOCK_H

// This header assumes Python.h (for the GIL API) and Cppyy.h (for the lock API)
// have already been included; both come in via CPyCppyy.h.
#include "Cppyy.h"


namespace CPyCppyy {

// RAII that holds the global interpreter lock (Cppyy::GetGlobalMutex) for the
// duration of a scope, releasing the GIL while it blocks waiting for it.
//
// The GIL used to serialize every PyROOT call into ROOT; on a free-threaded
// ("GIL-less") Python build it no longer does, so cppyy takes this lock around
// any work that drives the thread-hostile interpreter or touches state shared
// between threads. Dropping the GIL while waiting is essential: a thread parked
// on the lock must remain at a Python safe point, otherwise it would stall a
// stop-the-world pause (e.g. cyclic GC) that the lock-holding thread may need to
// complete before it can release the lock -- which would dead-lock the process.
class InterpreterLockGuard {
public:
    InterpreterLockGuard() {
#ifdef WITH_THREAD
        PyThreadState* save = PyEval_SaveThread();
        Cppyy::LockInterpreter();
        PyEval_RestoreThread(save);
#else
        Cppyy::LockInterpreter();
#endif
    }
    ~InterpreterLockGuard() { Cppyy::UnlockInterpreter(); }

    InterpreterLockGuard(const InterpreterLockGuard&) = delete;
    InterpreterLockGuard& operator=(const InterpreterLockGuard&) = delete;
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_INTERPRETERLOCK_H
