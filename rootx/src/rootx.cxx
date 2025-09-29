// @(#)root/rootx:$Id$
// Author: Fons Rademakers   19/02/98

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Rootx                                                                //
//                                                                      //
// Rootx is a small front-end program that starts the main ROOT module. //
// This program is called "root" in the $ROOTSYS/bin directory and the  //
// real ROOT executable is now called "root.exe" (formerly "root").     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"
#include "Rtypes.h"
#include "snprintf.h"
#include "rootCommandLineOptionsHelp.h"

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <csignal>
#include <sys/wait.h>
#include <sys/stat.h>
#include <cerrno>
#include <netdb.h>
#include <sys/socket.h>
#include <string>
#ifdef __APPLE__
#include <AvailabilityMacros.h>
#include <mach-o/dyld.h>
#endif

#ifdef R__FBSD
#include <sys/param.h>
#include <sys/user.h>
#include <sys/types.h>
#include <libutil.h>
#include <libprocstat.h>
#endif // R__FBSD

static int gChildpid;
static int GetErrno()
{
   return errno;
}

static void ResetErrno()
{
   errno = 0;
}

extern "C" {
   static void SigTerm(int);
}

static void SigTerm(int sig)
{
   // When we get terminated, terminate child, too.
   kill(gChildpid, sig);
}

static void WaitChild()
{
   // Wait till child (i.e. ROOT) is finished.

   int status;

   do {
      while (waitpid(gChildpid, &status, WUNTRACED) < 0) {
         if (GetErrno() != EINTR)
            break;
         ResetErrno();
      }

      if (WIFEXITED(status))
         exit(WEXITSTATUS(status));

      if (WIFSIGNALED(status))
         exit(WTERMSIG(status));

      if (WIFSTOPPED(status)) {         // child got ctlr-Z
         raise(SIGTSTP);                // stop also parent
         kill(gChildpid, SIGCONT);       // if parent wakes up, wake up child
      }
   } while (WIFSTOPPED(status));

   exit(0);
}

static void PrintUsage()
{
   fprintf(stderr, kCommandLineOptionsHelp);
}

int main(int argc, char **argv)
{
   // In batch mode don't show splash screen, idem for no logo mode,
   // in about mode show always splash screen
   int i;
   for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "-?") || !strncmp(argv[i], "-h", 2) ||
          !strncmp(argv[i], "--help", 6)) {
         PrintUsage();
         return 0;
      }
   }

   // Ignore SIGINT and SIGQUIT. Install handler for SIGUSR1.

   struct sigaction ignore, actTerm,
      saveintr, savequit, saveterm;

#if defined(__sun) && !defined(__i386) && !defined(__SVR4)
   ignore.sa_handler = (void (*)())SIG_IGN;
#elif defined(__sun) && defined(__SVR4)
   ignore.sa_handler = (void (*)(int))SIG_IGN;
#else
   ignore.sa_handler = SIG_IGN;
#endif
   sigemptyset(&ignore.sa_mask);
   ignore.sa_flags = 0;

   actTerm = ignore;
#if defined(__sun) && !defined(__i386) && !defined(__SVR4)
   actTerm.sa_handler = (void (*)())SigTerm;
#elif defined(__sun) && defined(__SVR4)
   actTerm.sa_handler = SigTerm;
#elif (defined(__sgi) && !defined(__KCC)) || defined(__Lynx__)
#   if defined(IRIX64) || (__GNUG__>=3)
   actTerm.sa_handler = SigTerm;
#   else
   actTerm.sa_handler = (void (*)(...))SigTerm;
#   endif
#else
   actTerm.sa_handler = SigTerm;
#endif
   sigaction(SIGINT,  &ignore, &saveintr);
   sigaction(SIGQUIT, &ignore, &savequit);
   sigaction(SIGTERM, &actTerm, &saveterm);

   // Create child...

   if ((gChildpid = fork()) < 0) {
      fprintf(stderr, "%s: error forking child\n", argv[0]);
      return 1;
   } else if (gChildpid > 0) {
      WaitChild();
   }

   // Continue with child...

   // Restore original signal actions
   sigaction(SIGINT,  &saveintr, 0);
   sigaction(SIGQUIT, &savequit, 0);
   sigaction(SIGTERM, &saveterm, 0);

   // Child is going to overlay itself with the actual ROOT module...

   // Build argv vector
   char  arg0[] = "root.exe\0";

   std::vector<char *> argvv(argc+1);
   argvv[0] = arg0;

   for (i = 1; i < argc; i++)
      argvv[i] = argv[i];
   argvv[i] = 0;

   // Execute actual ROOT module
   execv("root.exe", argvv.data());

   // Exec failed
   fprintf(stderr, "%s: can't start ROOT -- check that root.exe exists!\n",
           argv[0]);

   return 1;
}
