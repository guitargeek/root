import os
import sys
import subprocess

MYHOME = os.path.dirname(__file__)

def main():
    os.environ['LD_LIBRARY_PATH'] = os.path.join(MYHOME, 'lib')
    rootexe = os.path.join(MYHOME, 'bin', 'root.exe')

    args = [rootexe, *sys.argv[1:]]

    return subprocess.call(args)


if __name__ == "__main__":
    sys.exit(main())
