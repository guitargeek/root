#!/bin/bash

# Inspired by:
# https://github.com/root-project/root-conda-nightly/blob/main/test_conda_root.sh

source ~/.virtualenvs/pip-install-root/bin/activate

ROOTSYS=$(python -c "import pathlib, ROOT; print(pathlib.Path(ROOT.__file__).parent.resolve())")

ROOTTEST_SRC=~/tmp/pip-tests/roottest

cmake -DROOT_DIR=$ROOTSYS/cmake $ROOTTEST_SRC

cmake --build . -j$(nproc)


#--- run tests ---#
ctest -T test --no-compress-output -j$(nproc) || true  # ignore ctest exit code, we will parse the logs
