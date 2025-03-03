#!env bash
set -o pipefail
export C_INCLUDE_PATH="$CONDA_PREFIX/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH="$CONDA_PREFIX/include:$CPLUS_INCLUDE_PATH"
export LIBCLANG_PATH="$CONDA_PREFIX/lib:$LIBCLANG_PATH"
#export CFLAGS="-I$CONDA_PREFIX/include -I/usr/include"
export CFLAGS="-I $CONDA_PREFIX/include -MD"  # -MD: enables dependency generation
export CPPFLAGS="-I $CONDA_PREFIX/include"
export LDFLAGS="-L $CONDA_PREFIX/lib"
env | grep -E "LIBCLANG|CFLAGS|CPPFLAGS|LDFLAGS|INCLUDE_PATH"
pip install git+https://github.com/deeptools/deepTools.git@4.0.0
#exit 1
