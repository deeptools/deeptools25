#!/usr/bin/env bash
set -o pipefail

# Find the actual path to libclang.so
LIBCLANG_SO=$(find "$CONDA_PREFIX/lib" -name "libclang.so*" | head -n 1)
if [ -z "$LIBCLANG_SO" ]; then
    echo "ERROR: Could not find libclang.so in $CONDA_PREFIX"
    exit 1
fi
LIBCLANG_DIR=$(dirname "$LIBCLANG_SO")

# Set environment variables
export C_INCLUDE_PATH="$CONDA_PREFIX/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH="$CONDA_PREFIX/include:$CPLUS_INCLUDE_PATH"
export LIBCLANG_PATH="$LIBCLANG_DIR:$LIBCLANG_PATH"
export LLVM_CONFIG_PATH="$CONDA_PREFIX/bin/llvm-config"
export CFLAGS="-I$CONDA_PREFIX/include -MD"  # -MD: enables dependency generation
export CPPFLAGS="-I$CONDA_PREFIX/include"
export LDFLAGS="-L$CONDA_PREFIX/lib"

# Debug information
echo "=== Environment Variables ==="
env | grep -E "LIBCLANG|CFLAGS|CPPFLAGS|LDFLAGS|INCLUDE_PATH|LLVM"
echo "=== libclang location ==="
ls -l $LIBCLANG_DIR/libclang*

# Install dependencies
pip install git+https://github.com/deeptools/deepTools.git@4.0.0