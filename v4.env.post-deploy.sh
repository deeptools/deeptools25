#!/usr/bin/env bash
set -euxo pipefail

# Find the actual path to libclang.so
LIBCLANG_SO=$(find "$CONDA_PREFIX/lib" -name "libclang.so*" | head -n 1)
if [ -z "$LIBCLANG_SO" ]; then
    echo "ERROR: Could not find libclang.so in $CONDA_PREFIX"
    exit 1
fi
LIBCLANG_DIR=$(dirname "$LIBCLANG_SO")

# Ensure we have a proper libclang.so symlink that points to the right version
# Check if we have libclang-cpp.so.19.1
if [ -f "$LIBCLANG_DIR/libclang-cpp.so.19.1" ]; then
    # Create correct symlinks if they don't exist
    if [ ! -f "$LIBCLANG_DIR/libclang.so.19" ]; then
        ln -sf "$LIBCLANG_DIR/libclang.so.13" "$LIBCLANG_DIR/libclang.so.19"
    fi
    # Make sure the main symlink points to version 19
    ln -sf "$LIBCLANG_DIR/libclang.so.19" "$LIBCLANG_DIR/libclang.so"
fi

# Set environment variables
export C_INCLUDE_PATH="$CONDA_PREFIX/include:${C_INCLUDE_PATH:-}"
export CPLUS_INCLUDE_PATH="$CONDA_PREFIX/include:${CPLUS_INCLUDE_PATH:-}"
export LIBCLANG_PATH="$LIBCLANG_DIR"  # Don't append, replace
export LLVM_CONFIG_PATH="$CONDA_PREFIX/bin/llvm-config"
export CFLAGS="-I$CONDA_PREFIX/include -MD"
export CPPFLAGS="-I$CONDA_PREFIX/include"
export LDFLAGS="-L$CONDA_PREFIX/lib"

# For Rust's bindgen
export BINDGEN_EXTRA_CLANG_ARGS="-I$CONDA_PREFIX/include"
# Make rustc look in conda lib paths
export RUSTFLAGS="-L $CONDA_PREFIX/lib"

# Debug information
echo "=== Environment Variables ==="
env | grep -E "LIBCLANG|CFLAGS|CPPFLAGS|LDFLAGS|INCLUDE_PATH|LLVM|BINDGEN|RUST"
echo "=== libclang location ==="
ls -l $LIBCLANG_DIR/libclang*

# Install dependencies
pip install git+https://github.com/deeptools/deepTools.git@4.0.0
