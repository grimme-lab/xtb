#!/bin/bash

# Simple build script for xTB on macOS with Homebrew GCC
# Fixes the -lSystem linking issue

set -e

echo "Building xTB with macOS-compatible settings..."

# Export environment variables for macOS linking
export LDFLAGS="-Wl,-framework,CoreFoundation -Wl,-framework,Security"
export FC=gfortran-12
export CC=gcc-12
export CXX=g++-12

# Create build directory
rm -rf build_macos
mkdir build_macos
cd build_macos

# Try to use meson with custom compiler settings
echo "Configuring with meson..."
meson setup .. \
  --buildtype=release \
  -Dfortran_link_args="-Wl,-framework,CoreFoundation,-framework,Security" \
  --cross-file=../meson_macos.ini

echo "Building..."
meson compile

echo "Build completed!"
ls -la xtb*