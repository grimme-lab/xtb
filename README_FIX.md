# xTB v6.7.1 - Fortran Format String Bug Fix

This repository contains a fixed version of xTB v6.7.1 that resolves a critical Fortran format string bug that causes runtime crashes during geometry optimization.

## The Bug

**Location**: `src/optimizer.f90:853` (in original v6.7.1)

**Original Buggy Code**:
```fortran
write(env%unit,'(1x,"("f7.2"%)")')       (depred-echng)/echng*100
```

**Error**: Missing commas around the `f7.2` format descriptor, causing a Fortran runtime error:
```
"Missing comma between descriptors"
```

## The Fix

**Fixed Code**:
```fortran
write(env%unit,'(1x,"(",f7.2,"%)")')       (depred-echng)/(echng+1e-34_wp)*100
```

**Changes Made**:
1. Added missing commas around `f7.2`: `"("f7.2"%)"` → `"(",f7.2,"%)"` 
2. Improved division safety by adding small epsilon: `(depred-echng)/echng*100` → `(depred-echng)/(echng+1e-34_wp)*100`

## Impact

This bug affects:
- **Geometry optimization** (`--opt` flag)
- **Optimization output display** showing energy change predictions
- **Any workflow** that uses xTB optimization features

The bug manifests as a runtime crash when xTB tries to display percentage values during optimization cycles.

## Verification

The fix can be verified by running geometry optimization:

```bash
# Create test molecule
echo "5
Methane molecule
C    0.0000   0.0000   0.0000
H    1.0890   0.0000   0.0000
H   -0.3630   1.0276   0.0000
H   -0.3630  -0.5138   0.8900
H   -0.3630  -0.5138  -0.8900" > test.xyz

# Test optimization (should show percentage predictions without crash)
xtb test.xyz --opt
```

**Expected Output** (with fix):
```
predicted    -0.1224011E-03 ( -13.13%)
```

**Original Behavior** (without fix):
- Runtime crash with format string error

## Build Instructions

### macOS (with homebrew)

```bash
# Prerequisites
brew install cmake gfortran openblas

# Configure and build
export MACOSX_DEPLOYMENT_TARGET=14.0
export SDKROOT=$(xcrun --show-sdk-path)
cmake -B build -DCMAKE_Fortran_COMPILER=gfortran \
      -DBLAS_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib \
      -DLAPACK_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib

# Build
cmake --build build --parallel

# Test
./build/xtb --version
```

### Alternative macOS Build (using provided script)

```bash
./build_macos.sh
```

## Files Added for macOS Support

- `build_macos.sh`: Automated build script for macOS
- `meson_macos.ini`: Meson cross-compilation configuration for macOS
- This `README_FIX.md`: Documentation of the fix

## Original xTB Information

This is based on the official xTB v6.7.1 from the Grimme Lab:
- **Original Repository**: https://github.com/grimme-lab/xtb
- **Version**: 6.7.1
- **Commit**: 7618f60692ee3cb204b4dbac4961c15617aa4eb1

## Citation

If you use this fixed version, please cite the original xTB work:

```
C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht,
J. Seibert, S. Spicher, S. Grimme, WIREs Comput. Mol. Sci., 2020, 11,
e01493. DOI: 10.1002/wcms.1493
```

## License

This maintains the same LGPL v3 license as the original xTB project.

## Status

✅ **Fixed**: Fortran format string bug in optimizer.f90
✅ **Tested**: Geometry optimization works correctly
✅ **Compatible**: Drop-in replacement for xTB v6.7.1