#!/usr/bin/env bash
# One-time install of the NVIDIA HPC SDK (nvfortran + cuSolver/cuBLAS) in WSL.
# Run this in an interactive WSL terminal so sudo can prompt for your password:
#
#     bash /mnt/e/Prasanna/xTB/xtb/tools/gpu-dev/diagnostics/install_nvhpc.sh
#
# It lands in /opt/nvidia/hpc_sdk (inside the WSL distro, which is on E:).
set -euo pipefail

echo "==> Adding NVIDIA HPC SDK apt repository"
curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK \
  | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' \
  | sudo tee /etc/apt/sources.list.d/nvhpc.list >/dev/null

echo "==> apt-get update"
sudo apt-get update

echo "==> Selecting the latest nvhpc-XX-Y metapackage"
PKG=$(apt-cache search '^nvhpc-' \
      | awk '{print $1}' \
      | grep -E '^nvhpc-[0-9]+-[0-9]+$' \
      | sort -V | tail -1)
if [ -z "${PKG}" ]; then
  echo "ERROR: no nvhpc-XX-Y package found in the repo. 'apt-cache search nvhpc' output:"
  apt-cache search nvhpc || true
  exit 1
fi
echo "    -> installing ${PKG}  (a few GB download; ~10-20 min)"
sudo apt-get install -y "${PKG}"

echo
echo "==> Install finished. nvfortran location(s):"
ls -d /opt/nvidia/hpc_sdk/Linux_*/*/compilers/bin 2>/dev/null || true

echo
echo "==> Quick check:"
BIN=$(ls -d /opt/nvidia/hpc_sdk/Linux_*/*/compilers/bin 2>/dev/null | sort -V | tail -1 || true)
if [ -n "${BIN}" ]; then
  "${BIN}/nvfortran" --version 2>/dev/null | head -2 || echo "nvfortran present but --version failed"
  echo
  echo "ALL SET. Tell Claude it's done. (PATH for this shell: export PATH=${BIN}:\$PATH)"
else
  echo "Could not locate nvfortran under /opt/nvidia/hpc_sdk — check the install output above."
fi
