@echo off
rem ============================================================================
rem  xtbx  -  one unified xtb command (auto-picks the right engine).
rem
rem    xtbx mol.xyz --gfn 2 --opt        one molecule, full output
rem    xtbx <folder> --gfn 2 --gpu       persistent GPU batch pool (8 workers)
rem    xtbx <folder> --gfn 0 --gpu       native cross-molecule GPU batch
rem    xtbx <folder> --gfn 0 --gpu --gpu-batch
rem    xtbx <folder> --gfn 2 --gpu --gpu-batch
rem                                      native one-process xTB GPU batch
rem    xtbx *.xyz --gfn 0                parallel per-compound over the files
rem
rem  Path args = files/folders that exist (or globs); everything else (--gfn 0,
rem  --opt, --chrg, ...) is passed straight to xtb. Subsumes xtbg/xtbfolder/xtbgpu.
rem ============================================================================
setlocal
set "args=%*"
if defined args set "args=%args:\=/%"
set "envfwd="
if defined XTB_JOBS set "envfwd=XTB_JOBS=%XTB_JOBS%"
if defined XTB_OMP set "envfwd=%envfwd% XTB_OMP=%XTB_OMP%"
if defined XTB_GPU_JOBS set "envfwd=%envfwd% XTB_GPU_JOBS=%XTB_GPU_JOBS%"
if defined XTB_GPU_DEVICES set "envfwd=%envfwd% XTB_GPU_DEVICES=%XTB_GPU_DEVICES%"
if defined XTB_QUEUE_TRACE set "envfwd=%envfwd% XTB_QUEUE_TRACE=%XTB_QUEUE_TRACE%"
wsl.exe env %envfwd% bash /mnt/e/Prasanna/xTB/win/xtbx_run.sh %args%
endlocal
