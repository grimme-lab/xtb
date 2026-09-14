@echo off
rem ============================================================================
rem  One-time xtbx activation:
rem    1. Add this launcher folder to the Windows user PATH.
rem    2. Auto-detect CPU cores/RAM and write xtbg.conf.
rem
rem  After this, open a new PowerShell/CMD window and simply run:
rem      xtbx molecule.xyz --gpu
rem ============================================================================
setlocal
set "XTBX_BIN=%~dp0"
if "%XTBX_BIN:~-1%"=="\" set "XTBX_BIN=%XTBX_BIN:~0,-1%"

powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$bin=[IO.Path]::GetFullPath('%XTBX_BIN%');" ^
  "$user=[Environment]::GetEnvironmentVariable('Path','User');" ^
  "$parts=@($user -split ';' | Where-Object { $_ });" ^
  "if (-not ($parts | Where-Object { [IO.Path]::GetFullPath($_).TrimEnd('\') -ieq $bin.TrimEnd('\') })) {" ^
  "  [Environment]::SetEnvironmentVariable('Path',(($parts + $bin) -join ';'),'User');" ^
  "  Write-Host ('Added to user PATH: ' + $bin)" ^
  "} else { Write-Host ('Already on user PATH: ' + $bin) }"

if errorlevel 1 (
  echo ERROR: Could not update the Windows user PATH.
  exit /b 1
)

wsl.exe bash /mnt/e/Prasanna/xTB/win/xtb-setup.sh
if errorlevel 1 exit /b 1

echo.
echo xtbx activation complete.
echo Open a NEW terminal, then use:  xtbx molecule.xyz --gpu
endlocal
