@echo off
setlocal EnableDelayedExpansion

REM ================================================================
REM  Batch script to run coax_pack_cpu.exe simulation
REM  Similar format to md_ess.bat but adapted for 3D coaxial fill DEM
REM ================================================================

REM Record start time
for /f "tokens=1-4 delims=:." %%a in ("%TIME%") do (
    set /a START_H=%%a
    set /a START_M=%%b
    set /a START_S=%%c
    set /a START_MS=%%d
)

REM -------------------- User Parameters ---------------------------
set DEBUG_LEVEL=1
set NATOMS_MAX=25000
set DT=1e-6
set NITER=400000
set DUMP_INTERVAL=2500
set SEED=42

REM Geometry [m]
set RIN=22e-6
set ROUT=40e-6
set LENGTH=380e-6

REM Physics
set FLUX=60000
set GRAVITY=9.81
set SHAKE_FREQ=300
set SHAKE_AMP=0.1e-6

REM Process stages [s]
set FILL_TIME=8.0
set RAM_START=0.0
set RAM_DURATION=0.0
set RAM_SPEED=0.0

REM Volume Fraction Target
set VF=0.5
REM ----------------------------------------------------------------

cd /d "%~dp0"

echo ===============================================================
echo Cleaning previous output files...
del *.xyz >nul 2>&1
del *.dat >nul 2>&1
del *.vtk >nul 2>&1
del *.png >nul 2>&1
del *.mp4 >nul 2>&1
del *.csv >nul 2>&1
echo ===============================================================

echo Running coax_pack_cpu.exe ...
echo Parameters:
echo   NATOMS_MAX     = %NATOMS_MAX%
echo   DT             = %DT%
echo   NITER          = %NITER%
echo   DUMP_INTERVAL  = %DUMP_INTERVAL%
echo   DEBUG_LEVEL    = %DEBUG_LEVEL%
echo   SEED           = %SEED%
echo   Rin, Rout, L   = %RIN%, %ROUT%, %LENGTH%
echo   Flux           = %FLUX%
echo   Gravity        = %GRAVITY%
echo   Shake (Hz, Amp)= %SHAKE_FREQ%, %SHAKE_AMP%
echo   Fill/Ram       = %FILL_TIME%s / %RAM_START%s to %RAM_DURATION%s @ %RAM_SPEED%m/s
echo   VF Target      = %VF%
echo ===============================================================

REM -------------------- Run the executable ------------------------
..\src\coax_pack_cpu.exe ^
  %NATOMS_MAX% ^
  %DT% ^
  %NITER% ^
  %DUMP_INTERVAL% ^
  %DEBUG_LEVEL% ^
  %SEED% ^
  %RIN% ^
  %ROUT% ^
  %LENGTH% ^
  %FLUX% ^
  %GRAVITY% ^
  %SHAKE_FREQ% ^
  %SHAKE_AMP% ^
  %FILL_TIME% ^
  %RAM_START% ^
  %RAM_DURATION% ^
  %RAM_SPEED% ^
  %VF%

REM -------------------- Timing and Status -------------------------
for /f "tokens=1-4 delims=:." %%a in ("%TIME%") do (
    set /a END_H=%%a
    set /a END_M=%%b
    set /a END_S=%%c
    set /a END_MS=%%d
)

timeout /t 1 /nobreak >nul

set /a START_TOTAL_MS=(%START_H%*3600000)+(%START_M%*60000)+(%START_S%*1000)+%START_MS%
set /a END_TOTAL_MS=(%END_H%*3600000)+(%END_M%*60000)+(%END_S%*1000)+%END_MS%
set /a ELAPSED_MS=%END_TOTAL_MS% - %START_TOTAL_MS%
set /a ELAPSED_SEC=%ELAPSED_MS% / 1000
set /a ELAPSED_MS_REMAINDER=%ELAPSED_MS% %% 1000

echo ===============================================================
echo Simulation completed in %ELAPSED_SEC%.%ELAPSED_MS_REMAINDER% seconds.

if %ERRORLEVEL%==0 (
    echo Simulation completed successfully!
) else (
    echo Simulation failed! Error code: %ERRORLEVEL%
    pause
    exit /b %ERRORLEVEL%
)
echo ===============================================================

REM -------------------- Postprocessing ----------------------------
echo Renaming output files...
rename "output_step_*.xyz" "atoms_*.xyz" >nul 2>&1

echo Running Python plot script if available...
if exist "_atoms_plot_all.py" (
    py.exe _atoms_plot_all.py
) else (
    echo (No plot script found, skipping.)
)

echo ===============================================================
echo Process complete!
pause
