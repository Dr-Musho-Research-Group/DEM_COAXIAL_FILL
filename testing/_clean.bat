@echo off

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
