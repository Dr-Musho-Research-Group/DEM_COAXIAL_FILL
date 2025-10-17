@echo off
setlocal enabledelayedexpansion

echo Removing old xyz files...
del *.png >nul 2>&1
del *.mp4 >nul 2>&1

echo Running plot script...
REM Run the Python plot script
py.exe _atoms_plot_all.py

echo Process complete!
pause
