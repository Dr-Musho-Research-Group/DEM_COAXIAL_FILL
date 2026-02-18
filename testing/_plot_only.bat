@echo off
setlocal enabledelayedexpansion

echo Removing old xyz files...
del *.png >nul 2>&1
del *.mp4 >nul 2>&1

echo Running plot script...
REM Run the Python plot script
py _atoms_plot_all.py --pattern "atoms_*.xyz" --color element --fps 10

echo Process complete!
pause
