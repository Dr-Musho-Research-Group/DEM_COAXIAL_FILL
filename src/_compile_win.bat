@echo off
REM Set the path to your C++ compiler
set C_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC\14.39.33519\bin\Hostx64\x64
set PATH=%C_PATH%;%PATH%

echo Setting up environment...
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64

echo Compiling with CL...
cl /EHsc /O2 /openmp:llvm /std:c++17 /Fe:coax_pack_cpu.exe coax_pack_cpu.cpp

if %ERRORLEVEL%==0 (
    echo Compilation successful! Executable: coax_pack_cpu.exe 
) else (
    echo Compilation failed!
)

pause