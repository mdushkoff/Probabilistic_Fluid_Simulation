@echo off
rem 1) Call the native tools script to set up VC++/clang/nvcc environment
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

rem 2) Force the PATH to use new CMake version first
set PATH=C:\Program Files\CMake\bin;%PATH%