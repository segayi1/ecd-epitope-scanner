@echo off
chcp 65001 >nul
setlocal
title Stop ECD Epitope Workbench

set PORT=8765
echo.
echo   Stopping the workbench server listening on port %PORT% ...
set FOUND=0
for /f "tokens=5" %%a in ('netstat -ano ^| findstr /r /c:"LISTENING" ^| findstr /c:":%PORT% "') do (
    set FOUND=1
    echo   Terminating PID %%a
    taskkill /PID %%a /F >nul
)
if "%FOUND%"=="0" echo   Nothing was listening on port %PORT%.
echo.
echo   Done.
pause
