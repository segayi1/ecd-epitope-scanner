@echo off
chcp 65001 >nul
setlocal
title ECD Epitope Workbench

rem ============================================================
rem  Double-click launcher. This file lives INSIDE the package
rem  folder, so the whole folder is portable: move or copy it
rem  anywhere on this machine and double-click this file.
rem ============================================================

set "PKG_DIR=%~dp0"
rem `python -m ecd_epitope_scanner` resolves the package from the PARENT
rem directory, so run from there while the venv stays inside the folder.
cd /d "%~dp0.."

set "PYTHON_EXE=%PKG_DIR%.venv\Scripts\python.exe"
set "HOST=127.0.0.1"
set "PORT=8765"
set "URL=http://%HOST%:%PORT%"

echo.
echo   ECD Epitope Workbench
echo   =====================
echo   Interpreter : %PYTHON_EXE%
echo   Address     : %URL%
echo.

if not exist "%PYTHON_EXE%" goto novenv

rem ---- running already? open the browser and quit -----------------
netstat -ano | findstr /r /c:"LISTENING" | findstr /c:":%PORT% " >nul
if %errorlevel%==0 goto already

echo   Starting the server in a minimized window...
start "ECD Workbench Server" /min "%PYTHON_EXE%" -m ecd_epitope_scanner serve --host %HOST% --port %PORT%

echo   Waiting for the server to accept connections.
set /a tries=0
:wait
set /a tries+=1
timeout /t 1 >nul
netstat -ano | findstr /r /c:"LISTENING" | findstr /c:":%PORT% " >nul
if %errorlevel%==0 goto up
if %tries% GEQ 25 goto down
echo     ... %tries% s
goto wait

:up
echo   Server is up. Opening %URL%
start "" "%URL%"
echo.
echo   Keep the minimized "ECD Workbench Server" window running while
echo   you work. To stop it: run stop_ecd_epitope_workbench.bat
echo   or just close that window.
echo.
pause
exit /b 0

:already
echo   Port %PORT% is already serving - reusing the running instance.
start "" "%URL%"
echo.
pause
exit /b 0

:down
echo.
echo   [ERROR] The server did not start within 25 seconds.
echo   Open the minimized "ECD Workbench Server" window to read why.
echo   If dependencies are missing, recreate the venv once
echo   (run these from the folder that CONTAINS ecd_epitope_scanner):
echo     python -m venv ecd_epitope_scanner\.venv --system-site-packages
echo     ecd_epitope_scanner\.venv\Scripts\python.exe -m pip install -r ecd_epitope_scanner\requirements-dev.txt
echo.
pause
exit /b 1

:novenv
echo.
echo   [ERROR] Virtual environment not found:
echo          %PYTHON_EXE%
echo.
echo   Create it once (run these from the folder that CONTAINS
echo   ecd_epitope_scanner):
echo     python -m venv ecd_epitope_scanner\.venv --system-site-packages
echo     ecd_epitope_scanner\.venv\Scripts\python.exe -m pip install -r ecd_epitope_scanner\requirements-dev.txt
echo.
pause
exit /b 1
