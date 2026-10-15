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
rem Use a project-local temp directory when system TEMP is stale or unwritable.
set "TEMP=%PKG_DIR%workspace_data\runtime_tmp"
set "TMP=%TEMP%"
if not exist "%TEMP%" mkdir "%TEMP%"

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

"%PYTHON_EXE%" -c "import fastapi, uvicorn" >nul 2>&1
if errorlevel 1 goto dependencies

rem Model workers require WSL access from this Python process identity.
"%PYTHON_EXE%" -c "import subprocess, sys; sys.exit(subprocess.run(['wsl.exe','-d','Ubuntu-22.04','--exec','/bin/true'], timeout=30).returncode)"
if errorlevel 1 goto wslruntime

echo   Starting the server in the background...
if not exist "%PKG_DIR%workspace_data" mkdir "%PKG_DIR%workspace_data"
powershell -NoProfile -Command "$ErrorActionPreference='Stop'; Start-Process -FilePath $env:PYTHON_EXE -ArgumentList @('-m','ecd_epitope_scanner','serve','--host',$env:HOST,'--port',$env:PORT) -WorkingDirectory (Get-Location).Path -WindowStyle Hidden -RedirectStandardOutput (Join-Path $env:PKG_DIR 'workspace_data\server.stdout.log') -RedirectStandardError (Join-Path $env:PKG_DIR 'workspace_data\server.stderr.log')"
if errorlevel 1 goto down

echo   Waiting for the server to accept connections.
set /a tries=0
:wait
set /a tries+=1
timeout /t 1 >nul
netstat -ano | findstr /r /c:"LISTENING" | findstr /c:":%PORT% " >nul
if %errorlevel%==0 goto up
if %tries% GEQ 60 goto down
echo     ... %tries% s
goto wait

:up
echo   Server is up. Opening %URL%
start "" "%URL%"
echo.
echo   The server runs in the background.
echo   To stop it: run stop_ecd_epitope_workbench.bat
echo   Startup logs: %PKG_DIR%workspace_data\server.stderr.log
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
echo   [ERROR] The server did not start within 60 seconds.
echo   Read the startup log:
echo     %PKG_DIR%workspace_data\server.stderr.log
echo   If dependencies are missing, recreate the venv once
echo   (run these from the folder that CONTAINS ecd_epitope_scanner):
echo     python -m venv ecd_epitope_scanner\.venv
echo     ecd_epitope_scanner\.venv\Scripts\python.exe -m pip install -r ecd_epitope_scanner\requirements.lock
echo.
pause
exit /b 1

:wslruntime
echo.
echo   [ERROR] This Python process cannot start the WSL model workers.
echo   SPURS and DiscoTope require Ubuntu-22.04 to run successfully.
echo   Start this launcher from your normal Windows user session.
echo   See docs\full_runtime_restore_20261012.md for integrity-label repair.
echo.
pause
exit /b 1

:dependencies
echo.
echo   [ERROR] This Python environment cannot import FastAPI / Uvicorn.
echo   Restore the project dependencies before starting:
echo     uv pip install --python "%PYTHON_EXE%" -r "%PKG_DIR%requirements.lock"
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
echo     python -m venv ecd_epitope_scanner\.venv
echo     ecd_epitope_scanner\.venv\Scripts\python.exe -m pip install -r ecd_epitope_scanner\requirements.lock
echo.
pause
exit /b 1
