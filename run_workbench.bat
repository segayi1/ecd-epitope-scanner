@echo off
setlocal
set "APP_DIR=%~dp0"
set "PYTHON_EXE=%APP_DIR%.venv\Scripts\python.exe"
set "URL=http://127.0.0.1:8765"
if not exist "%PYTHON_EXE%" (
  echo Virtual environment not found: "%PYTHON_EXE%"
  echo Create it and install requirements-dev.txt first.
  pause
  exit /b 1
)
powershell -NoProfile -Command "try { Invoke-WebRequest -UseBasicParsing -TimeoutSec 1 '%URL%' ^| Out-Null; exit 0 } catch { exit 1 }"
if errorlevel 1 (
  start "ECD Epitope Workbench server" /min "%PYTHON_EXE%" -m ecd_epitope_scanner serve --host 127.0.0.1 --port 8765
  powershell -NoProfile -Command "Start-Sleep -Seconds 2"
)
start "" "%URL%"
endlocal
