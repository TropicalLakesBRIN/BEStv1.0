@echo off
:: Change directory to where this batch file lives
cd /d "%~dp0"

:: Execute Rscript directly on run_app.R
"portable\bin\Rscript.exe" "run_app.R"

pause