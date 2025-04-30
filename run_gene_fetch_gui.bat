@echo off
echo ==============================================
echo         Gene Fetch GUI Launcher
echo ==============================================
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in your PATH.
    echo Please install Python 3.9 or higher.
    pause
    exit /b 1
)

REM Run the application
echo Starting Gene Fetch GUI...
python gene_fetch_gui.py
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Error: Gene Fetch GUI failed to run properly.
    echo Please check the logs for more information.
    pause
    exit /b 1
)

echo Gene Fetch GUI completed successfully.
exit /b 0