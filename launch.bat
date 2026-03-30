@echo off
REM ============================================================================
REM Rangifer Assignment Tool — Windows Launcher
REM ============================================================================
REM Double-click this file to launch the app.
REM Requires R to be installed (https://cran.r-project.org/).
REM ============================================================================

cd /d "%~dp0"

REM Find Rscript.exe
where Rscript.exe >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set RSCRIPT=Rscript.exe
    goto :found
)

REM Check common R install locations
for /d %%d in ("C:\Program Files\R\R-*") do (
    if exist "%%d\bin\Rscript.exe" (
        set "RSCRIPT=%%d\bin\Rscript.exe"
        goto :found
    )
)

echo.
echo ERROR: R is not installed or not found in PATH.
echo Please install R from https://cran.r-project.org/
echo.
echo If R is installed, add its bin directory to your PATH:
echo   e.g., C:\Program Files\R\R-4.4.0\bin
echo.
pause
exit /b 1

:found
echo Using: %RSCRIPT%
"%RSCRIPT%" launch.R

echo.
echo App stopped.
pause
