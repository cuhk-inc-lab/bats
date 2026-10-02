@echo off
cd /d "%~dp0"
tools\tcc\tcc.exe -std=c11 -o bats_sim.exe field.c link.c source.c relay.c dest.c codec.c sim.c main.c psi_opt.c
if errorlevel 1 exit /b %ERRORLEVEL%
tools\tcc\tcc.exe -std=c11 -DPSI_OPT_MAIN -o psi_opt.exe psi_opt.c
exit /b %ERRORLEVEL%
