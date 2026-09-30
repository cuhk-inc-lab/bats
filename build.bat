@echo off
cd /d "%~dp0"
tools\tcc\tcc.exe -std=c11 -o bats_sim.exe field.c link.c source.c relay.c dest.c codec.c sim.c main.c
exit /b %ERRORLEVEL%
