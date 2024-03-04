:: Disable per-line exection output
@ECHO OFF
:: Get current working directory, i.e. directory this script is called from
SET WD=%cd%
:: Get this directory, i.e. directory this scripts is placed in, i.e. the SANS directory
SET PTH=%~dp0

:: Run SANS with check_N option
SET CALLVARS=%* -M %PTH%makefile
CALL %PTH%SANS.exe %CALLVARS%

:: If Ns okay, we're done. Re-compile and re-run if necessary, i.e. if exit code is 3
IF %ERRORLEVEL%==3 (
    :: Recompile SANS
    ECHO Re-compile with makefile_autoN in %PTH% ...
    CD %PTH%
	CALL make -f makefile_autoN
	
    :: Rerun SANS
    ECHO Re-run in working directory %WD% ... 
    CD %WD%
    CALL %PTH%SANS.exe %*
)
ECHO DONE
