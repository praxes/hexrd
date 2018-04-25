REM !!! need to replace for proper versioning under setuptools???
REM git describe --tags --dirty > %SRC_DIR%/__conda_version__.txt
REM %PYTHON% %RECIPE_DIR%/format_version.py %SRC_DIR%/__conda_version__.txt

rmdir build /s /q

%PYTHON% setup.py install --old-and-unmanageable
if errorlevel 1 exit 1

copy scripts\* %SCRIPTS%\
if errorlevel 1 exit 1
