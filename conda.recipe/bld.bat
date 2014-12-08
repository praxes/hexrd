rmdir build /s /q

%PYTHON% setup.py install --old-and-unmanageable
if errorlevel 1 exit 1

copy scripts\* %SCRIPTS%\
if errorlevel 1 exit 1
