rmdir build /s

%PYTHON% setup.py install
if errorlevel 1 exit 1

copy scripts\* %SCRIPTS%\
if errorlevel 1 exit 1
