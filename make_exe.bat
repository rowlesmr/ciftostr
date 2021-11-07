
set filename=ciftostr
set dest=exe
set src=%~dp0

REM make the ciftostr exe file and clean up
pyinstaller  --onefile src/ciftostr/__main__.py --name %filename%
del %dest%\%filename%.exe
copy dist\%filename%.exe %dest%
del dist\%filename%.exe
pause