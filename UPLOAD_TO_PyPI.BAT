python setup.py sdist bdist_wheel
python -m twine upload --skip-existing dist/*

rem ======== if problems ===============
rem python -m twine check dist/*
