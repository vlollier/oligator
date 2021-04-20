@echo off
if not exist venv virtualenv venv
call venv\Scripts\activate.bat
pip install -q -r .requirement.txt 
cd src
python main.py
