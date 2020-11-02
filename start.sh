#! /bin/bash
wd=`dirname $0`
target=$( readlink "$0" )

if [[ $target ]]
then
 wd=`dirname $target`
fi

if [[ ! -d $wd/.venv ]]
then
 virtualenv .venv
fi

source $wd/.venv/bin/activate

if [[ $? -eq 0 ]]
then
pip install -r $wd/.requirement.txt
cd $wd/src && python main.py
fi

