#! /usr/bin/env bash

mkdir tmp
cp .devcontainer/phyloModels-master.zip tmp
cd tmp
unzip phyloModels-master.zip
cd phyloModels-master
python3 -m pip install -e .
cd ../..
# rm -rf tmp
