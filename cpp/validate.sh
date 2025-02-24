#!/usr/bin/env bash

./cells $@
./reference $@
pushd validate

cat cells/* > cells.out
python validate.py cells.out reference.out
