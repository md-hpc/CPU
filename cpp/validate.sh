#!/usr/bin/env bash

name=$1
shift
./reference $@ > validate/reference.out
./$name $@ > validate/$name.out

cd validate
python validate.py reference.out $name.out
