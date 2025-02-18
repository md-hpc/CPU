#!/usr/bin/env bash

./cells $@ > cells.out
cp cells.out validate
python validate/validate-migration.py
