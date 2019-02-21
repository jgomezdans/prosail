#!/bin/bash

set -e # Abort on error

$PYTHON setup.py install --single-version-externally-managed --record=record.txt

