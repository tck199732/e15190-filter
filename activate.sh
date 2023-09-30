#!/bin/bash

export PROJECT_DIR=$(pwd)
# see ./env exists
if [ -d "./env" ]; then
    conda activate ./env
fi
