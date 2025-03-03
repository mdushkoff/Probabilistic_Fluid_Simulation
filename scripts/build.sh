#!/bin/bash

# Output directories
build_dir="`pwd`/build"

# Check if the output directory exists
if [ ! -d `pwd`/build ]; then
    mkdir `pwd`/build
fi

# Make
pushd $build_dir
cmake ..
make clean && make
popd