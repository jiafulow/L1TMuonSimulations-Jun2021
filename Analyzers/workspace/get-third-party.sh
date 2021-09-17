#!/bin/bash

URL="https://github.com/jiafulow/emtf-tree.git"
PKG="emtf-tree"
SRC="emtf_tree"
TMPDIR="tmp/${PKG}/"

# Create tmp directory
mkdir -p tmp/
rm -rf ${TMPDIR}

# Download via git
git clone ${URL} ${TMPDIR}

# Create directory
mkdir -p third_party/
touch third_party/__init__.py
mv ${TMPDIR}/${SRC} third_party/${SRC}

# Remove tmp directory
rm -rf ${TMPDIR}
rmdir tmp/
