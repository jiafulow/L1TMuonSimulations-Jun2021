#!/bin/bash

PKG_1="emtf-tree"
URL_1="https://github.com/jiafulow/emtf-tree.git"
SRC_1="emtf_tree"

PKG_2="emtf-nnet"
URL_2="https://github.com/jiafulow/emtf-nnet.git"
SRC_2="emtf_nnet"

# Create tmp directory
mkdir -p tmp/third_party/
rm -rf tmp/third_party/${PKG_1}/
rm -rf tmp/third_party/${PKG_2}/

# Download via git
git clone ${URL_1} tmp/third_party/${PKG_1}/
git clone ${URL_2} tmp/third_party/${PKG_2}/

# Create directory
mkdir -p third_party/
touch third_party/__init__.py
rsync -avu tmp/third_party/${PKG_1}/${SRC_1}/ third_party/${SRC_1}/
rsync -avu tmp/third_party/${PKG_2}/${SRC_2}/ third_party/${SRC_2}/

# Apply some hacks
sed -i "s/\<${SRC_2}\>/third_party.${SRC_2}/g" third_party/${SRC_2}/__init__.py

# Remove tmp directory
rm -rf tmp/third_party/${PKG_1}/
rm -rf tmp/third_party/${PKG_2}/
