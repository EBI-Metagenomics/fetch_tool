#!/bin/bash

set -e

INSTALL_DIR="${1:-$PWD/aspera-cli}"

export VERSION=3.9.6
export FILE=ibm-aspera-cli-3.9.6.1467.159c5b1-linux-64-release

wget -nc https://download.asperasoft.com/download/sw/cli/$VERSION/$FILE.sh

sed -i "s@INSTALL_DIR=\~\/.aspera@INSTALL_DIR=$INSTALL_DIR@" $FILE.sh

chmod +x $FILE.sh

# Run the installation
./$FILE.sh
