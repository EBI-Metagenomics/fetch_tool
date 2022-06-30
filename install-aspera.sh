#!/bin/bash
set -e

INSTALL_DIR="${1:-$PWD/aspera-cli}"
VERSION=3.9.6
FILE=ibm-aspera-cli-3.9.6.1467.159c5b1-linux-64-release

if [[ "$OSTYPE" == "darwin"* ]]; then
    FILE=ibm-aspera-cli-3.9.6.1467.159c5b1-mac-10.11-64-release
fi

curl https://download.asperasoft.com/download/sw/cli/$VERSION/$FILE.sh -o $FILE.sh

if [[ "$OSTYPE" == "darwin"* ]]; then
    LANG=C sed -i.bk -e "s@INSTALL_DIR=\"\$HOME\/Applications\"@INSTALL_DIR=$INSTALL_DIR@" $FILE.sh
else
    sed -i.bk "s@INSTALL_DIR=\~\/.aspera@INSTALL_DIR=$INSTALL_DIR@" $FILE.sh
fi

chmod +x $FILE.sh

./$FILE.sh
