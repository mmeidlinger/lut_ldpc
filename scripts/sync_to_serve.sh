#!/bin/sh

# Usage ./sync_to_serve [server] [folder]

# Specify your server here. More detailed ssh configureation
# should go into your ~/.ssh/config
SERVER_DEFAULT=gate
FOLDER_DEFAULT=params

BASEDIR='~/epfl-tuwien-ldpc/cpp'

SERVER=${1:-$SERVER_DEFAULT}
FOLDER=${2:-$FOLDER_DEFAULT}


rsync -avPzh --delete \
      --exclude '*.DS_Store' \
      --exclude '*.swp' \
      --exclude '*.swo' \
      $FOLDER/ $SERVER:$BASEDIR/$FOLDER/ 

   
