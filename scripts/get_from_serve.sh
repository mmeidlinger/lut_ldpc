#!/bin/sh

# Usage ./sync_to_serve [server] [folder]

# Specify your server here. More detailed ssh configureation
# should go into your ~/.ssh/config
SERVER_DEFAULT=gate
FOLDER_DEFAULT=results

BASEDIR='~/epfl-tuwien-ldpc/cpp'

SERVER=${1:-$SERVER_DEFAULT}
FOLDER=${2:-$FOLDER_DEFAULT}


rsync -avPzh \
      $SERVER:$BASEDIR/$FOLDER/ $FOLDER/  

   
