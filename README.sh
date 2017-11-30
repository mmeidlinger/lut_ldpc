#!/bin/sh
 
sed -e '/PARAMS_BER_INI/ {' -e 'r params/ber.ini.example' -e 'd' -e '}'  README.src.md > README.md

sed '/^#/ d'
