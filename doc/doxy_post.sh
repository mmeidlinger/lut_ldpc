#!/bin/sh

# Copy files into html directory a

rsync -a codes/ doc/html/
rsync -a ensembles/ doc/html/
rsync -a Makefile doc/html
rsync -aR prog doc/html
rsync -aR  scripts doc/html
rsync -aR  params/{ber,de}.ini.example doc/html
rsync -a trees/tree_file_example.ini doc/html
