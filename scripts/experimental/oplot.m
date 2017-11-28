#!/bin/octave-cli -qf

printf ("%s", program_name ());
arg_list = argv ();
for i = 1:nargin
  printf (" %s", arg_list{i});
endfor
printf ("\n");

addpath scripts;
analyze_results(arg_list);

% Wait before quit
pause()

