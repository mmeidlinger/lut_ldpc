# ber_sim

## Parameters

### -p \<relative/path/to/parameterfile>

### -b \<path/to/basedirectory>
The base directory is where the program reads and writes files that are specified by relative paths. 
Usually, this is the directory from which the program is called. However, when using Software like HTCondor, this needs to be specified
explicitly since HTCondor has its own base directory. There are templates in the `scripts` subfolder that automatically give the right directory to condor

### -s \<integer>
Set the random seed of the simulation

### -c \<_result_name_extension>
Extends the name of the result files and folders by the the specified string (do not use whitespaces)


# de_sim

## Parameters

### -p \<relative/path/to/parameterfile>

# Other programs
All other programs need changes in their source code or take positional arguments in order to configure them.
