# FUTURES FAQs
## General questions
### What does FUTURES stand for:
FUTure Urban-Regional Environment Simulation

## Running FUTURES

## Windows-specific
### Running r.futures.potential gives ``FileNotFoundError: [Errno 2] No such file or directory: 'Rscript'`` or ``WindowsError: The system cannot find the file specified``
Module r.futures.potential is calling Rscript binary and it can't find it.
First, check R is installed, if not, install it.
Then try to reinstall GRASS, during installation it automatically looks for the binary.
If it still doesn't run afterwards, look where your R is installed and add the path to the binary to your PATH variable.

### Running r.futures.potential gives ``IOError: Permission denied: "potential.csv"``
You are attempting to write output potential.csv file where you don't have permission, such as in C:\,
change working directory (in GUI menu in Settings) or change directory in terminal to where you can create files.

### During running r.futures.potential it can't find required R packages, although I have them installed.
Try to install them again from within GRASS terminal: launch R interactively and install the packages.

## MacOS-specific
### Running r.futures.potential gives ``OSError: [Errno 2] No such file or directory`
Module r.futures.potential is calling Rscript binary and it can't find it.
Providing you have R installed, open a new terminal and run:

    which R

In the GRASS terminal now append this path to PATH variable:

    export PATH="$PATH:/usr/local/bin/"
    
Reopen GUI with ``g.gui``, or run r.futures.potential from terminal.
