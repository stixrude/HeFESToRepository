# HeFESTo Linux Installation Guide

Written 2024-03-07 by Joseph Lewis-Merrill

Updated 2024-03-26 by David James

Special thanks to Matthew Bogumil and David James
***

## Introductions

Welcome to the unofficial installation guide for HeFESTo for a Linux type
environment. This guide was written especially for users of Debian based distros
like Ubuntu but should be generally applicable for most Linux distros. This
guide assumes you have already installed the basic compiler packages needed such
as gfortran and cmake. If not, you can generally install them using your package
manager, e.g.

```bash
 sudo apt update
 sudo apt install build-essential cmake libblas-dev liblapack-dev gfortran
```

HeFESTo is written in Fortran and must be compiled for each new environment. It
is recommended to back up your data before beginning this installation process
(as with all new code installations).

## Instruction

1. Download and install NLopt
   1. Download the latest version of NLopt at https://nlopt.readthedocs.io/en/latest/
   2. Follow the installation instructions for your distro at https://nlopt.readthedocs.io/en/latest/NLopt_Installation/. This generally looks like downloading the latest "nlopt-*.tar.gz" then extracting it to a directory using
		`tar -xf nlopt-*.tar.gz`
	and running
		`cmake . && make && sudo make install`
   3. Make note of the location to which "libnlopt.so" has been installed (this is usually something like "/usr/local/lib/libnlopt.so.0.11.1"). We'll call the path to this library [LIBNLOPT_PATH].

2. Download the HeFESTo code
   Navigate to the directory in which you'd like to install the HeFESTo source code
		`cd <installation path>`
   1. (Optional) Make a new directory to contain your HeFESTo environment
		`mkdir <HeFESTo directory name>`
   2. Download the latest version of the HeFESTo code from https://github.com/stixrude/HeFESToRepository using the following:
		`git clone https://github.com/stixrude/HeFESToRepository`
   3. Make note of this directory's path, which we'll call [HEFESTO_CODE_PATH]

3. Modify the HeFESTo code for Linux
   1. Open the file "makefile" under [HEFESTO_CODE_PATH] using your preferred text editor. In a gnome environment, this might look like
		`gedit <[HEFESTO_CODE_PATH]>/"makefile"`
   2. Replace the section
   ```
		LIB1 = -framework Accelerate
		LIB2 = -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

		$(COMMAND): $(MAIN) $(SUBS)
			$(LDR) $(LFLAGS) -o $(COMMAND) $(MAIN) $(SUBS) $(LIB1) $(LIB2)
   ```
   with
   ```
		LIB1 = -lGL -lblas -llapack <[LIBNLOPT_PATH]>

		$(COMMAND): $(MAIN) $(SUBS)
			$(LDR) $(LFLAGS) -o $(COMMAND) $(MAIN) $(SUBS) $(LIB1)
   ```
> Make sure to replace <[LIBNLOPT_PATH]> with the appropriate path from step 1 !!!
   3. Now, under [HEFESTO_CODE_PATH], run
		`make`
	This will compile HeFESTo in this directory.
   4. (Optional) Add HeFESTo to your $PATH environment so that you can run it from any directory by creating a symbolic link
		`sudo ln -s <ABSOLUTE PATH TO MAIN> /usr/local/bin/HeFESTo`
> Make sure to replace `ABSOLUTE PATH TO MAIN` with the full path to the newly created `main` file under `HEFESTO_CODE_PATH` (e.g. `/home/username/path/to/main`) !!!
	Note: you'll need to restart your terminal for this to take effect

4. Download your chosen HeFESTo parameters set
   1. Select a parameter set for HeFESTo under https://github.com/stixrude/. The latest at the time of this writing is https://github.com/stixrude/HeFESTo_Parameters_010123.git.
   2. Clone this parameter set to a directory with a short path length (e.g. close to your home directory like "/home/username/HeFESTo_env" using a command like
	`git clone https://github.com/stixrude/HeFESTo_Parameters_010123.git`
> Fortran generally does not handle symbolic links well so this parameter set needs to have a relatively short absolute path !!!
   3. Make note of this path, which we will call [PARAM_PATH]
   4. Under `PARAM_PATH`, change the name of the directory `phase` to `PHASE`, e.g.
   `mv phase PHASE`

5. Adjust your HeFESTo installation to reflect the location of `PARAM_PATH`
   1. Under `HEFESTO_CODE_PATH` find the directory named `BENCHMARK`
   2. Under `BENCHMARK` open the file `control` for editing, e.g.
   `gedit <[HEFESTO_CODE_PATH]>/BENCHMARK/control`
   3. **OPTIONAL** Make a backup copy of the control file before editing, e.g.
   `cp control .orig.control`
   4. Replace
   ```
		/Users/stixrude/qh/invert/par/inv010123
   ```
	with
    ```
		<[PARAM_PATH]>
    ```

6. Set up a clean workspace and test HeFESTo
   1. Copy `BENCHMARK` into a new directory outside of your HeFESTo code directory, e.g.
   `cp -r <[HEFESTO_CODE_PATH]>/BENCHMARK /home/username/HeFESTo_env/test_00`
   2. Navigate to this new directory. If you linked HeFESTo to your `$PATH` environment in section 3(iv), make sure you've restarted your terminal and simply run
		`HeFESTo`
	Otherwise, you'll need to invoke the full path to the HeFESTo executable, e.g.
		`<[HEFESTO_CODE_PATH]>/main`
   3. If no errors have occurred in the installation process, HeFESTo should run without issue. This may take a moment.
   4. Note: If you would like to save the output of running HeFESTo while still having it displayed, try:
		`HeFESTo | tee <log path>`
	If you would like to route the output only to a log file, try:
		`HeFESTo > <log path>`
