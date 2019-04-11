# FSPCVODE_MEX

Codes for building Matlab interface with Sundials's CVODES library using Mex.

Prerequisites:
- CMake.
- Sundials.
- Matlab.

To build (using command line):
- Create an empty subfolder within the FSPCVODE_MEX folder called "build". 
- Change directory (via, e.g, cd) to "build".
- Type "cmake ..".
- Type "make". The C source files are compiled into shared libraries.
- Open Matlab, run the script 'build_mex_files.m'.
- The Mex solver is ready to use. For more information, see the comments in "FspCVodeMex.c".
