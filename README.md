libWetHair
================
libWetHair is an open source project for the physical simulation of liquid and wet hairs. It has been compiled and tested on Mac OS X (with both Intel and Apple M1 chips), Ubuntu Linux, Windows, and licensed under the Mozilla Public License v. 2.0.

We would like to hear from you if you appreciate this work.

It is the original implementation of paper A Multi-Scale Model for Simulating Liquid-Hair Interactions (refer to our [project page](http://www.cs.columbia.edu/cg/liquidhair/) for more details). This code base contains the following parts:

 - A liquid simulator implementing the affine-particle-in-cell method.
 - A hair simulator implementing the elastic rods model.
 - A reduced-liquid simulator for the simulation of flow on hairs.
 - Cohesion effects between the hairs
 - Coupling between the hairs and liquid, including dragging, capturing and dripping effect.

Dependencies
--------------------
libWetHair depends on following libraries (some of them have been included in the code base, marked with an asterisk):

- Eigen* (http://eigen.tuxfamily.org/)
- RapidXML* (http://rapidxml.sourceforge.net)
- tclap* (http://tclap.sourceforge.net)
- AntTweakBar (http://anttweakbar.sourceforge.net/doc/)
- Intel TBB (https://www.threadingbuildingblocks.org)
- FreeGLUT (http://freeglut.sourceforge.net/)
- GLEW (http://glew.sourceforge.net)
- libPNG (https://libpng.sourceforge.io/)
- zlib (https://www.zlib.net/)

On Mac OS X or Linux-based systems, most of the dependencies are either included, or can be easily installed with Homebrew (https://brew.sh) or the APT package handling utility. For example, with Homebrew on Mac OS X, these external dependencies can be installed through
```
brew install anttweakbar tbb freeglut libpng zlib glew
```

On Windows you may need manually download and install some of them.

Compilation
-----------------
libWetHair has been tested with AppleClang (under Mac OS X), GCC 4.8+ (under Linux), and Microsoft Visual Studio (under Windows 10).

To compile libWetHair, you'll need CMake on Mac OS X or Linux, or CMake-GUI (https://cmake.org) on Windows.

On Mac OS X or Linux:
1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake .*
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

On Windows:
1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose your installed version of the Microsoft Visual Studio.
2. after configuration you may find several libraries not found, check the *Advanced* box and locate those missing libraries manually. Please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables to generate your Visual Studio solution.
4. open the solution and compile the code.

Run the Demo
--------------------
To run the demo of libWetHair, you may simply use the command line argument *-s [scene_file]* to specify the scene to be loaded. For example, you may type

./libWetHair -s assets/unit_tests/cylinder_2.xml

to run the simulation of the scene containing two hairs being cohesive to each other. 

All the parameters can be modified offline in the scene description XML files. Some can be changed online in the user interface provided by the demo program.

USAGE: 

   ./libWetHair  [-j <string>] [-o <integer>] [-d <boolean>] [-g <integer>] [-p <boolean>] [-s <string>] [--] [--version] [-h]

Where: 

   -j <string>,  --initfile <string>
     Binary file to load simulation state from for initialization

   -o <integer>,  --readableoutput <integer>
     readable file to save simulation state to; *after getting the output here, you may do surface reconstruction and rendering with the Houdini .HIP file provided* at http://www.cs.columbia.edu/cg/liquidhair/pseudo_dog.hipnc by renaming the HIP file into the name of the folder storing all output files (please check the python scripts in the HIP file for details).

   -d <boolean>,  --display <boolean>
     Run the simulation with display enabled if 1, without if 0

   -g <integer>,  --generate <integer>
     Generate PNG if 1, not if 0

   -p <boolean>,  --paused <boolean>
     Begin the simulation paused if 1, running if 0

   -s <string>,  --scene <string>
     Simulation to run; an xml scene file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

Contact
-----------
Please contact the author (yf2320@columbia.edu) for questions and bug report, or consulting for the usage of this code base.

BibTex Citation
----------------------
  @article{Fei:2017:liquidhair,
    title={A Multi-Scale Model for Simulating Liquid-Hair Interactions},
    author={Fei, Yun (Raymond) and Maia, Henrique Teles and Batty, Christopher and Zheng, Changxi and Grinspun, Eitan},
    journal={ACM Trans. Graph.},
    volume={36},
    number={4},
    year={2017},
    doi={10.1145/3072959.3073630},
  }
