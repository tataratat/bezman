# bezierManipulation
Light weight library in C++ to test out various composition methods

This small prototyping library is created to test out algorithms to create and modify small microstructures and analytically derive their control points. The longterm objective is to calculate analytical gradients of microstructures formed through functional composition between two splines


## Installation
The Project uses cmake, to facilitate integration in other projects.

Go to your destinated installation directory and create a build and install directory (can also be inside the bezierManipulation folder). Step into the build directory.
```
mkdir build install
cd build
```

Run CMake specifying an installation directory
```
cmake -DCMAKE_INSTALL_PREFIX=../install ..
```
Default build type is debug but if you aim for performace, choose Release instead. GTEST is inactive by default.

Now install the bezierManipulation by running the command
```
make install
```

That's it. There are no external dependencies, however the `c++17` standard must be supported and new compilers are recommended (I am running gcc 10.3.0)

If you want to run the unit test rerun cmake with 
```
cmake -DGOOGLETEST=ON ..
```
and build the tests using
```
make test
```
run the test using the `ctest` command (e.g. with the verbose option).

## Building an example
There is currently only one example provided (to be extended). To build it go inside the example directory. It is recommended to not build directly inside the directory itself, but to provide an additional build folder. To do so run:
```
mkdir build
cd build
```
Run CMake and build the executable
```
cmake ..
make
```
