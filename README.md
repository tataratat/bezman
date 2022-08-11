# BEZMAN
Light weight library in C++ for functional composition for Bézier splines.

This small prototyping library is created to test out algorithms to create and modify small microstructures and analytically derive their control points. The objective is the analytical calculation of gradients of microstructures formed through functional composition between splines.


## Installation
The Project uses cmake, to facilitate integration in other projects. *Project is tested on gcc `10.3.xx` and has also been tested on clang `11.1.xx` and `12.0.xx`.*

To build the project clone bezman and go into its directory. Here you create a build directory and step into it using the following command:
```
mkdir build
cd build
```
Now all that is left to do for installation is running the cmake command
```
cmake ..
```
... and you are ready to go. 

By default, cmake will create a new directory `install` in your bezman folder. If you want to specify another install directory, use the ` -DCMAKE_INSTALL_PREFIX=<link-to-your-install-directory>` option instead. Further, the default build type is `debug` (which sets debug flags and checks the code at run-time using assertions). If you aim for performace, choose Release instead, which sets compiler optimization flags and ignores these checks (set the type using `-DCMAKE_BUILD_TYPE=Release`. Google-test is also activated by default, however it is fetched from the internet, no need to have it preinstalled (can be turned off using `-DGOOGLETEST=OFF`).

Now install the library by running the command
```
make install
```

That's it. There are no external dependencies, however the `c++17` standard must be supported and new compilers are recommended. To run the unit tests, execute `ctest` (e.g. with the verbose option).

| :grey_question: Changing the standard to `c++20` :grey_question: |
|:---------------------------|
| `std::vector` types are used frequently to store all kinds of information (e.g. control points, Spline groups, etc.). The `c++17` standard library prohibits its use at compile time. The new standard allows for these operations|

## Building an example
There are also some simple examples provided, that show the usage of the library in a bit more detail. To build one of them go inside the example directory. It is recommended to not build directly inside the directory itself, but to provide an additional build folder. To do so run:
```
mkdir build
cd build
```
Run CMake and build the executable
```
cmake ..
make
```
The examples are also equiped with bezman's own logging functionality, which can be enabled using the cmake option `-DLOGGING=ON`, please check out the `CMakeLists.txt` file, if you want to use logging in your personal project.

## Building the documentation
If you want to build the documentation, you can do so using [doxygen](https://www.doxygen.nl/index.html), by running these commands:
```
# Go to the documentation directory
cd doc
# Build Documentation
doxygen Doxyfile
```
This will create a new folder named `doxydocs` in the current directory, in which you will find an `index.html` file. Open it with a browser of your choice.

## Python integration
Most of the functionality provided by this library is integrated into the software suite [gustaf](https://github.com/tataratat/gustaf). `Gustaf` is a user-friendly python package that offers all kinds of post- and preprocessing tools for vizualisation and much, much more.
