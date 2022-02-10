# MIRCO (MUSAM-IMCS Rough Contact cOde)

`MIRCO` is a Boundary element algorithm for simulating linear elastic frictionless normal contact between a rigid rough indentor and an elastic half-space. The research code is implemented throughout in object-oriented programming (C++) using modern software design and is parallelized with OpenMP for shared memory hardware architectures.

## Clone the Repository

You can clone the repository along with its submodules using:

```bash
cd <someBaseDir>
mkdir <sourceDir>
git clone --recursive https://github.com/imcs-compsim/MIRCO.git <sourceDir>
```
where `<someBaseDir>` is some directory in your machine and `<sourceDir>` will contain the `MIRCO` source code.

If you have already cloned the repository using:
```bash
git clone https://github.com/imcs-compsim/MIRCO.git <sourceDir>
```
you can pull the submodules using: 
```bash
cd <sourceDir>
git submodule update --init --recursive
```
To update the submodules, you can use the following command from your source directory:
```bash
git submodule update --recursive --remote
```

## Configure and Build code

To create an out-of-source build, first create a build directory using:
```bash
cd <someBaseDir>
mkdir <buildDir>
```
where `<buildDir>` is the build directory.

> Note: The exact location of `<buildDir>` is arbitrary, as long as it is _not_ a subdirectory of `<sourceDir>`.

Now, you have to navigate to the build directory and call the `do-configure` script in order to invoke `cmake`:
```bash
cd <buildDir>
<sourceDir>/do-configure
```
Build the `mirco` executable in the build directory using:
```bash
cd <buildDir>
make -j <numProc>
```
with `<numProc>` specifying the number of processes used for compilation.
The `mirco` executable will be created in the build directory.

## Run unit tests

You can run the unit tests from the the build directory using:
```bash
ctest
```

## Run code

To run the code with an input file, use the following command in your build directory:
```bash
./mirco <sourceDir>/Input/<someInputFile.json>
```
where `<someInputFile.json>` is any input file in the prescribed format.