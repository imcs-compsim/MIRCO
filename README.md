# BEM

It is a Boundary element algorithm for simulating linear elastic frictionless normal contact between rigid rough indentor and a half-space. The research code is implemented throughout in object-oriented programming (C++) using modern software design and is parallelized with OpenMP for shared memory hardware architectures.

## Clone the Repository

You can clone the repository along with its submodules using:

```bash
cd <someBaseDir>
mkdir <sourceDir>
git clone --recursive git@gitlab.com:compsim/codes/bem.git <sourceDir>
```
where `<someBaseDir>` is some directory in your machine and `<sourceDir>` will contain the BEM source code.

If you have already cloned the repository using:
```bash
git clone git@gitlab.com:compsim/codes/bem.git <sourceDir>
```
, you can pull the submodules using: 
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

Now, you have to navigate to the build directory and call `do-configure` using:
```bash
cd <buildDir>
<sourceDir>/do-configure
```
Build the `bem` executable in the build directory using:
```bash
<buildDir>/make
```
The `bem` executable will be created in the build directory.

## Run unit tests

You can run the unit tests from the the build directory using:
```bash
ctest
```

## Run code

To run the code with an input file, use the following command in your build directory:
```bash
./bem <sourceDir>/Input/<someInputFile.json>
```
where <someInputFile.json> is any input file in the prescribed format.