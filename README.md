# MIRCO (MUSAM-IMCS Rough Contact cOde)

`MIRCO` is a Boundary element algorithm for simulating linear elastic frictionless normal contact between a rigid rough indentor and an elastic half-space.
The research code is implemented throughout in object-oriented programming (C++) and is parallelized with OpenMP for shared memory hardware architectures.

## Getting up and running

### Prerequesites

MIRCO requires

- C++ compiler
- [CMake](www.cmake.org) (version: >= 3.14)
- OpenMP
- [Trilinos](https://github.com/trilinos/Trilinos) (version: >= 12.8)

### Clone the repository

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

### Configure and build the code

To create an out-of-source build, first create a build directory using:

```bash
cd <someBaseDir>
mkdir <buildDir>
```

where `<buildDir>` is the build directory.

> Note: The exact location of `<buildDir>` is arbitrary, as long as it is _not_ the same directory than `<sourceDir>`.

Now, you have to navigate to the build directory and invoke `CMake`. We use CMake Presets for defining build configurations.

```bash
cd <buildDir>
cmake --preset=<name_of_your_preset> <sourceDir>
```

> **IMPORTANT** Make sure to set `Trilinos_DIR` to point to you Trilinos installation.

Build the `mirco` executable in the build directory using:

```bash
cd <buildDir>
ninja -j <numProc>
```

with `<numProc>` specifying the number of processes used for compilation.
The `mirco` executable will be created in the build directory.

### Run all tests

You can run the tests from the build directory using:

```bash
ctest
```

### Run the code

To run the code with an input file, use the following command in your build directory:

```bash
./mirco <sourceDir>/Input/<someInputFile.xml>
```

where `<someInputFile.xml>` is any input file in the prescribed format.

### Developing MIRCO

To develop MIRCO,
make sure to setup a proper development environment by running

```
./create-mirco-python-venv.sh
```

in the <sourceDir>.
This will install all required packages such as a code formatter and commit hooks.

## How to cite MIRCO?

If you are using this code, please cite the following paper:

- J. Bonari, M. R. Marulli, N. Hagmeyer, M. Mayr, A. Popp, M. Paggi: **A multi-scale FEM-BEM formulation for contact mechanics between rough surfaces**, _Computational Mechanics_, 65(3):731-749, 2020, [DOI: 10.1007/s00466-019-01791-3](https://doi.org/10.1007/s00466-019-01791-3)

   ```
   @article{Bonari2020a,
      author = {Bonari, Jacopo and Marulli, Maria Rosaria and Hagmeyer, Nora and Mayr, Matthias and Popp, Alexander and Paggi, Marco},
      doi = {10.1007/s00466-019-01791-3},
      issue = {3},
      journal = {Computational Mechanics},
      pages = {731--749},
      title = {{A multi-scale FEM-BEM formulation for contact mechanics between rough surfaces}},
      url = {https://doi.org/10.1007/s00466-019-01791-3},
      volume = {65},
      year = {2020}}
   ```
