# Project work: ECRIS neutral density simulation

## Building
Create a directory for building. Let's call it `build_root` for now.

Install the libraries `GMP`, `MPFR` and `Boost` on your system via the package manager.
Boost >= 1.48 is recommended.
Also install the tools `gcc (>= 5)`, `cmake (>= 2.8.11)` and `git`. Note that a version of GCC
supporting the `C++14` standard is required. This should be 5, but hasn't
been tested.

### CGAL
The latest CGAL release 4.9 (as of writing this) has a bug relevant to this
code. Until a bugfix release (likely 4.9.1) is released, a checkout of the
`releases/CGAL-4.9-branch` branch should be used. The bug is fixed in that
branch.

Obtain [a checkout of the
branch](https://github.com/CGAL/cgal/archive/releases/CGAL-4.9-branch.zip) on
GitHub. Extract the zip file to `build_root`. There should now be a directory
named `cgal-releases-CGAL-4.9-branch`. Execute the following commands:
```sh
cd cgal-releases-CGAL-4.9-branch
# Create a directory for building
mkdir build
cd build

# Configure a minimal, static CGAL build
cmake -DBUILD_SHARED_LIBS=FALSE \
    -DWITH_CGAL_Core=FALSE \
    -DWITH_CGAL_ImageIO=FALSE \
    -DWITH_CGAL_Qt5=FALSE ..

# Build the library
make
```

### GSL
Download and extract [GSL 2.2.1](ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.1.tar.gz) to `build_root`. Go to folder `gsl-2.2.1` and run the usual
```sh
./configure
make
```
to build it.

### The library
Clone this repo to `build_root`:
```sh
git clone git://yousource.it.jyu.fi/erkka/samakapa-erkka.git
```
Then execute the following commands:
```sh
cd samakapa-erkka
# Run a custom script which runs CMake with the right flags
./run-cmake
make
```

### Apps / drivers
Simulations using the library are realized as frontends or "apps" which are
linked against the library. These reside in the `app` folder. To build an app,
ensure `libneutrals.a` is built and exists in the `samakapa-erkka` folder.

The simulations implemented for the research project are:

* `hiisi_electron_density`: an app that runs a single electron model simulation to obtain the electron distribution
* `hiisi_plasma_volume`: computes the volume in the HIISI plasma chamber which
  is bounded by field lines leading to the extraction hole
* `hiisi_neutral_density`: the main neutral density simulation, making use of
  the results of the two above simulations

As an example, building `hiisi_neutral_density` happens as follows:
```sh
cd app/hiisi_neutral_density
./run-cmake
make
```
An executable named `main` is generated.

The electron distribution files are already included in the folder
`app/hiisi_neutral_density/electron_data`. You can generate them yourself by
building and running `hiisi_electron_density` and copying the resulting CSV
files into `app/hiisi_neutral_density/electron_data`.

### Coding practices
If developing this code, please adhere to the following practices:

* Never ever even consider manually managing memory. You will eventually fail. Miserably. Elaboration follows in the next bullets.
* The use of `new` and `delete` is entirely prohibited. For arrays, use `std::vector`.
  For other heap allocated objects, use `std::make_unique` or
  `std::make_shared` depending on the desired object lifetime and resource
  sharing. If possible, the lifetime of objects shall be managed by scope.
* When indexing `std::vector`s, use the `at` method rather that the square brackets
  `[]`. E.g. use `vec.at(0) = 1`, not `vec[0] = 1`. The `at` method does bounds
  checking which saves you from memory corruption by throwing an exception when
  the vector index is out of bounds. If the profiler tells you
  that there's a bottle neck due to the usage of `at`, only remove them from the
  hot spots when sure that there's no out-of-bounds indexing bugs.
* Become familiar with the [C++ Core
  Guidelines](https://isocpp.github.io/CppCoreGuidelines/). Adopting these
  practices will save you a lot of headache.

