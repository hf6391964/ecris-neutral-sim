# Project work: ECRIS neutral density simulation

## Building
Create a directory for building. Let's call it `build_root` for now.

Install the libraries `GMP`, `MPFR` and `Boost` on your system via the package manager.
Boost >= 1.48 is recommended.
Also install the tools `gcc`, `cmake` and `git`.

### CGAL
Download [CGAL 4.9](https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.9/CGAL-4.9.tar.xz) to `build_root`. Extract it and run `cmake-gui .` in the folder where it was extracted (`CGAL-4.9`). Uncheck `BUILD_SHARED_LIBS` and all the `WITH_*` flags except `WITH_CGAL_CORE`. This should get you a pretty minimal build. Run `make` to make it.

### GSL
Download and extract [GSL 2.2.1](ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.1.tar.gz) to `build_root`. Go to folder `gsl-2.2.1` and run the usual
```
./configure
make
```
to build it.

### The library
Clone this repo to `build_root`:
`git clone git://yousource.it.jyu.fi/erkka/erkka.git`
Enter the folder `erkka`. If everything went well before, now the project should
build by simply running `make`.

