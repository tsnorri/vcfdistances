# vcfdistances

Calculate SMD and Hamming and Jaccard distances between each pair of samples in a set of variant files.

## Build/Runtime Requirements

On Linux the following libraries are required:

- [libdispatch](http://nickhutchinson.me/libdispatch/) (provided as a Git submodule)
- [libpthread_workqueue](https://github.com/mheily/libpwq) (provided as a Git submodule)
- [libkqueue](https://github.com/mheily/libkqueue)

## Build Requirements

- Reasonably new compilers for C and C++, e.g. GCC 7. C++17 support is required.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [Ragel State Machine Compiler](http://www.colm.net/open-source/ragel/) (tested with version 6.7)
- [CMake](http://cmake.org)
- [Boost](http://www.boost.org)

## Building

### Short version

1. `git clone https://github.com/tsnorri/vcfdistances.git`
2. `cd vcfdistances`
3. `git submodule update --init --recursive`
4. Edit local.mk
5. `make -j4`

### Long version

1. Clone the repository with `git clone https://github.com/tsnorri/vcfdistances.git`.
2. Change the working directory with `cd vcfdistances`.
3. Run `git submodule update --init --recursive`. This clones the missing submodules and updates their working tree.
4. Edit `local.mk` in the repository root to override build variables. Useful variables include `CC`, `CXX`, `RAGEL` and `GENGETOPT` for C and C++ compilers, gengetopt and Ragel respectively. `BOOST_INCLUDE` is used as preprocessor flags when Boost is required. `BOOST_LIBS` and `LIBDISPATCH_LIBS` are passed to the linker. See `common.mk` for additional variables.
5. Run make with a suitable numer of parallel jobs, e.g. `make -j4`

Useful make targets include:

<dl>
	<dt>all</dt>
	<dd>Build everything</dd>
	<dt>clean</dt>
	<dd>Remove build products except for dependencies (in the <code>lib</code> folder).</dd>
	<dt>clean-all</dt>
	<dd>Remove all build products.</dd>
</dl>


## Running

The tool takes one or more Variant Call Format files as its input. It reads the variants and makes pairwise comparisons between each sample. Remaining files are processed similarly. The requested distances are output as triangular matrices.

Please see `src/vcfdistances --help` for command line options.
