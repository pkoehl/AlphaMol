# AlphaMolWrap

# **1. Background**
This repository is a fork of [pkoehl/AlphaMol](https://github.com/pkoehl/AlphaMol) and wraps the C++ source code provided there using [JuliaInterop/CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl). 

# **2. License**

This code is distributed under the LGPL-LICENSE. See LGPL-LICENSE.txt for details.

# **3. Compiling the code**

This section is a copy of the "Compiling the C++ code" section from [JuliaInterop/CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl):

The recommended way to compile the C++ code is to use CMake to discover `libcxxwrap-julia` and the Julia libraries.
A full example is in the [`testlib` directory of `libcxxwrap-julia`](https://github.com/JuliaInterop/libcxxwrap-julia/tree/master/testlib-builder/src/testlib).
The following sequence of commands can be used to build:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/libcxxwrap-julia-prefix /path/to/sourcedirectory
cmake --build . --config Release
```

The path for `CMAKE_PREFIX_PATH` can be obtained from Julia using:

```julia
julia> using CxxWrap
julia> CxxWrap.prefix_path()
```
**4. Disclaimer**

The program is provided "as is". Contact koehl@cs.ucdavis.edu
if you have any problems with the C++ source code. Contact spirandelli@uni-potsdam.de if you have any questions regarding the wrapping. 
