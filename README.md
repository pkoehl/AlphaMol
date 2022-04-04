# AlphaMolWrap

# **1. Background**
This repository is a fork of [pkoehl/AlphaMol](https://github.com/pkoehl/AlphaMol) and wraps the C++ source code provided there using [JuliaInterop/CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl). There are some minor changes to the C++ source code but all functionality from the C++ side is exposed to the Julia side. 

# **2. License**

This code is distributed under the LGPL-LICENSE. See LGPL-LICENSE.txt for details.

# **3. Compiling the code**

This following is a copy of the "Compiling the C++ code" section from [JuliaInterop/CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl):

>The recommended way to compile the C++ code is to use CMake to discover `libcxxwrap-julia` and the Julia libraries.
>A full example is in the [`testlib` directory of `libcxxwrap-julia`](https://github.com/JuliaInterop/libcxxwrap-julia/tree/master/testlib-builder>/src/testlib).
>The following sequence of commands can be used to build:
>
>```bash
>mkdir build && cd build
>cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/libcxxwrap-julia-prefix /path/to/sourcedirectory
>cmake --build . --config Release
>```
>
>The path for `CMAKE_PREFIX_PATH` can be obtained from Julia using:
>
>```julia
>julia> using CxxWrap
>julia> CxxWrap.prefix_path()
>```

After building the Makefile using Cmake you can call ```make all``` to build the shared library. When the build process is completed you should see a shared library file called libAlphaMol.so
