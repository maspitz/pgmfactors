# Building with CMake

## Dependencies

The pgmfactors project requires Catch2, xtensor, and xtl.

## Build

You can build using cmake:

```sh
cmake -S . -B build
cmake --build build
```

## Install

Similarly:

```sh
cmake --install build
```

Here is the command for installing the release mode artifacts with a
multi-configuration generator, like the Visual Studio ones:

```sh
cmake --install build --config Release
```
