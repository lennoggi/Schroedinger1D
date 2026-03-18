# Evolution
This code evolves the Schroedinger equation.

## Minimal requirements
- A C++ compiler supporting the `C++17` standard
- The HDF5 library

## Usage
1. Tune the parameters in `Parameters.hh`. See some example parameter files under `Examples`.
2. Compile
   ```
   make -j3 options=OptionLists/<optionlist>
   ```
   Here `<optionlist>` is the list of compiler options for your machine.

3. Run
   ```
   ./Schroedinger1D_exe
   ```
4. To remove the executable and all object files, type
   ```
   make clean
   ```
5. To plot evolution snapshots and make movies, see `Utils/README.md` in the top directory of this repository.
