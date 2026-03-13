## Minimal requirements
- A C++ compiler supporting the `C++17` standard
- The HDF5 library

## Usage
1. Tune the parameters -- macros living in `Parameters.hh`
2. Compile
   ```
   make -j2 options=OptionLists/<optionlist>
   ```
   Here `<optionlist>` is the list of compiler options for your machine.

3. Edit the job submission script for your machine as needed, or create a new one if your machine is not listed there.
4. Run
   ```
   cd RunScripts
   ./<runscript>
   ```
5. To remove the executable and all object files, type
   ```
   make clean
   ```
