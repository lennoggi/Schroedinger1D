# DFT
This code computes the discrete Fourier transform (DFT) of a wave function and inverts it to get the original wave function back.

## Minimal requirements
- A C++ compiler supporting the `C++17` standard
- `cmake`
- The HDF5 library
- `python3` with `numpy`, `matplotlib`, and `h5py` to generate snapshots

## Usage
1. Tune the parameters -- macros living in `Parameters.hh`
2. Compile:
   ```
   ./build.sh
   ```
3. Run:
   ```
   ./install/bin/DFT_exe
   ```
