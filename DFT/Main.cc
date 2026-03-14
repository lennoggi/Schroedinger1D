#include <cassert>
#include <array>
#include <complex>
#include <iostream>

#include <hdf5.h>

#include "include/Check_parameters.hh"
#include "include/Declare_functions.hh"

#include "Parameters.hh"

using namespace std;
using namespace std::complex_literals;



int main() {
    // Print some info
    #if (WF == GAUSSIAN)
    cout << "Gaussian wave function" << endl
         << "Center:    " << X0      << endl
         << "Sigma:     " << SIGMA   << endl
         << "Momentum:  " << P0      << endl
         << endl;
    #elif (WF == BOX)
    cout << "Box wave function"    << endl
         << "Center:    " << X0    << endl
         << "Width:     " << SIGMA << endl
         << "Momentum:  " << P0    << endl
         << endl;
    #else
    #error "Invalid wave function"
    #endif


    // Build the wave function
    constexpr double nhalf = static_cast<double>(N)/2.;
    constexpr double dx    = L/static_cast<double>(N); 

    array<double, N> jj, x;
    array<complex<double>, N> wf;

    for (auto j = decltype(N){0}; j < N; ++j) {
        jj[j] = static_cast<double>(j);
         x[j] = j*dx;  // NOTE: only used for output

        #if (WF == GAUSSIAN)
        wf[j] = gaussian_wf(x[j]);
        #elif (WF == BOX)
        wf[j] = box_wf(x[j]);
        #else
        #error "Invalid wave function"
        #endif
    }


    /* ---------------------------------------------------------
     * Build the discrete Fourier transform of the wave function
     * --------------------------------------------------------- */
    const complex<double> cexpfac  = 2.i*M_PI/static_cast<double>(N);
    constexpr double      Fwf_norm = dx/sqrt(2.*M_PI*HBAR);
    constexpr double      dp       = 2.*M_PI*HBAR/L; 

    array<double, N> ks, p;
    array<complex<double>, N> Fwf, Fwf_exact;

    for (auto k = decltype(N){0}; k < N; ++k) {
        ks[k]  = static_cast<double>(k) - nhalf;
         p[k]  = ks[k]*dp;  // NOTE: only used for output
        Fwf[k] = 0.;

        for (auto j = decltype(N){0}; j < N; ++j) {
            Fwf[k] += exp(-cexpfac*jj[j]*ks[k])*wf[j];
        }

        Fwf[k] *= Fwf_norm;

        #if (WF == GAUSSIAN)
        Fwf_exact[k] = gaussian_Fwf(p[k]);
        #elif (WF == BOX)
        Fwf_exact[k] = box_Fwf(p[k]);
        #else
        #error "Invalid wave function"
        #endif
    }


    /* -----------------------------------------------------------------
     * Build the inverse discrete Fourier transform of the wave function
     * ----------------------------------------------------------------- */
    constexpr double iFwf_norm = sqrt(2.*M_PI*HBAR)/L;
    array<complex<double>, N> iFwf;

    for (auto j = decltype(N){0}; j < N; ++j) {
        iFwf[j] = 0.;

        for (auto k = decltype(N){0}; k < N; ++k) {
            iFwf[j] += exp(cexpfac*jj[j]*ks[k])*Fwf[k];
        }

        iFwf[j] *= iFwf_norm;
    }



    /* ------------------
     * Write data to file
     * ------------------ */
    const auto file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    constexpr hsize_t fdims = N;
    const auto fspace_id    = H5Screate_simple(1, &fdims, nullptr);
    assert(fspace_id >= 0);


    // ***** Write the position and momentum *****
    const auto x_dset_id = H5Dcreate(file_id, "Position", H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto p_dset_id = H5Dcreate(file_id, "Momentum", H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(x_dset_id >= 0);
    assert(p_dset_id >= 0);

    assert(H5Dwrite(x_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, x.data()) >= 0);
    assert(H5Dwrite(p_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, p.data()) >= 0);

    assert(H5Dclose(x_dset_id) >= 0);
    assert(H5Dclose(p_dset_id) >= 0);


    // **** Write the wave function and its direct and inverse Fourier transforms *****
    const auto       wf_group_id = H5Gcreate(file_id, "Wave function",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto      Fwf_group_id = H5Gcreate(file_id, "DFT of the wave function",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto true_Fwf_group_id = H5Gcreate(file_id, "Analytical DFT of the wave function",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto     iFwf_group_id = H5Gcreate(file_id, "Inverse DFT of the wave function",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(      wf_group_id >= 0);
    assert(     Fwf_group_id >= 0);
    assert(true_Fwf_group_id >= 0);
    assert(    iFwf_group_id >= 0);

    const auto        re_wf_dset_id = H5Dcreate(      wf_group_id, "Real part",      H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto        im_wf_dset_id = H5Dcreate(      wf_group_id, "Imaginary part", H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto       re_Fwf_dset_id = H5Dcreate(     Fwf_group_id, "Real part",      H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto       im_Fwf_dset_id = H5Dcreate(     Fwf_group_id, "Imaginary part", H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto  re_true_Fwf_dset_id = H5Dcreate(true_Fwf_group_id, "Real part",      H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto  im_true_Fwf_dset_id = H5Dcreate(true_Fwf_group_id, "Imaginary part", H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto      re_iFwf_dset_id = H5Dcreate(    iFwf_group_id, "Real part",      H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto      im_iFwf_dset_id = H5Dcreate(    iFwf_group_id, "Imaginary part", H5T_NATIVE_DOUBLE, fspace_id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(       re_wf_dset_id >= 0);
    assert(       im_wf_dset_id >= 0);
    assert(      re_Fwf_dset_id >= 0);
    assert(      im_Fwf_dset_id >= 0);
    assert( re_true_Fwf_dset_id >= 0);
    assert( im_true_Fwf_dset_id >= 0);
    assert(     re_iFwf_dset_id >= 0);
    assert(     im_iFwf_dset_id >= 0);

    constexpr hsize_t memdims = 2*N;
    const auto memspace_id    = H5Screate_simple(1, &memdims, nullptr);
    assert(memspace_id >= 0);

              hsize_t start  = 0;  // Get real elements
    constexpr hsize_t stride = 2;  //   (Complex array layout: { [re, im], [re, im], ... , [re, im] })
    constexpr hsize_t count  = N;
    constexpr hsize_t block  = 1;

    assert(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, &start, &stride, &count, &block) >= 0);

    assert(H5Dwrite(       re_wf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,   wf.data()) >= 0);
    assert(H5Dwrite(      re_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
    assert(H5Dwrite( re_true_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
    assert(H5Dwrite(     re_iFwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT, iFwf.data()) >= 0);

    assert(H5Dclose(       re_wf_dset_id) >= 0);
    assert(H5Dclose( re_true_Fwf_dset_id) >= 0);
    assert(H5Dclose(      re_Fwf_dset_id) >= 0);
    assert(H5Dclose(     re_iFwf_dset_id) >= 0);

    start = 1;  // Imaginary elements
    assert(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, &start, &stride, &count, &block) >= 0);

    assert(H5Dwrite(       im_wf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,   wf.data()) >= 0);
    assert(H5Dwrite(      im_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
    assert(H5Dwrite( im_true_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
    assert(H5Dwrite(     im_iFwf_dset_id, H5T_NATIVE_DOUBLE, memspace_id, fspace_id, H5P_DEFAULT, iFwf.data()) >= 0);

    assert(H5Dclose(       im_wf_dset_id) >= 0);
    assert(H5Dclose(      im_Fwf_dset_id) >= 0);
    assert(H5Dclose( im_true_Fwf_dset_id) >= 0);
    assert(H5Dclose(     im_iFwf_dset_id) >= 0);


    assert(H5Sclose(  fspace_id) >= 0);
    assert(H5Sclose(memspace_id) >= 0);

    assert(H5Gclose(      wf_group_id) >= 0);
    assert(H5Gclose(     Fwf_group_id) >= 0);
    assert(H5Gclose(true_Fwf_group_id) >= 0);
    assert(H5Gclose(    iFwf_group_id) >= 0);

    assert(H5Fclose(file_id) >= 0);

    cout << "Data written to '" << FILENAME << "'" << endl;

    return 0;
}
