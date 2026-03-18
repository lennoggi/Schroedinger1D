#include <cassert>
#include <cmath>
#include <array>
#include <iostream>

#include <hdf5.h>

#include "include/Check_parameters.hh"
#include "include/Declare_functions.hh"

using namespace std;


int main() {
    // Print some info
    #if (WF == GAUSSIAN)
    cout << "Gaussian wave function"       << endl
         << "Initial center:    " << X0    << endl
         << "Initial sigma:     " << SIGMA << endl
         << "Initial momentum:  " << P0    << endl
         << endl;
    #elif (WF == BOX)
    cout << "Box wave function"            << endl
         << "Initial center:    " << X0    << endl
         << "Initial width:     " << SIGMA << endl
         << "Initial momentum:  " << P0    << endl
         << endl;
    #else
    #error "Invalid wave function"
    #endif


    // Build the initial wave function, potential, and helpers for the evolution
    constexpr double      nhalf   = static_cast<double>(NX)/2.;
    const complex<double> cexpfac = 2.0i*M_PI/static_cast<double>(NX);

    constexpr double dx = L/static_cast<double>(NX);
    constexpr double dp = 2.*M_PI*HBAR/L; 

    constexpr double dt_over_2hbar  = DT/(2.*HBAR);
    constexpr double dt_over_2hbarM = dt_over_2hbar/M;

    constexpr double  Fwf_norm = dx/sqrt(2.*M_PI*HBAR);
    constexpr double iFwf_norm = sqrt(2.*M_PI*HBAR)/L;

    array<double, NX> jj, x, V;
    array<double, NX> ks, p;

    array<complex<double>, NX> wf;
    array<complex<double>, NX>  cexp1, iFcexp1;
    array<complex<double>, NX> Fcexp1;


    for (auto j = decltype(NX){0}; j < NX; ++j) {
        jj[j] = static_cast<double>(j);
         x[j] = j*dx;  // NOTE: only used for output

        #if (WF == GAUSSIAN)
        wf[j] = gaussian_wf(x[j]);
        #elif (WF == BOX)
        wf[j] = box_wf(x[j]);
        #else
        #error "Invalid wave function"
        #endif

        #if (POT == FREE_PROPAGATION)
        V[j] = 0.;
        #elif (POT == BARRIER_WELL)
        V[j] = barrier_well(x[j]);
        #elif (POT == STEP)
        V[j] = step(x[j]);
        #elif (POT == HARMONIC)
        V[j] = harmonic(x[j]);
        #else
        #error "Invalid potential"
        #endif

          cexp1[j] = exp(-1.0i*V[j]*dt_over_2hbar);
        iFcexp1[j] = iFwf_norm*cexp1[j];
    }


    for (auto k = decltype(NX){0}; k < NX; ++k) {
           ks[k]  = static_cast<double>(k) - nhalf;
            p[k]  = ks[k]*dp;  // NOTE: only used for output
        Fcexp1[k] = Fwf_norm*exp(-1.0i*p[k]*p[k]*dt_over_2hbarM);
    }


    constexpr hsize_t tdims = NT/OUT_EVERY;  // NOTE: integer division
    array<double, tdims> time;

    /* ------------------------
     * Initialize the HDF5 file
     * ------------------------ */
    const auto file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    constexpr hsize_t fdims = NX;
    const auto fspace_id    = H5Screate_simple(1, &fdims, nullptr);
    assert(fspace_id >= 0);

    constexpr hsize_t memdims = 2*NX;
    const auto memspace_re_id = H5Screate_simple(1, &memdims, nullptr);
    const auto memspace_im_id = H5Screate_simple(1, &memdims, nullptr);
    assert(memspace_re_id >= 0); 
    assert(memspace_im_id >= 0); 

              hsize_t start;
    constexpr hsize_t stride = 2;   // (Complex array layout: { [re, im], [re, im], ... , [re, im] })
    constexpr hsize_t count  = NX;
    constexpr hsize_t block  = 1;

    start = 0;  // Get real elements
    assert(H5Sselect_hyperslab(memspace_re_id, H5S_SELECT_SET, &start, &stride, &count, &block) >= 0);
    start = 1;  // Get imaginary elements
    assert(H5Sselect_hyperslab(memspace_im_id, H5S_SELECT_SET, &start, &stride, &count, &block) >= 0);


    // ***** Output frequency ****
    constexpr hsize_t one   = 1;
    const auto space_one_id = H5Screate_simple(1, &one, nullptr);
    assert(space_one_id >= 0);

    const auto out_freq_dset_id = H5Dcreate(file_id, "Output frequency", H5T_NATIVE_INT, space_one_id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(out_freq_dset_id >= 0);
    assert(H5Dwrite(out_freq_dset_id, H5T_NATIVE_INT, space_one_id, space_one_id, H5P_DEFAULT, &OUT_EVERY) >= 0);
    assert(H5Dclose(out_freq_dset_id) >= 0);
    assert(H5Sclose(space_one_id) >= 0);


    // ***** Position, momentum, and potential *****
    const auto x_dset_id = H5Dcreate(file_id, "Position",  H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto p_dset_id = H5Dcreate(file_id, "Momentum",  H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto V_dset_id = H5Dcreate(file_id, "Potential", H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(x_dset_id >= 0);
    assert(p_dset_id >= 0);
    assert(V_dset_id >= 0);

    assert(H5Dwrite(x_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, x.data()) >= 0);
    assert(H5Dwrite(p_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, p.data()) >= 0);
    assert(H5Dwrite(V_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, V.data()) >= 0);

    assert(H5Dclose(x_dset_id) >= 0);
    assert(H5Dclose(p_dset_id) >= 0);
    assert(H5Dclose(V_dset_id) >= 0);


    // **** Wave function and its direct and inverse Fourier transforms *****
    const auto   wf_group_id = H5Gcreate(file_id, "Wave function",
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto  Fwf_group_id = H5Gcreate(file_id, "DFT of the wave function",
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto iFwf_group_id = H5Gcreate(file_id, "Inverse DFT of the wave function",
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(  wf_group_id >= 0);
    assert( Fwf_group_id >= 0);
    assert(iFwf_group_id >= 0);

    const auto   re_wf_group_id = H5Gcreate(  wf_group_id, "Real part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto   im_wf_group_id = H5Gcreate(  wf_group_id, "Imaginary part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto  re_Fwf_group_id = H5Gcreate( Fwf_group_id, "Real part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto  im_Fwf_group_id = H5Gcreate( Fwf_group_id, "Imaginary part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto re_iFwf_group_id = H5Gcreate(iFwf_group_id, "Real part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const auto im_iFwf_group_id = H5Gcreate(iFwf_group_id, "Imaginary part",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(  re_wf_group_id >= 0); 
    assert(  im_wf_group_id >= 0); 
    assert( re_Fwf_group_id >= 0); 
    assert( im_Fwf_group_id >= 0); 
    assert(re_iFwf_group_id >= 0); 
    assert(im_iFwf_group_id >= 0);



    /* ====================
     * Begin time evolution
     * ==================== */
    // NOTE: Fwf and iFwf are only used for output/diagnostics purposes
    array<complex<double>, NX>         wf1;
    array<complex<double>, NX>  Fwf,  Fwf1;
    array<complex<double>, NX> iFwf, iFwf1;

    for (auto n = decltype(NT){0}; n < NT; ++n) {
        /* Write the wave function before updating it so that there's no delay
         * between the wave function and the its inverse DFT
         * NOTE: this implies that the very last update of the wave function is
         *       never written to file. Not a big deal.                         */
        static_assert(sizeof(std::size_t) <= sizeof(long long),
                      "size_t too big for long long");
        const auto [quot, rem] = lldiv(static_cast<long long>(n), static_cast<long long>(OUT_EVERY));

        if (rem == 0) {
            time[quot] = n*DT;

            ostringstream iteration_ss;
            iteration_ss << "Iteration " << n;

            const auto re_wf_dset_id = H5Dcreate(re_wf_group_id, iteration_ss.str().c_str(),
                                                 H5T_NATIVE_DOUBLE, fspace_id,
                                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            const auto im_wf_dset_id = H5Dcreate(im_wf_group_id, iteration_ss.str().c_str(),
                                                 H5T_NATIVE_DOUBLE, fspace_id,
                                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert(re_wf_dset_id >= 0);
            assert(im_wf_dset_id >= 0);

            assert(H5Dwrite(re_wf_dset_id, H5T_NATIVE_DOUBLE, memspace_re_id, fspace_id, H5P_DEFAULT, wf.data()) >= 0);
            assert(H5Dwrite(im_wf_dset_id, H5T_NATIVE_DOUBLE, memspace_im_id, fspace_id, H5P_DEFAULT, wf.data()) >= 0);

            assert(H5Dclose(re_wf_dset_id) >= 0);
            assert(H5Dclose(im_wf_dset_id) >= 0);
        }


        // Time evolution
        for (auto j = decltype(NX){0}; j < NX; ++j) {
            wf1[j] = cexp1[j]*wf[j];
        }

        Fwf.fill(0.);
        Fwf1.fill(0.);

        for (auto k = decltype(NX){0}; k < NX; ++k) {
            for (auto j = decltype(NX){0}; j < NX; ++j) {
                const auto cexp_jk = exp(-cexpfac*jj[j]*ks[k]);
                Fwf[k]  += cexp_jk*wf[j];
                Fwf1[k] += cexp_jk*wf1[j];
            }

            Fwf[k]  *= Fwf_norm;
            Fwf1[k] *= Fcexp1[k];
        }

        iFwf.fill(0.);
        iFwf1.fill(0.);

        for (auto j = decltype(NX){0}; j < NX; ++j) {
            for (auto k = decltype(NX){0}; k < NX; ++k) {
                const auto cexp_jk = exp(cexpfac*jj[j]*ks[k]);
                iFwf[j]  += cexp_jk*Fwf[k];
                iFwf1[j] += cexp_jk*Fwf1[k];
            }

            iFwf[j] *= iFwf_norm;
              wf[j]  = iFcexp1[j]*iFwf1[j];
        }


        // Write the DFT and its inverse to file
        if (rem == 0) {
            ostringstream iteration_ss;
            iteration_ss << "Iteration " << n;

            const auto  re_Fwf_dset_id = H5Dcreate( re_Fwf_group_id, iteration_ss.str().c_str(),
                                                   H5T_NATIVE_DOUBLE, fspace_id,
                                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            const auto  im_Fwf_dset_id = H5Dcreate( im_Fwf_group_id, iteration_ss.str().c_str(),
                                                   H5T_NATIVE_DOUBLE, fspace_id,
                                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            const auto re_iFwf_dset_id = H5Dcreate(re_iFwf_group_id, iteration_ss.str().c_str(),
                                                   H5T_NATIVE_DOUBLE, fspace_id,
                                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            const auto im_iFwf_dset_id = H5Dcreate(im_iFwf_group_id, iteration_ss.str().c_str(),
                                                   H5T_NATIVE_DOUBLE, fspace_id,
                                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert( re_Fwf_dset_id >= 0);
            assert( im_Fwf_dset_id >= 0);
            assert(re_iFwf_dset_id >= 0);
            assert(im_iFwf_dset_id >= 0);

            assert(H5Dwrite( re_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_re_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
            assert(H5Dwrite( im_Fwf_dset_id, H5T_NATIVE_DOUBLE, memspace_im_id, fspace_id, H5P_DEFAULT,  Fwf.data()) >= 0);
            assert(H5Dwrite(re_iFwf_dset_id, H5T_NATIVE_DOUBLE, memspace_re_id, fspace_id, H5P_DEFAULT, iFwf.data()) >= 0);
            assert(H5Dwrite(im_iFwf_dset_id, H5T_NATIVE_DOUBLE, memspace_im_id, fspace_id, H5P_DEFAULT, iFwf.data()) >= 0);

            assert(H5Dclose( re_Fwf_dset_id) >= 0);
            assert(H5Dclose( im_Fwf_dset_id) >= 0);
            assert(H5Dclose(re_iFwf_dset_id) >= 0);
            assert(H5Dclose(im_iFwf_dset_id) >= 0);

            if constexpr (VERBOSE) {
                cout << "Iteration " << n << ": data written to file" << endl;
            }
        }

        if constexpr (VERBOSE) {
            cout << "Iteration " << n << ": done" << endl;
        }
    }


    // Write time to file 
    const auto tspace_id = H5Screate_simple(1, &tdims, nullptr);
    assert(tspace_id >= 0);

    const auto t_dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, tspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(t_dset_id >= 0);

    assert(H5Dwrite(t_dset_id, H5T_NATIVE_DOUBLE, tspace_id, tspace_id, H5P_DEFAULT, time.data()) >= 0);

    assert(H5Dclose(t_dset_id) >= 0);
    assert(H5Sclose(tspace_id) >= 0);


    // Close up
    assert(H5Sclose(     fspace_id) >= 0);
    assert(H5Sclose(memspace_re_id) >= 0);
    assert(H5Sclose(memspace_im_id) >= 0);

    assert(H5Gclose(  re_wf_group_id) >= 0);
    assert(H5Gclose( re_Fwf_group_id) >= 0);
    assert(H5Gclose(re_iFwf_group_id) >= 0);

    assert(H5Gclose(  im_wf_group_id) >= 0);
    assert(H5Gclose( im_Fwf_group_id) >= 0);
    assert(H5Gclose(im_iFwf_group_id) >= 0);

    assert(H5Gclose(  wf_group_id) >= 0);
    assert(H5Gclose( Fwf_group_id) >= 0);
    assert(H5Gclose(iFwf_group_id) >= 0);

    assert(H5Fclose(file_id) >= 0);

    cout << "Data written to '" << FILENAME << "'" << endl;

    return 0;
}
