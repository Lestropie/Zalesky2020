The data and script provided in this repository are a response to the 'toy example' phantom presented in the article:

    Andrew Zalesky, Tabinda Sarwar and Kotagiri Ramamohanarao
    A cautionary note on the use of SIFT in pathological connectomes
    Magnetic Resonance in Medicine 2020;83(3):791-794

The results are presented in:

    Robert E. Smith, Fernando Calamante, Alan Connelly
    Notes on “A cautionary note on the use of SIFT in pathological connectomes”
    Magnetic Resonance in Medicine 2020;84(5):2303-2307

There are two experiments encapsulated within this repository:

-   Demonstration of the operation of the SIFT and SIFT2 methods on data for this phantom if those streamlines failing to traverse the pathological voxel in bundle 1 were not erroneously removed prior to the application of these methods. This demonstration is performed by executing file "`script`". This requires a working installation of the *MRtrix3* software (https://github.com/MRtrix3/mrtrix3).

-   Reproduction of the results of the analytic solution presented in the aforementioned manuscript, but with additional crucial details demonstrated. This was written for Octave - though it therefore may additionally work if performed in Matlab - and is executed via file "`process.m`". The additional files "`calc_cf.m`", "`calc_mu.m`" and "`calc_n1.m`" must all appear within the Octave path.
