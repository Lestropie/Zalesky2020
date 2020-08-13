The data and script provided in this repository accompany results presented in the [article](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28266):

    Robert E. Smith, Fernando Calamante, Alan Connelly
    Notes on “A cautionary note on the use of SIFT in pathological connectomes”
    Magnetic Resonance in Medicine 2020;84(5):2303-2307

These are a response to the 'toy example' phantom presented in the [article](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.28037):

    Andrew Zalesky, Tabinda Sarwar and Kotagiri Ramamohanarao
    A cautionary note on the use of SIFT in pathological connectomes
    Magnetic Resonance in Medicine 2020;83(3):791-794

There are two experiments encapsulated within this repository:

-   Demonstration of the operation of the SIFT and SIFT2 methods on data for this phantom if those streamlines failing to traverse the pathological voxel in bundle 1 were not erroneously removed prior to the application of these methods. This demonstration is performed by executing file "`script`". This requires a working installation of the *MRtrix3* software (https://github.com/MRtrix3/mrtrix3).

-   Reproduction of the results of the analytic solution presented in the aforementioned manuscript, but with additional crucial details demonstrated. This was written for Octave - though it therefore may additionally work if performed in Matlab - and is executed via file "`process.m`". The additional files "`calc_cf.m`", "`calc_mu.m`" and "`calc_n1.m`" must all appear within the Octave path.

-----

This repository is the second of three demonstrations related to a series of articles in the Magnetic Resonance in Medicine journal. For demonstrations related to the prior and subsequent articles see:

http://github.com/Lestropie/Sarwar2019

http://github.com/Lestropie/Zalesky2020b
