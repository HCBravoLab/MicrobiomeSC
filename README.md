# MicrobiomeSC

This code compute Presence/Absence, Differential Abundance, and Alpha Diversity testing 
over OTU data in a secure manner using garbled circuits.

In order to use this code:

1) Clone this repository

2) open two terminal windows

3) in the first terminal, run ./x_GenChiPGP.sh where x can be pc (Pre-compute), pf (Per Feature), or pfs (Per Feature Sparse)

4) in the second terminal, run ./x_EvaChiPGP.sh where x can be pc (Pre-compute), pf (Per Feature), or pfs (Per Feature Sparse)

5) This produced a verification of the correctness of the algorithm and a circuit size

6) In order to run the analysis for real, modify the Config.conf line #MODE to REAL.
