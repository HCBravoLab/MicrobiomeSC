# MicrobiomeSC

This project is the first module of secureSeq, a suite of secure sequencing analysis tools that is in development.

The paper for this project is available at http://biorxiv.org/content/early/2015/09/02/025999. 

This code computes Presence/Absence, Differential Abundance, and Alpha Diversity testing over OTU data using garbled circuits.

Dependencies: FlexSC (https://github.com/wangxiao1254/FlexSC) included as jar in lib/; Java 8 (jdk-1.8, jre-1.8)

In order to use this code:

1) Clone this repository

2) open two terminal windows

3) in the first terminal, run ./x_GenChiPGP.sh where x can be pc (Pre-compute), naive (Per Feature Naive), or sparse (Per Feature Sparse)

4) in the second terminal, run ./x_EvaChiPGP.sh where x can be pc (Pre-compute), naive_(Per Feature Naive), or sparse (Per Feature Sparse)

5) This produces an output of the program for verification of the correctness of the algorithm and circuit size

6) In order to run the analysis, modify the Config.conf by commenting line 4 MODE VERIFY and uncommenting line 5 MODE REAL.  Before running an analysis in REAL mode, it is recommended to check the running time against the results in accuracy/
