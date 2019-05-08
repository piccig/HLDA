
HLDA allows deriving low dimensional CVs out of a large set of descriptors in the form of a normalized linear combination to be easily used with the COMBINE command in PLUMED.
The repository allows reproducing two examples for a 2-class problem (the Diels-Alder reaction) and a 3-class problem (a symmetric SN2 reaction with three identical metastable states).

Both 2- and 3-class directories contain a *_class_HLDA.py script and 1/, 2/, and/or 3/ directories containing the COLVAR files for the unbiased fluctuations in reference metastable states. 
Run the script:

python *_class_HLDA.py

this will produce two output files with eigenvalues and first eigenvectors coefficients to be used as convergence check.
The screen output displays the eigenvalues and the related eigenvectors coefficients.
These can serve as HLDA CV.

Examples of PLUMED inputs (plumed.dat) can be found in the meta/ directories with relative CP2K inputs.

The python scripts are rather straightforward to extend to more than 3 classes and an internal PLUMED tool will be released asap.
