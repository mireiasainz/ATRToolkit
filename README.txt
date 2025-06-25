This repository contains the code used to perform the calibration procedure and determine the permittivity presented in the article "Permittivity determination of inorganic nanoparticles by ATR spectroscopy: comparison with reflectance techniques and limitations" by Sainz-Menchón et al.

The following files are included:

polarization_fraction.ipynb: Jupyter notebook for determining the polarization fraction induced by the spectrometer and ATR setup.

permittivity_determination.ipynb: Jupyter notebook for determining the permittivity of inorganic nanoparticles. The code is demonstrated using the delta-Al₂O₃ sample data presented in the paper.

module.py: Python module containing all the functions used in both Jupyter notebooks.

ATR_absorbance_H2O.dat: Experimental ATR absorbance data for distilled water.

ATR_absorbance_delta_Al2O3_NPs.dat: Experimental ATR absorbance data for delta-Al₂O₃ nanoparticles.

ATR_absorbance_gamma_Al2O3_NPs.dat: Experimental ATR absorbance data for gamma-Al₂O₃ nanoparticles.

Reference_permittivity_gammadelta_Al2O3_spectrochimica.csv: Reference permittivity values for delta- and gamma-Al₂O₃ from González de Arrieta et al. (DOI: 10.1016/j.saa.2023.122795).

requirements.txt: List of all the exact packages and their versions needed to reproduce the environment and run the project. You can easily install the packages using pip.
