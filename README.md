## Physiological synchrony in EEG, electrodermal activity and heart rate reflects shared selective auditory attention

This repository contains the MATLAB scripts reproducing the results as presented in (Stuldreher et al., 2020). Please read carefully the information below. The presented analysis is dependent on lots of externally produced scripts, that are referred to here. For questions, please contact [Ivo Stuldreher](ivo.stuldreher@tno.nl).

### Dependencies
The code in this repository is dependent on the following add-ons:
- The physiological data corresponding to this paper, that is publicly available through https://osf.io/8kh36/. The data files should be put in the directory 'data'. Note: you should not put this data in the pre-assigned folders 'eeg_processed' or 'aut_processed'.
- EEGLAB toolbox (https://sccn.ucsd.edu/eeglab/index.php)
- MARA plugin for EEGLAB
- Ledalab toolbox (http://www.ledalab.de)
- corrca.m from the Correlated Component Analysis package (https://www.parralab.org/corrca/)
- sigstar.m (https://nl.mathworks.com/matlabcentral/fileexchange/39696-raacampbell-sigstar)
- Pan-Tompkin algorithm (for instance: https://nl.mathworks.com/matlabcentral/fileexchange/45840-complete-pan-tompkins-implementation-ecg-qrs-detector)
- fdr_bh.m (https://nl.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh)
- swtest.m (https://nl.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests)

### REFERENCES
Stuldreher, I.V., Thammasan, N., Van Erp, J.B.F., Brouwer, A.M. (2020). *Journal of Neural Engineering* [doi:10.1088/1741-2552/aba87d](https://iopscience.iop.org/article/10.1088/1741-2552/aba87d)
