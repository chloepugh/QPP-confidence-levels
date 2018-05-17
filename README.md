# QPP-confidence-levels
find_cl.py and find_cl_rebin.py calculate confidence levels for power spectra with a power-law dependence. They can be used 
to test the significance of quasi-periodic pulsation (QPP) signatures in flare time series data.

The time series data used to calculate the power spectrum should be evenly time-spaced. Fourier power spectra or periodograms calculated at the Fourier frequencies may be used with these codes.

A power-law or broken power-law model should be fitted to the power spectrum. It is usually easier to work in log space for this. Uncertainties on the fitted power-law model should be estimated for each frequency index, and written to an ascii file named 'pow_fit_err.dat' so that it can be read in by either find_cl.py or find_cl_rebin.py. The uncertainties could be taken from the least-squares fitting, or a random sampling method could be used.

find_cl_rebin.py should be used with a rebinned power spectrum, where the powers in every n frequency bins are summed together (usually n=2 or n=3 are most suitable).

find_cl_example.pro and find_cl_rebin_example.pro show examples of how find_cl.py and find_cl_rebin.py (respectively) can be
used with IDL.

Used in the following publications:
- Pugh, C. E., Broomhall, A.-M., Nakariakov, V. M. 2017, Astronomy & Astrophysics, 602, A47
- Pugh, C. E., Nakariakov, V. M., Broomhall, A.-M., Bogomolov, A. V., Myagkova, I. N. 2017, Astronomy & Astrophysics, 608, A101
- Kolotkov, D. Y., Pugh, C. E., Broomhall, A.-M., Nakariakov, V. M. 2018, The Astrophysical Journal Letters, 858, L3
