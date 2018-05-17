# *SeismoSoil*

#### Software package for 1D site response analysis and engineering seismology tools

-----



*SeismoSoil* provides the following site response analysis routines:

- Linear visco-elastic analysis in both time and frequency domain
- Equivalent linear analysis
  - Original algorithm (Seed & Idriss, 1970)
  - Frequency-dependent algorithm ([Asimaki & Kausel, 2002](http://www.sciencedirect.com/science/article/pii/S0267726102001203))
- Nonlinear analysis (total stress) in time domain
  - Masing rule with modified Kondner-Zelasko (MKZ) constitutive model
  - Non-masing rule with MKZ model
  - Non-masing rule with Hybrid Hyperbolic (HH) constitutive model ([Shi & Asimaki, 2017](http://resolver.caltech.edu/CaltechAUTHORS:20170404-150827374))

*SeismoSoil* also provides the following tools for manipulating earthquake time series:

- Band pass and high/low pass filtering
- Baseline correction of acceleration time histories
- Fourier spectra and elastic response spectra of acceleration time histories
- Goodness-of-fit calculation of two time series
- Parsing of PEER and SMC formatted ground motion records



### System requirements

To run from source codes within MATLAB:

- MATLAB R2013b+
- Required MATLAB toolboxes:
  - Curve Fitting Toolbox
  - Global Optimization Toolbox
  - Parallel Computing Toolbox
  - Signal Processing Toolbox
  - System Identification Toolbox
  - Wavelet Toolbox

To run compiled executables:

- Windows 7/8/8.1/10, or macOS (Mavericks and newer)
- Memory > 2 GB



### Installation

To run from source:

- No installation necessary
- Just download the source codes and put them in a folder

To run from compiled executable:

- Download executable files and extract to hard drive
- Before running *SeismoSoil* for the first time, **make sure to set up MATLAB Runtime**:
  - Go to http://www.mathworks.com/products/compiler/mcr/
  - Download MATLAB Runtime (64-bit, R2017b version) for your operating system
  - Install MATLAB Runtime to your computer



### How to run *SeismoSoil*

To run from source:

- Run `SeismoSoil.m` from the folder that contains it
- You can also add the folder that contains `SeismoSoil.m` into your MATLAB search path (by running `pathtool` from the command window)

To run from compiled executable:

- On Windows
  - Double-click `SeismoSoil.exe`
- On macOS
  - Place `SeismoSoil` in the Applications folder
  - Open Terminal and execute: `Applications/SeismoSoil.app/Contents/MacOS/applauncher`



### User manual

The user manual, which contains a tutorial and technical details, is [here](https://github.com/jsh9/SeismoSoil-manual/blob/master/SeismoSoil_manual.pdf).



### Copyright

(c) 2014-2018, [GeoQuake Research Group](http://asimaki.caltech.edu/), California Institute of Technology



### Contributors

Domniki Asimaki, Jian Shi, Wei Li, Peyman Ayoubi



### Terms and conditions

By downloading and using this software, you agree to the terms and conditions as listed in the `LICENSE.txt` file.


