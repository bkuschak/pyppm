# pyppm/software/geomagnetic

This directory contains software for using the PyPPM as a geomagnetic observatory, to measure the Earth's field over time.

In this role, PyPPM makes continuous measurements at a fixed rate, and logs raw data to a set of files.  Scripts then analyze the data and calculate Earth field strength and other metrics.  These scripts can generate optional plots for inclusion on a website.

### Solver

The Free Induction Decay (FID) frequency is solved using harmonic inversion by the Filter Diagonalization Method (FDM).  We use the [pharminv](https://github.com/aaren/pharminv) Python interface to the [harminv](https://github.com/stevengj/harminv) solver by Steven G. Johnson.  The advantage of this technique is a much better frequency resolution than what is possible by using an FFT.  This is possible because we assume a specific form for the signal, namely a small number of sinusoids in a relatively narrow band of interest.  

From the harminv README:

>In general, [FDM] is useful when you know on physical grounds that your system consists of a small number of decaying & oscillating modes in the bandwidth of interest, plus a limited amount of noise, and is not appropriate to analyze completely arbitrary waveforms.

FDM calculates a set of frequencies, decay constants, amplitudes, phases, and quality factors.  For each of these sets, it also calculates a figure of merit that indicates the quality of the fit.

### Data Acquisition

To support offline and possibly iterative analysis of collected data, the scripts are split into two components: data acquisition and analysis.

There is a data acquisition script that communicates with the PyPPM hardware, runs pulse programs, and collects raw data.  This data is saved to a file along with some metadata.  At a minimum, the metadata includes timestamps and actual sampling frequency (including any adjustments for ADC clock error).  It may also include data such as polarization time / current, sensor orientation, location, or other relevant information.  To improve the analysis, two data collections are done, one prior to the polarization, and one immediately afterwards.  This provides a background measurement and a measurement of the FID signal.  The data files will be standalone and contain everything that is useful in (re)intepreting them, perhaps at a much later time, or on another machine.

Logger scripts require that numpy be installed.

### Analysis

The analysis script(s) run the solver and compute the FID parameters and the Earth field strength.  They might also perform filtering, plotting, and audio file generation.  It is possible they might also influence the operation of the data acquisition, such as adjusting polarization time or ADC sampling.

Analysis scripts require a few prerequisites be installed: numpy, scipy, matplotlib, harminv, and pharminv

Currently there is analysis.py, that performs a basic analysis of the data using FDM and generates some plots, such as [ppm_plot.png.](./ppm_plot.png)  While the analysis seems to work during times of quiet background, when there is high background noise, it has trouble picking the right frequency that corresponds to the FID decay.  Much work is still to be done.

### Hardware Notes

This was initially tested with a [slightly modified version](../../designs/ppm-1.3/v1_redlined.pdf) of the PyPPM-1.3 hardware, and a coil very similar to the coil-2.0.  

I learned a few things during construction that I'll just briefly summarize.  The FID signal is very low amplitude and many things must be done correctly to get any detectable signal at all.

- Move the sensor coils as far away from cars, buildings, mains, (and the PyPPM hardware?) as practical, and closely balance the inductance of the coils.
- Use quality shielded twisted pair wire cable.  Shield wire should connect to power supply ground, not to the BNC barrel.  (I used a different 3-pin connector rather than the BNC.) I'm using 150 ft of 18AWG.
- Adding resonating capacitance (across INA_IN+/-) helps to increase FID amplitude and adds some amount of bandpass filtering.
- If your ADC is saturating, reduce background noise or reduce circuit gain.
- Use a clean power supply, and adjust for sufficiently high coil current (~2A seems good.)
- Make some measurements in the dead of night, when things tend to be quieter.

### Miscellaneous Notes

I've successfully run the capture script on MacOS Yosemite and Debian Jessie running on a Beaglebone Black.

On MacOS (Yosemite at least), the OS seems to truncate the serial number of the USB device and add an index. I had to create this link so PyPPM can find the device:

    sudo ln -s /dev/cu.usbmodemPyPP1 cu.usbmodemPyPPM

