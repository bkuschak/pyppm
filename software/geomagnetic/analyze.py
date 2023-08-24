# Analyze the data files collected by the proton precession magnetometer.
# Using the FDM method of harmonic inversion, calculate the FID precession
# frequency, and field strength.
# 
# 1/13/18 B. Kuschak <bkuschak@yahoo.com> 
#

import pyppm            
import harminv              # https://github.com/aaren/pharminv
from pykalman import KalmanFilter # https://github.com/pykalman/pykalman

import sys
import time
import pytz
import glob, os
import pickle
import getopt
import datetime as dt
import scipy.io.wavfile
import numpy as np
import numpy.lib.recfunctions as rfn

# Allow MPL to work with no display attached.
import matplotlib
matplotlib.use('Agg') # set the backend before importing pyplot

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from matplotlib.font_manager import FontProperties
from scipy import signal

location = 'Bend, Oregon'

# Import iirpeak from a recent scipy, since these are not present in scipy-0.14
def iirpeak(w0, Q, fs=2.0):
    """
    Design second-order IIR peak (resonant) digital filter.
    A peak filter is a band-pass filter with a narrow bandwidth
    (high quality factor). It rejects components outside a narrow
    frequency band.
    Parameters
    ----------
    w0 : float
        Frequency to be retained in a signal. If `fs` is specified, this is in
        the same units as `fs`. By default, it is a normalized scalar that must
        satisfy  ``0 < w0 < 1``, with ``w0 = 1`` corresponding to half of the
        sampling frequency.
    Q : float
        Quality factor. Dimensionless parameter that characterizes
        peak filter -3 dB bandwidth ``bw`` relative to its center
        frequency, ``Q = w0/bw``.
    fs : float, optional
        The sampling frequency of the digital system.
        .. versionadded:: 1.2.0
    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (``b``) and denominator (``a``) polynomials
        of the IIR filter.
    See Also
    --------
    iirnotch
    Notes
    -----
    .. versionadded:: 0.19.0
    References
    ----------
    .. [1] Sophocles J. Orfanidis, "Introduction To Signal Processing",
           Prentice-Hall, 1996
    Examples
    --------
    Design and plot filter to remove the frequencies other than the 300 Hz
    component from a signal sampled at 1000 Hz, using a quality factor Q = 30
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    >>> fs = 1000.0  # Sample frequency (Hz)
    >>> f0 = 300.0  # Frequency to be retained (Hz)
    >>> Q = 30.0  # Quality factor
    >>> # Design peak filter
    >>> b, a = signal.iirpeak(f0, Q, fs)
    >>> # Frequency response
    >>> freq, h = signal.freqz(b, a, fs=fs)
    >>> # Plot
    >>> fig, ax = plt.subplots(2, 1, figsize=(8, 6))
    >>> ax[0].plot(freq, 20*np.log10(np.maximum(abs(h), 1e-5)), color='blue')
    >>> ax[0].set_title("Frequency Response")
    >>> ax[0].set_ylabel("Amplitude (dB)", color='blue')
    >>> ax[0].set_xlim([0, 500])
    >>> ax[0].set_ylim([-50, 10])
    >>> ax[0].grid(True)
    >>> ax[1].plot(freq, np.unwrap(np.angle(h))*180/np.pi, color='green')
    >>> ax[1].set_ylabel("Angle (degrees)", color='green')
    >>> ax[1].set_xlabel("Frequency (Hz)")
    >>> ax[1].set_xlim([0, 500])
    >>> ax[1].set_yticks([-90, -60, -30, 0, 30, 60, 90])
    >>> ax[1].set_ylim([-90, 90])
    >>> ax[1].grid(True)
    >>> plt.show()
    """

    return _design_notch_peak_filter(w0, Q, "peak", fs)


def _design_notch_peak_filter(w0, Q, ftype, fs=2.0):
    """
    Design notch or peak digital filter.
    Parameters
    ----------
    w0 : float
        Normalized frequency to remove from a signal. If `fs` is specified,
        this is in the same units as `fs`. By default, it is a normalized
        scalar that must satisfy  ``0 < w0 < 1``, with ``w0 = 1``
        corresponding to half of the sampling frequency.
    Q : float
        Quality factor. Dimensionless parameter that characterizes
        notch filter -3 dB bandwidth ``bw`` relative to its center
        frequency, ``Q = w0/bw``.
    ftype : str
        The type of IIR filter to design:
            - notch filter : ``notch``
            - peak filter  : ``peak``
    fs : float, optional
        The sampling frequency of the digital system.
        .. versionadded:: 1.2.0:
    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (``b``) and denominator (``a``) polynomials
        of the IIR filter.
    """

    # Guarantee that the inputs are floats
    w0 = float(w0)
    Q = float(Q)
    w0 = 2*w0/fs

    # Checks if w0 is within the range
    if w0 > 1.0 or w0 < 0.0:
        raise ValueError("w0 should be such that 0 < w0 < 1")

    # Get bandwidth
    bw = w0/Q

    # Normalize inputs
    bw = bw*np.pi
    w0 = w0*np.pi

    # Compute -3dB attenuation
    gb = 1/np.sqrt(2)

    if ftype == "notch":
        # Compute beta: formula 11.3.4 (p.575) from reference [1]
        beta = (np.sqrt(1.0-gb**2.0)/gb)*np.tan(bw/2.0)
    elif ftype == "peak":
        # Compute beta: formula 11.3.19 (p.579) from reference [1]
        beta = (gb/np.sqrt(1.0-gb**2.0))*np.tan(bw/2.0)
    else:
        raise ValueError("Unknown ftype.")

    # Compute gain: formula 11.3.6 (p.575) from reference [1]
    gain = 1.0/(1.0+beta)

    # Compute numerator b and denominator a
    # formulas 11.3.7 (p.575) and 11.3.21 (p.579)
    # from reference [1]
    if ftype == "notch":
        b = gain*np.array([1.0, -2.0*np.cos(w0), 1.0])
    else:
        b = (1.0-gain)*np.array([1.0, 0.0, -1.0])
    a = np.array([1.0, -2.0*gain*np.cos(w0), (2.0*gain-1.0)])

    return b, a





class ppm_analysis:
    gp = 42.576         # actual gyromagnetic ratio of proton in earth field

    def __init__(self, log_fname, verbose):
        self.t0 = []
        self.a0 = []
        self.t1 = []
        self.a1 = []
        self.starttime = 0
        self.fs = 20000.0
        self.set_expected_range()
        self.best_candidate = None
        self.log_fname = log_fname
        self.verbose = verbose
        #self.interactive = True
        self.interactive = False
        self.intermag_recent = []
        self.intermag_name = None
        self.geometrics_recent = []
        self.geometrics_name = None

    def set_local_offset(self, offset):
        self.local_offset = offset
    
    def debug(self, level, msg):
        if self.debug_level >= level:
            print msg
    
    # Read raw ADC samples from the data file
    # (t0,a0) is the background. (t1,a1) is the real measurement.
    # If freq_error=None, then error value is taken from the data file 
    def load_data_file(self, fname, freq_error=None):
        if self.verbose > 0:
            print 'Loading file:', fname

        file_data = pickle.load(open(fname, 'rb'))
        self.a0 = file_data['a0']               # background
        self.a1 = file_data['a1']               # FID measurement
        self.fs = file_data['fs']               # nominal sample rate
        if freq_error != None:
            self.freq_error = freq_error            # override sample rate clock error
        else:
            self.freq_error = file_data['freq_error']   # error in the sample rate, from datafile
        self.start_time = file_data['start_time']       # Unix epoch of first (background) measurement
        self.polarize = file_data['polarize_time']      # polarization time
        self.comments = file_data['comments']
        self.data_length = len(self.a0)
        self.freq_error = float(self.freq_error)

        if self.verbose > 1:
            print   'Data File: Fs = %u KHz (%+.1f ppm), %u samples. Polarization time = %.3f sec.\n' \
                'Comments: %s' \
                % (self.fs, self.freq_error*1e6, self.data_length, self.polarize, self.comments)

        # Create time arrays
        ts = 1.0 / (self.fs * (1.0 + self.freq_error))
        self.t0 = tuple(np.arange(0, len(self.a0) * ts, ts))    # pyppm requires tuples
        self.t1 = tuple(np.arange(0, len(self.a1) * ts, ts))

    # Given a directory containing files downloaded from Intermagnet, load them for plotting
    def load_intermagnet_data(self, path, name=None):
        self.intermag_recent = []
        self.intermag_name = name
        
        for f in glob.glob(path):
            if self.verbose > 0:
                print 'scanning intermagnet data file', f
            # Filter out comments and invalid data
            for line in open(f, 'r'):
                if '99999.00' not in line:
                    # UTC time
                    # DATE       TIME         DOY     H      D      Z      F
                    data = [i for i in line.split()]
                    ts = data[0] + ' ' + data[1]
                    try:
                        ts = mdates.date2num(dt.datetime.strptime(ts, "%Y-%m-%d %H:%M:%S.%f"))
                        field = float(data[6])
                        self.intermag_recent.append((ts, field))
                    except ValueError:
                        if self.verbose > 1:
                            print 'rejecting line', line
                        continue

    # Given a Geometrics data file, load it for plotting
    def load_geometrics_data(self, path, name=None):
        self.geometrics_recent = []
        self.geometrics_name = name
        
        for f in glob.glob(path):
            if self.verbose > 0:
                print 'scanning geometrics data file', f
            # Filter out comments and invalid data
            for line in open(f, 'r'):
                # Local time
                # $ FIELD(nT),unknown   DATE       TIME 
                    data = [i for i in line.split()]
                    try:
                        ts = data[2] + ' ' + data[3]
                        ts = dt.datetime.utcfromtimestamp(time.mktime(dt.datetime.strptime(ts, "%m/%d/%y %H:%M:%S.%f").timetuple()))
                        ts = mdates.date2num(ts)
                        #ts = mdates.date2num(dt.datetime.strptime(ts, "%m/%d/%y %H:%M:%S.%f"))
                        field_data = [i for i in data[1].split(',')]
                        field = float(field_data[0])
                        self.geometrics_recent.append((ts, field))
                    except:
                        if self.verbose > 1:
                            print 'rejecting line', line
                        continue

    # Larmor frequency conversion
    def field_to_freq(self, field_nt):
        return field_nt / 1000.0 * self.gp

    def freq_to_field(self, freq_hz):
        return freq_hz * 1000.0 / self.gp

    # Set some bounds so we can filter appropriately.  Values in nanotesla.
    def set_expected_range(self, nt_low=47000, nt_high=49000):
        self.expected_field_low = nt_low
        self.expected_field_high = nt_high
        self.expected_field = (nt_high + nt_low) / 2.0
        self.expected_freq_low = self.field_to_freq(self.expected_field_low)
        self.expected_freq_high = self.field_to_freq(self.expected_field_high)
        self.expected_freq = self.field_to_freq(self.expected_field)

    # Analyze the data by running the FDM
    def analyze(self, plot_fname=None, wavfilename=None):
        if self.verbose > 0:
            print '\nAnalyzing...'
        # narrowband filter the data around the expected peak.  
        fnaught = self.field_to_freq(self.expected_field)
        Q = 40.0
        w0 = fnaught/(self.fs/2)    # normalized

        #b, a = signal.iirpeak(w0, Q)   # FIXME - not available in scipy-0.14
        #self.a1f = np.array(self.a1)
        #b, a = iirpeak(w0, Q)       # use local - not available in scipy-0.14
        b, a = signal.butter(4, [self.expected_freq_low / (self.fs/2), self.expected_freq_high / (self.fs/2)], btype='bandpass', output='ba')
        self.a1f = signal.filtfilt(b, a, self.a1)   

        # Filter Frequency response
        f, h = signal.freqz(b, a, 2048)
        self.filt_freq = f * self.fs / (2 * np.pi)
        #self.filt_mag = 20*np.log10(np.maximum(abs(h), 1e-5))
        self.filt_mag = np.maximum(abs(h), 1e-5)

        # FIXME - we should truncate the beginning and end that are subject to filter startup, so they don't throw off
        # the harmonic inversion decay calculation

        if wavfilename != None:
            self.save_wavfile(wavfilename)

        # Use FDM harmonic inversion to calculate a set of frequencies in this band.
        self.signals = harminv.invert(self.a1f, fmin=self.expected_freq_low, fmax=self.expected_freq_high, dt=1.0/self.fs, nf=30)

        # Next compute FFTs of background, measurement, and filtered measurement data.
        # We will use this to help select the correct signal from the list of FDM candidates.
        #self.a1f = tuple(self.a1f)
        (self.fft_f0, self.fft_a0) = pyppm.fft(self.t0, self.a0)
        (self.fft_f1, self.fft_a1) = pyppm.fft(self.t1, self.a1)
        (self.fft_f1f, self.fft_a1f) = pyppm.fft(self.t1, tuple(self.a1f))

        # For each of the candidate FDM signals, compute a narrowband SNR which is defined as the
        # difference between the measurement and the background, at that frequency.
        self.compute_narrowband_snrs()

        # Select best candidate signal, or None if we couldn't find one.
        self.fdm = self.select_candidate(self.signals)
        if self.fdm == None:
            return None, None, 0.0, 0.0, 0.0

        self.fdm_time_constant = 1.0 / self.fdm.decay
        self.field = self.freq_to_field(self.fdm.frequency)

        # Compute wideband SNR across the expected frequency band.
        self.compute_wideband_snr()

        # GellerLabs calls the FDM error the Figure of Merit. Lower is better.
        # Shouldn't be converted to nT, as harminv says 'error is not really error bars on frequency'
        self.fom = self.freq_to_field(self.fdm.error * self.fdm.frequency)
        #if self.verbose > 1:
            #print "FDM wideband SNR (dB)", self.wb_snr, "FOM", self.fdm.error

        return self.field, self.fdm, self.fom, self.fdm.nb_snr, self.wb_snr

#    # Given a real frequency, interpolate the FFT magnitude from the nearest points.
#    def fft_magnitude(self, freqs, amplitudes, frequency):
#        if frequency > freqs[-1]:
#            print 'Frequency out of range. Clamping to max.'
#            frequency = freqs[-1]
#        if frequency < freqs[0]:
#            print 'Frequency out of range. Clamping to min.'
#            frequency = freqs[0]
#
#        # Find bounding indexes for this frequency.
#        f_higher = np.argmax(freqs > frequency)
#        f_lower = np.min(f_higher - 1, 0)
#
#        # Interpolate the magnitude.
#        diff_amplitude = amplitudes[f_higher] - amplitudes[f_lower]
#        fractional_f = (frequency - freqs[f_lower]) / (freqs[f_higher] - freqs[f_lower])
#        interpolated_amplitude = amplitudes[f_lower] + fractional_f * diff_amplitude
#        print 'Found frequency %.3f between indexes %d and %d (%.3f and %.3f). ' \
#            'amplitudes %.3f and %.3f. Fractional freq: %.3f. Interpolated amplitude %.3f' \
#            % (frequency, f_lower, f_higher, freqs[f_lower], freqs[f_higher], amplitudes[f_lower], amplitudes[f_higher], fractional_f, interpolated_amplitude)
#        return interpolated_amplitude

    # Given length 2 lists for x and y, and xp value somewhere between x[0] and x[1], return an interpolated value y at xp.
    def interpolate(self, x, y, xp):
            diff_y = y[1] - y[0]
            diff_x = x[1] - x[0]
            fractional_x = (xp - x[0]) / (x[1] - x[0])
            return y[0] + fractional_x * diff_y

    # For each of the FDM signals detected, compute a narrowband SNR, defined as the ratio between the measurement and the background.
    def compute_narrowband_snrs(self):
        # TODO Sort the FDM by frequency if not already done.

        # Add new field to structured array, to hold the narrowband SNR. Default to zero SNR.
        # For some reason asrecarray=True does not work, so explicitly convert the structured array to recarray,
        # so we can access by member names.
        self.signals = rfn.append_fields(self.signals, names='nb_snr', data=np.zeros(self.signals.shape), usemask=False)
        self.signals = np.rec.array(self.signals) 

        # Iterate over the FDM signals by increasing frequency.
        # This code assumes frequencies fft_f0 and fft_f1 are the same, which they are since both background and measurement are the same length and sample rate.
        f_higher_idx = 0
        for (i, s) in enumerate(self.signals):

            # Find the fft frequency bucket just higher than this FDM frequency.
            while f_higher_idx < len(self.fft_f0) and self.fft_f0[f_higher_idx] <= s.frequency:
                f_higher_idx += 1

            if f_higher_idx >= len(self.fft_f0):
                # We have reached the end of the FFT. Cannot go any further.
                break

            # Interpolate the magnitude between FFT buckets.
            f_lower_idx = np.min(f_higher_idx - 1, 0)
            background_ampl = self.interpolate(self.fft_f0[f_lower_idx:f_higher_idx+1], self.fft_a0[f_lower_idx:f_higher_idx+1], s.frequency)
            measurement_ampl = self.interpolate(self.fft_f1[f_lower_idx:f_higher_idx+1], self.fft_a1[f_lower_idx:f_higher_idx+1], s.frequency)

            # Compute a narrowband SNR at this frequency.
            nb_snr = 20 * np.log10(measurement_ampl / background_ampl)
            self.signals[i].nb_snr = nb_snr
            if self.verbose:
                print 'narrowband SNR: %.3f: %.1f dB' % (s.frequency, s.nb_snr)

        print self.signals
        return
        #if self.verbose:
            #for i,s in enumerate(self.signals):


    # See http://www.gellerlabs.com/PMAG%20Docs.htm
    # Like GellerLabs, compute an SNR which is 20 log (FID amplitude / RMS_Sum(other inband amplitudes))
    # We call this the wideband SNR since it is computed over the full bandwidth of the expected frequency range.
    def compute_wideband_snr(self):
        rms_sum = 0
        for s in self.signals:              
            if  s.frequency >= self.expected_freq_low and s.frequency <= self.expected_freq_high and \
                s.frequency != self.fdm.frequency:
                rms_sum += (s.amplitude ** 2)
        if rms_sum > 0:
            self.wb_snr = 20.0 * np.log10(self.fdm.amplitude / np.sqrt(rms_sum))
        else:
            self.wb_snr = 0
        return


    # given the signals resulting from harmonic inversion, select the most likely FID signal
    # FIXME - use the background capture data to improve selection
    def select_candidate(self, signals):

        # mode frequencies, absolute amplitudes, phase shift, decay rates (reciprocal of time constant), Q factor, error
        # TODO - print 1/decay instead of decay. This is tau.
        if self.verbose > 0:
            print 'freq, amplitude, phase, decay, q, error'
            print signals

        # Exclude any candidates that have values out of the expected range.
        # The proton signal should have a small tau2 (0.3 to 2.0), high Q (over 5000?), low error, largish amplitude
        idx = 0
        while idx < len(signals):
            s = signals[idx]
            if s.frequency < self.expected_freq_low or s.frequency > self.expected_freq_high:
                if self.verbose > 1: print "discarding bad freq", s
                signals = np.delete(signals, idx)
            elif s.decay < 1/0.8 or s.decay > 1/0.2:    # tau2 is 1/decay
                if self.verbose > 1: print "discarding bad decay", s
                signals = np.delete(signals, idx)
            elif abs(s.Q) < 2000:
                if self.verbose > 1: print "discarding bad Q", s
                signals = np.delete(signals, idx)
            elif s.error > 5e-7:
                if self.verbose > 1: print "discarding bad error", s
                signals = np.delete(signals, idx)
            elif s.amplitude < 0.03 or s.amplitude > 0.50:  # RMS amplitude
                if self.verbose > 1: print "discarding bad amplitude", s
                signals = np.delete(signals, idx)
            else:
                idx += 1


#        # select the one that has the highest delta FFT amplitude between the background and measurement.
#        # TODO - fft_magnitude should be computed on the entire list of frequencies rahter than each one individually 
#        # to make the computation more efficient.
#        #delta_fft = []
#        narrowband_snrs = []
#        for s in signals:
#            background = self.fft_magnitude(self.fft_f0, self.fft_a0, s.frequency)
#            measurement = self.fft_magnitude(self.fft_f1, self.fft_a1, s.frequency)
#            delta = measurement - background
#            narrowband_snr = measurement / background
#            print 'Narrowband SNR at %f Hz = %.3f' % (s.frequency, narrowband_snr)
#            #print 'delta_fft at %f Hz = %.3f' % (s.frequency, measurement - background)
#            #delta_fft.append(delta)
#            narrowband_snrs.append(narrowband_snr)
#        #best_delta_fft = signals[np.argmax(delta_fft)]
#        best_narrowband_snr = signals[np.argmax(narrowband_snrs)]
#        #print 'best delta_fft: %f Hz = %.3f' % (best_delta_fft.frequency, delta_fft[np.argmax(delta_fft)])
#        print 'best narrowband_snr: %f Hz = %.3f' % (best_narrowband_snr.frequency, narrowband_snrs[np.argmax(narrowband_snrs)])

        if self.verbose > 0:
            print 'candidates: freq, amplitude, phase, decay, q, error'
            print signals

        # We may not find any that match our criteria...
        if len(signals) == 0:
            if self.verbose > 0:
                print 'None'
            return None

        # Now, just select the one with the highest narrowband SNR.
        idx = np.argmax(signals.nb_snr)
        return signals[idx]


    # Output the filtered proton signal as a wave file
    def save_wavfile(self, wavfilename):
        scipy.io.wavfile.write(wavfilename, self.fs, self.a1f)

    # Create some nice plots
    def plot_results(self):
        # 60 Hz harmonics 
        # FIXME - line rate isn't always exactly 60 Hz...
        powerline = np.asarray(range(60,100*60,60)) 

        # Compute FFTs
        self.a1f = tuple(self.a1f)
        (f0, A0) = pyppm.fft(self.t0, self.a0)
        (f1, A1) = pyppm.fft(self.t1, self.a1)
        (f1f, A1f) = pyppm.fft(self.t1, self.a1f)

        # Plot time domain
        plt.plot(np.array(self.t0), np.array(self.a0), 'b', label='background')
        plt.plot(np.array(self.t1), np.array(self.a1), 'r', label='measurement')
        plt.title('Proton Precession Magnetometer')
        plt.ylabel('ADC Voltage (V)')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show()

        # Plot filtered time domain
        plt.plot(np.array(self.t1), np.array(self.a1f), 'g', label='post-filtered')
        plt.title('Proton Precession Magnetometer')
        plt.ylabel('ADC Voltage (V)')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show()

        # Plot frequency domain
        l, = plt.plot(np.array(f0), np.array(A0), 'b', label='background')
        l, = plt.plot(np.array(f1), np.array(A1), 'r', label='measurement')
        l, = plt.plot(np.array(f1f), np.array(A1f), 'g', label='post-filtered')
        plt.axvline(x=self.expected_freq_low, color='red', linestyle='--', linewidth=2, label='expected signal range')
        plt.axvline(x=self.expected_freq_high, color='red', linestyle='--', linewidth=2)
        plt.axvline(x=self.fdm.frequency, color='red', linestyle=':', linewidth=3, label='detected signal', alpha=0.5)
        label = '60 Hz harmonic'
        for harmonic in powerline:
            plt.axvline(harmonic, color='grey', linestyle='--', linewidth=1, label=label)
            label = None
        plt.xlim(1970,2120)
        plt.title('Proton Precession Magnetometer')
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Amplitude')
        font = FontProperties().copy()
        font.set_weight('bold')
        plt.grid()
        plt.legend()
        plt.show()

    # Create some nice plots.
    # 
    # A large page multiplot for webpage: time domain filtered FID, filtered envelope, FFT FID
    # history of our data.  Maybe add: intermagnet station for reference? filtered dB/dt?
    #
    def multiplot_results(self, plot_fname=None):
        # 60 Hz harmonics 
        powerline = np.asarray(range(60,100*60,60)) 

        # Gather recent processed data from the log file
        t_recent = []
        f_recent = []
        #t_earliest = time.time() - (4*24*60*60)     # Limit to most recent 4 days
        t_earliest = time.time() - (1*24*60*60)     # Limit to most recent 1 days
        with open(self.log_fname, 'r') as f:
            for line in f:
                # ts, frequency, amplitude, decay, Q, error, fom, wb_snr))
                data = [float(i) for i in line.split()]
                ts = data[0]
                if ts >= t_earliest:
                    # timestamp in a format suitable for plotting
                    t_recent.append(mdates.date2num(dt.datetime.utcfromtimestamp(ts)))  # UTC time
                    f_recent.append(data[1] - self.local_offset)        # offset corrected

        # Multiple plots of different sizes
        fig = plt.figure(1)
        mpl.rcParams['font.size'] =  12
        fig.suptitle('PyPPM Proton Precession Magnetometer. %s. Updated %s UTC.' % 
            (location, dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")))
            #(dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")))
        mpl.rcParams['font.size'] = 9
        gs = gridspec.GridSpec(8, 3)         # 8 rows, 3 columns. TODO wspace

        # Plot time domain
        ax = plt.subplot(gs[0:2, 0])        # top rows, first column
        plt.plot(np.array(self.t0), np.array(self.a0), 'b', label='background')
        plt.plot(np.array(self.t1), np.array(self.a1), 'r', label='measurement')
        plt.title('%d ADC samples at %.1f kHz' % (len(self.t0), self.fs/1000.0))
        plt.ylabel('Voltage (V)')
        plt.xlabel('Time (s)')
        plt.legend(loc='lower right', framealpha=0.7)

        # Plot filtered time domain and FDM decay fit
        ax = plt.subplot(gs[0:2, 1])        # top rows, second column
        plt.plot(np.array(self.t1), np.array(self.a1f), 'g')
        fdm_decay = 1.414 * self.fdm.amplitude * np.exp(-self.fdm.decay * np.array(self.t1))
        plt.plot(np.array(self.t1), fdm_decay, color='black', linestyle='--')
        plt.plot(np.array(self.t1), -fdm_decay, color='black', linestyle='--')
        plt.title('Filtered ADC and FDM fit')
        plt.ylabel('Voltage (V)')
        plt.xlabel('Time (s)')

        # Plot text
        ax = plt.subplot(gs[0:2, 2])        # top rows, third column
        ax.axis('off')
        ax.text(0.0, 1.00, 'F(scalar):', ha='left', va='top', fontsize=12)
        ax.text(0.5, 1.00, '%.2f nT' % (self.field), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.85, 'Frequency:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.85, '%.3f Hz' % (self.fdm.frequency), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.70, 'Amplitude:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.70, '%.1f mVpp' % (2*1.414*self.fdm.amplitude*1000), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.55, 'Tau2:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.55, '%.2f sec' % (1/self.fdm.decay), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.40, 'FOM:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.40, '%.1g' % (self.fdm.error), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.25, 'NB SNR:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.25, '%.1f dB' % (self.fdm.nb_snr), ha='left', va='top', fontsize=12)
        ax.text(0.0, 0.10, 'WB SNR:', ha='left', va='top', fontsize=12)
        ax.text(0.5, 0.10, '%.1f dB' % (self.wb_snr), ha='left', va='top', fontsize=12)

        # Plot frequency domain
        ax = plt.subplot(gs[2:4, :])        # next rows, all columns
        plt.plot(np.array(self.fft_f0), np.array(self.fft_a0), 'b', label='background')
        plt.plot(np.array(self.fft_f1), np.array(self.fft_a1), 'r', label='measurement')
        plt.plot(np.array(self.fft_f1f), np.array(self.fft_a1f), 'g', label='post-filtered')
        plt.axvspan(self.expected_freq_low-100, self.expected_freq_low, alpha=0.2, color='grey')
        plt.axvspan(self.expected_freq_high, self.expected_freq_high+100, alpha=0.2, color='grey')
        plt.axvline(x=self.fdm.frequency, color='red', linestyle=':', linewidth=3, label='detected signal', alpha=0.5)
        label = '60 Hz harmonic'
        for harmonic in powerline:
            plt.axvline(harmonic, color='grey', linestyle='--', linewidth=1, label=label)
            label = None
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Amplitude')
        font = FontProperties().copy()
        font.set_weight('bold')
        plt.grid()
        ax.legend(loc='upper left', framealpha=0.7, fontsize=10)
        ax2 = ax.twinx()
        ax2.plot(self.filt_freq, self.filt_mag, '-', label='filter', color='orange')
        ax2.legend(loc='upper right', framealpha=0.7)
        plt.xlim(self.expected_freq_low-100, self.expected_freq_high+100)

        # Recent field strength plot
        ax = plt.subplot(gs[4:8, :])        # next rows, all columns
        timefmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S\n')      # workaround TZ issue
        #timefmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S\n%Z')   
        #ax.xaxis.set_major_locator(mdates.DayLocator()) 
        #ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6)) 
        ax.xaxis.set_major_formatter(timefmt)
        plt.plot([x[0] for x in self.intermag_recent], [x[1] for x in self.intermag_recent], '-', color='red', label=self.intermag_name)
        plt.plot([x[0] for x in self.geometrics_recent], [x[1] for x in self.geometrics_recent], '-', color='orange', label=self.geometrics_name)
        plt.plot(t_recent, f_recent, '.', markersize=1, label='PyPPM', color='purple')  # offset corrected
        plt.ylim(self.expected_field_low+1000, self.expected_field_high-925)
        plt.ylabel('Fscalar (nT)')
        plt.grid()
        plt.legend(loc='upper left', framealpha=0.7)
        ax.xaxis.grid(which='minor', linestyle='--')
    
        # Textual data
        # Annotate the plot: fscalar, precession frequency, fid amplitude, q, polarization time, tau2, FOM, 
        # actual sample rate, number samples, NB SNR, total meas, success rate any metadata text
        #t = 'F(scalar): 47423 nT     Frequency: 2032.423 Hz   Tau2: 1.82 sec    Amplitude: 0.123 V   FOM: 3.2e-7    Q: 10234        \n' \
        #    'Polarization: 2.0 sec   Fs: 10002.35 Hz          Nsamples: 16384                        Total: 24232   Rejected: 234   \n' \
        #    '01-01-1923 00:01:02 UTC                                                                                                \n' 
        #text = fig.text(0.5, 0.00, t, ha='center', va='bottom', size=10, family='monospace')

        # show and save
        if self.interactive:
            plt.show()
        #fig.set_size_inches(w=11,h=8.5)
        fig.set_size_inches(11, 8.5)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])       # leave some space for suptitle
        fig.savefig(plot_fname, bbox_inches='tight')        # FIXME should we add timestamp to filename and create a link to latest?


    # Generate some known signals for testing.
    def generate_test_data(self, fs=None, length=16384):
        if fs == None:
            fs = self.fs

        t0 = tuple(np.arange(0, float(length)/fs, 1.0/fs))  # pyppm requires tuples
        t1 = t0
        n0 = np.random.normal(0, 0.005, len(t0))        # background noise - std dev of 5mV 
        # FIXME - add BP filtered background noise

        # 1 Vrms amplitude signal centered at expected frequency
        f = field_to_freq(self.expected_field_nt)
        a1 = np.sqrt(2) * np.sin(2 * np.pi * f * t1)        # 1 Vrms signal
        n1 = np.random.normal(0, 0.005, len(t1))        # noise - std dev of 5mV 
        a1 += n1
        # FIXME add 60 Hz harmonics
        self.a0 = a0
        self.t0 = t0
        self.a1 = a1
        self.t1 = t1
        self.fs = fs


def usage(s):
    print \
'''
usage: python %s [OPTION]... 

Analyze the data files produced by the PyPPM capture script.
This script will scan all files with the given base filename (including
all the different timestamps). It will process them and record results
to the output log file, if it hasn't already been done.

  -f, --filename=FILENAME    Base filename, appended with timestamp
  -o, --output=STRING        Output log file name
  -p, --plot=FILENAME        Generate a plot to the specified file.
  -v, --verbose              Add this flag (multiple times) to increase verbosity
  -w, --wavfile=FILENAME     Create an audio file of the filtered proton signal
  -a, --audio                Play the audio file after each analysis
  -t, --time=SECONDS         Only process newer files. (The unix epoch timestamp if >0, else number of seconds into the past.)
  -i, --intermagnet=PATH     Path to Intermagnet data files (such as /directory/*.min)
  -g, --geometrics=PATH      Path to Geometrics data files (such as /directory/*.dat)
  -O, --offset=FLOAT         Set our local magnetic offset in nanoTesla
''' % (s)
    sys.exit(1)

def main():
    # Defaults
    plot = False
    testdata = False
    verbose = 0
    play_audio = False
    output_fname = 'analysis.out'
    input_fname = 'data/ppm_measurements.dat'
    wavfname = None
    earliest_filetime = 0
    intermagnet_path = None 
    geometrics_path = None 
    #local_offset = -369.0
    local_offset = 0
    plot_decimation = 1

    # Parse command line options
    try:
        opts,args = getopt.getopt (sys.argv[1:],
                'f:o:vp:d:w:at:i:g:hO:',
                ['filename=', 'output=', 'verbose', 'plot', 'plot_decimation=', 'wavfile=', 'audio', 'time=',
                 'intermagnet=', 'geometrics=', 'offset=', 'help'])
    except getopt.GetoptError:
        usage (sys.argv[0])
        sys.exit (1)
    for o,a in opts:
        if o in ('-h', '--help'):
            usage (sys.argv[0])
        elif o in ('-o', '--output'):
            output_fname = a
        elif o in ('-v', '--verbose'):
            if a != None and len(a) > 0:
                verbose = int(a)
            else:
                verbose += 1
        elif o in ('-t', '--time'):
            earliest_filetime = float(a)
            if earliest_filetime < 0:
                earliest_filetime += time.time()
        elif o in ('-p', '--plot'):
            plot = True
            if a:
                plot_fname = a
            else:
                plot_fname = 'ppm_plot.png'
        elif o in ('-d', '--plot_decimation'):
            plot_decimation = int(a)
        elif o in ('-f', '--filename'):
            input_fname = a
        #elif o in ('-t', '--testdata'):
            #testdata = True
        elif o in ('-w', '--wavfile'):
            wavfname = a
        elif o in ('-a', '--audio'):
            play_audio = True
        elif o in ('-i', '--intermagnet'):
            intermagnet_path = a
        elif o in ('-g', '--geometrics'):
            geometrics_path = a
        elif o in ('-O', '--offset'):
            local_offset = float(a)


    # Process all pending data:
    if verbose > 0:
        print 'Preparing file list'

    # Read previously recorded and timestamped files
    counter = 0
    for f in glob.glob(input_fname + '.*'):

        # Time-based file filtering
        if os.path.getctime(f) < earliest_filetime:
            continue    # skip

        fname, fext = os.path.splitext(f);
        ts = int(fext[1:])

        # Skip the input file if it already exists in the output file
        Skip = False
        try:
            for l in open(output_fname, 'r'):
                if str(ts) in l:
                    Skip = True
                    break
        except IOError:
            pass

        if Skip == False:
            ppm = ppm_analysis(output_fname, verbose)
            ppm.set_local_offset(local_offset)
            try:
                ppm.load_data_file(f, freq_error=+60e-6)
            except:
                print "Failed loading file", f
                continue

            #ppm.set_expected_range(47000, 49000)
            #ppm.set_expected_range(50000, 52700)
            #ppm.set_expected_range(51100, 51700)
            #ppm.set_expected_range(50000, 52600)
            ppm.set_expected_range(51960-1000, 51960+1000)

            # FIXME - wasteful to load each time...
            if intermagnet_path:
                ppm.load_intermagnet_data(intermagnet_path, 'Intermagnet')

            if geometrics_path:
                ppm.load_geometrics_data(geometrics_path, 'Geometrics')

            # try to identify the FID signal
            field, fdm, fom, nb_snr, wb_snr = ppm.analyze(plot_fname=plot_fname, wavfilename=wavfname)
            if fdm == None:
                print 'Failed to find a signal in file %s' % (f)
                continue

            # append analysis result to logfile
            print   ts, 'Fscalar:', field, 'Freq:', fdm.frequency, 'Amplitude:', fdm.amplitude, 'Decay:', fdm.decay, \
                'Q:', fdm.Q, 'FOM:', fdm.error, 'NB SNR:', nb_snr, 'WB SNR:', wb_snr

            with open(output_fname, 'a') as outf:
                outf.write('%u %f %f %f %f %f %g %f %f\n' % \
                       (ts, field, fdm.frequency, fdm.amplitude, fdm.decay, fdm.Q, fdm.error, fom, wb_snr))

            # Plot only every Nth analysis, since plotting takes a bit of time.
            if plot_fname and counter % plot_decimation == 0:
                if verbose:
                    print 'Plotting...'
                ppm.multiplot_results(plot_fname)
                if verbose:
                    print 'Plotting complete.'

            if play_audio == True:
                os.system('afplay protons.wav')     # MacOS only

            counter += 1

if __name__ == "__main__":
    main()
