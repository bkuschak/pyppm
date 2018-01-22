# Analyze the data files collected by the proton precession magnetometer.
# Using the FDM method of harmonic inversion, calculate the FID precession
# frequency, and field strength.
# 
# Status: works reasonably well for clean data, but in the presence of large
# background noise it frequently selects the wrong signal as the FID signal.
# The calculated SNR is pretty meaningless at this point.
#
# 1/13/18 B. Kuschak <bkuschak@yahoo.com> 
#

import pyppm			
import harminv				# https://github.com/aaren/pharminv

import sys
import time
import pytz
import glob, os
import pickle
import getopt
import datetime as dt
import scipy.io.wavfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from matplotlib.font_manager import FontProperties
from scipy import signal

class ppm_analysis:
	#gp = 42.57747892		# gyromagnetic ratio of proton
	gp = 42.576			# actual gyromagnetic ratio of proton in earth field

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
		self.interactive = True
		self.intermag_recent = []

	def set_local_offset(self, offset):
		self.local_offset = offset
	
	def debug(self, level, msg):
		if self.debug_level >= level:
			print msg
	
	# Read raw ADC samples from the data file
	# (t0,a0) is the background. (t1,a1) is the real measurement.
	# If freq_error=None, then error value is taken from the data file 
	def load_data_file(self, fname, freq_error=None):
		file_data = pickle.load(open(fname, 'rb'))
		self.a0 = file_data['a0']				# background
		self.a1 = file_data['a1']				# FID measurement
		self.fs = file_data['fs']				# nominal sample rate
		if freq_error != None:
			self.freq_error = freq_error			# override sample rate clock error
		else:
			self.freq_error = file_data['freq_error']	# error in the sample rate, from datafile
		self.start_time = file_data['start_time']		# Unix epoch of first (background) measurement
		self.polarize = file_data['polarize_time']		# polarization time
		self.comments = file_data['comments']
		self.data_length = len(self.a0)
		self.freq_error = float(self.freq_error)

		if self.verbose > 1:
			print   'Data File: Fs = %u KHz (%+.1f ppm), %u samples. Polarization time = %.3f sec.\n' \
				'Comments: %s' \
				% (self.fs, self.freq_error*1e6, self.data_length, self.polarize, self.comments)

		# Create time arrays
		ts = 1.0 / (self.fs * (1.0 + self.freq_error))
		self.t0 = tuple(np.arange(0, len(self.a0) * ts, ts))	# pyppm requires tuples
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
	def analyze(self, plot=False, wavfilename=None):
		# narrowband filter the data around the expected peak.  
		fnaught = self.field_to_freq(self.expected_field)
		Q = 60.0
		w0 = fnaught/(self.fs/2)	# normalized
		b, a = signal.iirpeak(w0, Q)
		self.a1f = signal.filtfilt(b, a, self.a1) 

		# FIXME - we should truncate the beginning and end that are subject to filter startup, so they don't throw off
		# the harmonic inversion decay calculation

		if wavfilename != None:
			self.save_wavfile(wavfilename)

		# Use FDM harmonic inversion to calculate a set of frequencies in this band
		self.signals = harminv.invert(self.a1f, fmin=1000, fmax=3000, dt=1.0/self.fs)

		# Select best candidate signal, or None if we couldn't find one
		self.fdm = self.select_candidate(self.signals)
		if self.fdm == None:
			return None, None, 0.0, 0.0

		self.fdm_time_constant = 1.0 / self.fdm.decay
		self.field = self.freq_to_field(self.fdm.frequency)

		# See http://www.gellerlabs.com/PMAG%20Docs.htm
		# Like GellerLabs, compute an SNR which is 20 log (FID amplitude / RMS_Sum(other inband amplitudes))
		# This doesn't appear to be very meaningful, based on the way FDM works.  We should compute this differently.
		rms_sum = 0
		for s in self.signals:				
			# restrict to narrow band
			if  s.frequency >= self.expected_freq_low and s.frequency <= self.expected_freq_high and \
			    s.frequency != self.fdm.frequency:
				rms_sum += (s.amplitude ** 2)
		self.fdm_snr = 20.0 * np.log10(self.fdm.amplitude / np.sqrt(rms_sum))

		# GellerLabs calls the FDM error the Figure of Merit. 
		# Shouldn't be converted to nT, as harminv says 'error is not really error bars on frequency'
		self.fom_nt = self.freq_to_field(self.fdm.error * self.fdm.frequency)
		if self.verbose > 0:
			print "FDM SNR (dB)", self.fdm_snr, "FOM", self.fdm.error

		if plot:
			#self.plot_results()
			self.multiplot_results()

		return self.field, self.fdm, self.fom_nt, self.fdm_snr

	# given the signals resulting from harmonic inversion, select the most likely FID signal
	# FIXME - use the background capture data to improve selection
	def select_candidate(self, signals):

		# mode frequencies, absolute amplitudes, phase shift, decay rates (reciprocal of time constant), Q factor, error
		if self.verbose > 0:
			print 'freq, amplitude, phase, decay, q, error'
			print signals

		# find the right one. just picking the highest peak isn't enough.
		# the proton signal should have a small decay parameter (0.3 to 2.0?), high Q (over 5000?), low error, largish amplitude
		idx = 0
		while idx < len(signals):
			s = signals[idx]
			if s.frequency < self.expected_freq_low or s.frequency > self.expected_freq_high:
				signals = np.delete(signals, idx)
				if self.verbose > 1:
					print "discarding bad freq ", s
			elif s.decay < 0.05 or s.decay > 3.0:
				signals = np.delete(signals, idx)
				if self.verbose > 1:
					print "discarding bad decay", s
			elif abs(s.Q) < 2000.0:
				signals = np.delete(signals, idx)
				if self.verbose > 1:
					print "discarding bad q    ", s
			elif s.error > 1e-5:
				signals = np.delete(signals, idx)
				if self.verbose > 1:
					print "discarding bad error", s
			elif s.amplitude < 0.02 or s.amplitude > 0.500:
				signals = np.delete(signals, idx)
				if self.verbose > 1:
					print "discarding bad amplitude", s
			else:
				idx += 1

		if self.verbose > 0:
			print 'candidates: freq, amplitude, phase, decay, q, error'
			print signals

		# We may not find any that match our criteria...
		if len(signals) == 0:
			if self.verbose > 0:
				print 'None'
			return None

		# Prefer higher Q, higher amplitude, lower error
		best_q = signals[np.argmax(abs(signals.Q))]
		best_ampl = signals[np.argmax(signals.amplitude)]
		best_error = signals[np.argmin(signals.error)]

		# Select the one that gets at least 2 out of 3
		if best_q == best_ampl:
			if self.verbose > 1:
				print 'best q, ampl', best_q
			return best_q
		elif best_q == best_error:
			if self.verbose > 1:
				print 'best q, error', best_q
			return best_q
		elif best_ampl == best_error:
			if self.verbose > 1:
				print 'best ampl, error', best_ampl
			return best_ampl
		else:
			return None		# shouldn't happen

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
	def multiplot_results(self):
		# 60 Hz harmonics 
		powerline = np.asarray(range(60,100*60,60)) 

		# Compute FFTs
		self.a1f = tuple(self.a1f)
		(f0, A0) = pyppm.fft(self.t0, self.a0)
		(f1, A1) = pyppm.fft(self.t1, self.a1)
		(f1f, A1f) = pyppm.fft(self.t1, self.a1f)

		# Gather recent processed data from the log file
		t_recent = []
		f_recent = []
		t_earliest = time.time() - (3*24*60*60) 	# Limit to most recent 3 days
		with open(self.log_fname, 'r') as f:
			for line in f:
				# ts, frequency, amplitude, decay, Q, error, fom_nt, fdm_snr))
				data = [float(i) for i in line.split()]
				ts = data[0]
				if ts >= t_earliest:
					# timestamp in a format suitable for plotting
					t_recent.append(mdates.date2num(dt.datetime.utcfromtimestamp(ts)))	# UTC time
					f_recent.append(data[1] - self.local_offset)		# offset corrected

		# Multiple plots of different sizes
		fig = plt.figure(1)
		mpl.rcParams['font.size'] =  12
		fig.suptitle('PyPPM Proton Precession Magnetometer')
		mpl.rcParams['font.size'] = 9
		gs = gridspec.GridSpec(8,2)

		# Plot time domain
		ax = plt.subplot(gs[0:2, 0])		# top rows, first column
		plt.plot(np.array(self.t0), np.array(self.a0), 'b', label='background')
		plt.plot(np.array(self.t1), np.array(self.a1), 'r', label='measurement')
		plt.title('Raw ADC measurements')
		plt.ylabel('Voltage (V)')
		plt.xlabel('Time (s)')
		plt.legend(loc='lower left')

		# Plot filtered time domain and FDM decay fit
		ax = plt.subplot(gs[0:2, 1])		# top rows, second column
		plt.plot(np.array(self.t1), np.array(self.a1f), 'g')
		fdm_decay = 1.414 * self.fdm.amplitude * np.exp(-self.fdm.decay * np.array(self.t1))
		plt.plot(np.array(self.t1), fdm_decay, color='grey', linestyle='--')
		plt.plot(np.array(self.t1), -fdm_decay, color='grey', linestyle='--')
		plt.title('Filtered ADC measurements and FDM fit')
		plt.ylabel('Voltage (V)')
		plt.xlabel('Time (s)')

		# Plot frequency domain
		ax = plt.subplot(gs[2:4, :])		# next rows, all columns
		plt.plot(np.array(f0), np.array(A0), 'b', label='background')
		plt.plot(np.array(f1), np.array(A1), 'r', label='measurement')
		plt.plot(np.array(f1f), np.array(A1f), 'g', label='post-filtered')
		plt.axvline(x=self.expected_freq_low, color='red', linestyle='--', linewidth=2, label='expected signal range')
		plt.axvline(x=self.expected_freq_high, color='red', linestyle='--', linewidth=2)
		plt.axvline(x=self.fdm.frequency, color='red', linestyle=':', linewidth=3, label='detected signal', alpha=0.5)
		label = '60 Hz harmonic'
		for harmonic in powerline:
			plt.axvline(harmonic, color='grey', linestyle='--', linewidth=1, label=label)
			label = None
		plt.xlim(self.expected_freq_low-100, self.expected_freq_high+100)
		plt.xlabel('Freq (Hz)')
		plt.ylabel('Amplitude')
		font = FontProperties().copy()
		font.set_weight('bold')
		plt.grid()
		plt.legend(loc='upper left')

		# Recent field strength plot
		#ax = plt.subplot(gs[4:7, :])		# next rows, all columns
		ax = plt.subplot(gs[4:8, :])		# next rows, all columns
		timefmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S\n%Z')	
		#ax.xaxis.set_major_locator(mdates.DayLocator()) 
		#ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6)) 
		ax.xaxis.set_major_formatter(timefmt)
		plt.plot(t_recent, f_recent, '.', markersize=1, label='PyPPM')	# offset corrected
		plt.plot([x[0] for x in self.intermag_recent], [x[1] for x in self.intermag_recent], '-', color='red', label=self.intermag_name)
		plt.ylabel('Fscalar (nT)')
		plt.grid()
		plt.legend(loc='upper left')
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
		fig.set_size_inches(w=11,h=8.5)
		fig.tight_layout(rect=[0, 0.03, 1, 0.95])		# leave some space for suptitle
		fig.savefig('ppm_plot.png', bbox_inches='tight')	# FIXME should we add timestamp to filename and create a link to latest?


	# Generate some known signals for testing.
	def generate_test_data(self, fs=None, length=16384):
		if fs == None:
			fs = self.fs

		t0 = tuple(np.arange(0, float(length)/fs, 1.0/fs))	# pyppm requires tuples
		t1 = t0
		n0 = np.random.normal(0, 0.005, len(t0))		# background noise - std dev of 5mV 
		# FIXME - add BP filtered background noise

		# 1 Vrms amplitude signal centered at expected frequency
		f = field_to_freq(self.expected_field_nt)
		a1 = np.sqrt(2) * np.sin(2 * np.pi * f * t1)		# 1 Vrms signal
		n1 = np.random.normal(0, 0.005, len(t1))		# noise - std dev of 5mV 
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
  -p, --plot                 Generate a plot file ppm_plot.png
  -v, --verbose              Add this flag (multiple times) to increase verbosity
  -w, --wavfile=FILENAME     Create an audio file of the filtered proton signal
  -a, --audio                Play the audio file after each analysis
  -t, --time=SECONDS         Only process newer files. (The unix epoch timestamp if >0, else number of seconds into the past.)
  -i, --intermagnet=PATH     Path to Intermagnet data files (such as /directory/*.min)
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

	# Parse command line options
	try:
		opts,args = getopt.getopt (sys.argv[1:],
				'f:o:vpw:at:i:h',
				['filename=', 'output=', 'verbose', 'plot', 'wavfile=', 'audio', 'time=',
				 'intermagnet=', 'help'])
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

	# Process all pending data:
	# Read previously recorded and timestamped files
	for f in glob.glob(input_fname + '.*'):

		# Time-based file filtering
		if os.path.getctime(f) < earliest_filetime:
			continue	# skip

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
			ppm.set_local_offset(-369.0)
			ppm.load_data_file(f, freq_error=+60e-6)
			ppm.set_expected_range(47000, 49000)

			# FIXME - wasteful to load each time...
			if intermagnet_path:
				ppm.load_intermagnet_data(intermagnet_path, 'Intermagnet')

			# try to identify the FID signal
			field, fdm, fom_nt, fdm_snr = ppm.analyze(plot=plot, wavfilename=wavfname)
			if fdm == None:
				print 'Failed to find a signal in file %s' % (f)
				continue

			# append analysis result to logfile
			print 	ts, 'Fscalar:', field, 'Freq:', fdm.frequency, 'Amplitude:', fdm.amplitude, 'Decay:', fdm.decay, \
				'Q:', fdm.Q, 'FOM:', fdm.error, 'SNR:', fdm_snr

			with open(output_fname, 'a') as outf:
				outf.write('%u %f %f %f %f %f %g %f %f\n' % \
					   (ts, field, fdm.frequency, fdm.amplitude, fdm.decay, fdm.Q, fdm.error, fom_nt, fdm_snr))

			if play_audio == True:
				os.system('afplay protons.wav')		# MacOS only

if __name__ == "__main__":
	main()
