# Run the PPM and capture raw data continuously, logging to files.
#
# B. Kuschak <bkuschak@yahoo.com> 1/13/18
#
import pyppm
import time
import pickle
import getopt
import sys

class ppm_logger():

	def __init__(self, sample_rate_khz, data_length, polarization_time, delay, output_fname, comments, freq_error=0.0, protocol=0):
		self.sample_rate_khz = sample_rate_khz
		self.frequency_error = freq_error
		self.data_length = data_length
		self.polarization_time = polarization_time
		self.delay = delay
		self.output_fname = output_fname
		self.comments = comments
		self.freq_error = float(freq_error)
		self.protocol = protocol

	def connect(self):
		# connect to the device.
		self.dev = pyppm.PPM()

		# download a pulse program.
		self.dev.pulprog = [
		  [pyppm.ACQUIRE, self.data_length, self.sample_rate_khz],
		  [pyppm.RELAY, True],
		  [pyppm.DEADTIME, 40],				# relay settling time
		  [pyppm.POLARIZE, True],
		  [pyppm.DELAY, self.polarization_time],
		  [pyppm.POLARIZE, False],
		  [pyppm.DEADTIME, 2],				# allow full decay before switching relay
		  [pyppm.RELAY, False],
		  [pyppm.DEADTIME, 20], 			# relay takes about ~3 msec to switch, followed by some settling
		  [pyppm.ACQUIRE, self.data_length, self.sample_rate_khz],
		  [pyppm.END]
		]

	# run the pulse program and save data to a file.  Pickled format is not the most space efficient,
	# but for development purposes it's easy to use 
	def acquire(self, comments=None, protocol=None):
		ts = time.time()

		# execute the pulse program.
		(t, a) = self.dev.execute()

		# slice out the first acquisition.
		#t0 = t[0 : 16384]
		#a0 = a[0 : 16384]
		t0 = t[0 : self.data_length]
		a0 = a[0 : self.data_length]

		# slice out the second acquisition.
		#t1 = t[16384 : 32768]
		#a1 = a[16384 : 32768]
		t1 = t[self.data_length : 2*self.data_length]
		a1 = a[self.data_length : 2*self.data_length]

		if self.output_fname: 
			# Append Unix epoch timestamp to end of filename
			fname = self.output_fname + '.' + str(int(ts))

			if comments == None:
				comments = self.comments
				
			if protocol == None:
				protocol = self.protocol
				
			# What about other things like board temperature, fluid temperature, polarization current?
			# Should include timestamps of each acquisition
			file_data = { 	'start_time' : ts,			# time of first acquisition
					'comments': self.comments, 
					'fs' : self.sample_rate_khz * 1000.0,	# nominal
					'freq_error' : self.frequency_error,	# example: +20e-6 for clock running 20ppm fast
					'polarize_time' : self.polarization_time,
					#'t0' : t0,				# first acquisition (background)
					'a0' : a0,
					#'t1' : t1,				# second acquisition (FID)
					'a1' : a1 }

			# File format is Python pickle format
			with open(fname, 'wb') as data_file:
				pickle.dump(file_data, data_file, protocol)

def usage(s):

  print '''
usage: python %s [OPTION]... 

Run the pulse program periodically, logging collected data to files.

  -o, --output=FILENAME      base filename, appended with timestamp
  -f, --fs=KHZ               sampling frequency in kHz (FIXME resolution?)
  -e, --freqerror=FLOAT      sampling frequency error multiplier (-10e-6 for clock running 10ppm slow, default 0)
  -l, --length=SAMPLES       number of samples for each acquisition (default 16384)
  -p, --polarization=SEC     polarization time (default 2)
  -d, --delay=SEC            delay a number of seconds between collections (default 8)
  -c, --comment=STRING       comment to include in the data file
  -P, --protocol=0,1,2       Python pickle protocol (default 2)
  -q, --quiet                Don't print anything
''' % (s)

def main():
	# defaults
	sample_rate_khz = 20
	freq_error = 0.0
	data_length = 16384		# each acquisition
	polarization_time = 2		
	delay = 8
	output_fname = 'data/ppm_measurements.dat'
	protocol = 2			# pickle protocol
	comments = None
	quiet = False

	# Parse command line options
	try:
		opts,args = getopt.getopt (sys.argv[1:],
				'ho:l:f:e:p:d:c:qP:',
				['help', 'output=', 'length=', 'fs=', 'freqerror=', 
				 'polarization=', 'delay=', 'count=', 'quiet=',
				 'protocol='])
	except getopt.GetoptError:
		usage(sys.argv[0])
		sys.exit (1)
	for o,a in opts:
		if o in ('-h', '--help'):
			usage(sys.argv[0])
			sys.exit(1)
		elif o in ('-o', '--output'):
			output_fname = a
		elif o in ('-l', '--length'):
			data_length = int(a)
		elif o in ('-f', '--fs'):
			sample_rate_khz = int(a)
		elif o in ('-e', '--freqerror'):
			freq_error = float(a)
		elif o in ('-p', '--polarization'):
			polarization_time = float(a)
		elif o in ('-d', '--delay'):
			delay = float(a)
		elif o in ('-q', '--quiet'):
			quiet = True
		elif o in ('-P', '--protocol'):
			protocol = int(a)
		elif o in ('-c', '--comment'):
			comments = a

	if not quiet:
		print 	'Acquiring: Fs = %u KHz (%+.1f ppm), %u samples. Polarization time = %.3f sec, cycle delay = %.3f sec.\n' \
			'Base filename: "%s". Protocol = %d\nComments: %s' \
			% (sample_rate_khz, freq_error*1e6, data_length, polarization_time, delay, output_fname, protocol, comments)

	# connect to the device.
	ppm = ppm_logger(sample_rate_khz, data_length, polarization_time, delay, output_fname, comments, protocol=protocol)
	ppm.connect()

	# Run forever
	while True:
		ppm.acquire()
		time.sleep(delay)


if __name__ == "__main__":
	main()

