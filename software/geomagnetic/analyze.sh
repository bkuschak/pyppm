#!/bin/bash
# This script is run from the pyppm-analyze systemd service. 
# Generate plots and upload to webserver.

PNG=/var/www/html/ppm_plot.png
#WAV=/var/www/html/audio.wav
WEBSITE=webserver_host:~/pyppm/

#python analyze.py -t -30 -f /data/ppm_measurements.dat -v -w audio.wav 
#python3 analyze.py -t -30 -f /data/ppm_measurements.dat -v -w audio.wav -p
#python analyze.py -t -12600 -f /data/ppm_measurements.dat -v -v
#python analyze.py -t -300 -f /data/ppm_measurements.dat -p
#python analyze.py -t -3600 -f /data/ppm_measurements.dat
#python analyze.py -t -20 -f /data/ppm_measurements.dat  -w audio.wav
#python analyze.py -t -21600 -f /data/ppm_measurements.dat
#python analyze.py -t -20 -f /data/ppm_measurements.dat -p ${PNG} -w ${WAV}
#python analyze.py -t -21600 -f /data/ppm_measurements.dat -p ${PNG} -w ${WAV}
#python analyze.py -t -86400 -f /data/ppm_measurements.dat -p ${PNG} -w ${WAV}
#python analyze.py -t -259200 -f /data/ppm_measurements.dat -p ${PNG}
#python analyze.py -t -20 -f /data/ppm_measurements.dat -p ${PNG} -v
#python analyze.py -t -1800 -f /data/ppm_measurements.dat -p ${PNG}
#sudo cp audio.wav /var/www/html
#sudo cp ppm_plot.png /var/www/html

while true; do
	python analyze.py -t -43240 -f /data/ppm_measurements.dat -p ${PNG} --plot_decimation 30 -vv
    scp $(PNG) $(WEBSITE)
    sleep 10
done
