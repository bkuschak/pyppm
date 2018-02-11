#!/bin/sh
#
# Run the capture script in background at a higher priority.  Need to be root to adjust priority, 
# but don't run the capture script as root.

USER=$(whoami)
sudo -b nice -n -5 sudo -b -u $USER sh -c 'nohup python pyppm_capture_raw.py -p4 -l32768 -P2 -q > /dev/null 2>&1  &'

