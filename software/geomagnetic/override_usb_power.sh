#!/bin/bash
#
# This is a hack to get around around reported insufficient power on the USB bus
# when using PyPPM on an embedded system.  
#
# Be careful, this could be dangerous!
#

ID_BP=$(dmesg |grep "insufficient available bus power" |tail -1 |sed 's?.*usb ??g;s?:.*??g')
ID_SERNUM=$(dmesg |grep "SerialNumber: PyPPM" |tail -1 |sed 's?.*usb ??g;s?:.*??g')

# Disable power safety mechanism
if [ "$ID_BP" == "$ID_SERNUM" -a "$ID_BP" != "" ]; then

	echo Disabling USB current limit on PyPPM USB \"$ID_BP\"
	sh -c "echo 1 > /sys/bus/usb/devices/$ID_BP/bConfigurationValue"
	sleep 2
fi
