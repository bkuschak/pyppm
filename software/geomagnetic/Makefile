INSTALL_BIN	:= /usr/local/bin
INSTALL_SERVICE	:= /etc/systemd/system/

all: 

install: install-bin install-systemd-service

install-bin:
	chmod +x override_usb_power.sh
	cp pyppm_capture_raw.py $(INSTALL_BIN)
	cp override_usb_power.sh $(INSTALL_BIN)

# Run as sudo make install
install-systemd-service:
	-systemctl stop pyppm.service
	-systemctl stop pyppm-analyze.service
	cp pyppm.service pyppm-analyze.service $(INSTALL_SERVICE)
	systemctl daemon-reload
	systemctl enable pyppm.service
	systemctl enable pyppm-analyze.service
	systemctl start pyppm.service
	systemctl start pyppm-analyze.service
	systemctl status pyppm.service
	systemctl status pyppm-analyze.service

.PHONY: install install-bin install-systemd-service
