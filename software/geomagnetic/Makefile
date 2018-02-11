INSTALL_BIN	:= /usr/local/bin
INSTALL_SERVICE	:= /etc/systemd/system/

all: 

install: install-bin install-systemd-service

install-bin:
	cp pyppm_capture_raw.py $(INSTALL_BIN)

# Run as sudo make install
install-systemd-service:
	-systemctl disable pyppm.service
	cp pyppm.service $(INSTALL_SERVICE)
	systemctl enable pyppm.service
	systemctl start pyppm.service
	systemctl status pyppm.service

.PHONY: install install-bin install-systemd-service
