#!/bin/bash

echo "Script for exporting files starts!"
for files in 1 2
do
	echo "file $files"
	sudo tshark -a duration:3600 -T fields -e wlan.fc.type_subtype -e frame.len -e wlan_radio.data_rate -e wlan_radio.duration -e wlan.fcs.status -E header=y -E separator=/t > library$files.txt
done
