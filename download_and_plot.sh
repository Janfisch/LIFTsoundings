#!/bin/bash

#this gives me name of the latest file on SFTP server
#fileNameVS=$(echo "ls -1rt *0.bfr" | sftp -P 24 sftpuser@thredds.atmohub.kit.edu:/home/sftpuser/sftp-upload/RASO/VS/ | tail -1)
#fileNameVS=$(echo "ls -1rt *0.bfr" | sftp -P 24 sftpuser@thredds.atmohub.kit.edu:/home/sftpuser/sftp-upload/RASO/KOE/ | tail -1)
fileNameVS=$((cd /home/jfischer/KITCUBEmount/sftp-upload/RASO/VS/ && ls -1rt bufr309057*0.bfr ) | tail -1)
fileNameKO=$((cd /home/jfischer/KITCUBEmount/sftp-upload/RASO/KOE/ && ls -1rt bufr309057*2.bfr ) | tail -1)

source /home/jfischer/mambaforge/etc/profile.d/conda.sh
conda activate LIFTsoundings

# check if this file already exists, otherwise copy and create plot
name=${fileNameVS%.bfr}.png
if [ -f "$name" ]; then
	echo "VS sounding is up to date"
else
	echo $name	
	python LIFTsounding.py /home/jfischer/KITCUBEmount/sftp-upload/RASO/VS/$fileNameVS Villingen
	echo “Soundings attached” | mailx –s “Sounding delivery” –a ./$fileNameVS cuckoo_janni@live.de	
fi
