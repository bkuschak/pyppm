#!/bin/bash
# Raw data files take up a considerable amount of disk space. Run this script
# to purge old files from the disk.

# Delete files older than this:
DAYS=30
DATADIR=/mnt/ssd/pyppm

TMPDIR=/tmp/cleanup
rm -rf ${TMPDIR}
mkdir -p ${TMPDIR}

# Entire list of files.
echo "Finding files older than ${DAYS} days..."
sudo find ${DATADIR} -type f -mtime +${DAYS} > ${TMPDIR}/list
echo "Found $(wc -l ${TMPDIR}/list) files"

# Break into 50 file chunks.
cd ${TMPDIR}
split -l 50 list split_

# Delete the files.
for i in $(ls split_*); do 
    echo $i
    files=$(cat $i|xargs)
    cd ${DATADIR}
    rm -f $files
    cd -
    rm $i
done

echo "Results:"
df -h ${DATADIR}
sync
