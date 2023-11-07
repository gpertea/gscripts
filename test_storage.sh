#!/bin/bash
SCRATCH_PATH_1="$HOME/myscratch_test"
SCRATCH_PATH_2="$HOME/dcs04_test"
SCRATCH_PATH_3="$HOME/dcs05_test"
## file size in megabytes
SMALL_FILE_SIZE=1
LARGE_FILE_SIZE=1024
SMALL_FILES_COUNT=1000
LARGE_FILES_COUNT=2

# Function to prepare and test write speed with random data
test_write_speed() {
 local path=$1
 local bs=$2
 local count=$3
 local size=$(($bs * $count))
 echo "Writing $count files of size $bs to $path"
 local start=$(date +%s.%N)
 for i in $(seq 1 $count); do
   dd if=/dev/urandom of=${path}/testfile_${i} bs=${bs}M count=1 iflag=fullblock oflag=sync status=none
 done
 local end=$(date +%s.%N)
 local duration=$(echo "$end - $start" | bc)
 local speed=$(echo "scale=2; $size / $duration" | bc)
 echo "Write speed for $path: $speed MB/s"
}

# Function to test read speed
test_read_speed() {
  local path=$1
  local bs=$2
  local count=$3
  local size=$((bs*count))
  echo "Reading $count files of size $bs (Total: $size MB) from $path"
  local start=$(date +%s.%N)
  for i in $(seq 1 $count); do
    dd if=${path}/testfile_${i} of=/dev/null bs=${bs}M count=1 iflag=dsync 2> /dev/null
  done
  local end=$(date +%s.%N)
  local duration=$(echo "$end - $start" | bc)
  local speed=$(echo "scale=2; $size / $duration" | bc)
  echo "Speed for reading from $path: $speed MB/s"
}

# Testing small files
test_write_speed $SCRATCH_PATH_1 $SMALL_FILE_SIZE $SMALL_FILES_COUNT
test_write_speed $SCRATCH_PATH_2 $SMALL_FILE_SIZE $SMALL_FILES_COUNT
test_write_speed $SCRATCH_PATH_2 $SMALL_FILE_SIZE $SMALL_FILES_COUNT

test_read_speed $SCRATCH_PATH_1 $SMALL_FILE_SIZE $SMALL_FILES_COUNT
test_read_speed $SCRATCH_PATH_2 $SMALL_FILE_SIZE $SMALL_FILES_COUNT
test_read_speed $SCRATCH_PATH_3 $SMALL_FILE_SIZE $SMALL_FILES_COUNT

# Testing large files
test_write_speed $SCRATCH_PATH_1 $LARGE_FILE_SIZE $LARGE_FILES_COUNT
test_write_speed $SCRATCH_PATH_2 $LARGE_FILE_SIZE $LARGE_FILES_COUNT
test_write_speed $SCRATCH_PATH_3 $LARGE_FILE_SIZE $LARGE_FILES_COUNT

test_read_speed $SCRATCH_PATH_1 $LARGE_FILE_SIZE $LARGE_FILES_COUNT
test_read_speed $SCRATCH_PATH_2 $LARGE_FILE_SIZE $LARGE_FILES_COUNT
test_read_speed $SCRATCH_PATH_3 $LARGE_FILE_SIZE $LARGE_FILES_COUNT

# Cleanup
rm -f ${SCRATCH_PATH_1}/testfile_*
rm -f ${SCRATCH_PATH_2}/testfile_*
rm -f ${SCRATCH_PATH_3}/testfile_*
