#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <text_file> <num_tests>"
    exit 1
fi

text_file=$1
num_tests=$2
index_file="${text_file}.lidx"
if [[ -z $num_tests ]]; then
num_tests=20
fi

# Create the index file if it doesn't exist
if [ ! -f "$index_file" ]; then
    linix "$text_file"
fi

# Get the total number of lines in the text file
total_lines=$(wc -l < "$text_file")

# Run the specified number of tests
for ((i = 0; i < num_tests; i++)); do
    # Generate a random line number
    line_num=$(( RANDOM % total_lines + 1 ))

    # Get the line using the linix script
    #echo ".. checking line $line_num"
    linix_line=$(linix "$text_file" "$line_num")

    # Get the line using sed
    sed_line=$(sed -n "${line_num}p" "$text_file")

    # Compare the outputs
    if [ "$linix_line" == "$sed_line" ]; then
        echo "Test $((i + 1)): Pass"
    else
        echo "Test $((i + 1)): FAIL (line $line_num)"
        echo "linix: $linix_line"
        echo "sed:   $sed_line"
    fi
done
