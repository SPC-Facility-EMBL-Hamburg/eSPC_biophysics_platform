#!/bin/bash

# Directory to iterate over (you can change this to your specific directory)
directory="."

# Loop through all files in the directory
for file in "$directory"/*; do
    # Check if the file contains the pattern 'helix' or 'HELX' (case-insensitive)
    if ! grep -qiE 'helix|HEL' "$file"; then
        # If neither pattern is found, print the filename
        echo "$file"
    fi
done
