#!/bin/bash

# Loop over all directories starting with "export"
for dir in $(find . -type d -name "export*"); do
    # Check if .tsv file exists in the directory
    if [ -e "$dir/*.tsv" ]; then
        # Go inside the directory
        cd "$dir" || exit
        # Remove the first line from each .tsv file and remove "#" character from the new header
        for file in *.tsv; do
            awk 'NR > 1 {sub(/^#/, ""); print}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
        done
        # Return to the original directory
        cd - || exit
    fi
done
