#!/bin/bash

# Specify the input file
input_file="strains_metadata.csv"

# Specify the source directory where the files are located
source_directory="all_assemblies"

# Specify the destination directory for removed files
filtered_assemblies="filtered_assemblies"

# Initialize an empty list (array)
elements=()

# Read each line of the file
while IFS= read -r line; do
    # Use awk to process each line and extract the non-empty value from columns 4, 5, or 6 from the end
    value=$(echo "$line" | awk -F, '
    {
        # Calculate positions from the end
        col1 = NF - 3;
        col2 = NF - 4;
        col3 = NF - 5;

        # Remove leading and trailing quotes
        gsub(/^"|"$/, "", $(col1));
        gsub(/^"|"$/, "", $(col2));
        gsub(/^"|"$/, "", $(col3));

        # Return the first non-empty value from columns
        if ($(col1) != "") print $(col1);
        else if ($(col2) != "") print $(col2);
        else if ($(col3) != "") print $(col3);
    }')

    # If a value was found, construct the full path and add it to the list
    if [ -n "$value" ]; then
        elements+=("all_assemblies/$value.fna")
    fi
done < "$input_file"

echo ${elements[1]}
# Move all matched files to the removed directory
for filepath in "${elements[@]}"; do
    
    # Check if the file exists and move it
    if [ -f "$filepath" ]; then
        #cp "$filepath" "$filtered_assemblies"
        rm -f "$filepath" 
    fi
done

echo "Selected files have been moved to '$filtered_assemblies'."
