#!/bin/bash

# Output file to store the combined 7th columns
output_file="combined_columns.txt"

# Create a temporary directory to store individual columns
temp_dir=$(mktemp -d)

for i in {1..282}; do
  file="counts_$i"
  # Check if the file exists and is not empty
  if [ -s "$file" ]; then
    # Use awk to extract the 7th column, using a tab as the field separator and ignoring lines starting with "#"
    column=$(awk -F'\t' '!/^#/ {print $7}' "$file")

    # Save the extracted column in a temporary file
    echo -e "$column" > "$temp_dir/column_$i.txt"
  fi
done

# Use paste to combine the columns horizontally and save the result in the output file
paste "$temp_dir"/* > "$output_file"

# Remove the temporary directory
rm -r "$temp_dir"

echo "Combined 7th columns from counts_1 to counts_282 into $output_file"
