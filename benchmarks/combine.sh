#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_csv1> <input_csv2> <output_csv>"
  exit 1
fi

input_csv1="$1"
input_csv2="$2"
output_csv="$3"

# Function to find the index of the 'querytime' column
find_querytime_index() {
  local header="$1"
  IFS=',' read -r -a columns <<< "$header"
  for i in "${!columns[@]}"; do
    if [ "${columns[$i]}" == "querytime" ]; then
      echo "$i"
      return
    fi
  done
  echo "-1"
}

# Extract the headers and find the column index for 'querytime'
header1=$(head -n 1 "$input_csv1")
header2=$(head -n 1 "$input_csv2")

querytime_index1=$(find_querytime_index "$header1")
querytime_index2=$(find_querytime_index "$header2")

if [ "$querytime_index1" -eq -1 ]; then
  echo "Error: 'querytime' column not found in $input_csv1"
  exit 1
fi

if [ "$querytime_index2" -eq -1 ]; then
  echo "Error: 'querytime' column not found in $input_csv2"
  exit 1
fi

# Create the output CSV file and write the header
echo "file,ourtime,pibiritime" > "$output_csv"

# Function to extract the required columns from a CSV file
extract_columns() {
  local input_csv="$1"
  local querytime_index="$2"
  tail -n +2 "$input_csv" | while IFS=',' read -r -a row; do
    echo "${row[0]},${row[$querytime_index]}"
  done
}

# Extract columns from both input CSV files and join them
paste -d ',' <(extract_columns "$input_csv1" "$querytime_index1") <(extract_columns "$input_csv2" "$querytime_index2" | cut -d ',' -f2) >> "$output_csv"

echo "Combined CSV file created: $output_csv"