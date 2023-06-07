#!/bin/bash

# file_tuples=(
#   ("total_cell_mass_result.txt" "results\total_cell_mass_result.txt")
#   ("total_pressure_result.txt" "results\total_pressure_result.txt")
#   ("total_sie_result.txt" "results\total_sie_result.txt")
#   ("velocity_result.txt" "results\velocity_result.txt")
# )

compare_files() {
  local file1="$1"
  local file2="$2"
  
  if cmp -s "$file1" "$file2"; then
    echo "Files '$file1' and '$file2' have the same content."
    return 0 
  else
    echo "Files '$file1' and '$file2' have different content."
    return 1 
  fi
}

# for ((i = 0; i < ${#file_tuples[@]}; i++)); do
#   file1="${file_tuples[i][0]}"
#   file2="${file_tuples[i][1]}"
  
#   compare_files "$file1" "$file2"
# done

compare_files "total_cell_mass_result.txt" "results\total_cell_mass_result.txt"
compare_files "total_pressure_result.txt" "results\total_pressure_result.txt"
compare_files "total_sie_result.txt" "results\total_sie_result.txt"
compare_files "velocity_result.txt" "results\velocity_result.txt"