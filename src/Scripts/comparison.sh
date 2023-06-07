#!/bin/bash

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



compare_files "total_cell_mass_result.txt" "results/total_cell_mass_result.txt"
compare_files "total_pressure_result.txt" "results/total_pressure_result.txt"
compare_files "total_sie_result.txt" "results/total_sie_result.txt"
compare_files "velocity_result.txt" "results/velocity_result.txt"
