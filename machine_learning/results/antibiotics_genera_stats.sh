output="combined_genera_summaryClasses.csv"
echo "antibiotic,<your_header_line>" > "$output"

for dir in /mnt/raid1b/kdan_data/Paper/Machine_Learning/results_new/*; do
    antibiotic=$(basename "$dir")
    file=$(find "$dir/models" -maxdepth 1 -name '*genera_summary.csv' 2>/dev/null)

    if [ -f "$file" ]; then
       tail -n +2 "$file" | sed "s/^/$antibiotic,/" >> "$output"
    fi
done

