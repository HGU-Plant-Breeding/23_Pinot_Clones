#!/bin/bash

# Script Name: process_te_divergence.sh
# Description: Parses EDTA GFF3 annotations to calculate TE divergence (age).
#              1. Generates a master BED file with divergence encoded in the name.
#              2. Splits TEs into "Young" and "Old" categories for specific genomic partitions.
#
# Logic: Divergence = 100 - (Identity_Fraction * 100)
#        (e.g., 0.98 identity -> 2% divergence)

set -euo pipefail

# --- Help Function ---
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -m, --main <file>       Main TE annotation GFF3 (from EDTA)"
    echo "  --hetero <file>         GFF3 of TEs in Heterozygous regions"
    echo "  --homo <file>           GFF3 of TEs in Homozygous regions"
    echo "  --hemi <file>           GFF3 of TEs in Hemizygous regions"
    echo "  --young <float>         Threshold for Young TEs (default: 2.0)"
    echo "  --old <float>           Threshold for Old TEs (default: 20.0)"
    echo "  -h, --help              Display this help"
    exit 1
}

# --- Defaults ---
MAIN_GFF=""
HETERO_GFF=""
HOMO_GFF=""
HEMI_GFF=""
YOUNG_THRESH=2
OLD_THRESH=20

# --- Parse Arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -m|--main) MAIN_GFF="$2"; shift ;;
        --hetero) HETERO_GFF="$2"; shift ;;
        --homo) HOMO_GFF="$2"; shift ;;
        --hemi) HEMI_GFF="$2"; shift ;;
        --young) YOUNG_THRESH="$2"; shift ;;
        --old) OLD_THRESH="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$MAIN_GFF" || -z "$HETERO_GFF" || -z "$HOMO_GFF" || -z "$HEMI_GFF" ]]; then
    echo "Error: All GFF input files must be specified."
    usage
fi

echo "--- Starting TE Divergence Analysis ---"
echo "Main GFF: $MAIN_GFF"
echo "Thresholds: Young < $YOUNG_THRESH% | Old > $OLD_THRESH%"

# --- Task 1: Master BED File ---
echo -e "\n[Task 1] Generating master BED file with divergence info..."

# Awk logic: Extract identity, calc divergence, print BED
# Name format: Chr__Start__End__Divergence
awk -F'\t' -v OFS='\t' '
$9 !~ /Parent=/ && $3 != "repeat_region" {
    n=split($9,a,";"); id_val=-1; ltr_id_val=-1;
    for(i=1;i<=n;i++){
        if(a[i]~/ltr_identity=/){sub(/ltr_identity=/,"",a[i]); ltr_id_val=a[i]}
        if(a[i]~/identity=/){sub(/identity=/,"",a[i]); id_val=a[i]}
    }
    div=-1;
    if(ltr_id_val!=-1){div=100-(ltr_id_val*100)} else if(id_val!=-1){div=100-(id_val*100)}
    
    if(div!=-1){
        unique_id = $1"__"$4"__"$5"__"div;
        print $1, $4, $5, unique_id, "0", $7;
    }
}' "$MAIN_GFF" | sort -k1,1 -k2,2n | uniq > all_tes_with_divergence.UNIQUE.bed

echo " -> Created: all_tes_with_divergence.UNIQUE.bed"

# --- Task 2: Partition Analysis ---
echo -e "\n[Task 2] Splitting Young/Old TEs by partition..."

# Function to process a partition file
process_partition() {
    local input_gff="$1"
    local partition_name="$2"
    
    echo " -> Processing $partition_name..."
    
    # Process Young TEs
    awk -F'\t' -v OFS='\t' -v thresh="$YOUNG_THRESH" '
    {
        n=split($9,a,";"); id=-1; ltr_id=-1; 
        for(i=1;i<=n;i++){
            if(a[i]~/ltr_identity=/){sub(/ltr_identity=/,"",a[i]); ltr_id=a[i]} 
            if(a[i]~/identity=/){sub(/identity=/,"",a[i]); id=a[i]}
        } 
        div=-1;
        if(ltr_id!=-1){div=100-(ltr_id*100)} else if(id!=-1){div=100-(id*100)}
        
        if(div!=-1 && div < thresh){ 
            print $1, $4, $5, "young_te_"NR, "0", $7 
        }
    }' "$input_gff" > "Diploid_TE_Young_${partition_name}.bed"

    # Process Old TEs
    awk -F'\t' -v OFS='\t' -v thresh="$OLD_THRESH" '
    {
        n=split($9,a,";"); id=-1; ltr_id=-1; 
        for(i=1;i<=n;i++){
            if(a[i]~/ltr_identity=/){sub(/ltr_identity=/,"",a[i]); ltr_id=a[i]} 
            if(a[i]~/identity=/){sub(/identity=/,"",a[i]); id=a[i]}
        } 
        div=-1;
        if(ltr_id!=-1){div=100-(ltr_id*100)} else if(id!=-1){div=100-(id*100)}
        
        if(div!=-1 && div > thresh){ 
            print $1, $4, $5, "old_te_"NR, "0", $7 
        }
    }' "$input_gff" > "Diploid_TE_Old_${partition_name}.bed"
}

# Run for each partition
process_partition "$HETERO_GFF" "heterozygous"
process_partition "$HOMO_GFF"   "homozygous"
process_partition "$HEMI_GFF"   "hemizygous"

echo -e "\n--- Processing Complete ---"
