import pandas as pd
import pysam
import os

############ INPUTS ############
## paths are resolved from the repository root so this runs from any directory
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mutations_file = os.path.join(repo_root, "results", "close_muts_181.tsv") ## filtered maf from R
bam_dir = os.path.join(repo_root, "inputs", "BAMs")   ## directory containing BAM files
output_file = os.path.join(repo_root, "results", "phasing_181_indels.csv")
################################

# Load mutation table
mutations = pd.read_csv(mutations_file, sep= "\t")
print(mutations)


# Initialize results
all_results = []
check_pos = 7675069 ## germline 181 mutation position

# Loop over rows in mutation file
for idx, row in mutations.iterrows():
    var_type = row["Variant_Type"]
    bam_folder = f'{row["Tumor_Sample_Barcode"]}'
    bam_name = f'{row["Tumor_Sample_Barcode"]}.consensus.bam'
    chrom = row["Chromosome"]
    mut_pos = int(row["Start_Position"])
    #mut_pos_end = int(row["End_Position"])
    ref_allele = row["Reference_Allele"]
    mut_base = row["Tumor_Seq_Allele2"]
    ref_len = len(ref_allele) if ref_allele != "-" else 0
    mut_len = len(mut_base) if mut_base != "-" else 0

    # mutation type (distinguish between indel types)
    if ref_len > 0 and mut_len == 0:
        var_type = "DEL"
    elif ref_len == 0 and mut_len > 0:
        var_type = "INS"
    elif ref_len > 0 and mut_len > 0 and ref_len != mut_len:
        var_type = "DELINS"
    else:
        var_type = "SNV" if ref_len == 1 and mut_len == 1 else "ONP"


    ## does .bam exist?
    bam_path = os.path.join(bam_dir, f"{bam_folder}/{bam_name}")
    if not os.path.exists(bam_path):
        print(f"WARNING: BAM file not found: {bam_path}")
        continue
    
    ## does .bai exist and is samtools happy with it being up to date?
    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path) or os.path.getmtime(bai_path) < os.path.getmtime(bam_path):
        print(f"Rebuilding index for {bam_path}...")
        pysam.index(bam_path)   # rebuild .bai

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    results = []

    for read in bamfile.fetch(chrom, mut_pos - 1, mut_pos):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        has_mut = False
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

        ############ Handle deletions ############
        if var_type == "DEL":
            deleted_positions = range(mut_pos - 1, mut_pos - 1 + ref_len)  # shift to 0-based

            # Show all aligned pairs, highlight deletions
            gaps = []
            for query_pos, ref_pos in aligned_pairs:
                if query_pos is None and ref_pos is not None:
                    gaps.append(ref_pos)

            # Check if expected deleted positions are all in gaps
            deletion_detected = all(pos in gaps for pos in deleted_positions)

            if deletion_detected:
                has_mut = True

        ############ Handle insertions ############
        elif var_type == "INS":
            # Collect query-only positions (inserted bases)
            insertion_bases = []
            insertion_pos = None

            for query_pos, ref_pos in aligned_pairs:
                if ref_pos == mut_pos - 1:
                    # Start recording after this position
                    insertion_pos = True
                elif insertion_pos and ref_pos is None and query_pos is not None:
                    insertion_bases.append(read.query_sequence[query_pos])
                elif insertion_pos and ref_pos is not None:
                    break

            # Check if inserted bases match expected mutation sequence
            if "".join(insertion_bases) == mut_base:
                has_mut = True

        ############ Handle complex indels ############
        elif var_type == "DELINS":

            # check deletion
            deleted_positions = range(mut_pos - 1, mut_pos - 1 + ref_len)
            gaps = [ref_pos for query_pos, ref_pos in aligned_pairs if query_pos is None and ref_pos is not None]

            # check for overlap
            overlap = len(set(deleted_positions) & set(gaps))

            # may need more generalizability, but few mutations in dataset permits manual validation 
            deletion_detected = overlap > 0 
            if deletion_detected:
                has_mut = True


        ############ Handle SNV and ONVs ############
        else:
            has_mut = True
            for i in range(mut_len):
                target_pos = mut_pos - 1 + i
                expected_base = mut_base[i]
                found_match = False
                for qpos, rpos in aligned_pairs:
                    if rpos == target_pos and qpos is not None:
                        base = read.query_sequence[qpos]
                        if base == expected_base:
                            found_match = True
                        break
                if not found_match:
                    has_mut = False
                    break

        if not has_mut:
            continue

        ############ Check if reads cover germline 181 ############
        start_1based = read.reference_start + 1
        end_1based = read.reference_end
        covers_check = False
        base_at_check = None

        if var_type == "DELINS":
            print("deletion checking for 181\n")

        if read.reference_start <= check_pos - 1 < read.reference_end:
            for query_pos, ref_pos in aligned_pairs:
                if ref_pos == check_pos and query_pos is not None:
                    base_at_check = read.query_sequence[query_pos]
                    covers_check = True
                    break

        results.append({
            "bam_file": bam_name,
            "chrom": chrom,
            "mut_pos": mut_pos,
            "ref_base": ref_allele,
            "mut_base": mut_base,
            "check_pos": check_pos,
            "covers_check": covers_check,
            "base_at_check": base_at_check,
            "read_name": read.query_name,
            "read_start": start_1based,
            "read_end": end_1based
        })

    bamfile.close()
    all_results.extend(results)

# Convert to dataframe and save
if not all_results:
    raise SystemExit(
        f"No phasing results generated - no BAM files were read from {bam_dir}\n"
        "Consensus BAMs are not included in this repository; place them at "
        "inputs/BAMs/<sample>/<sample>.consensus.bam."
    )

df_out = pd.DataFrame(all_results)
df_out.to_csv(output_file, index=False)
print(f"Results written to {output_file}")