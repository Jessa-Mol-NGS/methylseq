#!/usr/bin/env python3

import os
import platform
import pysam
import pandas as pd

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def read_overlap(forward_read: pysam.AlignedSegment, reverse_read: pysam.AlignedSegment) -> tuple[bool, int]:
    """
    Check if two reads overlap on the reference genome.

    Args:
        forward_read: AlignedSegment object containing the forward read of a read pair
        reverse_read: AlignedSegment object containing the reverse read of a read pair

    Returns:
        Tuple containing a boolean that indicates whether or not the reads overlap, and
        an integer that represents the number of bases that overlap.    
    """
    
    # Case 1: The fragment is shorter than the individual reads
    if forward_read.reference_start > reverse_read.reference_start:
        overlap_length = reverse_read.reference_end - forward_read.reference_start
    
    # Case 2: The fragment is longer than the individual reads
    else:
        overlap_length = forward_read.reference_end - reverse_read.reference_start

    return (overlap_length > 0, max(overlap_length, 0))

def combine_overlapping_read_sequences(forward_read: pysam.AlignedSegment, reverse_read: pysam.AlignedSegment, overlap_length: int, use_reference: bool = False) -> str:
    """
    Reconstructs the complete fragment from paired-end reads that overlap. Overlapping bases are taken from read 1, given that
    they typically have higher quality than read 2.

    Args:
        forward_read: AlignedSegment object containing the forward read of a read pair
        reverse_read: AlignedSegment object containing the reverse read of a read pair
        overlap_length: Integer equal to the number of bases that overlap in the reads
        use_reference: Boolean indicating whether or not to use the reference sequence to which the reads map
    
    Returns:
        String containing the sequence of the original fragment
    """
    if not use_reference:

        # Keep all bases from the forward read
        if forward_read.is_read1:
            fragment_seq = forward_read.query_sequence + reverse_read.query_sequence[overlap_length:]

        # Keep all bases from the reverse read
        else:
            fragment_seq = forward_read.query_sequence[:-overlap_length] + reverse_read.query_sequence
    
    # Reconstruct the fragment from the reference genome
    else:
         fragment_seq = forward_read.get_reference_sequence() + reverse_read.get_reference_sequence()[overlap_length:]
    
    return fragment_seq

def extract_fragment_data(bam_file_path: str, include_quality_scores: bool, include_reference_sequence: bool, reconstruct_inserts: bool) -> pd.DataFrame:
    """
    Extract insert sizes and other fragment information from a name-sorted BAM file.
    
    Args:
        bam_file_path: Path to the BAM file (must be sorted by read name)
        include_quality_scores: Boolean indicating whether base quality scores should be included in the output
        include_reference_sequence: Boolean indicating whether reference sequences should be included in the output
        reconstruct_inserts: Boolean indicating whether inserts should be reconstructed from overlapping reads
   
    Returns:
        DataFrame with insert sizes indexed by read name
    """
    # Valid flag combinations for paired-end reads
    VALID_FLAGS = {(99, 147), (147, 99), (83, 163), (163, 83)}
    
    # Initialize lists to store fragment information
    read_names, chrom, pos = [], [], []
    insert_sizes, overlap_bool, forward_R1 = [], [], []
    forward_read_seqs, reverse_read_seqs = [], []

    if include_quality_scores:
        forward_qual, reverse_qual = [], []
    if include_reference_sequence:
        forward_ref_seqs, reverse_ref_seqs = [], []
    if reconstruct_inserts:
        fragment_read_seqs = []
        if include_reference_sequence:
            fragment_ref_seqs = []
    
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        it = iter(bam_file)
        prev = None
        try:
            while True:
                # Get first non-secondary, non-supplementary read
                while True:
                    if prev is None:
                        first_read = next(it)
                    else:
                        first_read = prev
                        prev = None
                    if first_read.is_secondary or first_read.is_supplementary:
                        # skip these and continue looking for a valid first_read
                        continue
                    break

                # Get second non-secondary, non-supplementary read
                while True:
                    second_read = next(it)
                    if second_read.is_secondary or second_read.is_supplementary:
                        continue
                    break
                
                # Verify that the reads have valid flag combinations
                if (first_read.flag, second_read.flag) not in VALID_FLAGS:  
                    print(f"Warning: Unexpected flags: {first_read.flag}, {second_read.flag}. Skipping this fragment.")
                    continue

                # Verify that the reads are a pair
                if first_read.query_name != second_read.query_name:
                    print(f"Warning: singleton or out-of-order read: {first_read.query_name}")
                    prev = second_read
                    continue

                # Identify forward and reverse reads
                forward_read = first_read if first_read.is_forward else second_read
                reverse_read = second_read if first_read.is_forward else first_read

                # Define boolean indicating whether the forward read is R1
                forward_R1.append(forward_read.is_read1)

                # Compute insert size
                insert_size = reverse_read.reference_end - forward_read.reference_start
                insert_sizes.append(insert_size)

                # Check for overlapping reads
                overlap = read_overlap(forward_read, reverse_read)
                overlap_bool.append(overlap[0])

                # Get fragment start position
                chrom.append(forward_read.reference_name)
                pos.append(forward_read.reference_start)

                # Extract read sequences
                forward_read_seqs.append(forward_read.query_sequence)
                reverse_read_seqs.append(reverse_read.query_sequence)

                # Extract reference sequences
                if include_reference_sequence:
                    forward_ref_seqs.append(forward_read.get_reference_sequence())
                    reverse_ref_seqs.append(reverse_read.get_reference_sequence())

                # Extract base quality scores
                if include_quality_scores:
                    forward_qual.append(forward_read.query_qualities_str if hasattr(forward_read, "query_qualities_str") else "")
                    reverse_qual.append(reverse_read.query_qualities_str if hasattr(reverse_read, "query_qualities_str") else "")

                # Combine reads into fragment if they overlap
                if reconstruct_inserts:
                    if overlap[0]:
                        fragment_read_seq = combine_overlapping_read_sequences(forward_read, reverse_read, overlap[1])
                        if include_reference_sequence:
                            fragment_ref_seq = combine_overlapping_read_sequences(forward_read, reverse_read, overlap[1], True)
                    else:
                        fragment_read_seq = ""
                        if include_reference_sequence:
                            fragment_ref_seq = ""
                    fragment_read_seqs.append(fragment_read_seq)
                    if include_reference_sequence:
                        fragment_ref_seqs.append(fragment_ref_seq)

                # Store read name information
                read_names.append(first_read.query_name)

        except StopIteration:
            # if prev holds a leftover unpaired read, warn about singleton at EOF
            if prev is not None:
                print(f"Warning: singleton read at EOF: {prev.query_name}")
            pass
    
    # Construct and return DataFrame
    data = {
        "chrom": chrom, "pos": pos, "insert_size": insert_sizes,
        "reads_overlap": overlap_bool, "forward_read_R1": forward_R1,
        "forward_read": forward_read_seqs, "reverse_read": reverse_read_seqs
    }
    result = pd.DataFrame(
        {
            "chrom": chrom, "pos": pos, "insert_size": insert_sizes,
            "reads_overlap": overlap_bool, "forward_read_R1": forward_R1,
            "forward_read": forward_read_seqs, "reverse_read": reverse_read_seqs
        }, index = read_names)
    
    if include_reference_sequence:
        data.update({"forward_ref": forward_ref_seqs, "reverse_ref": reverse_ref_seqs})    
    if include_quality_scores:
        data.update({"forward_qual": forward_qual, "reverse_qual": reverse_qual})
    if reconstruct_inserts:
        if include_reference_sequence:
            data.update({"fragment_read_sequence": fragment_read_seqs, "fragment_ref_sequence": fragment_ref_seqs})
        else:
            data.update({"fragment_read_sequence": fragment_read_seqs})
    
    result = pd.DataFrame(data, index = read_names)
    return result

if __name__ == "__main__":
    
    # Get parameter values
    include_quality_scores = True if "$params.include_quality_scores" == "true" else False
    include_reference_sequence = True if "$params.include_reference_sequence" == "true" else False
    reconstruct_inserts = True if "$params.reconstruct_inserts" == "true" else False
    sample_id = "$meta.id"
   
    # Extract insert sizes and other relevant data
    results_df = extract_fragment_data("$bam", include_quality_scores, include_reference_sequence, reconstruct_inserts)
    results_df.to_csv(f"{sample_id}.insert.sizes.csv.gz", sep = "\t", index_label = "qname", compression = "gzip")

    # Write the versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
