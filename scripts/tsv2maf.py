#!/usr/bin/env python3

import argparse
import csv
import sys


options = argparse.ArgumentParser()
options.add_argument("input", type=argparse.FileType("r"), default="-")
options.add_argument("--output", "-o", type=argparse.FileType("w"), default="-")


csv.register_dialect("tsv",
                     delimiter="\t",
                     lineterminator="\n")
csv.field_size_limit(sys.maxsize)


def norm(pos, ref, alt):
    ref_length = len(ref)
    alt_length = len(alt)
    
    # Remove any prefixed reference bps from all alleles, using "-" for simple indels
    while ref and alt and ref[0] == alt[0] and ref != alt:
        ref = ref[1:] or "-"
        alt = alt[1:] or "-"
        pos += 1
        ref_length -= 1
        alt_length -= 1

    if ref_length == alt_length:
        var_type = "ONP" if ref_length > 3 else ["SNP", "DNP", "TNP"][ref_length - 1]
        return pos, pos + alt_length - 1, ref, alt, var_type
    elif ref_length < alt_length:
        if ref == "-":
            return pos - 1, pos, ref, alt, "INS"
        else:
            return pos, pos + ref_length - 1, ref, alt, "INS"
    else:
        return pos, pos + ref_length - 1, ref, alt, "DEL"


def main(args):
    output_field_names = [
        "Hugo_Symbol",
        "NCBI_Build",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Strand",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode"
    ]
    
    reader = csv.DictReader(args.input, dialect="tsv")
    writer = csv.DictWriter(args.output, output_field_names, dialect="tsv")
    writer.writeheader()
    
    for variant in reader:
        start_position, end_position, ref, alt, var_type = norm(
            float(variant["start_position"]), variant["ref"], variant["alt"])
        writer.writerow({
            "Hugo_Symbol": (
                variant["gene_symbol"] if variant["gene_symbol"] != ""
                else "Unknown"),
            "NCBI_Build": "GRCh37",
            "Chromosome": variant["chrom"],
            "Start_Position": start_position,
            "End_Position": end_position,
            "Strand": "+",
            "Variant_Classification": (
                variant["variant_classification"] if variant["variant_classification"] != ""
                else "IGR"),
            "Variant_Type": var_type,
            "Reference_Allele": ref,
            "Tumor_Seq_Allele2": alt,
            "Tumor_Sample_Barcode": variant["TCGA_ID"]
        })


if __name__ == "__main__":
    main(options.parse_args())
