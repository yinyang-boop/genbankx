import argparse
import sys
from .parser import parse_genbank
from .translation import translate, reverse_complement
from .io import write_fasta

def extract_cds_sequence(genome_seq: str, cds) -> str:
    """
    Extract CDS sequence from genome and handle strand orientation.
    
    Args:
        genome_seq: Complete genome sequence
        cds: CDSFeature object
        
    Returns:
        Oriented CDS sequence
    """
    fragments = []
    for interval in cds:
        start_idx = interval.start - 1
        end_idx = interval.end
        fragment = genome_seq[start_idx:end_idx]
        fragments.append(fragment)
    
    sequence = "".join(fragments)
    
    if cds.strand.value == -1:
        sequence = reverse_complement(sequence)
        
    return sequence

def main(args=None):
    """Command line interface for genbankx."""
    parser = argparse.ArgumentParser(
        description="Extract CDS sequences from GenBank files"
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("-o", "--output", help="Output FASTA file")
    parser.add_argument("-p", "--proteins", action="store_true",
                       help="Output protein sequences instead of DNA")
    parser.add_argument("--prefix", default="", help="Prefix for sequence headers")
    
    if args is None:
        args = parser.parse_args()
    
    try:
        with open(args.input, 'r') as f:
            records = list(parse_genbank(f))
    except Exception as e:
        print(f"Error reading GenBank file: {e}", file=sys.stderr)
        return 1
    
    if not records:
        print(f"No records found in {args.input}", file=sys.stderr)
        return 1
    
    output_records = {}
    
    for record_idx, record in enumerate(records):
        for cds_idx, cds in enumerate(record.cds, 1):
            try:
                # Extract DNA sequence
                dna_seq = extract_cds_sequence(record.sequence, cds)
                
                # Create header
                gene_name = cds.gene or f"CDS{cds_idx}"
                header_parts = [args.prefix + gene_name, record.locus]
                if cds.product:
                    header_parts.append(cds.product.replace(" ", "_"))
                
                header = "|".join(header_parts)
                
                if args.proteins:
                    # Translate to protein
                    protein_seq = translate(dna_seq, to_stop=True)
                    output_records[header] = protein_seq
                else:
                    output_records[header] = dna_seq
                    
            except Exception as e:
                print(f"Warning: Failed to process CDS {cds_idx} in record {record_idx}: {e}",
                      file=sys.stderr)
                continue
    
    # Write output
    if args.output:
        write_fasta(args.output, output_records)
        print(f"Successfully wrote {len(output_records)} sequences to {args.output}")
    else:
        # Print to stdout
        for header, seq in output_records.items():
            print(f">{header}")
            for i in range(0, len(seq), 60):
                print(seq[i:i+60])
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
