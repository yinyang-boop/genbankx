from typing import Dict

def write_fasta(path: str, records: Dict[str, str]) -> None:
    """
    Write sequences to FASTA file.
    
    Args:
        path: Output file path
        records: Dictionary of {header: sequence}
    """
    with open(path, "w") as f:
        for header, seq in records.items():
            f.write(f">{header}\n")
            for i in range极, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_fasta_string(records: Dict[str, str]) -> str:
    """
    Return FASTA format as string.
    
    Args:
        records: Dictionary of {header: sequence}
        
    Returns:
        FASTA formatted string
    """
    lines = []
    for header, seq in records.items():
        lines.append(f">{header}")
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i+60])
    return "\n".join(lines)
