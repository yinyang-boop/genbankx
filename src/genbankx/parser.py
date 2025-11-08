from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
import re
import sys
from .interval import Interval, IntervalList, Strand
from .features import CDSFeature, Feature

@dataclass
class Record:
    """GenBank record containing sequence and features."""
    locus: str
    length: int
    definition: str
    accession: Optional[str]
    sequence: str
    features: List[Feature]
    cds: List[CDSFeature]

def _parse_location(loc_str: str) -> Tuple[List[Interval], Strand]:
    """
    Parse GenBank location string.
    
    Supports:
    - Simple: 1..100
    - Complement: complement(1..100)
    - Join: join(1..100,200..300)
    - Complement join: complement(join(1..100,200..300))
    """
    strand = Strand.PLUS
    loc_str = loc_str.strip()
    
    # Handle complement
    if loc_str.startswith("complement(") and loc_str.endswith(")"):
        strand = Strand.MINUS
        loc_str = loc_str[11:-1]  # Remove complement(...)
    
    # Handle join
    if loc_str.startswith("join(") and loc_str.endswith(")"):
        loc_str = loc_str[5:-1]
    
    intervals = []
    
    # Split by commas and parse each interval
    for part in loc_str.split(','):
        part = part.strip()
        if not part:
            continue
            
        # Remove partial markers < >
        part = part.replace('<', '').replace('>', '')
        
        if '..' in part:
            try:
                start_str, end_str = part.split('..', 1)
                start = int(start_str.strip())
                end = int(end_str.strip())
                intervals.append(Interval(start, end))
            except ValueError:
                continue
        elif part.isdigit():
            # Single base position
            pos = int(part)
            intervals.append(Interval(pos, pos))
    
    return intervals, strand

def _flush_feature(record: dict, key: Optional[str], location_lines: List[str], qualifiers: Dict[str, List[str]]):
    """Convert parsed feature data into Feature/CDSFeature object and add to record."""
    if key:
        full_location_str = "".join(location_lines).replace(" ", "")
        try:
            intervals, strand = _parse_location(full_location_str)
            if intervals:
                # Create base Feature object
                feature = Feature(key=key, location=IntervalList(intervals, strand=strand), qualifiers=qualifiers)
                record.setdefault("features", []).append(feature)

                # If it's a CDS, also create and store as CDSFeature
                if key == "CDS":
                     cds_feature = CDSFeature(intervals=intervals, strand=strand, qualifiers=qualifiers)
                     record.setdefault("cds", []).append(cds_feature)

        except Exception as e:
            print(f"Warning: Failed to parse feature '{key}' location '{full_location_str}': {e}", file=sys.stderr)

def _create_record(record_data: dict, sequence_lines: list) -> Record:
    """Create Record object from parsed data."""
    sequence = "".join(sequence_lines).upper()
    return Record(
        locus=record_data.get("locus", ""),
        length=record_data.get("length", 0),
        definition=record_data.get("definition", ""),
        accession=record_data.get("accession"),
        sequence=sequence,
        features=record_data.get("features", []),
        cds=record_data.get("cds", [])
    )

def parse_genbank(handle: Iterable[str]) -> Iterator[Record]:
    """
    Parse GenBank file and yield Record objects.
    
    Args:
        handle: Iterable of lines from GenBank file
        
    Yields:
        Record objects containing sequence and features
    """
    current_record = {}
    in_features = False
    in_origin = False
    current_feature_key = None
    current_feature_location_lines = []
    current_feature_qualifiers = {}
    sequence_lines = []
    
    for line in handle:
        line = line.rstrip('\n\r')
        
        # LOCUS line
        if line.startswith("LOCUS"):
            if current_record and 'locus' in current_record:
                _flush_feature(current_record, current_feature_key, current_feature_location_lines, current_feature_qualifiers)
                yield _create_record(current_record, sequence_lines)
                current_record = {}
                sequence_lines = []
                current_feature_key = None
                current_feature_location_lines = []
                current_feature_qualifiers = {}
                
            parts = line.split()
            if len(parts) >= 2:
                current_record["locus"] = parts[1]
                if len(parts) >= 3 and parts[2].isdigit():
                    current_record["length"] = int(parts[2])
                current_record.setdefault("features", [])
                current_record.setdefault("cds", [])
        
        # DEFINITION line
        elif line.startswith("DEFINITION"):
            current_record["definition"] = line[12:].strip()
        
        # ACCESSION line
        elif line.startswith("ACCESSION"):
            accession_part = line[10:].strip()
            current_record["accession"] = accession_part.split()[0] if accession_part else None
        
        # FEATURES section start
        elif line.startswith("FEATURES"):
            in_features = True
            in_origin = False
            current_feature_key = None
            current_feature_location_lines = []
            current_feature_qualifiers = {}
        
        # ORIGIN section start
        elif line.startswith("ORIGIN"):
            _flush_feature(current_record, current_feature_key, current_feature_location_lines, current_feature_qualifiers)
            in_features = False
            in_origin = True
            current_feature_key = None
            current_feature_location_lines = []
            current_feature_qualifiers = {}
        
        # End of record
        elif line.startswith("//"):
            _flush_feature(current_record, current_feature_key, current_feature_location_lines, current_feature_qualifiers)
            if current_record and 'locus' in current_record:
                yield _create_record(current_record, sequence_lines)
                current_record = {}
                sequence_lines = []
                current_feature_key = None
                current_feature_location_lines = []
                current_fe极_qualifiers = {}
            
            in_features = False
            in_origin = False
        
        # Parse features and locations
        elif in_features:
            if line.startswith(" " * 5) and not line.startswith(" " * 21):
                # New feature line
                _flush_feature(current_record, current极ature_key, current_feature_location_lines, current_feature_qualifiers)
                
                current_feature_key = line[5:15].strip()
                current_feature_location_lines = [line[21:].strip()]
                current_feature_qualifiers = {}
            elif line.startswith(" " * 21) and current_feature_key:
                # Continuation of location or qualifier
                trimmed_line = line[21:].strip()
                if trimmed_line.startswith("/"):
                     # Qualifier line
                     qual_line = trimmed_line[1:]
                     if "=" in qual_line:
                        key, value = qual_line.split("=", 1)
                        value = value.strip('"')
                        current_feature_qualifiers.setdefault(key, []).append(value)
                     else:
                        # Boolean qualifiers
                        current_feature_qualifiers.setdefault(qual_line, []).append("")
                else:
                    # Continuation of location string
                    current_feature_location_lines.append(trimmed_line)
        
        # Parse sequence
        elif in_origin and line and not line.startswith("//"):
            if any(c in 'acgtnACGTN' for c in line):
                # Extract sequence part
                seq_part = "".join(re.findall('[a-zA-Z]', line))
                sequence_lines.append(seq_part)
    
    # Yield final record if any
    _flush_feature(current_record, current_feature_key, current_feature_location_lines, current_feature_qualifiers)
    if current_record and 'locus' in current_record:
        yield _create_record(current_record, sequence_lines)
