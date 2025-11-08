from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from .interval import Interval, IntervalList, Strand

@dataclass
class Feature:
    """Base feature class for GenBank features."""
    key: str
    location: IntervalList
    qualifiers: Dict[str, List[str]] = field(default_factory=dict)

    def get(self, name: str, default: Optional[str] = None) -> Optional[str]:
        """Get first value of a qualifier."""
        vals = self.qualifiers.get(name)
        return vals[0] if vals else default

class CDSFeature(IntervalList):
    """CDS feature with intervals, strand, and qualifiers."""
    
    def __init__(self, intervals: List[Interval], strand: Strand,
                 qualifiers: Optional[Dict[str, List[str]]] = None):
        super().__init__(intervals, strand=strand)
        self.qualifiers = qualifiers or {}
    
    @property
    def gene(self) -> Optional[str]:
        """Get gene name from qualifiers."""
        return self.get_qualifier("gene")
    
    @property 
    def product(self) -> Optional[str]:
        """Get product description from qualifiers."""
        return self.get_qualifier("product")
    
    def get_qualifier(self, name: str) -> Optional[str]:
        """Get first value of a qualifier."""
        vals = self.qualifiers.get(name)
        return vals[0] if vals else None
        
    @property
    def translation(self) -> Optional[str]:
        """Get translation from qualifiers."""
        return self.get_qualifier("translation")
