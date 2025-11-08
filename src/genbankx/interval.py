from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Iterable, List

class Strand(Enum):
    PLUS = 1
    MINUS = -1
    UNKNOWN = 0

@dataclass(frozen=True, order=True)
class Interval:
    """Genomic interval [start, end], 1-based inclusive for GenBank."""
    start: int
    end: int

    def length(self) -> int:
        return self.end - self.start + 1

    def __post_init__(self):
        if self.start < 1 or self.end < self.start:
            raise ValueError(f"Invalid interval: {self.start}..{self.end}")

class IntervalList(List[Interval]):
    """Ordered list of intervals preserving biological order (5'->3')."""
    
    def __init__(self, intervals: Iterable[Interval] = (), strand: Strand = Strand.PLUS):
        super().__init__(intervals)
        self.strand = strand

    def total_length(self) -> int:
        return sum(iv.length() for iv in self)
