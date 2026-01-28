"""
XLEAP-SBS Color Balance Validator

Validates index combinations for Illumina XLEAP-SBS chemistry
(NovaSeq X Plus, NextSeq 1000/2000).

XLEAP-SBS Color Channels:
- A: Blue only (single channel)
- C: Blue + Green (dual channel)
- G: Dark (no signal)
- T: Green only (single channel)

For successful sequencing:
- Each position needs at least one blue signal (A or C)
- Each position needs at least one green signal (C or T)
- Avoid all G's at any position (causes registration failure)
"""

from dataclasses import dataclass
from typing import Optional


# XLEAP-SBS color profiles
XLEAP_COLORS = {
    'A': {'blue': True,  'green': False, 'label': 'Blue only'},
    'C': {'blue': True,  'green': True,  'label': 'Both (dual)'},
    'G': {'blue': False, 'green': False, 'label': 'Dark (no signal)'},
    'T': {'blue': False, 'green': True,  'label': 'Green only'},
}


@dataclass
class PositionResult:
    """Result for a single index position."""
    position: int
    bases: list[str]
    has_blue: bool
    has_green: bool
    is_balanced: bool
    unique_bases: int
    recommendation: Optional[str] = None


@dataclass
class ValidationResult:
    """Overall validation result."""
    is_valid: bool
    total_positions: int
    positions: list[PositionResult]
    failed_positions: list[int]
    warnings: list[str]
    summary: str


def parse_indexes(input_text: str) -> list[str]:
    """
    Parse user input into a list of index sequences.

    Accepts:
    - One index per line
    - Comma-separated indexes
    - Tab-separated indexes

    Returns uppercase, validated sequences.
    """
    # Handle different separators
    text = input_text.strip()

    # Try newlines first
    if '\n' in text:
        raw_indexes = text.split('\n')
    elif ',' in text:
        raw_indexes = text.split(',')
    elif '\t' in text:
        raw_indexes = text.split('\t')
    else:
        raw_indexes = [text]

    # Clean and validate
    indexes = []
    for idx in raw_indexes:
        # Remove whitespace and convert to uppercase
        cleaned = idx.strip().upper()

        # Skip empty lines
        if not cleaned:
            continue

        # Validate characters (only A, C, G, T allowed)
        if all(base in 'ACGT' for base in cleaned):
            indexes.append(cleaned)

    return indexes


def get_color_profile(base: str) -> dict:
    """
    Get the XLEAP color channel profile for a base.

    Args:
        base: Single nucleotide (A, C, G, or T)

    Returns:
        Dict with 'blue', 'green', and 'label' keys
    """
    return XLEAP_COLORS.get(base.upper(), {'blue': False, 'green': False, 'label': 'Unknown'})


def check_position_balance(position: int, bases: list[str]) -> PositionResult:
    """
    Check color balance at a single index position.

    Args:
        position: The position number (0-indexed internally, displayed as 1-indexed)
        bases: List of bases at this position across all indexes

    Returns:
        PositionResult with balance analysis
    """
    # Check for blue and green signals
    has_blue = any(XLEAP_COLORS[b]['blue'] for b in bases)
    has_green = any(XLEAP_COLORS[b]['green'] for b in bases)

    # Count unique bases for diversity
    unique_bases = len(set(bases))

    # Determine if balanced
    is_balanced = has_blue and has_green

    # Generate recommendation if not balanced
    recommendation = None
    if not is_balanced:
        if not has_blue and not has_green:
            recommendation = f"Position {position + 1}: All bases are G (dark). Add indexes with A, C, or T."
        elif not has_blue:
            recommendation = f"Position {position + 1}: No blue signal. Add indexes with A or C."
        elif not has_green:
            recommendation = f"Position {position + 1}: No green signal. Add indexes with C or T."
    elif unique_bases < 2:
        recommendation = f"Position {position + 1}: Low diversity ({unique_bases} unique base). Consider more varied indexes."

    return PositionResult(
        position=position + 1,  # 1-indexed for display
        bases=bases,
        has_blue=has_blue,
        has_green=has_green,
        is_balanced=is_balanced,
        unique_bases=unique_bases,
        recommendation=recommendation
    )


def validate_indexes(indexes: list[str]) -> ValidationResult:
    """
    Validate a set of indexes for XLEAP-SBS color balance.

    Args:
        indexes: List of index sequences (must be same length)

    Returns:
        ValidationResult with full analysis
    """
    # Handle empty input
    if not indexes:
        return ValidationResult(
            is_valid=False,
            total_positions=0,
            positions=[],
            failed_positions=[],
            warnings=["No indexes provided."],
            summary="No indexes to validate."
        )

    # Check all indexes are same length
    lengths = set(len(idx) for idx in indexes)
    if len(lengths) > 1:
        return ValidationResult(
            is_valid=False,
            total_positions=0,
            positions=[],
            failed_positions=[],
            warnings=[f"Indexes have different lengths: {lengths}. All indexes must be the same length."],
            summary="Index length mismatch."
        )

    index_length = len(indexes[0])

    # Check each position
    positions = []
    failed_positions = []
    warnings = []

    for pos in range(index_length):
        # Get all bases at this position
        bases_at_pos = [idx[pos] for idx in indexes]

        # Check balance
        result = check_position_balance(pos, bases_at_pos)
        positions.append(result)

        if not result.is_balanced:
            failed_positions.append(result.position)

        if result.recommendation:
            warnings.append(result.recommendation)

    # Additional warnings
    if len(indexes) < 2:
        warnings.insert(0, "Only 1 index provided. Color balance requires multiple pooled samples.")

    # Determine overall validity
    is_valid = len(failed_positions) == 0

    # Generate summary
    if is_valid:
        if warnings:
            summary = f"PASS with warnings. All {index_length} positions are color-balanced."
        else:
            summary = f"PASS. All {index_length} positions are color-balanced."
    else:
        summary = f"FAIL. {len(failed_positions)} of {index_length} positions lack color balance."

    return ValidationResult(
        is_valid=is_valid,
        total_positions=index_length,
        positions=positions,
        failed_positions=failed_positions,
        warnings=warnings,
        summary=summary
    )


def validate_from_text(input_text: str) -> ValidationResult:
    """
    Convenience function: parse and validate indexes from raw text.

    Args:
        input_text: Raw text with indexes (one per line, comma-separated, etc.)

    Returns:
        ValidationResult with full analysis
    """
    indexes = parse_indexes(input_text)
    return validate_indexes(indexes)


def get_color_balance_table(indexes: list[str]) -> list[dict]:
    """
    Generate a position-by-position color balance table for display.

    Returns list of dicts with position info for easy rendering.
    """
    if not indexes:
        return []

    table = []
    index_length = len(indexes[0])

    for pos in range(index_length):
        bases = [idx[pos] for idx in indexes]

        blue_count = sum(1 for b in bases if XLEAP_COLORS[b]['blue'])
        green_count = sum(1 for b in bases if XLEAP_COLORS[b]['green'])

        has_blue = blue_count > 0
        has_green = green_count > 0

        table.append({
            'position': pos + 1,
            'bases': bases,
            'base_counts': {b: bases.count(b) for b in 'ACGT'},
            'blue_count': blue_count,
            'green_count': green_count,
            'has_blue': has_blue,
            'has_green': has_green,
            'is_balanced': has_blue and has_green,
            'status': 'pass' if (has_blue and has_green) else 'fail'
        })

    return table
