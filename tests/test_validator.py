"""
Unit tests for XLEAP-SBS Color Balance Validator.
"""

import pytest
from core.validator import (
    parse_indexes,
    get_color_profile,
    check_position_balance,
    validate_indexes,
    validate_from_text,
    get_color_balance_table,
    XLEAP_COLORS,
)


class TestXLEAPColors:
    """Test the XLEAP color definitions."""

    def test_a_is_blue_only(self):
        assert XLEAP_COLORS['A']['blue'] is True
        assert XLEAP_COLORS['A']['green'] is False

    def test_c_is_dual_channel(self):
        assert XLEAP_COLORS['C']['blue'] is True
        assert XLEAP_COLORS['C']['green'] is True

    def test_g_is_dark(self):
        assert XLEAP_COLORS['G']['blue'] is False
        assert XLEAP_COLORS['G']['green'] is False

    def test_t_is_green_only(self):
        assert XLEAP_COLORS['T']['blue'] is False
        assert XLEAP_COLORS['T']['green'] is True


class TestParseIndexes:
    """Test index parsing."""

    def test_parse_newline_separated(self):
        text = "ATCG\nGCTA\nTTAA"
        result = parse_indexes(text)
        assert result == ['ATCG', 'GCTA', 'TTAA']

    def test_parse_comma_separated(self):
        text = "ATCG, GCTA, TTAA"
        result = parse_indexes(text)
        assert result == ['ATCG', 'GCTA', 'TTAA']

    def test_parse_lowercase_converted(self):
        text = "atcg\ngcta"
        result = parse_indexes(text)
        assert result == ['ATCG', 'GCTA']

    def test_parse_skips_empty_lines(self):
        text = "ATCG\n\nGCTA\n"
        result = parse_indexes(text)
        assert result == ['ATCG', 'GCTA']

    def test_parse_skips_invalid_characters(self):
        text = "ATCG\nXYZZ\nGCTA"
        result = parse_indexes(text)
        assert result == ['ATCG', 'GCTA']

    def test_parse_single_index(self):
        text = "ATCGATCG"
        result = parse_indexes(text)
        assert result == ['ATCGATCG']


class TestColorProfile:
    """Test color profile retrieval."""

    def test_get_color_profile_a(self):
        profile = get_color_profile('A')
        assert profile['blue'] is True
        assert profile['green'] is False

    def test_get_color_profile_lowercase(self):
        profile = get_color_profile('a')
        assert profile['blue'] is True
        assert profile['green'] is False

    def test_get_color_profile_unknown(self):
        profile = get_color_profile('X')
        assert profile['blue'] is False
        assert profile['green'] is False


class TestCheckPositionBalance:
    """Test single position balance checking."""

    def test_balanced_with_c(self):
        # C provides both blue and green
        result = check_position_balance(0, ['C', 'C', 'C'])
        assert result.is_balanced is True
        assert result.has_blue is True
        assert result.has_green is True

    def test_balanced_with_a_and_t(self):
        # A provides blue, T provides green
        result = check_position_balance(0, ['A', 'T'])
        assert result.is_balanced is True

    def test_unbalanced_all_g(self):
        # G provides no signal
        result = check_position_balance(0, ['G', 'G', 'G'])
        assert result.is_balanced is False
        assert result.has_blue is False
        assert result.has_green is False
        assert 'dark' in result.recommendation.lower()

    def test_unbalanced_no_blue(self):
        # T only provides green
        result = check_position_balance(0, ['T', 'T', 'G'])
        assert result.is_balanced is False
        assert result.has_blue is False
        assert result.has_green is True

    def test_unbalanced_no_green(self):
        # A only provides blue
        result = check_position_balance(0, ['A', 'A', 'G'])
        assert result.is_balanced is False
        assert result.has_blue is True
        assert result.has_green is False


class TestValidateIndexes:
    """Test full index validation."""

    def test_balanced_indexes_pass(self):
        # Well-balanced set
        indexes = ['ATCG', 'CGTA', 'TGCA', 'GACT']
        result = validate_indexes(indexes)
        assert result.is_valid is True
        assert len(result.failed_positions) == 0

    def test_all_g_position_fails(self):
        # Position 1 is all G
        indexes = ['GATC', 'GATC', 'GATC', 'GATC']
        result = validate_indexes(indexes)
        assert result.is_valid is False
        assert 1 in result.failed_positions

    def test_empty_input(self):
        result = validate_indexes([])
        assert result.is_valid is False
        assert 'No indexes' in result.warnings[0]

    def test_mismatched_lengths(self):
        indexes = ['ATCG', 'ATCGATCG']
        result = validate_indexes(indexes)
        assert result.is_valid is False
        assert 'length' in result.warnings[0].lower()

    def test_single_index_warning(self):
        indexes = ['ATCGATCG']
        result = validate_indexes(indexes)
        # Single index with all different bases is technically balanced
        # but should warn about single sample
        assert 'Only 1 index' in result.warnings[0]


class TestValidateFromText:
    """Test convenience function."""

    def test_validate_from_text_balanced(self):
        text = """ATCGATCG
GCTAGCTA
TTAACCGG
CGATCGAT"""
        result = validate_from_text(text)
        assert result.is_valid is True

    def test_validate_from_text_unbalanced(self):
        text = """GGGGGGGG
GGGGGGGG"""
        result = validate_from_text(text)
        assert result.is_valid is False


class TestColorBalanceTable:
    """Test table generation for display."""

    def test_table_generation(self):
        indexes = ['AT', 'CG']
        table = get_color_balance_table(indexes)

        assert len(table) == 2

        # Position 1: A and C (both have blue)
        assert table[0]['position'] == 1
        assert table[0]['has_blue'] is True
        assert table[0]['has_green'] is True  # C has green

        # Position 2: T and G
        assert table[1]['position'] == 2
        assert table[1]['has_green'] is True  # T has green
        assert table[1]['has_blue'] is False  # Neither T nor G has blue

    def test_empty_indexes(self):
        table = get_color_balance_table([])
        assert table == []


# Real-world test cases
class TestRealWorldScenarios:
    """Test with realistic index combinations."""

    def test_illumina_udp_style_balanced(self):
        """Simulated balanced UDI set."""
        indexes = [
            'GCGATCTA',
            'AACAGGTT',
            'TGGTCACA',
            'CTCAGAGC',
        ]
        result = validate_indexes(indexes)
        assert result.is_valid is True

    def test_problematic_low_plex(self):
        """Low-plex pool that might fail."""
        indexes = [
            'GGGGAAAA',
            'GGGGAAAA',
        ]
        result = validate_indexes(indexes)
        # Position 1-4 are all G (no signal), positions 5-8 are all A (no green)
        assert result.is_valid is False
        assert len(result.failed_positions) == 8  # All positions fail
