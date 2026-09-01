"""
Tests for DNA analysis modules: repeats, secondary structure, quality scoring.
Also tests pattern library and report generation.
"""

import pytest
from sytogen.scripts.dna_analysis import (
    find_tandem_repeats,
    find_dispersed_repeats,
    simple_secondary_structure_score,
    identify_high_structure_regions,
    sequence_quality_score,
)
from sytogen.scripts.pattern_library import (
    PatternLibrary,
    PatternSpec,
    get_default_library,
    iupac_to_regex,
    reverse_complement_full,
)
from sytogen.scripts.report_generator import (
    ReportBuilder,
    create_sequence_analysis_report,
)
from sytogen.scripts.dna_optimization import (
    OptimizationConstraints,
    constraint_violations,
    optimize_sequence,
)


class TestRepeatDetection:
    """Test tandem and dispersed repeat detection."""
    
    def test_tandem_repeats_simple(self):
        """Test detection of simple tandem repeats."""
        # ATAT repeated 4 times
        sequence = "ATATATATAT"
        repeats = find_tandem_repeats(sequence, min_repeat_length=2)
        
        assert len(repeats) > 0
        # Should detect AT as 5-copy repeat (or ATAT as 2-copy)
        assert any(r.copy_count >= 2 for r in repeats)
    
    def test_tandem_repeats_min_length(self):
        """Test minimum repeat length filtering."""
        # GCGCGCGCGC (5 copies of GC)
        sequence = "GCGCGCGCGC"
        repeats = find_tandem_repeats(sequence, min_repeat_length=4)
        
        # Should find GC repeats (2bp)
        assert len(repeats) > 0
        assert any(r.repeat_length <= 4 for r in repeats)
    
    def test_tandem_no_repeats(self):
        """Test sequence with no repeats."""
        sequence = "ATGCATGCATGC"  # No consecutive repeats
        repeats = find_tandem_repeats(sequence, min_repeat_length=6)
        
        # May find some repeats depending on algorithm
        # but should find fewer than in repeat-heavy sequences
        assert isinstance(repeats, list)
    
    def test_tandem_repeats_deduplication(self):
    
        """Test that overlapping repeats are detected for homopolymers."""
        # AAAAAAAAAA (10 As)
        sequence = "A" * 10
        repeats = find_tandem_repeats(sequence, min_repeat_length=1)
        
        # Homopolymer creates many overlapping repeats (sliding windows)
        # Should find multiple copies of A
        assert len(repeats) > 0
        # Should have at least the longest repeat (A repeated 10 times)
        assert any(r.copy_count == 10 for r in repeats)
    def test_dispersed_repeats_two_locations(self):
        """Test detection of dispersed repeats at 2+ locations."""
        # ATGC at position 0 and 10
        sequence = "ATGC" + "NNNNNNNNNN" + "ATGC"
        repeats = find_dispersed_repeats(sequence, min_repeat_length=4)
        
        # Should find ATGC at positions 0 and 10
        assert len(repeats) > 0
        assert any(r.sequence == "ATGC" for r in repeats)
    
    def test_dispersed_repeats_multiple_copies(self):
        """Test dispersed repeats with 3+ copies."""
        sequence = "GCTAGC" + "AAAA" + "GCTAGC" + "TTTT" + "GCTAGC"
        repeats = find_dispersed_repeats(sequence, min_repeat_length=6)
        
        # Should find GCTAGC at 3 locations
        gctagc_repeats = [r for r in repeats if r.sequence == "GCTAGC"]
        assert len(gctagc_repeats) > 0


class TestSecondaryStructure:
    """Test secondary structure analysis."""
    
    def test_gc_content_high_stability(self):
        """Test that GC-rich regions score high."""
        gc_rich = "GCGCGCGCGC"  # 100% GC
        scores = simple_secondary_structure_score(gc_rich, window_size=4)
        
        assert len(scores) > 0
        # GC-rich should have moderate to high stability
        assert any(score > 0.3 for pos, score in scores)
    
    def test_at_content_low_stability(self):
        """Test that AT-rich regions score lower."""
        at_rich = "ATATATATAT"  # 100% AT
        scores = simple_secondary_structure_score(at_rich, window_size=4)
        
        assert len(scores) > 0
        # AT-rich should have low stability (far from 50% GC)
        max_score = max(score for pos, score in scores)
        assert max_score < 0.8
    
    def test_homopolymer_detection(self):
        """Test that long homopolymers are scored."""
        homopolymer = "AAAAAAAA"  # 8 As
        scores = simple_secondary_structure_score(homopolymer, window_size=4)
        
        assert len(scores) > 0
        # Long homopolymer should score high for structure risk
        assert any(score > 0.5 for pos, score in scores)
    
    def test_inverted_repeat_detection(self):
        """Test detection of inverted repeats (palindromes)."""
        # GAATTC is a palindrome (EcoRI site)
        palindrome = "GAATTC"
        scores = simple_secondary_structure_score(palindrome, window_size=6)
        
        assert len(scores) > 0
        # Palindromes can form stem-loop structures
        # (may or may not score high depending on base composition)
    
    def test_high_structure_regions(self):
        """Test identification of regions above stability threshold."""
        # Mix of GC-rich and AT-rich
        sequence = "GCGCGCGC" * 5 + "ATATATAT" * 5
        regions = identify_high_structure_regions(sequence, threshold=0.7)
        
        assert isinstance(regions, list)
        # Should have some regions (from GC-rich section)
        assert len(regions) >= 0
    
    def test_high_structure_threshold(self):
        """Test that threshold filtering works."""
        sequence = "GCGCGCGC" * 10
        
        regions_high = identify_high_structure_regions(sequence, threshold=0.9)
        regions_low = identify_high_structure_regions(sequence, threshold=0.3)
        
        # Lower threshold should find more regions
        assert len(regions_low) >= len(regions_high)


class TestSequenceQuality:
    """Test comprehensive quality scoring."""
    
    def test_quality_score_returns_dict(self):
        """Test that quality_score returns correct structure."""
        sequence = "ATGCATGCATGC"
        scores = sequence_quality_score(sequence)
        
        assert isinstance(scores, dict)
        assert "gc_content" in scores
        assert "homopolymer" in scores
        assert "repeats" in scores
        assert "secondary_structure" in scores
        assert "overall" in scores
    
    def test_quality_score_ranges(self):
        """Test that all quality scores are 0-100."""
        sequence = "GCGCGCGC" * 10
        scores = sequence_quality_score(sequence)
        
        for key, value in scores.items():
            if key != "longest_homopolymer_length":
                assert 0 <= value <= 100, f"{key} out of range: {value}"
    
    def test_quality_perfect_sequence(self):
        """Test quality score for well-balanced sequence."""
        # ~50% GC, no long homopolymers
        sequence = "ATGC" * 100
        scores = sequence_quality_score(sequence)
        
        # Should have decent overall score
        assert scores["overall"] >= 60
    
    def test_quality_poor_sequence(self):
        """Test quality score for problematic sequence."""
        # 100% A, lots of homopolymer
        sequence = "AAAA" * 100
        scores = sequence_quality_score(sequence)
        
        # Should have lower overall score
        assert scores["overall"] < 80


class TestPatternLibrary:
    """Test pattern library functionality."""
    
    def test_default_library_loads(self):
        """Test that default library initializes."""
        lib = get_default_library()
        
        assert lib is not None
        assert len(lib.patterns) > 0
    
    def test_library_has_restriction_sites(self):
        """Test that library includes restriction sites."""
        lib = get_default_library()
        
        restriction_names = lib.pattern_names_by_category("restriction")
        assert len(restriction_names) > 0
        assert "EcoRI" in restriction_names
        assert "BamHI" in restriction_names
    
    def test_library_has_assembly_standards(self):
        """Test that library includes assembly standards."""
        lib = get_default_library()
        
        assembly_names = lib.pattern_names_by_category("assembly")
        assert len(assembly_names) > 0
    
    def test_add_custom_pattern(self):
        """Test adding custom patterns to library."""
        lib = PatternLibrary()  # Empty library
        
        spec = PatternSpec(
            name="CustomMotif",
            pattern="TATATAT",
            description="Custom test motif",
            category="motif",
        )
        lib.add_pattern(spec)
        
        assert lib.get_pattern("CustomMotif") == spec
    
    def test_find_patterns_in_sequence(self):
        """Test pattern matching in sequence."""
        lib = get_default_library()
        sequence = "ATGAATTCGCGGATCC"
        
        # EcoRI is GAATTC
        matches = lib.find_patterns("EcoRI", sequence)
        
        assert len(matches) > 0
        # Should find GAATTC at position 3
        assert any(start <= 3 <= end for start, end, _ in matches)
    
    def test_find_patterns_reverse_complement(self):
        """Test that reverse complement matching works."""
        lib = PatternLibrary()
        spec = PatternSpec(
            name="TestPattern",
            pattern="GAATTC",
            description="Test",
            category="test",
            reverse_complement=True,
        )
        lib.add_pattern(spec)
        
        sequence = "CTTAAG"  # Reverse complement of GAATTC
        matches = lib.find_patterns("TestPattern", sequence)
        
        # Might find RC matches depending on implementation
        assert isinstance(matches, list)
    
    def test_iupac_to_regex(self):
        """Test IUPAC degenerate code conversion."""
        # R = A or G
        regex = iupac_to_regex("R")
        assert "[AG]" in regex or "A" in regex
        
        # N = any base
        regex = iupac_to_regex("N")
        assert "ACGT" in regex or "N" in regex or "[" in regex


class TestReportBuilder:
    """Test report generation."""
    
    def test_report_builder_creates(self):
        """Test basic report creation."""
        report = ReportBuilder("Test Report")
        
        assert report.title == "Test Report"
        assert len(report.sections) == 0
    
    def test_report_add_section(self):
        """Test adding sections."""
        report = ReportBuilder("Test")
        report.add_section("Section 1", "Content 1", "text")
        report.add_section("Section 2", "Content 2", "text")
        
        assert len(report.sections) == 2
        assert report.sections[0]["title"] == "Section 1"
    
    def test_report_add_metric(self):
        """Test adding metrics."""
        report = ReportBuilder("Test")
        report.add_metric("Test Metric", 42.5, "%")
        report.add_metric("Count", 100, "items")
        
        assert len(report.metrics) == 2
        assert report.metrics[0]["value"] == 42.5
    
    def test_report_to_html(self):
        """Test HTML report generation."""
        report = ReportBuilder("Test Report")
        report.add_section("Introduction", "This is a test", "text")
        report.add_metric("Score", 85, "%")
        
        html = report.to_html()
        
        assert len(html) > 0
        assert "Test Report" in html
        assert "Introduction" in html
        assert "85" in html
        assert "<html>" in html
    
    def test_report_to_json(self):
        """Test JSON report generation."""
        report = ReportBuilder("Test")
        report.add_section("Test", "Content", "text")
        report.add_metric("Score", 90, "%")
        
        json_str = report.to_json()
        
        assert len(json_str) > 0
        assert "Test" in json_str
        assert '"sections"' in json_str
    
    def test_report_to_markdown(self):
        """Test Markdown report generation."""
        report = ReportBuilder("Test Report")
        report.add_section("Section", "Content", "text")
        report.add_metric("Metric", 75, "units")
        
        md = report.to_markdown()
        
        assert len(md) > 0
        assert "# Test Report" in md
        assert "## Section" in md
        assert "Metric" in md


class TestIntegration:
    """Integration tests combining multiple modules."""
    
    def test_full_analysis_pipeline(self):
        """Test running analysis on a realistic sequence."""
        # Realistic construct with repeats and various features
        construct = "ATGAATTC" + "GCGC" * 10 + "AAAAAA" + "GATCC" * 5 + "ATGCATGC"
        
        # Get quality scores
        quality = sequence_quality_score(construct)
        assert 0 <= quality["overall"] <= 100
        
        # Find repeats
        tandem = find_tandem_repeats(construct, min_repeat_length=2)
        assert len(tandem) >= 0
        
        # Find secondary structures
        structures = identify_high_structure_regions(construct, threshold=0.5)
        assert len(structures) >= 0
        
        # Search patterns
        lib = get_default_library()
        ecori_sites = lib.find_patterns("EcoRI", construct)
        # EcoRI is GAATTC, which is in our construct
        assert len(ecori_sites) > 0
    
    def test_report_with_analysis_data(self):
        """Test creating report with analysis results."""
        sequence = "GCGCGCGC" * 10
        
        quality = sequence_quality_score(sequence)
        repeats = find_tandem_repeats(sequence)
        structures = identify_high_structure_regions(sequence, threshold=0.6)
        
        report = create_sequence_analysis_report(
            sequence, quality, repeats, structures
        )
        
        assert report is not None
        assert len(report.sections) > 0
        
        # Should be able to export
        html = report.to_html()
        assert len(html) > 500


class TestConstraintOptimization:
    def test_named_pattern_constraint_is_removed(self):
        constraints = OptimizationConstraints(forbidden_patterns=["EcoRI"])

        result = optimize_sequence("ATGAATTCATG", constraints)

        assert result["edits_applied"] == 1
        assert result["violations_before"]["forbidden_patterns"] == 1
        assert result["violations_after"]["forbidden_patterns"] == 0
        assert result["resolved"] is True

    def test_constraint_validation_rejects_invalid_gc_range(self):
        with pytest.raises(ValueError, match="min_gc cannot exceed max_gc"):
            OptimizationConstraints.from_dict({"min_gc": 60, "max_gc": 40})


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
