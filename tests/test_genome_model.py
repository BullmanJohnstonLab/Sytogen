from sytogen.scripts.genome_model import GenomeModel, Motif, Mutation


def test_neutral_candidates_ignore_edits_that_only_change_wildcard_motif_bases():
    sequence = ["A"] * 1200
    # 0-based positions 903..914 for ACCNNNNNRTGT
    motif_site = "ACCGTTTTAGTGT"
    sequence[903:903 + len(motif_site)] = list(motif_site)
    sequence = "".join(sequence)

    target = Motif("ACCNNNNNRTGT", 903, 914, "+")
    # Overlapping concrete motif that mutation at 906 can destroy.
    overlapping = Motif("CGTT", 905, 908, "+")

    genome = GenomeModel(
        sequence=sequence,
        topology="linear",
        genes=[],
        motifs=[target, overlapping],
        protected_regions=[],
        codon_usage={},
    )

    mutation = Mutation(position=906, old="G", new="A")

    result = genome.evaluate_mutation(mutation)
    assert result["valid"] is True
    assert result["destroyed"] >= 1
    assert genome._mutation_destroys_motif_occurrence(mutation, target) is False

    candidates = genome.generate_neutral_candidates(target)
    assert all(candidate.mutation.position != 906 for candidate in candidates)


def test_overlap_priority_prefers_shared_concrete_bases_over_wildcard_only_slots():
    sequence = ["A"] * 1200
    motif_site = "ACCGTTTTAGTGT"
    sequence[903:903 + len(motif_site)] = list(motif_site)
    sequence = "".join(sequence)

    target = Motif("ACCNNNNNRTGT", 903, 914, "+")
    overlapping = Motif("CGTT", 905, 908, "+")

    genome = GenomeModel(
        sequence=sequence,
        topology="linear",
        genes=[],
        motifs=[target, overlapping],
        protected_regions=[],
        codon_usage={},
    )

    shared_concrete = genome._mutation_overlap_priority(Mutation(position=905, old="C", new="A"), [target, overlapping])
    wildcard_only = genome._mutation_overlap_priority(Mutation(position=906, old="G", new="A"), [target, overlapping])

    assert shared_concrete > wildcard_only
    assert shared_concrete == 2
    assert wildcard_only == 1
