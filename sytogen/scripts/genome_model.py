# build the genome model from the genbank file and codon usage table
# Some oop to encapsulate the genome and its features, and provide methods for editing and scoring

from sytogen.scripts.legacy_sytogen import reverse_complement
class GenomeModel:
    def __init__(
        self,
        topology="circular",
        sequence=None,
        genes=None,
        motifs=None,
        protected_regions=None,
        codon_usage=None
    ):
        self.sequence = sequence
        self.length = len(sequence)
        self.topology = topology  # "circular" or "linear"

        self.genes = genes or []
        self.motifs = motifs or []
        self.protected_regions = protected_regions or []
        self.codon_usage = codon_usage or {}

        self.topology_engine = (
            CircularTopology(sequence)
            if topology == "circular"
            else LinearTopology(sequence))

class Gene:
    def __init__(self, id, start, end, strand):
        self.id = id
        self.start = start
        self.end = end
        self.strand = strand

class Motif:
    def __init__(self, motif, start, end, strand):
        self.motif = motif
        self.start = start
        self.end = end
        self.strand = strand
        self.regex = compile_iupac(motif)
        self.length = len(motif)

# Define circular and linear topology classes to handle sequence indexing and motif counting
class CircularTopology:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)
    
    def get_window(self, center, radius):
        seq = self.sequence
        L = self.length
        start = (center - radius) % L
        end = (center + radius) % L
        if start < end:
            return seq[start:end]
        else:
            return seq[start:] + seq[:end]
    
class LinearTopology:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)
    
    def get_window(self, center, radius):
        seq = self.sequence
        start = max(0, center - radius)
        end = min(len(seq), center + radius)
        return seq[start:end]

class Mutation:
    def __init__(self, position, old, new):
        self.position = position
        self.old = old
        self.new = new

def evaluate_mutation(genome, mutation, window_radius=25):

    # --- STEP 1: REGION CHECK ---
    region = get_region(genome, mutation.position)

    if is_in_protected_region(genome, mutation.position):
        return {"valid": False, "reason": "Protected region"}

    # --- STEP 2: CDS CONSTRAINT ---
    if region == "CDS":
        if not is_synonymous(genome, mutation):
            return {"valid": False, "reason": "Not a synonymous change"}

    # --- STEP 3: APPLY MUTATION ---
    seq_list = list(genome.sequence)
    seq_list[mutation.position:mutation.position+len(mutation.old)] = list(mutation.new)
    mutated_seq = "".join(seq_list)

    # --- STEP 4: GET WINDOWS ---
    original_window = genome.topology_engine.get_window(
        mutation.position,
        window_radius)

    # Important: build *temporary topology* for mutated sequence
    if genome.topology == "circular":
        mutated_topology = CircularTopology(mutated_seq)
    else:
        mutated_topology = LinearTopology(mutated_seq)

    mutated_window = mutated_topology.get_window(
        mutation.position,
        window_radius)

    # --- STEP 5: REVERSE COMPLEMENTS ---
    original_rc = reverse_complement(original_window)
    mutated_rc = reverse_complement(mutated_window)

    # --- STEP 6: MOTIF COMPARISON ---
    destroyed = 0
    created = 0

    for motif in genome.motifs:

        before = (
            motif.regex.search(original_window) or
            motif.regex.search(original_rc))

        after = (
            motif.regex.search(mutated_window) or
            motif.regex.search(mutated_rc))

        if before and not after:
            destroyed += 1
        elif not before and after:
            created += 1

        # Early rejection
        if created > 0:
            return {"valid": False, "reason": "Creates new motif"}

    # --- STEP 7: RETURN RESULT ---
    return {"valid": True,
        "destroyed": destroyed,
        "created": created,
        "edits": 1
        }       
