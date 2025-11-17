# DNA Replication & Sequencing Modules

Comprehensive Lua modules for DNA sequence analysis, replication, and sequencing operations.

## Modules

### dna_replication.lua

Implements DNA replication and complementary base pairing operations.

#### Functions

- **`is_valid_dna(sequence)`** - Validates if a string is a valid DNA sequence
  - Returns: `boolean`
  - Example: `is_valid_dna("ATGC")` → `true`

- **`get_complement(base)`** - Gets the complementary base (A↔T, G↔C)
  - Returns: `string` or `nil`
  - Example: `get_complement("A")` → `"T"`

- **`complement_strand(strand)`** - Generates complementary DNA strand
  - Returns: `string` or `nil`
  - Example: `complement_strand("ATGC")` → `"TACG"`

- **`replicate(template)`** - Replicates DNA (creates complementary strand)
  - Returns: `table` with `original`, `complement`, and `length` fields
  - Example: `replicate("ATGC")` → `{original="ATGC", complement="TACG", length=4}`

- **`reverse_complement(strand)`** - Creates reverse complement of strand
  - Returns: `string` or `nil`
  - Example: `reverse_complement("ATGC")` → `"GCAT"`

- **`transcribe_to_rna(strand)`** - Transcribes DNA to RNA (T→U)
  - Returns: `string` or `nil`
  - Example: `transcribe_to_rna("ATGC")` → `"AUGC"`

- **`gc_content(strand)`** - Calculates GC content percentage
  - Returns: `number` (percentage)
  - Example: `gc_content("ATGC")` → `50.0`

- **`split_into_codons(strand)`** - Splits DNA into triplet codons
  - Returns: `table` of codon strings
  - Example: `split_into_codons("ATGCGA")` → `{"ATG", "CGA"}`

### dna_sequencing.lua

Implements DNA sequence analysis and pattern matching operations.

#### Functions

- **`count_bases(sequence)`** - Counts occurrences of each base
  - Returns: `table` with A, T, G, C counts
  - Example: `count_bases("ATGC")` → `{A=1, T=1, G=1, C=1}`

- **`find_pattern(sequence, pattern)`** - Finds all occurrences of a pattern
  - Returns: `table` of positions (1-indexed)
  - Example: `find_pattern("ATGATG", "ATG")` → `{1, 4}`

- **`hamming_distance(seq1, seq2)`** - Calculates Hamming distance
  - Returns: `number` or `nil`
  - Example: `hamming_distance("ATGC", "TTGC")` → `1`

- **`most_frequent_kmer(sequence, k)`** - Finds most frequent k-mer
  - Returns: `string`, `number` (k-mer and count)
  - Example: `most_frequent_kmer("ATGATG", 3)` → `"ATG", 2`

- **`find_kmers_by_frequency(sequence, k, frequency)`** - Finds k-mers with specific frequency
  - Returns: `table` of k-mer strings
  - Example: `find_kmers_by_frequency("ATGATG", 3, 2)` → `{"ATG"}`

- **`sequence_similarity(seq1, seq2)`** - Calculates sequence similarity percentage
  - Returns: `number` (percentage)
  - Example: `sequence_similarity("ATGC", "ATGC")` → `100.0`

- **`longest_common_substring(seq1, seq2)`** - Finds longest common substring
  - Returns: `string`
  - Example: `longest_common_substring("ATGCG", "TGCGA")` → `"TGCG"`

- **`sequence_stats(sequence)`** - Generates comprehensive sequence statistics
  - Returns: `table` with length, base_counts, gc_content, at_content, composition
  - Example: `sequence_stats("ATGC")` → detailed statistics table

- **`is_palindrome(sequence)`** - Checks if sequence is a DNA palindrome (reverse complement)
  - Returns: `boolean`
  - Example: `is_palindrome("ATCGAT")` → `true`

- **`find_restriction_sites(sequence, min_len, max_len)`** - Finds palindromic restriction sites
  - Returns: `table` of sites with position, sequence, and length
  - Example: `find_restriction_sites("ATGCGCAT", 4, 8)` → `{{position=3, sequence="GCGC", length=4}}`

## Usage Examples

### Basic Usage

```lua
local dna_replication = require("dna_replication")
local dna_sequencing = require("dna_sequencing")

-- Validate and replicate DNA
local strand = "ATGCGATCG"
if dna_replication.is_valid_dna(strand) then
    local result = dna_replication.replicate(strand)
    print("Original:   " .. result.original)
    print("Complement: " .. result.complement)
end

-- Analyze sequence
local stats = dna_sequencing.sequence_stats(strand)
print("GC Content: " .. stats.gc_content .. "%")

-- Find patterns
local positions = dna_sequencing.find_pattern(strand, "GC")
print("Pattern 'GC' found at: " .. table.concat(positions, ", "))
```

### Advanced Operations

```lua
-- Find restriction sites
local sites = dna_sequencing.find_restriction_sites("ATGCGCAT", 4, 8)
for _, site in ipairs(sites) do
    print("Site at position " .. site.position .. ": " .. site.sequence)
end

-- Calculate sequence similarity
local seq1 = "ATGCGATCG"
local seq2 = "ATGCAATCG"
local similarity = dna_sequencing.sequence_similarity(seq1, seq2)
print("Similarity: " .. similarity .. "%")

-- Transcribe to RNA
local rna = dna_replication.transcribe_to_rna(strand)
print("RNA: " .. rna)
```

## Running Examples and Tests

### Run All Examples
```bash
lua5.3 dna_examples.lua
```

### Run Tests
```bash
lua5.3 dna_tests.lua
```

## Technical Notes

### DNA Palindromes

The modules implement **biological DNA palindromes**, which are sequences that read the same as their reverse complement. For example:
- `ATCGAT` is a palindrome because its reverse complement is also `ATCGAT`
- This is different from simple string palindromes

### Base Pairing

Standard Watson-Crick base pairing rules:
- Adenine (A) pairs with Thymine (T)
- Guanine (G) pairs with Cytosine (C)

### Coordinate System

All position values use 1-based indexing (Lua convention).

## Requirements

- Lua 5.3 or higher

## License

Part of the SampSquadCheckpointLua repository.
