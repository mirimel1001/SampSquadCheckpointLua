-- DNA Replication and Sequencing Examples
-- This file demonstrates how to use the DNA modules

local dna_replication = require("dna_replication")
local dna_sequencing = require("dna_sequencing")

print("=== DNA Replication Examples ===\n")

-- Example 1: Basic DNA validation
local dna_strand = "ATGCGATCGTAGC"
print("1. Validating DNA strand: " .. dna_strand)
print("   Valid: " .. tostring(dna_replication.is_valid_dna(dna_strand)))
print()

-- Example 2: Complement strand
print("2. Creating complementary strand:")
print("   Original:     " .. dna_strand)
local complement = dna_replication.complement_strand(dna_strand)
print("   Complement:   " .. complement)
print()

-- Example 3: DNA replication
print("3. DNA Replication:")
local replication = dna_replication.replicate(dna_strand)
print("   Original:     " .. replication.original)
print("   Complement:   " .. replication.complement)
print("   Length:       " .. replication.length .. " base pairs")
print()

-- Example 4: Reverse complement
print("4. Reverse Complement:")
local rev_comp = dna_replication.reverse_complement(dna_strand)
print("   Original:     " .. dna_strand)
print("   Rev Comp:     " .. rev_comp)
print()

-- Example 5: Transcription to RNA
print("5. DNA Transcription to RNA:")
local rna = dna_replication.transcribe_to_rna(dna_strand)
print("   DNA:  " .. dna_strand)
print("   RNA:  " .. rna)
print()

-- Example 6: GC content
print("6. GC Content Analysis:")
local gc = dna_replication.gc_content(dna_strand)
print("   Sequence: " .. dna_strand)
print("   GC%:      " .. string.format("%.2f%%", gc))
print()

-- Example 7: Codons
print("7. Splitting into Codons:")
local codons = dna_replication.split_into_codons(dna_strand)
print("   Sequence: " .. dna_strand)
print("   Codons:   " .. table.concat(codons, " "))
print()

print("\n=== DNA Sequencing Examples ===\n")

-- Example 8: Base counting
print("8. Base Counting:")
local counts = dna_sequencing.count_bases(dna_strand)
print("   Sequence: " .. dna_strand)
print("   A: " .. counts.A .. ", T: " .. counts.T .. ", G: " .. counts.G .. ", C: " .. counts.C)
print()

-- Example 9: Pattern finding
print("9. Pattern Finding:")
local pattern = "GC"
local positions = dna_sequencing.find_pattern(dna_strand, pattern)
print("   Sequence: " .. dna_strand)
print("   Pattern:  " .. pattern)
print("   Found at positions: " .. table.concat(positions, ", "))
print()

-- Example 10: Hamming distance
print("10. Hamming Distance:")
local seq1 = "ATGCGATCG"
local seq2 = "ATGCAATCG"
local distance = dna_sequencing.hamming_distance(seq1, seq2)
print("   Seq1:     " .. seq1)
print("   Seq2:     " .. seq2)
print("   Distance: " .. distance .. " base(s) differ")
print()

-- Example 11: Most frequent k-mer
print("11. Most Frequent K-mer:")
local test_seq = "ATGATGATGATGCCC"
local kmer, count = dna_sequencing.most_frequent_kmer(test_seq, 3)
print("   Sequence:         " .. test_seq)
print("   Most frequent 3-mer: " .. kmer .. " (appears " .. count .. " times)")
print()

-- Example 12: Sequence similarity
print("12. Sequence Similarity:")
local similarity = dna_sequencing.sequence_similarity(seq1, seq2)
print("   Seq1:       " .. seq1)
print("   Seq2:       " .. seq2)
print("   Similarity: " .. string.format("%.2f%%", similarity))
print()

-- Example 13: Longest common substring
print("13. Longest Common Substring:")
local lcs = dna_sequencing.longest_common_substring(seq1, seq2)
print("   Seq1: " .. seq1)
print("   Seq2: " .. seq2)
print("   LCS:  " .. lcs)
print()

-- Example 14: Sequence statistics
print("14. Sequence Statistics:")
local stats = dna_sequencing.sequence_stats(dna_strand)
print("   Sequence:   " .. dna_strand)
print("   Length:     " .. stats.length)
print("   GC Content: " .. string.format("%.2f%%", stats.gc_content))
print("   AT Content: " .. string.format("%.2f%%", stats.at_content))
print("   Composition:")
print("     A: " .. string.format("%.2f%%", stats.composition.A))
print("     T: " .. string.format("%.2f%%", stats.composition.T))
print("     G: " .. string.format("%.2f%%", stats.composition.G))
print("     C: " .. string.format("%.2f%%", stats.composition.C))
print()

-- Example 15: Palindrome detection
print("15. Palindrome Detection:")
local palindrome = "ATCGAT"
print("   Sequence:    " .. palindrome)
print("   Palindrome:  " .. tostring(dna_sequencing.is_palindrome(palindrome)))
print()

-- Example 16: Restriction sites
print("16. Finding Restriction Sites:")
local restriction_seq = "ATGCGCATTAGCAT"
local sites = dna_sequencing.find_restriction_sites(restriction_seq, 4, 6)
print("   Sequence: " .. restriction_seq)
print("   Found " .. #sites .. " restriction site(s):")
for i, site in ipairs(sites) do
    print("     Position " .. site.position .. ": " .. site.sequence .. " (length " .. site.length .. ")")
end
print()

print("=== Examples Complete ===")
