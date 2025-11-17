-- DNA Module Tests
-- Comprehensive tests for DNA replication and sequencing modules

local dna_replication = require("dna_replication")
local dna_sequencing = require("dna_sequencing")

local tests_passed = 0
local tests_failed = 0

-- Test helper function
local function assert_equal(actual, expected, test_name)
    if actual == expected then
        print("✓ PASS: " .. test_name)
        tests_passed = tests_passed + 1
        return true
    else
        print("✗ FAIL: " .. test_name)
        print("  Expected: " .. tostring(expected))
        print("  Got:      " .. tostring(actual))
        tests_failed = tests_failed + 1
        return false
    end
end

local function assert_true(condition, test_name)
    return assert_equal(condition, true, test_name)
end

local function assert_false(condition, test_name)
    return assert_equal(condition, false, test_name)
end

print("=== DNA Replication Tests ===\n")

-- Test 1: Valid DNA validation
assert_true(dna_replication.is_valid_dna("ATGC"), "Valid DNA uppercase")
assert_true(dna_replication.is_valid_dna("atgc"), "Valid DNA lowercase")
assert_true(dna_replication.is_valid_dna("ATGCatgc"), "Valid DNA mixed case")
assert_false(dna_replication.is_valid_dna("ATGCX"), "Invalid DNA with X")
assert_false(dna_replication.is_valid_dna(""), "Empty string")
assert_false(dna_replication.is_valid_dna(nil), "Nil input")

-- Test 2: Complement generation
assert_equal(dna_replication.get_complement("A"), "T", "Complement of A")
assert_equal(dna_replication.get_complement("T"), "A", "Complement of T")
assert_equal(dna_replication.get_complement("G"), "C", "Complement of G")
assert_equal(dna_replication.get_complement("C"), "G", "Complement of C")

-- Test 3: Complement strand
assert_equal(dna_replication.complement_strand("ATGC"), "TACG", "Complement strand ATGC")
assert_equal(dna_replication.complement_strand("AAAA"), "TTTT", "Complement strand AAAA")
assert_equal(dna_replication.complement_strand("GCGC"), "CGCG", "Complement strand GCGC")

-- Test 4: Replication
local rep = dna_replication.replicate("ATGC")
assert_equal(rep.original, "ATGC", "Replication original")
assert_equal(rep.complement, "TACG", "Replication complement")
assert_equal(rep.length, 4, "Replication length")

-- Test 5: Reverse complement
assert_equal(dna_replication.reverse_complement("ATGC"), "GCAT", "Reverse complement ATGC")
assert_equal(dna_replication.reverse_complement("AAAA"), "TTTT", "Reverse complement AAAA")

-- Test 6: Transcription
assert_equal(dna_replication.transcribe_to_rna("ATGC"), "AUGC", "Transcription ATGC")
assert_equal(dna_replication.transcribe_to_rna("TTTT"), "UUUU", "Transcription TTTT")
assert_equal(dna_replication.transcribe_to_rna("atgc"), "augc", "Transcription lowercase")

-- Test 7: GC content
assert_equal(dna_replication.gc_content("ATGC"), 50, "GC content 50%")
assert_equal(dna_replication.gc_content("AAAA"), 0, "GC content 0%")
assert_equal(dna_replication.gc_content("GGCC"), 100, "GC content 100%")

-- Test 8: Codons
local codons = dna_replication.split_into_codons("ATGCGATAA")
assert_equal(#codons, 3, "Number of codons")
assert_equal(codons[1], "ATG", "First codon")
assert_equal(codons[2], "CGA", "Second codon")
assert_equal(codons[3], "TAA", "Third codon")

print("\n=== DNA Sequencing Tests ===\n")

-- Test 9: Base counting
local counts = dna_sequencing.count_bases("ATGCATGC")
assert_equal(counts.A, 2, "Count A bases")
assert_equal(counts.T, 2, "Count T bases")
assert_equal(counts.G, 2, "Count G bases")
assert_equal(counts.C, 2, "Count C bases")

-- Test 10: Pattern finding
local positions = dna_sequencing.find_pattern("ATGATGATG", "ATG")
assert_equal(#positions, 3, "Number of pattern occurrences")
assert_equal(positions[1], 1, "First pattern position")
assert_equal(positions[2], 4, "Second pattern position")
assert_equal(positions[3], 7, "Third pattern position")

-- Test 11: Hamming distance
assert_equal(dna_sequencing.hamming_distance("ATGC", "ATGC"), 0, "Hamming distance identical")
assert_equal(dna_sequencing.hamming_distance("ATGC", "TTGC"), 1, "Hamming distance 1")
assert_equal(dna_sequencing.hamming_distance("ATGC", "CGTA"), 4, "Hamming distance 4")

-- Test 12: Most frequent k-mer
local kmer, count = dna_sequencing.most_frequent_kmer("ATGATGATG", 3)
assert_equal(kmer, "ATG", "Most frequent k-mer")
assert_equal(count, 3, "K-mer frequency")

-- Test 13: K-mers by frequency
local kmers = dna_sequencing.find_kmers_by_frequency("ATGATGATG", 3, 3)
assert_equal(#kmers, 1, "Number of k-mers with frequency 3")
assert_equal(kmers[1], "ATG", "K-mer with frequency 3")

-- Test 14: Sequence similarity
assert_equal(dna_sequencing.sequence_similarity("ATGC", "ATGC"), 100, "Identical sequences")
assert_equal(dna_sequencing.sequence_similarity("ATGC", "TTTT"), 25, "25% similarity")

-- Test 15: Longest common substring
assert_equal(dna_sequencing.longest_common_substring("ATGCG", "TGCGA"), "TGCG", "Longest common substring")
assert_equal(dna_sequencing.longest_common_substring("AAAA", "TTTT"), "", "No common substring")

-- Test 16: Sequence statistics
local stats = dna_sequencing.sequence_stats("ATGC")
assert_equal(stats.length, 4, "Stats length")
assert_equal(stats.gc_content, 50, "Stats GC content")
assert_equal(stats.at_content, 50, "Stats AT content")

-- Test 17: Palindrome detection (DNA palindrome = reverse complement)
assert_true(dna_sequencing.is_palindrome("ATCGAT"), "Palindrome ATCGAT")
assert_true(dna_sequencing.is_palindrome("GAATTC"), "Palindrome GAATTC (EcoRI)")
assert_false(dna_sequencing.is_palindrome("AAAA"), "Not palindrome AAAA")
assert_false(dna_sequencing.is_palindrome("ATGC"), "Not palindrome ATGC")

-- Test 18: Restriction sites
local sites = dna_sequencing.find_restriction_sites("ATCGATAA", 4, 8)
assert_true(#sites >= 1, "Find restriction sites")
-- Check if ATCGAT is found (it's a palindrome)
local found = false
for _, site in ipairs(sites) do
    if site.sequence == "ATCGAT" then
        found = true
        break
    end
end
assert_true(found, "Find palindromic restriction site")

-- Print summary
print("\n=== Test Summary ===")
print("Tests Passed: " .. tests_passed)
print("Tests Failed: " .. tests_failed)
print("Total Tests:  " .. (tests_passed + tests_failed))

if tests_failed == 0 then
    print("\n✓ All tests passed!")
    os.exit(0)
else
    print("\n✗ Some tests failed!")
    os.exit(1)
end
