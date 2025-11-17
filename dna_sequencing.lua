-- DNA Sequencing Module
-- Implements DNA sequence analysis and pattern matching

local dna_sequencing = {}

-- Count occurrences of each base in a sequence
-- @param sequence: DNA sequence
-- @return: table with counts of each base
function dna_sequencing.count_bases(sequence)
    if not sequence or type(sequence) ~= "string" then
        return nil, "Invalid sequence"
    end
    
    local counts = {A = 0, T = 0, G = 0, C = 0}
    local upper_seq = sequence:upper()
    
    for i = 1, #upper_seq do
        local base = upper_seq:sub(i, i)
        if counts[base] then
            counts[base] = counts[base] + 1
        end
    end
    
    return counts
end

-- Find all occurrences of a pattern in a sequence
-- @param sequence: DNA sequence to search in
-- @param pattern: pattern to search for
-- @return: table of starting positions (1-indexed)
function dna_sequencing.find_pattern(sequence, pattern)
    if not sequence or not pattern then
        return nil, "Invalid input"
    end
    
    local positions = {}
    local upper_seq = sequence:upper()
    local upper_pat = pattern:upper()
    
    local start_pos = 1
    while true do
        local pos = upper_seq:find(upper_pat, start_pos, true)
        if not pos then
            break
        end
        table.insert(positions, pos)
        start_pos = pos + 1
    end
    
    return positions
end

-- Calculate Hamming distance between two sequences
-- (number of positions at which bases differ)
-- @param seq1: first DNA sequence
-- @param seq2: second DNA sequence
-- @return: Hamming distance or nil if lengths differ
function dna_sequencing.hamming_distance(seq1, seq2)
    if not seq1 or not seq2 then
        return nil, "Invalid sequences"
    end
    
    if #seq1 ~= #seq2 then
        return nil, "Sequences must be of equal length"
    end
    
    local distance = 0
    local upper_seq1 = seq1:upper()
    local upper_seq2 = seq2:upper()
    
    for i = 1, #upper_seq1 do
        if upper_seq1:sub(i, i) ~= upper_seq2:sub(i, i) then
            distance = distance + 1
        end
    end
    
    return distance
end

-- Find the most frequent k-mer (substring of length k)
-- @param sequence: DNA sequence
-- @param k: length of k-mer
-- @return: most frequent k-mer and its count
function dna_sequencing.most_frequent_kmer(sequence, k)
    if not sequence or not k or k < 1 or k > #sequence then
        return nil, "Invalid input"
    end
    
    local kmer_counts = {}
    local upper_seq = sequence:upper()
    
    -- Count all k-mers
    for i = 1, #upper_seq - k + 1 do
        local kmer = upper_seq:sub(i, i + k - 1)
        kmer_counts[kmer] = (kmer_counts[kmer] or 0) + 1
    end
    
    -- Find most frequent
    local max_count = 0
    local max_kmer = nil
    
    for kmer, count in pairs(kmer_counts) do
        if count > max_count then
            max_count = count
            max_kmer = kmer
        end
    end
    
    return max_kmer, max_count
end

-- Find all k-mers with a given frequency
-- @param sequence: DNA sequence
-- @param k: length of k-mer
-- @param frequency: target frequency
-- @return: table of k-mers with the specified frequency
function dna_sequencing.find_kmers_by_frequency(sequence, k, frequency)
    if not sequence or not k or not frequency then
        return nil, "Invalid input"
    end
    
    local kmer_counts = {}
    local upper_seq = sequence:upper()
    
    -- Count all k-mers
    for i = 1, #upper_seq - k + 1 do
        local kmer = upper_seq:sub(i, i + k - 1)
        kmer_counts[kmer] = (kmer_counts[kmer] or 0) + 1
    end
    
    -- Find k-mers with target frequency
    local result = {}
    for kmer, count in pairs(kmer_counts) do
        if count == frequency then
            table.insert(result, kmer)
        end
    end
    
    return result
end

-- Calculate sequence similarity (percentage of matching bases)
-- @param seq1: first DNA sequence
-- @param seq2: second DNA sequence
-- @return: similarity percentage
function dna_sequencing.sequence_similarity(seq1, seq2)
    if not seq1 or not seq2 then
        return nil, "Invalid sequences"
    end
    
    local min_length = math.min(#seq1, #seq2)
    local matches = 0
    local upper_seq1 = seq1:upper()
    local upper_seq2 = seq2:upper()
    
    for i = 1, min_length do
        if upper_seq1:sub(i, i) == upper_seq2:sub(i, i) then
            matches = matches + 1
        end
    end
    
    return (matches / min_length) * 100
end

-- Find longest common substring between two sequences
-- @param seq1: first DNA sequence
-- @param seq2: second DNA sequence
-- @return: longest common substring
function dna_sequencing.longest_common_substring(seq1, seq2)
    if not seq1 or not seq2 then
        return nil, "Invalid sequences"
    end
    
    local upper_seq1 = seq1:upper()
    local upper_seq2 = seq2:upper()
    local longest = ""
    
    for i = 1, #upper_seq1 do
        for j = i, #upper_seq1 do
            local substring = upper_seq1:sub(i, j)
            if upper_seq2:find(substring, 1, true) and #substring > #longest then
                longest = substring
            end
        end
    end
    
    return longest
end

-- Generate sequence statistics
-- @param sequence: DNA sequence
-- @return: table with various statistics
function dna_sequencing.sequence_stats(sequence)
    if not sequence or type(sequence) ~= "string" then
        return nil, "Invalid sequence"
    end
    
    local counts = dna_sequencing.count_bases(sequence)
    local total = counts.A + counts.T + counts.G + counts.C
    
    return {
        length = total,
        base_counts = counts,
        gc_content = ((counts.G + counts.C) / total) * 100,
        at_content = ((counts.A + counts.T) / total) * 100,
        composition = {
            A = (counts.A / total) * 100,
            T = (counts.T / total) * 100,
            G = (counts.G / total) * 100,
            C = (counts.C / total) * 100
        }
    }
end

-- Check if a sequence is a DNA palindrome (reverse complement equals original)
-- This is the biological definition used for restriction sites
-- @param sequence: DNA sequence
-- @return: boolean indicating if palindromic
function dna_sequencing.is_palindrome(sequence)
    if not sequence then
        return false
    end
    
    local upper_seq = sequence:upper()
    
    -- Create reverse complement
    local complement_map = {
        A = "T",
        T = "A",
        G = "C",
        C = "G"
    }
    
    local rev_comp = {}
    for i = #upper_seq, 1, -1 do
        local base = upper_seq:sub(i, i)
        local comp = complement_map[base]
        if not comp then
            return false
        end
        table.insert(rev_comp, comp)
    end
    
    return upper_seq == table.concat(rev_comp)
end

-- Find restriction sites (palindromic sequences of length 4-12)
-- @param sequence: DNA sequence
-- @param min_len: minimum length (default 4)
-- @param max_len: maximum length (default 12)
-- @return: table of restriction sites with positions
function dna_sequencing.find_restriction_sites(sequence, min_len, max_len)
    if not sequence then
        return nil, "Invalid sequence"
    end
    
    min_len = min_len or 4
    max_len = max_len or 12
    local sites = {}
    local upper_seq = sequence:upper()
    
    for len = min_len, math.min(max_len, #upper_seq) do
        for i = 1, #upper_seq - len + 1 do
            local substr = upper_seq:sub(i, i + len - 1)
            if dna_sequencing.is_palindrome(substr) then
                table.insert(sites, {
                    position = i,
                    sequence = substr,
                    length = len
                })
            end
        end
    end
    
    return sites
end

return dna_sequencing
