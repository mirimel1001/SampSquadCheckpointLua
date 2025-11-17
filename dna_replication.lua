-- DNA Replication Module
-- Implements DNA replication and complementary base pairing

local dna_replication = {}

-- Base pairing rules (A-T, G-C)
local complement_map = {
    A = "T",
    T = "A",
    G = "C",
    C = "G",
    a = "t",
    t = "a",
    g = "c",
    c = "g"
}

-- Validate if a string is a valid DNA sequence
-- @param sequence: string containing DNA bases
-- @return: boolean indicating if sequence is valid
function dna_replication.is_valid_dna(sequence)
    if not sequence or type(sequence) ~= "string" or #sequence == 0 then
        return false
    end
    
    -- Check if all characters are valid DNA bases
    for i = 1, #sequence do
        local base = sequence:sub(i, i)
        if not (base == "A" or base == "T" or base == "G" or base == "C" or
                base == "a" or base == "t" or base == "g" or base == "c") then
            return false
        end
    end
    
    return true
end

-- Get the complementary base for a given base
-- @param base: single character DNA base
-- @return: complementary base or nil if invalid
function dna_replication.get_complement(base)
    return complement_map[base]
end

-- Generate complementary DNA strand (3' to 5' direction)
-- @param strand: original DNA strand (5' to 3')
-- @return: complementary strand or nil if invalid
function dna_replication.complement_strand(strand)
    if not dna_replication.is_valid_dna(strand) then
        return nil, "Invalid DNA sequence"
    end
    
    local complement = {}
    for i = 1, #strand do
        local base = strand:sub(i, i)
        table.insert(complement, complement_map[base])
    end
    
    return table.concat(complement)
end

-- Replicate DNA (create complementary strand and reverse it)
-- This simulates the biological process of DNA replication
-- @param template: template DNA strand
-- @return: table with original and complementary strands
function dna_replication.replicate(template)
    if not dna_replication.is_valid_dna(template) then
        return nil, "Invalid DNA sequence"
    end
    
    local complement = dna_replication.complement_strand(template)
    
    return {
        original = template,
        complement = complement,
        length = #template
    }
end

-- Reverse complement (for replication fork simulation)
-- @param strand: DNA strand
-- @return: reverse complement of the strand
function dna_replication.reverse_complement(strand)
    if not dna_replication.is_valid_dna(strand) then
        return nil, "Invalid DNA sequence"
    end
    
    local complement = dna_replication.complement_strand(strand)
    local reversed = {}
    
    for i = #complement, 1, -1 do
        table.insert(reversed, complement:sub(i, i))
    end
    
    return table.concat(reversed)
end

-- Transcribe DNA to RNA (replace T with U)
-- @param strand: DNA strand
-- @return: RNA sequence
function dna_replication.transcribe_to_rna(strand)
    if not dna_replication.is_valid_dna(strand) then
        return nil, "Invalid DNA sequence"
    end
    
    local rna = strand:gsub("[Tt]", function(t)
        return t == "T" and "U" or "u"
    end)
    
    return rna
end

-- Calculate GC content percentage
-- @param strand: DNA strand
-- @return: percentage of G and C bases
function dna_replication.gc_content(strand)
    if not dna_replication.is_valid_dna(strand) then
        return nil, "Invalid DNA sequence"
    end
    
    local gc_count = 0
    local upper_strand = strand:upper()
    
    for i = 1, #upper_strand do
        local base = upper_strand:sub(i, i)
        if base == "G" or base == "C" then
            gc_count = gc_count + 1
        end
    end
    
    return (gc_count / #strand) * 100
end

-- Split DNA into codons (triplets)
-- @param strand: DNA strand
-- @return: table of codons
function dna_replication.split_into_codons(strand)
    if not dna_replication.is_valid_dna(strand) then
        return nil, "Invalid DNA sequence"
    end
    
    local codons = {}
    for i = 1, #strand, 3 do
        local codon = strand:sub(i, i + 2)
        if #codon == 3 then
            table.insert(codons, codon)
        end
    end
    
    return codons
end

return dna_replication
