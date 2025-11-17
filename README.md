# Samp Squad Checkpoint++
A lua modification dedicated to the Horizong Gaming San Andreas Multiplayer server by the discord community of SA:MP Squad (https://discord.gg/mkMy3Cd).

This is a Moonloader-LUA modification that allows you to set a checkpoint on a SA:MP server. Developers of the modification are players from HZG so almost every location built in the mod that you can find in the GUI is based on the HZG servers hotspots or most important location for players to find. There are four different custom locations that can be set by players, accessed one at a time and can be re-adjusted as well. It is a very adjustable and user friendly modification.


# cpganglist.ini

This has the updated gang list and coordinates.


# version

This includes the version name and cpganglist.ini link.


# DNA Replication & Sequencing Modules

This repository also includes comprehensive Lua modules for DNA sequence analysis:

- **dna_replication.lua** - DNA replication, complementary base pairing, transcription
- **dna_sequencing.lua** - Sequence analysis, pattern matching, restriction sites

See [DNA_README.md](DNA_README.md) for detailed documentation and usage examples.

### Quick Start

```lua
local dna_replication = require("dna_replication")
local dna_sequencing = require("dna_sequencing")

-- Replicate DNA
local result = dna_replication.replicate("ATGCGATCG")
print("Original:   " .. result.original)
print("Complement: " .. result.complement)

-- Analyze sequence
local stats = dna_sequencing.sequence_stats("ATGCGATCG")
print("GC Content: " .. stats.gc_content .. "%")
```

Run examples: `lua5.3 dna_examples.lua`  
Run tests: `lua5.3 dna_tests.lua`
