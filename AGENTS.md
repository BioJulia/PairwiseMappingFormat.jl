# AGENTS.md
This file contains instructions for AI agents working with this repository.

Purpose: PairwiseMappingFormat.jl provides parsing for the bioinformatics PAF format.

Design goals:
* Be correct, even in edge cases. Unrealistic edge cases like 2 GiB lines or 4 GiB references may not be optimized for, but should be handled.
* Be performant

Scope:
* Writing or mutating records is intentionally not supported yet.

Guidelines for developing
* After making changes, format with `runic -i .`. If Runic is not installed, alert the user and do not attempt to install it.
