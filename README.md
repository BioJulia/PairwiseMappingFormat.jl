# <img src="./sticker.svg" width="30%" align="right" /> PairwiseMappingFormat

[![Latest Release](https://img.shields.io/github/release/BioJulia/PairwiseMappingFormat.jl.svg)](https://github.com/BioJulia/PairwiseMappingFormat.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/PairwiseMappingFormat.jl/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/PairwiseMappingFormat.jl/stable)

PairwiseMappingFormat.jl provide a parser for Pairwise Mapping Format (PAF) files.
PAF is a simple, tab-delimited format created by programs such as minimap2.

To learn how to use the package, [read the documentation](https://biojulia.github.io/PairwiseMappingFormat.jl/stable/)

## Example
```julia
PAFReader(open("file.paf")) do reader
    for record in reader
        if aln_identity(record) > 0.95 && record.alnlen > 2_000
            println(record.qname, " aligns well to ", record.tname)
        end
    end
end
```

## Installation
Install PairwiseMappingFormat.jl from the julia
REPL. Press `]` to enter Pkg mode, and enter the following:

```julia
add PairwiseMappingFormat
```

## Contributing
Get in touch with the BioJulia community over at the [Julia Slack](https://julialang.org/slack/) or Zulip servers.
