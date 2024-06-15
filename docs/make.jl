using Documenter, PairwiseMappingFormat

meta = quote
    using PairwiseMappingFormat
    path_to_paf = joinpath(pkgdir(PairwiseMappingFormat), "test", "example.paf")
    record = PAFReader(first, open(path_to_paf))
    const PAF = PairwiseMappingFormat
    using PairwiseMappingFormat: try_parse
end

DocMeta.setdocmeta!(
    PairwiseMappingFormat,
    :DocTestSetup,
    meta,
    recursive=true
)

makedocs(
    sitename = "PairwiseMappingFormat.jl",
    modules = [PairwiseMappingFormat],
    pages = [
        "Home" => "index.md",
    ],
    authors = "Jakob Nybo Nissen",
    checkdocs = :public,
    remotes=nothing, # TODO: Remove
)

# TODO: Call deploydocs
