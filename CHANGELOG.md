# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0
### Changed
- **Breaking:** Allow XAMAuxData 0.2. Since `aux_data` returns a `XAMAuxData.SAM.Auxiliary`,
  the breaking changes of XAMAuxData 0.2 apply to auxiliary data of PAF records.
  For example, `H` values with lowercase hex digits are now loaded as an `Error`,
  and `length` and `haskey` throw on malformed auxiliary data.
  See the migration guide in the documentation of XAMAuxData.jl.
- Allow MemoryViews 0.5.
