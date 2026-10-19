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
- **Breaking:** `PAFReader` now reads from a `BufferIO.AbstractBufReader`, and its type
  parameter is the type of the buffered reader. An `IO` passed to `PAFReader` is
  wrapped in a `BufferIO.BufReader`, and so becomes `PAFReader{BufReader{T}}`.
  Reading from an `AbstractBufReader` avoids copying the data into the reader's own buffer.
- **Breaking:** A `PAFRecord` may now be a view into the memory it was parsed from, instead
  of owning a copy of the data:
  - With `copy=false`, the record iterated by `PAFReader` and returned by `try_next!`
    is a view into the reader's buffer, and is invalidated when the reader is advanced.
  - `try_parse` and `parse` return a view into the input, so mutating the input
    also mutates the record.
  - `copy(record)` returns a record that owns its data.
- **Breaking:** `PAFRecord` is now immutable. With `copy=false`, iterating a `PAFReader`
  no longer returns the same object on each iteration, but a new record each time.
  On Julia versions before 1.14, returning a `PAFRecord` from a non-inlined function,
  such as `try_next!` or `iterate`, allocates.
- **Breaking:** `aux_data` now returns an immutable `Auxiliary{ImmutableMemoryView{UInt8}}`,
  so the auxiliary data of records can no longer be modified.
- **Breaking:** Reading from a closed `PAFReader` now throws, instead of
  behaving as if the reader is empty.
- **Breaking:** Accessing any property other than `qname` and `qlen` of an unmapped
  record now throws the new exception `PairwiseMappingFormat.UnmappedError`, as does calling `query_coverage`, `target_coverage`
  or `aln_identity` on it. Previously, `tname` and `is_rc` returned `nothing`, and the
  remaining properties returned arbitrary values. As a result, `tname` is now always
  a `StringView`, and `is_rc` always a `Bool`.
- Sequence lengths and positions up to `typemax(UInt32)` (about 4.3 Gbp) are now supported,
  up from `typemax(Int32)`. Larger integers are rejected with `IntegerOverflow`.
- **Breaking:** Fields 6-12 of unmapped records are now validated like those of mapped
  records, except that the integer fields may be zero. Previously, their content was ignored.

### Removed
- **Breaking:** The `PAFRecord(buffer_size::Int)` constructor.

### Fixed
- `isempty(reader)` no longer consumes a record from a `PAFReader`. It may invalidate
  the most recently returned record if that record is a view into the reader's buffer,
  i.e. if it was returned by `try_next!`, or iterated with `copy=false`.
  `PAFReader` now implements `Base.isdone`, which never invalidates records, and returns
  `missing` when it can't determine whether the reader is done without doing so.
