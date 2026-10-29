# Analytics Contract

Archived gnomAD lookup gaps can be supplemented by an explicit
`--gnomad-reannotation` snapshot. Version selection, source binding, cache
isolation, and coverage semantics are owned by
[the reannotation contract](../reannotation/README.md). Analytics remains a
read-only consumer and never acquires provider data to fill those gaps.

Read this document when changing analytics input compatibility, cache identity,
scientific derivations, or report interpretation. Use
`report_generation.md` only for the operator launch procedure.

## Ownership And Input Boundary

`analytics/` is a read-only consumer of one or more completed GAPH runs. It
does not repair pipeline evidence, fetch missing pipeline-owned annotations, or
write below a source run.

Each source run must contain:

- canonical Stage 1 target, selected-ortholog, feature, and taxonomy evidence;
- partitioned Stage 2 summaries, segments, compact events, and exact support;
- the Stage 3 partitioned ClinVar/gnomAD/VEP variant dataset and
  event-to-variant map;
- a complete root run manifest and its bound evidence inventory.

Multi-run analyses require the same current pipeline contracts, target and
reference identities, alignment strategy configuration, ClinVar and gnomAD
contracts, VEP backend/release, and variant columns. Accepted Gene IDs must be
disjoint. Missing evidence, incompatible runs, overlap, obsolete schemas, or an
incomplete shard set fail explicitly.

Alignment manifest Gene IDs must exactly match the query genes represented in
`fetch/orthologs.selected.tsv.gz`, all of which must belong to the fetched target
genes. A fetched target without selected orthologs remains in the target cohort
but has no alignment result. Combined alignment counts use the aligned gene
union, not the fetched target count.

There is no fallback to old pipeline aggregates, a synthetic cohort, or a
separate bulk-VEP artifact.

## Immutability And Verification

A completed source run is immutable. Before using a source cache, analytics:

1. validates the root manifest and its binding to `evidence_inventory.json`;
2. verifies every inventoried evidence file by content;
3. records a verification marker only under the external analytics root.

Later runs may reuse that marker while the exact file set and filesystem
metadata remain unchanged. A metadata change triggers full verification; a
content or membership mismatch fails.

Large variant shards stay in source runs and are queried as one virtual table.
Analytics materializes only reusable derivations and genuine run-set results.
It does not concatenate source evidence, create target-FASTA symlink trees, or
modify source metadata.

An optional Google Drive archive provider restores and verifies one complete
source at a time. Its prepared-source receipt binds the original root manifest
and evidence identity to hashed analytics outputs. Prepared sources retain the
same `source-id`; archive location, transport, and cache retention do not alter
scientific identity. This is a derived analytics input, not a synthetic pipeline
run or an alternative pipeline evidence schema.

Only the archive provider's own restored directories can be removed during
preparation. Local `--run-dir` trees remain read-only. `keep` retains reusable
prepared data and calculations. `discard` places new reusable data in a
separate work namespace, retains it on failure for resume, and removes it after
success. Existing persistent caches are not removed by this policy.

Candidate conservation operates on prepared genomic blocks with all contexts
of an allele together. Score frequencies are merged globally before exact
linear quantiles, Freedman-Diaconis histogram edges, and box summaries are
calculated. Medians and quantiles are not averaged across source runs or blocks.
Variant summaries merge allele-disjoint blocks with exact allele memberships
and AF frequencies. In discard mode, filtering consumes block-local relations
and retains only context/consequence histograms, not the full filter-score store.
Ordinary ClinVar validation stores only universe allele memberships; enabling
the target-space null still requires the full observed-variant store.

## Workspace And Identity

Every command requires an external `--analytics-root`:

```text
<analytics-root>/
  .strategy_report.lock
  cache/<source-id>/
    evidence_inventory.verified.json
    calculations/<calculation-id>/
      alignment_aggregates/
      annotation_support/
      taxonomy_summary/
  cache/calculations/<calculation-id>/clinvar_conditions/
  cache/reference_identity/
  analyses/<analysis-id>/
    manifest.json
    derived/
    reports/
    performance/
  slurm/
```

- `source-id` depends only on completed pipeline provenance.
- `calculation-id` depends on explicit algorithm/cache versions and the
  Python, R, and computation-library versions actually loaded.
- `analysis-id` depends on the unordered source-ID set, scientific options,
  and calculation identity.
- A report name selects an HTML file within one analysis; it does not create a
  second scientific cache.

Every scientific calculation change must bump its owning cache version, or the
top-level analytics semantics version when no narrower owner exists.
Documentation and HTML-only changes do not require a scientific cache bump.

One process owns an analytics root at a time. A concurrent writer exits rather
than racing shared caches. The lock file is only a rendezvous point; process
exit releases ownership even if the file remains.

A local phyloP BigWig is identified by streamed SHA-256 content. Its cached
filesystem identity avoids repeated hashing while the file is unchanged.
Moving identical content preserves analysis identity; replacing it invalidates
affected artifacts.

## Reusable Derivations

Analytics derives and fingerprints:

- strategy summary and feature coverage from Stage 2 summaries and segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus Stage 3 event lineage;
- `analyses/<analysis-id>/derived/pathogenic_clinvar_hits.tsv.gz`, the complete
  pathogenic/likely-pathogenic allele detail and support table;
- run-set scientific tables used by report sections.

Completed continuous Firth result and distribution tables are cached by their
exact cohort, observed memberships, eligible genes, strategy order, model
definition, and calculation runtime. Rebuilding presentation or retrying a
later failed report stage therefore does not fit identical models again.

Annotation-support preparation is partition-local. Each worker reads only the
variant shards for its evidence partition, uses one DuckDB thread, and
atomically checkpoints its small derived tables. An identical rerun reuses
complete checkpoints and rebuilds missing ones. Partial checkpoint directories
are invalid and completed checkpoints are removed after the final per-source
cache is assembled.

Pipeline manifests may expose integrity counts, but analytics derives scientific
counts from row-level evidence. Report-specific thresholds, bins, taxonomic
scopes, ranks, and plotting tables never become pipeline source data.

## External Evidence Semantics

ClinVar and gnomAD record fields are allele-level. Repeated gene-context rows
are reconciled by canonical `variant_key`: a successful non-empty value is
shared, while conflicting non-empty values or a non-numeric gnomAD allele
frequency fail. The pipeline publishes AF with six significant digits;
supplements retain full precision. One precise AF and its exact six-digit
round-trip are compatible, and reconciliation retains that precise value.
Distinct precise values still conflict even if they round alike; there is no
arbitrary tolerance or choice between them. This reconciliation reads existing
prepared packs unchanged. Previously successful reconciliations are unchanged,
so their cache identities remain valid. VEP consequences, target context, and
lookup outcome remain gene-context evidence.

The indexed ClinVar VCF used for validation must match the portable content
identity recorded by every source run. It may be located at a different path
from the pipeline launch.

Rows with non-`ok` pipeline VEP status remain explicit and appear as
`Not annotated` where a complete denominator is required. The optional
target-space null may issue new VEP and gnomAD requests for generated controls;
it uses the pipeline-pinned VEP release, and its backend/release are part of the
cache identity. Local VEP is probed before an analysis workspace is created.

## Report Semantics

### Strategy Concordance

Candidate Profile retains the pairwise Jaccard matrix and adds one interactive
area-proportional two-set diagram. A selection compares two strategies or one
strategy with the union of all other strategies. The matrix can select a pair.
Both views count distinct normalized alleles across gene contexts and include
all candidate variant types. “Only” refers to the selected pair, not necessarily
to all strategies. Directional coverage is shared / opposite-set size; an empty
denominator is undefined and displayed as `n/a`.

The diagram reuses summary counts: for strategy A, its intersection with the
other-strategy union is |A| minus its exclusive count, and that union's size is
the complete union minus A's exclusive count. Pairwise unions must never be
summed to derive the other-strategy union. No new evidence scan or scientific
cache version is required for this presentation.

### Firth Model Diagnostics

Each continuous model retains its R warnings in the result reason. Models with
nonconvergence warnings or invalid effects, confidence limits, or p-values are
`not_estimable`; their inferential values are missing and do not enter BH
correction. Other warnings remain visible as `estimated_warning` alongside the
usable estimate. Missing or duplicate model results are an error rather than
silently leaving a pending result.

### Pathogenic ClinVar Hits

Pathogenic and likely-pathogenic subtypes remain distinct in detail tables but
are combined in plots. `By strategy` preserves separate strategy groups.
`All strategies · Union` counts each normalized allele once and offers two
independent minimum-support filters, initially 6 orthologs and 2 strategies.
Enabled filters combine with AND; variant type applies to every chart and QC.
Ortholog support is the maximum in one gene–strategy context, never a sum.
Strategy support counts distinct calling strategies. The joint rule requires
six orthologs in at least one context and calls from at least two strategies,
not six orthologs in each strategy. These are exploratory stability filters,
not independent biological replications or calibrated benignity thresholds.

The compact P/LP payload retains allele–gene consequence rows with their own
strategy memberships and allele–condition memberships. Client-side filtering
updates counts, consequences, support and conditions together, without reading
candidate archives. Consequences count allele–gene observations, including in
Union; other counts use distinct alleles. The summary cache version changes
when this retained relation changes; incompatible cached summaries are rebuilt.

QC owns the complete paginated `Variant details` table. Its selection is shared
with this tab, with a visible filter description, reset and return link. The
TSV download always contains the complete unfiltered P/LP set. Alleles without
named conditions remain in denominators; an empty selection has undefined
fractions, not zero fractions.

Exact-support distributions count one SNV per strategy. When an allele has
several target-gene contexts, the row with maximum ALT-supporting ortholog count
is selected and Gene ID breaks ties. Violin densities use log10 counts while
axes and individual-point hover show original counts; sparse distributions use
points and a box.
To keep large cohorts legible, individual points are shown only for groups of at
most 250 observations; the violin, box, mean line, and group size still reflect
every observation.
phyloP availability does not restrict this tab.

Conditions compare unique pathogenic/likely-pathogenic alleles with either
ClinVar alleles in eligible target genes or the complete pinned GRCh38 VCF.
Union uses the union of the original eligible genes; support filters never
shrink this background to the genes of retained hits.
Both arms use the selected variant type. Denominators include alleles without a
named condition, and one allele may contribute several conditions, so displayed
fractions need not sum to one. Disease identity prefers MedGen, then MONDO, then
OMIM; names are used only when identifiers are absent.
The paired condition chart can order its displayed conditions by the GAPH
fraction, the selected ClinVar-background fraction, or the absolute difference
between them.

The whole-VCF background is streamed locally and cached per calculation
identity. It counts distinct VCF alleles, merges repeated records, and excludes
mixed benign/likely-benign and pathogenic/likely-pathogenic classifications.
It describes variants and associated conditions in the VCF, not all clinical
assertion records.

### Basic Filtering

Variant type, target context, candidate VEP consequence, filter, and strategy
selection apply to both retention and gnomAD. Basic Filtering is descriptive:
it does not calculate ClinVar associations, odds ratios, or threshold p-values.

`Compare all strategies` uses each strategy's own candidate denominator.
`Union (any strategy)` counts each normalized allele once and accepts it when
at least one calling strategy/context passes.

Context is assigned from canonical target features at the first affected
normalized allele position, using the same context priority as the variant
summary. Consequence uses all `vep_consequence_terms` of a successful candidate
VEP annotation, not ClinVar molecular-consequence terms. Memberships can
overlap: synonymous includes `synonymous_variant` and `stop_retained_variant`
and does not exclude splice-region or other transcript effects. Failed or empty
VEP annotations remain in `All consequences` and `Not annotated`.
Both context and consequence must match the same gene-context row before its
support is eligible. Only then are support scores reduced per allele/strategy
and union, preventing support from another gene or consequence from passing a
selected threshold. Strategy support counts all calling strategies for the
allele; it is not recomputed within a consequence selection.

Exact ALT support and strategy support use minimum thresholds. Supporting
families count distinct known family IDs; an unresolved family is not a new
family. Site-aligned ortholog filters support both maximum and minimum
thresholds. Family and site-depth filters are SNV-only. Within an eligible
allele/strategy selection, maximum support is used for minimum thresholds and
minimum site depth for maximum thresholds. Zero known families is a valid
called-SNV score. Invalid required support fails validation.

Retention's fixed denominator is all candidate alleles in the selected cohort
for that strategy. Its numerator includes failed gnomAD lookups. The gnomAD
fraction is found / (found + not found) among retained alleles; failed lookups
are excluded and shown separately on hover. An empty eligible denominator
produces a missing fraction, not zero. Presence in gnomAD is population overlap,
not a clinical benignity label or a calibrated measure of calling accuracy.

Both plots default to a logarithmic threshold axis on `N + 1`, keeping threshold
zero visible. Tick labels and hover report the original threshold `N`; a shared
linear-scale option is available. The retained fraction and gnomAD fraction
remain on linear percentage axes.

## Code Owners And Checks

- Input resolution, evidence verification, workspace identity, and artifact
  contracts: `analytics/io/`
- Evidence-derived reusable tables: `analytics/derivations/`,
  `analytics/io/alignment_aggregates.py`, and
  `analytics/io/annotation_support.py`
- Scientific calculations: `analytics/analyses/`
- Presentation only: `analytics/reporting/`
- Orchestration and CLI: `analytics/strategy_report.py`

Run the test matching the changed owner. Cross-boundary changes normally need
`tests/test_analysis_inputs.py`, the relevant cache/analysis test, and
`tests/test_strategy_report.py`. The full environment command is documented
in `run_validation.md`.
