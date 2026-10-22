from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses.basic_filtering import (
    _filter_score_connection,
    _filter_score_histograms,
    _resolve_filter_sources,
    _stream_filter_data,
    build_or_load_filter_score_store,
    candidate_curves_from_histograms,
    read_filter_score_histograms,
)


def _write_gzip(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_filter_score_store_and_candidate_curves_are_allele_level(tmp_path: Path, references) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    annotation_fields = [
        "variant_key",
        "gene_id",
        "event_type",
        "lookup_status",
        "strategies",
        "gnomad_af",
        "vep_status",
        "vep_primary_consequence",
        "vep_consequence_terms",
    ]
    _write_gzip(
        annotations,
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "lookup_status": "ok",
                "strategies": "s1,s2",
                "gnomad_af": "0.01",
                "vep_status": "ok",
                "vep_primary_consequence": "missense_variant",
                "vep_consequence_terms": "missense_variant",
            },
            {
                "variant_key": "1:10000020:C>T",
                "gene_id": "1",
                "event_type": "snv",
                "lookup_status": "ok",
                "strategies": "s1",
                "gnomad_af": "",
                "vep_status": "ok",
                "vep_primary_consequence": "missense_variant",
                "vep_consequence_terms": "missense_variant",
            },
        ],
        annotation_fields,
    )
    with gzip.open(annotations, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.append({**rows[0], "gene_id": "2"})
    _write_gzip(annotations, rows, annotation_fields)
    support = tmp_path / "variant_strategy_support.tsv.gz"
    support_fields = [
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "alt_support_family_count",
        "site_aligned_ortholog_count",
    ]
    _write_gzip(
        support,
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 3,
                "alt_support_family_count": 2,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s2",
                "alt_support_ortholog_count": 2,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:10:A>G",
                "gene_id": "2",
                "strategy": "s1",
                "alt_support_ortholog_count": 2,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 20,
            },
            {
                "variant_key": "1:10000020:C>T",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 1,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 10,
            },
        ],
        support_fields,
    )
    failures = tmp_path / "failures.tsv.gz"
    _write_gzip(
        failures,
        [],
        ["source", "scope", "chrom", "start", "end", "failure_type", "message"],
    )

    score_path, _manifest, cache_hit = build_or_load_filter_score_store(
        variant_annotations_source=annotations,
        variant_strategy_support_tsv=support,
        annotation_failures_tsv=failures,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1", "s2"],
        **references,
    )
    all_curves = candidate_curves_from_histograms(read_filter_score_histograms(score_path))
    curves = all_curves[all_curves.target_context.eq("all") & all_curves.consequence.eq("all")]

    assert cache_hit is False
    strategy_curve = curves[
        curves["strategy"].eq("s1")
        & curves["variant_type"].eq("snv")
        & curves["filter_key"].eq("strategy")
    ].set_index("threshold")
    assert strategy_curve.loc[1, "retained_variant_count"] == 2
    assert strategy_curve.loc[2, "retained_variant_count"] == 1
    assert strategy_curve.loc[2, "gnomad_found_fraction"] == 1.0
    union = curves[curves["strategy"].eq("union") & curves["filter_key"].eq("ortholog")]
    assert union.set_index("threshold").loc[1, "total_variant_count"] == 2
    assert union.set_index("threshold").loc[3, "retained_variant_count"] == 1
    assert union.set_index("threshold").loc[4, "retained_variant_count"] == 0
    minimum = curves[
        curves["strategy"].eq("s1") & curves["filter_key"].eq("aligned_min")
    ].set_index("threshold")
    maximum = curves[
        curves["strategy"].eq("s1") & curves["filter_key"].eq("aligned_max")
    ].set_index("threshold")
    assert minimum.loc[11, "retained_variant_count"] == 1
    assert maximum.loc[10, "retained_variant_count"] == 2

    source, columns = _resolve_filter_sources(annotations, support)
    with _filter_score_connection(
        source=source, support_path=support, support_columns=columns,
        annotation_failures_tsv=failures, temp_dir=tmp_path, **references,
    ) as connection:
        streamed_curves = candidate_curves_from_histograms(
            _filter_score_histograms(connection, "filter_scores")
        )
    pd.testing.assert_frame_equal(streamed_curves, all_curves)
    from analytics.io.archive_sources import prepare_parquet_dataset
    from analytics.io.variant_source import resolve_variant_table_source
    annotation_pack, support_pack = tmp_path / "annotation_pack", tmp_path / "support_pack"
    for input_path, output in ((annotations, annotation_pack), (support, support_pack)):
        prepare_parquet_dataset(
            resolve_variant_table_source(input_path, required_columns={"variant_key"}), output, input_path.name,
        )
    histograms = _stream_filter_data(
        variant_annotations_source=annotation_pack / "manifest.json",
        variant_strategy_support_tsv=support, filter_support_source=support_pack / "manifest.json",
        annotation_failures_tsv=failures, analytics_dir=tmp_path / "streamed", **references,
    )
    pd.testing.assert_frame_equal(candidate_curves_from_histograms(histograms), all_curves)


@pytest.fixture
def references(tmp_path):
    genes, features = tmp_path / "genes.tsv.gz", tmp_path / "features.tsv.gz"
    _write_gzip(genes, [{"gene_id": str(gene), "begin": 1} for gene in (1, 2)], ["gene_id", "begin"])
    _write_gzip(features, [
        {"gene_id": str(gene), "feature_type": kind, "target_start0": start, "target_end0": end}
        for gene in (1, 2)
        for kind, start, end in [("gene", 0, 20000000),
                                 ("cds", 0 if gene == 1 else 100, 100 if gene == 1 else 300),
                                 ("intron", 100 if gene == 1 else 0, 20000000 if gene == 1 else 100)]
    ], ["gene_id", "feature_type", "target_start0", "target_end0"])
    return dict(genes_tsv=genes, target_features_tsv=features)


def test_scoped_support_does_not_borrow_other_gene_or_consequence(tmp_path, references):
    annotations, support, failures = [tmp_path / name for name in ("annotations.tsv.gz", "support.tsv.gz", "failures.tsv.gz")]
    _write_gzip(annotations, [
        dict(variant_key=key, gene_id=gene, event_type="snv", lookup_status="ok", gnomad_af=af,
             vep_status=status, vep_consequence_terms=terms)
        for key, gene, af, status, terms in [
            ("1:10:A>G", "1", "0.01", "ok", "synonymous_variant|splice_region_variant"),
            ("1:10:A>G", "2", "0.01", "ok", "missense_variant"),
            ("1:200:A>G", "1", "", "ok", "synonymous_variant"),
            ("1:20:A>G", "1", "", "failed", "synonymous_variant"),
        ]
    ], ["variant_key", "gene_id", "event_type", "lookup_status", "gnomad_af", "vep_status", "vep_consequence_terms"])
    _write_gzip(support, [
        dict(variant_key=key, gene_id=gene, strategy=strategy, alt_support_ortholog_count=score,
             alt_support_family_count=0, site_aligned_ortholog_count=100)
        for key, gene, strategy, score in [
            ("1:10:A>G", "1", "s1", 2), ("1:10:A>G", "2", "s1", 50),
            ("1:10:A>G", "1", "s2", 4), ("1:200:A>G", "1", "s1", 10),
            ("1:20:A>G", "1", "s1", 1),
        ]
    ], ["variant_key", "gene_id", "strategy", "alt_support_ortholog_count", "alt_support_family_count", "site_aligned_ortholog_count"])
    _write_gzip(failures, [{"source": "gnomad", "scope": "region", "chrom": "1", "start": 20, "end": 20}],
                ["source", "scope", "chrom", "start", "end"])
    from analytics.analyses.basic_filtering import build_basic_filtering_analysis
    kwargs = dict(variant_annotations_source=annotations, variant_strategy_support_tsv=support,
                  annotation_failures_tsv=failures, analytics_dir=tmp_path / "analytics", strategies=["s1", "s2"], **references)
    analysis = build_basic_filtering_analysis(**kwargs)
    assert not analysis.cache_hit
    assert build_basic_filtering_analysis(**kwargs).cache_hit
    discarded = build_basic_filtering_analysis(**kwargs, cache_policy="discard")
    pd.testing.assert_frame_equal(discarded.candidate_curves, analysis.candidate_curves)
    assert not list((tmp_path / "analytics").rglob("*clinvar*"))
    curves = analysis.candidate_curves
    def curve(strategy="s1", context="all", consequence="synonymous", metric="ortholog"):
        return curves[curves.strategy.eq(strategy) & curves.target_context.eq(context)
                      & curves.consequence.eq(consequence) & curves.filter_key.eq(metric)].set_index("threshold")
    assert curve().loc[1, "total_variant_count"] == 2
    assert curve().loc[3, "gnomad_found_fraction"] == 0
    assert curve().loc[3, "retained_variant_count"] == 1
    assert curve(context="cds").loc[3, "retained_variant_count"] == 0
    assert curve(context="intron").loc[3, "retained_variant_count"] == 1
    assert curve(strategy="union", context="cds").loc[4, "retained_variant_count"] == 1
    assert curve(strategy="union", context="cds").loc[5, "retained_variant_count"] == 0
    assert curve(consequence="splice_region").loc[3, "retained_variant_count"] == 0
    assert curve(context="cds", consequence="missense").empty
    assert curve(context="intron", consequence="missense").loc[50, "retained_variant_count"] == 1
    assert curve(consequence="all").loc[50, "retained_variant_count"] == 1
    assert curve(consequence="all").loc[1, "total_variant_count"] == 3
    assert curve(consequence="all").loc[1, "gnomad_eligible_count"] == 2
    assert curve(consequence="not_annotated").loc[1, "gnomad_lookup_failed_count"] == 1
    assert curve(metric="family").loc[0, "retained_variant_count"] == 2
    assert curve(metric="aligned_max").loc[0, "retained_variant_count"] == 0
    # VEP failure must not leak its stale consequence into the synonymous cohort.
    assert curve(consequence="not_annotated").loc[1, "total_variant_count"] == 1
    assert pd.isna(curve(consequence="not_annotated").loc[1, "gnomad_found_fraction"])


def test_indel_context_uses_normalized_first_affected_base(tmp_path, references):
    from analytics.analyses.basic_filtering import build_basic_filtering_analysis
    annotations, support, failures = [tmp_path / name for name in ("a.tsv.gz", "s.tsv.gz", "f.tsv.gz")]
    keys = ["1:100:A>AT", "1:100:AT>A"]
    _write_gzip(annotations, [dict(variant_key=key, gene_id="1", event_type=kind,
                lookup_status="ok", gnomad_af="", vep_status="ok", vep_consequence_terms="intron_variant")
                for key, kind in zip(keys, ["ins", "del"])],
                ["variant_key", "gene_id", "event_type", "lookup_status", "gnomad_af", "vep_status", "vep_consequence_terms"])
    _write_gzip(support, [dict(variant_key=key, gene_id="1", strategy="s1", alt_support_ortholog_count=2,
                alt_support_family_count="", site_aligned_ortholog_count="") for key in keys],
                ["variant_key", "gene_id", "strategy", "alt_support_ortholog_count", "alt_support_family_count", "site_aligned_ortholog_count"])
    _write_gzip(failures, [], ["source", "scope", "chrom", "start", "end"])
    result = build_basic_filtering_analysis(variant_annotations_source=annotations,
                variant_strategy_support_tsv=support, annotation_failures_tsv=failures,
                analytics_dir=tmp_path / "analytics", strategies=["s1"], **references)
    curves = result.candidate_curves
    selected = curves[curves.strategy.eq("s1") & curves.consequence.eq("intronic")
                      & curves.filter_key.eq("ortholog") & curves.threshold.eq(1)]
    assert set(selected.variant_type) == {"indel"}
    assert set(selected.target_context) == {"all", "intron"}
    assert set(selected.retained_variant_count) == {2}
