"""Simple support-threshold analyses for candidate variant filtering."""

from __future__ import annotations

import csv
import gzip
import json
import tempfile
from collections.abc import Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic
from analytics.io.duckdb import available_cpu_count, configure_duckdb_memory
from analytics.io.variant_source import (
    resolve_variant_table_source,
    sql_string,
    variant_source_sql,
    variant_source_blocks,
)
from analytics.vep.consequences import (
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
    VALIDATION_CONSEQUENCE_TERMS,
)
from .variant_summary_aggregation import _common_prefix_sql, _register_reference_tables

CANDIDATE_CONSEQUENCE_OPTIONS = (*CONSEQUENCE_OPTIONS, ("not_annotated", "Not annotated"))


FILTER_SCORE_SCHEMA_VERSION = 3
UNION_STRATEGY = "union"
FILTER_SCORE_COLUMNS = [
    "variant_key",
    "strategy",
    "variant_type",
    "gnomad_status",
    "target_context",
    "consequence",
    "ortholog_support",
    "strategy_support",
    "family_support",
    "site_aligned_min",
    "site_aligned_max",
]
FILTER_OPTIONS = [
    ("ortholog", "Exact-ALT ortholog support", "ortholog_support"),
    ("strategy", "Strategy support", "strategy_support"),
    ("family", "Supporting families", "family_support"),
    ("aligned_max", "Site-aligned orthologs (at most)", "site_aligned_min"),
    ("aligned_min", "Site-aligned orthologs (at least)", "site_aligned_max"),
]


@dataclass(frozen=True)
class BasicFilteringAnalysis:
    score_path: Path | None
    score_manifest_path: Path | None
    candidate_curves: pd.DataFrame
    cache_hit: bool = False


def build_basic_filtering_analysis(
    *,
    variant_annotations_source: Path | Sequence[Path],
    variant_strategy_support_tsv: Path | Sequence[Path],
    annotation_failures_tsv: Path | Sequence[Path],
    analytics_dir: Path,
    genes_tsv: Path | Sequence[Path],
    target_features_tsv: Path | Sequence[Path],
    strategies: list[str],
    cache_policy: str = "keep",
    filter_support_source: Path | Sequence[Path] | None = None,
) -> BasicFilteringAnalysis:
    """Build context-specific candidate-retention and gnomAD threshold curves."""

    if cache_policy == "keep":
        score_path, manifest_path, cache_hit = build_or_load_filter_score_store(
            variant_annotations_source=variant_annotations_source,
            variant_strategy_support_tsv=variant_strategy_support_tsv,
            annotation_failures_tsv=annotation_failures_tsv,
            analytics_dir=analytics_dir,
            strategies=strategies,
            genes_tsv=genes_tsv, target_features_tsv=target_features_tsv,
        )
        histograms = read_filter_score_histograms(score_path)
    elif cache_policy == "discard":
        histograms = _stream_filter_data(
            variant_annotations_source=variant_annotations_source,
            variant_strategy_support_tsv=variant_strategy_support_tsv,
            filter_support_source=filter_support_source, annotation_failures_tsv=annotation_failures_tsv,
            analytics_dir=analytics_dir, genes_tsv=genes_tsv, target_features_tsv=target_features_tsv,
        )
        score_path = manifest_path = None
        cache_hit = False
    else:
        raise ValueError(f"Unknown cache policy: {cache_policy}")
    candidate_curves = candidate_curves_from_histograms(histograms)
    return BasicFilteringAnalysis(
        score_path,
        manifest_path,
        candidate_curves,
        cache_hit,
    )


def _stream_filter_data(*, variant_annotations_source, variant_strategy_support_tsv,
                        filter_support_source, annotation_failures_tsv, analytics_dir, genes_tsv, target_features_tsv):
    source, support_columns = _resolve_filter_sources(
        variant_annotations_source, variant_strategy_support_tsv,
    )
    prepared_support = (
        resolve_variant_table_source(filter_support_source, required_columns=set(support_columns))
        if filter_support_source is not None else None
    )
    if prepared_support is not None and source.format != "parquet":
        raise ValueError("Prepared support requires prepared variant inputs")
    blocks = variant_source_blocks(source) if prepared_support is not None else (source,)
    support_blocks = {
        block.identity["block"]: block for block in variant_source_blocks(prepared_support)
    } if prepared_support is not None else {}
    if prepared_support is not None and set(support_blocks) - {block.identity["block"] for block in blocks}:
        raise ValueError("Prepared support contains blocks absent from variant annotations")
    analytics_dir.mkdir(parents=True, exist_ok=True)
    histograms = []
    with tempfile.TemporaryDirectory(prefix=".filter_scores.", dir=analytics_dir) as temporary:
        for index, block in enumerate(blocks, 1):
            support = support_blocks.get(block.identity["block"]) if prepared_support is not None else None
            if prepared_support is not None and support is None:
                continue
            print(f"Filter summary block {index}/{len(blocks)}", flush=True)
            with _filter_score_connection(
                source=block, support_path=variant_strategy_support_tsv,
                support_columns=support_columns, annotation_failures_tsv=annotation_failures_tsv,
                temp_dir=Path(temporary), prepared_support=support,
                genes_tsv=genes_tsv, target_features_tsv=target_features_tsv,
            ) as connection:
                histograms.append(_filter_score_histograms(connection, "filter_scores"))
    if not histograms:
        return pd.DataFrame()
    histogram = pd.concat(histograms, ignore_index=True)
    columns = [column for column in histogram.columns if column != "variant_count"]
    return histogram.groupby(columns, as_index=False, sort=True)["variant_count"].sum()


def build_or_load_filter_score_store(
    *,
    variant_annotations_source: Path | Sequence[Path],
    variant_strategy_support_tsv: Path | Sequence[Path],
    annotation_failures_tsv: Path | Sequence[Path],
    analytics_dir: Path,
    strategies: list[str],
    genes_tsv: Path | Sequence[Path],
    target_features_tsv: Path | Sequence[Path],
) -> tuple[Path, Path, bool]:
    """Materialize one bounded-width allele/strategy table for all filter views."""

    source, support_columns = _resolve_filter_sources(
        variant_annotations_source, variant_strategy_support_tsv,
    )
    support_paths = _paths(variant_strategy_support_tsv)
    outdir = analytics_dir / "basic_filtering"
    score_path = outdir / "filter_scores.parquet"
    manifest_path = outdir / "manifest.json"
    expected_inputs = {
        "schema_version": FILTER_SCORE_SCHEMA_VERSION,
        "variant_source": source.identity,
        "strategy_support": [path_metadata(path) for path in support_paths],
        "annotation_failures": [
            path_metadata(path) for path in _paths(annotation_failures_tsv)
        ],
        "strategies": sorted(str(value) for value in strategies),
        "genes": [path_metadata(path) for path in _paths(genes_tsv)],
        "target_features": [path_metadata(path) for path in _paths(target_features_tsv)],
    }
    if score_path.exists() and manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            if (
                manifest.get("status") == "complete"
                and manifest.get("inputs") == expected_inputs
                and manifest.get("output") == path_metadata(score_path)
            ):
                _validate_score_store(score_path, int(manifest["row_count"]))
                return score_path, manifest_path, True
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
            pass

    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".filter_scores.", dir=outdir) as temporary:
        row_count = _build_filter_score_store(
            source=source,
            support_path=variant_strategy_support_tsv,
            support_columns=support_columns,
            annotation_failures_tsv=annotation_failures_tsv,
            score_path=score_path,
            temp_dir=Path(temporary),
            genes_tsv=genes_tsv, target_features_tsv=target_features_tsv,
        )
    manifest = {
        "status": "complete",
        "inputs": expected_inputs,
        "row_count": row_count,
        "output": path_metadata(score_path),
    }
    write_json_atomic(manifest_path, manifest)
    return score_path, manifest_path, False


def _resolve_filter_sources(variant_annotations_source, variant_strategy_support_tsv):
    source = resolve_variant_table_source(
        variant_annotations_source,
        required_columns={"variant_key", "gene_id", "event_type", "lookup_status", "gnomad_af",
                          "vep_status", "vep_consequence_terms"},
    )
    support_paths = _paths(variant_strategy_support_tsv)
    support_columns = _read_header(support_paths[0])
    if any(_read_header(path) != support_columns for path in support_paths[1:]):
        raise ValueError("Variant strategy support columns differ across source runs")
    required_support = {
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "alt_support_family_count",
        "site_aligned_ortholog_count",
    }
    missing = required_support - set(support_columns)
    if missing:
        raise ValueError(
            f"Variant strategy support {variant_strategy_support_tsv} missing columns: "
            + ", ".join(sorted(missing))
        )
    return source, support_columns


def _build_filter_score_store(
    *,
    source,
    support_path: Path | Sequence[Path],
    support_columns: list[str],
    annotation_failures_tsv: Path | Sequence[Path],
    score_path: Path,
    temp_dir: Path,
    genes_tsv: Path | Sequence[Path],
    target_features_tsv: Path | Sequence[Path],
) -> int:
    with tempfile.NamedTemporaryFile(
        dir=score_path.parent,
        prefix=f".{score_path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)
    temporary_path.unlink(missing_ok=True)
    try:
        with _filter_score_connection(
            source=source, support_path=support_path, support_columns=support_columns,
            annotation_failures_tsv=annotation_failures_tsv, temp_dir=temp_dir,
            genes_tsv=genes_tsv, target_features_tsv=target_features_tsv,
        ) as connection:
            row_count = int(connection.execute("SELECT count(*) FROM filter_scores").fetchone()[0])
            connection.execute(
                f"COPY filter_scores TO {sql_string(temporary_path)} "
                "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
            )
        _validate_score_store(temporary_path, row_count)
        temporary_path.chmod(0o644)
        temporary_path.replace(score_path)
        return row_count
    finally:
        temporary_path.unlink(missing_ok=True)


@contextmanager
def _filter_score_connection(*, source, support_path, support_columns, annotation_failures_tsv,
                             temp_dir, genes_tsv, target_features_tsv, prepared_support=None):
    repaired_status = "gnomad_lookup_status" if "gnomad_lookup_status" in source.columns else "''"
    thread_count = available_cpu_count()
    with duckdb.connect() as connection:
        connection.execute(f"SET threads={thread_count}")
        connection.execute("SET preserve_insertion_order=false")
        connection.execute(f"SET temp_directory={sql_string(temp_dir)}")
        configure_duckdb_memory(connection, thread_count)
        _register_reference_tables(
            connection, genes_path=genes_tsv, target_features_path=target_features_tsv,
            annotation_failures_path=annotation_failures_tsv,
        )
        connection.execute(
            f"CREATE VIEW variant_source_rows AS SELECT * FROM {variant_source_sql(source)}"
        )
        support_sql = (
            variant_source_sql(prepared_support) if prepared_support is not None
            else _support_source_sql(support_path, support_columns)
        )
        connection.execute(f"CREATE VIEW strategy_support_rows AS SELECT * FROM {support_sql}")

        failed = (
            "EXISTS (SELECT 1 FROM gnomad_failures f WHERE f.chrom = key_chrom "
            "AND key_pos BETWEEN f.start1 AND f.end1)"
        )
        connection.execute(
            "CREATE TEMP TABLE global_annotations AS WITH parsed AS ("
            "SELECT variant_key, lower(event_type) AS event_type, lookup_status, gnomad_af, "
            f"{repaired_status} AS gnomad_lookup_status, "
            "regexp_replace(split_part(variant_key, ':', 1), '^chr', '') AS key_chrom, "
            "try_cast(split_part(variant_key, ':', 2) AS BIGINT) AS key_pos "
            "FROM variant_source_rows WHERE variant_key <> ''"
            "), classified AS (SELECT *, "
            "CASE WHEN gnomad_lookup_status IN ('found','not_found') THEN gnomad_lookup_status "
            "WHEN try_cast(nullif(gnomad_af, '') AS DOUBLE) IS NOT NULL THEN 'found' "
            "WHEN coalesce(lookup_status, '') NOT IN ('', 'ok') OR key_pos IS NULL "
            "THEN 'lookup_failed' "
            f"WHEN {failed} THEN 'lookup_failed' ELSE 'not_found' END AS row_gnomad_status "
            "FROM parsed) SELECT variant_key, "
            "CASE WHEN bool_or(event_type = 'snv') THEN 'snv' ELSE 'indel' END AS variant_type, "
            "CASE WHEN bool_or(row_gnomad_status = 'found') THEN 'found' "
            "WHEN bool_or(row_gnomad_status = 'lookup_failed') THEN 'lookup_failed' "
            "ELSE 'not_found' END AS gnomad_status "
            "FROM classified GROUP BY variant_key"
        )
        connection.execute(
            "CREATE TEMP VIEW typed_support AS SELECT variant_key, gene_id, trim(strategy) AS strategy, "
            "try_cast(alt_support_ortholog_count AS BIGINT) AS ortholog_support, "
            "try_cast(nullif(alt_support_family_count, '') AS BIGINT) AS family_support, "
            "try_cast(nullif(site_aligned_ortholog_count, '') AS BIGINT) AS site_aligned "
            "FROM strategy_support_rows"
        )
        invalid = int(
            connection.execute(
                "SELECT count(*) FROM typed_support s JOIN global_annotations a USING (variant_key) "
                "WHERE coalesce(s.strategy, '') IN ('', 'union') OR ortholog_support IS NULL OR ortholog_support < 1 "
                "OR family_support < 0 OR family_support > ortholog_support "
                "OR (variant_type = 'snv' AND (family_support IS NULL OR site_aligned IS NULL "
                "OR site_aligned < ortholog_support))"
            ).fetchone()[0]
        )
        if invalid:
            raise ValueError(f"Variant strategy support contains {invalid} invalid filter scores")
        missing_annotations = connection.execute(
            "SELECT count(*) FROM typed_support s LEFT JOIN variant_source_rows a "
            "USING (variant_key, gene_id) WHERE a.variant_key IS NULL"
        ).fetchone()[0]
        if missing_annotations:
            raise ValueError(f"Variant strategy support contains {missing_annotations} contexts absent from annotations")
        missing_genes = connection.execute(
            "SELECT count(*) FROM variant_source_rows a LEFT JOIN gene_begins g USING (gene_id) "
            "WHERE g.gene_begin IS NULL"
        ).fetchone()[0]
        if missing_genes:
            raise ValueError(f"Variant annotations contain {missing_genes} contexts absent from target genes")
        ref = "split_part(split_part(a.variant_key, ':', 3), '>', 1)"
        alt = "split_part(split_part(a.variant_key, ':', 3), '>', 2)"
        position = f"try_cast(split_part(a.variant_key, ':', 2) AS BIGINT) - g.gene_begin + ({_common_prefix_sql(ref, alt)})"
        terms = "string_split(coalesce(a.vep_consequence_terms, ''), '|')"
        memberships = [
            f"({sql_string(key)}, a.vep_status = 'ok' AND list_has_any({terms}, "
            f"[{','.join(sql_string(term) for term in sorted(accepted))}]))"
            for key, accepted in VALIDATION_CONSEQUENCE_TERMS.items()
        ]
        recognized = set().union(*VALIDATION_CONSEQUENCE_TERMS.values())
        recognized_sql = '[' + ','.join(sql_string(term) for term in sorted(recognized)) + ']'
        annotated = "(a.vep_status = 'ok' AND coalesce(a.vep_consequence_terms, '') <> '')"
        memberships += [
            "('all', true)",
            f"('other', {annotated} AND len(list_filter({terms}, t -> NOT list_contains({recognized_sql}, t))) > 0)",
            f"('not_annotated', NOT coalesce({annotated}, false))",
        ]
        # Scope before reducing support: a synonymous call cannot borrow depth from
        # an intronic or different-gene row of the same genomic allele.
        connection.execute(
            "CREATE TEMP TABLE scoped_support AS SELECT s.*, c.target_context, m.consequence "
            "FROM typed_support s JOIN variant_source_rows a USING (variant_key, gene_id) "
            "JOIN gene_begins g ON g.gene_id = a.gene_id "
            f"LEFT JOIN target_contexts c ON c.gene_id = a.gene_id AND ({position}) >= c.start0 "
            f"AND ({position}) < c.end0, LATERAL (VALUES {','.join(memberships)}) m(consequence, included) "
            "WHERE m.included"
        )
        connection.execute(
            "CREATE TEMP TABLE compact_support AS SELECT variant_key, strategy, "
            "coalesce(scope.target_context, 'other') AS target_context, consequence, "
            "max(ortholog_support) AS ortholog_support, max(family_support) AS family_support, "
            "min(site_aligned) AS site_aligned_min, max(site_aligned) AS site_aligned_max "
            "FROM scoped_support, LATERAL (VALUES ('all'), (target_context)) scope(target_context) "
            "GROUP BY variant_key, strategy, scope.target_context, consequence"
        )
        connection.execute(
            "CREATE TEMP TABLE filter_scores AS WITH strategy_counts AS ("
            "SELECT variant_key, count(DISTINCT strategy) AS strategy_support FROM typed_support GROUP BY variant_key) "
            "SELECT s.variant_key, s.strategy, a.variant_type, a.gnomad_status, "
            "s.target_context, s.consequence, s.ortholog_support, n.strategy_support, "
            "s.family_support, s.site_aligned_min, s.site_aligned_max "
            "FROM compact_support s JOIN global_annotations a USING (variant_key) "
            "JOIN strategy_counts n USING (variant_key)"
        )
        yield connection


def read_filter_score_histograms(score_path: Path) -> pd.DataFrame:
    with duckdb.connect() as connection:
        thread_count = available_cpu_count()
        connection.execute(f"SET threads={thread_count}")
        configure_duckdb_memory(connection, thread_count)
        return _filter_score_histograms(connection, f"read_parquet({sql_string(score_path)})")


def _filter_score_histograms(connection, relation: str) -> pd.DataFrame:
    # Expand metrics only after deduplicating the union, so no allele is counted twice.
    metrics = ", ".join(f"({sql_string(key)}, {column})" for key, _label, column in FILTER_OPTIONS)
    return connection.execute(
            f"WITH individual AS MATERIALIZED (SELECT * FROM {relation}), "
            "combined AS (SELECT variant_key, 'union' AS strategy, variant_type, gnomad_status, target_context, consequence, "
            "max(ortholog_support) AS ortholog_support, max(strategy_support) AS strategy_support, "
            "max(family_support) AS family_support, min(site_aligned_min) AS site_aligned_min, "
            "max(site_aligned_max) AS site_aligned_max FROM individual "
            "GROUP BY variant_key, variant_type, gnomad_status, target_context, consequence), "
            "scores AS (SELECT * FROM individual UNION ALL SELECT * FROM combined) "
            "SELECT strategy, variant_type, target_context, consequence, filter_key, score, gnomad_status, count(*) AS variant_count "
            f"FROM scores, LATERAL (VALUES {metrics}) AS metric(filter_key, score) "
            "WHERE score IS NOT NULL GROUP BY ALL",
        ).fetchdf()


def candidate_curves_from_histograms(histograms: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "strategy",
        "variant_type",
        "filter_key",
        "target_context",
        "consequence",
        "threshold",
        "retained_variant_count",
        "total_variant_count",
        "retained_fraction",
        "gnomad_found_count",
        "gnomad_eligible_count",
        "gnomad_lookup_failed_count",
        "gnomad_found_fraction",
    ]
    if histograms.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, object]] = []
    group_columns = ["strategy", "variant_type", "filter_key", "target_context", "consequence"]
    for keys, subset in histograms.groupby(group_columns, sort=True):
        strategy, variant_type, filter_key, target_context, consequence = keys
        subset = subset.copy()
        subset["score"] = pd.to_numeric(subset["score"], errors="raise").astype(int)
        subset["variant_count"] = pd.to_numeric(
            subset["variant_count"], errors="raise"
        ).astype(int)
        max_score = int(subset["score"].max())
        totals = {
            status: np.bincount(
                subset.loc[subset["gnomad_status"].eq(status), "score"],
                weights=subset.loc[subset["gnomad_status"].eq(status), "variant_count"],
                minlength=max_score + 2,
            )
            for status in ("found", "not_found", "lookup_failed")
        }
        cumulative = {
            status: (
                np.cumsum(values) if filter_key == "aligned_max" else np.cumsum(values[::-1])[::-1]
            )
            for status, values in totals.items()
        }
        total = int(subset["variant_count"].sum())
        for threshold in range(0 if filter_key in {"family", "aligned_max"} else 1, max_score + 2):
            found = int(cumulative["found"][threshold])
            not_found = int(cumulative["not_found"][threshold])
            failed = int(cumulative["lookup_failed"][threshold])
            retained = found + not_found + failed
            eligible = found + not_found
            rows.append(
                {
                    "strategy": str(strategy),
                    "variant_type": str(variant_type),
                    "filter_key": str(filter_key),
                    "target_context": str(target_context),
                    "consequence": str(consequence),
                    "threshold": threshold,
                    "retained_variant_count": retained,
                    "total_variant_count": total,
                    "retained_fraction": retained / total if total else float("nan"),
                    "gnomad_found_count": found,
                    "gnomad_eligible_count": eligible,
                    "gnomad_lookup_failed_count": failed,
                    "gnomad_found_fraction": found / eligible if eligible else float("nan"),
                }
            )
    return pd.DataFrame(rows, columns=columns)


def _support_source_sql(
    path: Path | Sequence[Path],
    columns: list[str],
) -> str:
    schema = "{" + ",".join(f"{sql_string(column)}:'VARCHAR'" for column in columns) + "}"
    paths = "[" + ",".join(sql_string(item) for item in _paths(path)) + "]"
    return (
        f"read_csv({paths}, delim='\\t', header=true, columns={schema}, "
        "auto_detect=false, compression='auto', parallel=true, "
        "nullstr='__GAPH_NULL_SENTINEL__')"
    )


def _read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def _paths(value: Path | Sequence[Path]) -> tuple[Path, ...]:
    paths = (value,) if isinstance(value, Path) else tuple(value)
    if not paths:
        raise ValueError("At least one source table is required")
    return paths


def _validate_score_store(path: Path, expected_rows: int) -> None:
    with duckdb.connect() as connection:
        thread_count = available_cpu_count()
        connection.execute(f"SET threads={thread_count}")
        configure_duckdb_memory(connection, thread_count)
        relation = connection.read_parquet(str(path))
        if relation.columns != FILTER_SCORE_COLUMNS:
            raise ValueError(
                f"Basic-filter score columns changed in {path}: "
                + ", ".join(relation.columns)
            )
        observed = int(relation.count("*").fetchone()[0])
        if observed != expected_rows:
            raise ValueError(
                f"Basic-filter score row count changed: {observed} != {expected_rows}"
            )
