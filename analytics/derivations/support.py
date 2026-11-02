"""Report support aggregates derived from event lineage and exact supporters."""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path



VARIANT_STRATEGY_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "alt_support_row_count",
    "alt_support_ortholog_count",
    "human_lineage_depth",
]
EVENT_ORTHOLOG_SUPPORT_FIELDS = [
    "event_group_id",
    "ortholog_gene_id",
    "tax_id",
    "mapq",
    "native_alignment_type",
    "support_row_count",
]

@dataclass(slots=True)
class StrategySupport:
    row_count: int = 0
    ortholog_count: int = 0


class ExactSupportSpool:
    """Write exact supporters as narrow integer edges for bounded aggregation."""

    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("w", buffering=1024 * 1024)
        self.strategy_ids: dict[str, int] = {}
        self.strategy_names: list[str] = []
        self.ortholog_ids: dict[tuple[str, str], int] = {}
        self.input_edge_count = 0

    def __enter__(self) -> "ExactSupportSpool":
        return self

    def __exit__(self, *_args) -> None:
        self.close()

    def close(self) -> None:
        if not self.handle.closed:
            self.handle.close()

    def _strategy_id(self, strategy: str) -> int:
        strategy_id = self.strategy_ids.get(strategy)
        if strategy_id is None:
            strategy_id = len(self.strategy_names) + 1
            self.strategy_ids[strategy] = strategy_id
            self.strategy_names.append(strategy)
        return strategy_id

    def _ortholog_id(self, gene_id: str, ortholog_gene_id: str) -> int:
        key = (gene_id, ortholog_gene_id)
        ortholog_id = self.ortholog_ids.get(key)
        if ortholog_id is None:
            ortholog_id = len(self.ortholog_ids) + 1
            self.ortholog_ids[key] = ortholog_id
        return ortholog_id

    def add_group(
        self,
        *,
        variant_context_id: int,
        gene_id: str,
        strategy: str,
        support_rows: Iterable[Mapping[str, str]],
    ) -> None:
        if variant_context_id < 1:
            raise ValueError("Exact support requires a positive variant_context_id")
        if not gene_id:
            raise ValueError("Exact support requires gene_id")
        if not strategy:
            raise ValueError("Exact support requires one alignment strategy")
        strategy_id = self._strategy_id(strategy)
        for row in support_rows:
            ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
            if not ortholog_gene_id:
                raise ValueError("Exact event support requires ortholog_gene_id")
            support_row_count = _positive_int(
                row.get("support_row_count"),
                "exact support_row_count",
            )
            ortholog_id = self._ortholog_id(gene_id, ortholog_gene_id)
            self.handle.write(
                f"{variant_context_id}\t{strategy_id}\t{ortholog_id}\t"
                f"{support_row_count}\n"
            )
            self.input_edge_count += 1


class EventOrthologSupportStream:
    """Read one partition's exact supporters one event group at a time."""

    def __init__(self, path: Path):
        self.path = path
        self.handle = None
        self.reader: Iterator[dict[str, str]] | None = None
        self.current: dict[str, str] | None = None
        self.row_count = 0

    def __enter__(self) -> "EventOrthologSupportStream":
        self.handle = gzip.open(self.path, "rt", newline="")
        reader = csv.DictReader(self.handle, delimiter="\t")
        if reader.fieldnames != EVENT_ORTHOLOG_SUPPORT_FIELDS:
            self.handle.close()
            raise ValueError(
                "Event ortholog support table has an invalid schema: "
                f"expected {EVENT_ORTHOLOG_SUPPORT_FIELDS}, observed {reader.fieldnames}"
            )
        self.reader = iter(reader)
        self.current = next(self.reader, None)
        return self

    def __exit__(self, *_args) -> None:
        if self.handle is not None:
            self.handle.close()

    @staticmethod
    def group_id(row: Mapping[str, str]) -> int:
        try:
            group_id = int(row["event_group_id"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("event_group_id must be a positive integer") from exc
        if group_id < 1:
            raise ValueError("event_group_id must be a positive integer")
        return group_id

    def take(self, expected_group_id: int) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if self.current is None:
            return rows
        observed_group_id = self.group_id(self.current)
        if observed_group_id < expected_group_id:
            raise ValueError(
                "Event ortholog support is out of order or has no matching event: "
                f"event_group_id={observed_group_id}"
            )
        if observed_group_id > expected_group_id:
            return rows
        while self.current is not None and self.group_id(self.current) == expected_group_id:
            rows.append(self.current)
            self.row_count += 1
            if self.reader is None:
                raise RuntimeError("Event ortholog support stream is not open")
            self.current = next(self.reader, None)
        return rows

    def finish(self) -> None:
        if self.current is not None:
            raise ValueError(
                "Event ortholog support has no matching compact event: "
                f"event_group_id={self.group_id(self.current)}"
            )


def aggregate_exact_support(
    connection,
    spool: ExactSupportSpool,
    aggregates_by_id: list[dict | None],
) -> int:
    """Collapse exact edges and attach scalar per-strategy counts to variants."""

    spool.close()
    connection.execute("DROP TABLE IF EXISTS exact_support")
    try:
        if spool.input_edge_count == 0:
            if len(aggregates_by_id) > 1:
                raise ValueError("Canonical variants have no exact ortholog support")
            return 0
        connection.execute(
            f"""
            CREATE TEMP TABLE exact_support AS
            SELECT
                variant_context_id,
                strategy_id,
                ortholog_id,
                CAST(SUM(support_row_count) AS UBIGINT) AS support_row_count
            FROM read_csv(
                {_sql_string(spool.path)},
                delim = '\t',
                header = false,
                auto_detect = false,
                columns = {{
                    'variant_context_id': 'UBIGINT',
                    'strategy_id': 'UINTEGER',
                    'ortholog_id': 'UBIGINT',
                    'support_row_count': 'UBIGINT'
                }}
            )
            GROUP BY variant_context_id, strategy_id, ortholog_id
            """
        )
        cursor = connection.execute(
            """
            SELECT
                variant_context_id,
                strategy_id,
                COUNT(*) AS ortholog_count,
                SUM(support_row_count) AS row_count
            FROM exact_support
            GROUP BY variant_context_id, strategy_id
            ORDER BY variant_context_id, strategy_id
            """
        )
        while batch := cursor.fetchmany(100_000):
            for variant_context_id, strategy_id, ortholog_count, row_count in batch:
                variant_index = int(variant_context_id)
                if variant_index >= len(aggregates_by_id):
                    raise ValueError(f"Unknown variant_context_id: {variant_index}")
                aggregate = aggregates_by_id[variant_index]
                if aggregate is None:
                    raise ValueError(f"Unknown variant_context_id: {variant_index}")
                strategy_index = int(strategy_id) - 1
                if strategy_index < 0 or strategy_index >= len(spool.strategy_names):
                    raise ValueError(
                        f"Unknown exact-support strategy_id: {strategy_id}"
                    )
                strategy = spool.strategy_names[strategy_index]
                support_by_strategy = aggregate["_support_by_strategy"]
                if strategy in support_by_strategy:
                    raise ValueError(
                        "Duplicate exact-support aggregate for "
                        f"variant_context_id={variant_index}, strategy={strategy}"
                    )
                support_by_strategy[strategy] = StrategySupport(
                    row_count=int(row_count),
                    ortholog_count=int(ortholog_count),
                )

        for variant_context_id, aggregate in enumerate(aggregates_by_id[1:], start=1):
            if aggregate is None:
                raise ValueError(f"Missing variant aggregate: {variant_context_id}")
            if not aggregate["_support_by_strategy"]:
                raise ValueError(
                    "Canonical variant has no exact ortholog support: "
                    f"variant_context_id={variant_context_id}"
                )
        return int(
            connection.execute("SELECT COUNT(*) FROM exact_support").fetchone()[0]
        )
    finally:
        connection.execute("DROP TABLE IF EXISTS exact_support")
        spool.path.unlink(missing_ok=True)


def build_variant_strategy_support(aggregates: Iterable[dict]) -> tuple[list[dict[str, object]], int]:
    rows = []
    missing_key_count = 0
    for aggregate in aggregates:
        support_by_strategy = aggregate["_support_by_strategy"]
        if not aggregate.get("variant_key"):
            missing_key_count += len(support_by_strategy)
            continue
        for strategy, support in support_by_strategy.items():
            rows.append({
                "variant_key": aggregate["variant_key"],
                "gene_id": aggregate.get("gene_id", ""),
                "strategy": strategy,
                "alt_support_row_count": support.row_count,
                "alt_support_ortholog_count": support.ortholog_count,
                "human_lineage_depth": aggregate["_taxonomy_depth_by_strategy"][strategy],
            })
    rows.sort(key=lambda row: (_int_or_default(row["gene_id"], 10**18), row["variant_key"], row["strategy"]))
    return rows, missing_key_count


def _int_or_default(value: object, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _positive_int(value: object, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer: {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{label} must be a positive integer: {parsed}")
    return parsed


def _sql_string(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"
