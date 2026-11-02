"""Taxonomic scopes and grouping semantics for ortholog evidence."""

from __future__ import annotations

import csv
import gzip
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from genomics.taxonomy import TAXONOMY_FIELDS


# NCBI taxon 9606 parents, from Eukaryota to Homo (2026-11-02).
# https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/9606/dataset_report
# The target is fixed; memberships always use the archived Stage 1 lineage.
HUMAN_LINEAGE = (
    ("", "All selected"),
    ("2759", "Eukaryota"), ("33154", "Opisthokonta"),
    ("33208", "Metazoa"), ("6072", "Eumetazoa"),
    ("33213", "Bilateria"), ("33511", "Deuterostomia"),
    ("7711", "Chordata"), ("89593", "Craniata"),
    ("7742", "Vertebrata"), ("7776", "Gnathostomata"),
    ("117570", "Teleostomi"), ("117571", "Euteleostomi"),
    ("8287", "Sarcopterygii"), ("1338369", "Dipnotetrapodomorpha"),
    ("32523", "Tetrapoda"), ("32524", "Amniota"),
    ("40674", "Mammalia"), ("32525", "Theria"),
    ("9347", "Eutheria"), ("1437010", "Boreoeutheria"),
    ("314146", "Euarchontoglires"), ("9443", "Primates"),
    ("376913", "Haplorrhini"), ("314293", "Simiiformes"),
    ("9526", "Catarrhini"), ("314295", "Hominoidea"),
    ("9604", "Hominidae"), ("207598", "Homininae"), ("9605", "Homo"),
)

SCOPE_ANCESTORS = {
    "all": "",
    "eukaryota": "2759",
    "metazoa": "33208",
    "vertebrata": "7742",
    "tetrapoda": "32523",
    "amniota": "32524",
    "mammalia": "40674",
    "primates": "9443",
}
SCOPE_ORDER = tuple(SCOPE_ANCESTORS)
_SCOPE_ANCESTOR_IDS = tuple(SCOPE_ANCESTORS.values())
UNIT_ORDER = ("ortholog", "species", "genus", "family", "order")
COUNT_KEYS = tuple(
    f"{scope}__{unit}"
    for scope in SCOPE_ORDER
    for unit in UNIT_ORDER
)
TAXONOMY_SUMMARY_FIELDS = [
    "taxonomic_scope",
    "evidence_unit",
    "gene_count",
    "ortholog_count",
    "taxon_count",
    "unit_count",
    "orthologs_per_gene_min",
    "orthologs_per_gene_median",
    "orthologs_per_gene_mean",
    "orthologs_per_gene_max",
    "units_per_gene_min",
    "units_per_gene_median",
    "units_per_gene_mean",
    "units_per_gene_max",
]


@dataclass(frozen=True)
class TaxonomyProfile:
    tax_id: str
    ancestor_ids: frozenset[str]
    species_id: str
    genus_id: str
    family_id: str
    order_id: str
    scope_mask: int = field(init=False, repr=False)
    rank_unit_ids: tuple[object, ...] = field(init=False, repr=False)
    human_lineage_depth: int = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "human_lineage_depth", max(
            (index for index, (tax_id, _name) in enumerate(HUMAN_LINEAGE)
             if tax_id in self.ancestor_ids), default=0,
        ))
        scope_mask = 1
        for scope_index, ancestor_id in enumerate(
            _SCOPE_ANCESTOR_IDS[1:],
            start=1,
        ):
            if ancestor_id in self.ancestor_ids:
                scope_mask |= 1 << scope_index
        fallback_id = ("taxon", self.tax_id)
        object.__setattr__(self, "scope_mask", scope_mask)
        object.__setattr__(
            self,
            "rank_unit_ids",
            (
                self.species_id or fallback_id,
                self.genus_id or fallback_id,
                self.family_id or fallback_id,
                self.order_id or fallback_id,
            ),
        )

    def scopes(self) -> tuple[str, ...]:
        return tuple(
            scope
            for scope_index, scope in enumerate(SCOPE_ORDER)
            if self.scope_mask & (1 << scope_index)
        )

    def unit_id(self, unit: str, ortholog_gene_id: str) -> str:
        if unit == "ortholog":
            return f"ortholog:{ortholog_gene_id}"
        value = getattr(self, f"{unit}_id")
        return f"{unit}:{value}" if value else f"{unit}:taxon:{self.tax_id}"


def _split_ids(value: str) -> frozenset[str]:
    return frozenset(
        item.strip()
        for item in value.replace(";", ",").split(",")
        if item.strip()
    )


def load_taxonomy_profiles(path: Path) -> dict[str, TaxonomyProfile]:
    profiles: dict[str, TaxonomyProfile] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if fieldnames != TAXONOMY_FIELDS:
            raise ValueError(
                f"Taxonomy table {path} must use exact canonical fields: "
                + ", ".join(TAXONOMY_FIELDS)
            )
        for row_number, row in enumerate(reader, start=2):
            tax_id = str(row.get("tax_id") or "").strip()
            if not tax_id:
                raise ValueError(f"Taxonomy table {path} has empty tax_id at row {row_number}")
            if tax_id in profiles:
                raise ValueError(f"Duplicate taxonomy tax_id in {path}: {tax_id}")
            status = str(row.get("taxonomy_status") or "").strip()
            if status not in {"resolved", "not_returned"}:
                raise ValueError(
                    f"Taxonomy table {path} has invalid taxonomy_status for tax_id {tax_id}: "
                    f"{status!r}"
                )
            lineage_value = str(row.get("lineage_tax_ids") or "")
            profiles[tax_id] = TaxonomyProfile(
                tax_id=tax_id,
                ancestor_ids=frozenset(
                    {*_split_ids(lineage_value), tax_id}
                ),
                species_id=str(row.get("species_id") or ""),
                genus_id=str(row.get("genus_id") or ""),
                family_id=str(row.get("family_id") or ""),
                order_id=str(row.get("order_id") or ""),
            )
    return profiles


def lineage_ortholog_counts(
    taxonomy_paths: Sequence[Path], ortholog_paths: Sequence[Path],
) -> list[dict[str, object]]:
    """Count unique selected NCBI ortholog Gene IDs, never target/strategy rows."""
    if not taxonomy_paths or len(taxonomy_paths) != len(ortholog_paths):
        raise ValueError("Lineage counts require paired taxonomy and selected-ortholog inputs")
    profiles: dict[str, TaxonomyProfile] = {}
    for path in taxonomy_paths:
        for tax_id, profile in load_taxonomy_profiles(path).items():
            if profiles.setdefault(tax_id, profile) != profile:
                raise ValueError(f"Source taxonomy tables disagree for tax_id {tax_id}: {path}")
    orthologs: dict[str, str] = {}
    for path in ortholog_paths:
        with gzip.open(path, "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene, tax_id = row["ortholog_gene_id"], row["tax_id"]
                if not gene or tax_id not in profiles:
                    raise ValueError(f"Invalid selected ortholog or missing taxonomy: {gene!r}, {tax_id!r}")
                if orthologs.setdefault(gene, tax_id) != tax_id:
                    raise ValueError(f"Selected ortholog {gene} has conflicting tax IDs")
    depths = [0] * len(HUMAN_LINEAGE)
    for tax_id in orthologs.values():
        depths[profiles[tax_id].human_lineage_depth] += 1
    counts = [sum(depths[index:]) for index in range(len(depths))]
    return [dict(threshold=index, tax_id=tax_id, taxon=name, ortholog_count=counts[index])
            for index, (tax_id, name) in enumerate(HUMAN_LINEAGE)]


def member_group_keys(
    ortholog_gene_id: str,
    tax_id: str,
    profiles: dict[str, TaxonomyProfile],
) -> tuple[str, ...]:
    if not ortholog_gene_id:
        return ()
    profile = profiles.get(tax_id)
    if profile is None:
        raise ValueError(
            "Ortholog evidence references tax_id absent from canonical taxonomy: "
            f"ortholog_gene_id={ortholog_gene_id!r}, tax_id={tax_id!r}"
        )
    return tuple(
        f"{scope}__{unit}={profile.unit_id(unit, ortholog_gene_id)}"
        for scope in profile.scopes()
        for unit in UNIT_ORDER
    )


def count_member_groups(
    members: Iterable[tuple[str, str]],
    profiles: dict[str, TaxonomyProfile],
) -> dict[str, int]:
    unit_scope_masks: tuple[dict[object, int], ...] = tuple(
        {} for _unit in UNIT_ORDER
    )
    for ortholog_gene_id, tax_id in members:
        if not ortholog_gene_id:
            continue
        profile = profiles.get(tax_id)
        if profile is None:
            raise ValueError(
                "Ortholog evidence references tax_id absent from canonical taxonomy: "
                f"ortholog_gene_id={ortholog_gene_id!r}, tax_id={tax_id!r}"
            )
        unit_ids = (ortholog_gene_id, *profile.rank_unit_ids)
        for groups, unit_id in zip(unit_scope_masks, unit_ids, strict=True):
            groups[unit_id] = groups.get(unit_id, 0) | profile.scope_mask

    counts = [0] * len(COUNT_KEYS)
    unit_count = len(UNIT_ORDER)
    for unit_index, groups in enumerate(unit_scope_masks):
        mask_counts: dict[int, int] = {}
        for scope_mask in groups.values():
            mask_counts[scope_mask] = mask_counts.get(scope_mask, 0) + 1
        # Canonical lineages are scope prefixes, so one suffix sum replaces
        # repeated set insertions. The mask path preserves arbitrary profiles.
        if all(scope_mask & (scope_mask + 1) == 0 for scope_mask in mask_counts):
            depth_counts = [0] * len(SCOPE_ORDER)
            for scope_mask, group_count in mask_counts.items():
                depth_counts[scope_mask.bit_length() - 1] += group_count
            running_count = 0
            for scope_index in range(len(SCOPE_ORDER) - 1, -1, -1):
                running_count += depth_counts[scope_index]
                counts[scope_index * unit_count + unit_index] = running_count
            continue
        for scope_mask, group_count in mask_counts.items():
            remaining_mask = scope_mask
            while remaining_mask:
                scope_bit = remaining_mask & -remaining_mask
                scope_index = scope_bit.bit_length() - 1
                counts[scope_index * unit_count + unit_index] += group_count
                remaining_mask ^= scope_bit
    return dict(zip(COUNT_KEYS, counts, strict=True))


def _format_stat(value: float) -> str:
    return f"{value:.1f}"


def build_taxonomy_summary_rows(
    ortholog_rows: Iterable[dict[str, str]],
    profiles: dict[str, TaxonomyProfile],
) -> list[dict[str, object]]:
    members_by_gene: dict[str, list[tuple[str, str]]] = {}
    selected_tax_ids: set[str] = set()
    for row_number, row in enumerate(ortholog_rows, start=2):
        gene_id = str(row.get("query_gene_id") or "")
        ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
        tax_id = str(row.get("tax_id") or "").strip()
        if not tax_id:
            raise ValueError(
                "Selected ortholog table has empty tax_id at row "
                f"{row_number}: ortholog_gene_id={ortholog_gene_id!r}"
            )
        selected_tax_ids.add(tax_id)
        if gene_id and ortholog_gene_id:
            members_by_gene.setdefault(gene_id, []).append((ortholog_gene_id, tax_id))

    taxonomy_tax_ids = set(profiles)
    missing = selected_tax_ids - taxonomy_tax_ids
    unexpected = taxonomy_tax_ids - selected_tax_ids
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing taxonomy tax_id(s): " + ", ".join(sorted(missing)))
        if unexpected:
            details.append("unexpected taxonomy tax_id(s): " + ", ".join(sorted(unexpected)))
        raise ValueError("Canonical taxonomy coverage mismatch; " + "; ".join(details))

    rows: list[dict[str, object]] = []
    for scope in SCOPE_ORDER:
        for unit in UNIT_ORDER:
            ortholog_counts: list[int] = []
            unit_counts: list[int] = []
            all_taxa: set[str] = set()
            all_units: set[str] = set()
            for members in members_by_gene.values():
                selected_orthologs: set[str] = set()
                selected_units: set[str] = set()
                for ortholog_gene_id, tax_id in members:
                    profile = profiles[tax_id]
                    if scope not in profile.scopes():
                        continue
                    selected_orthologs.add(ortholog_gene_id)
                    if tax_id:
                        all_taxa.add(tax_id)
                    group_id = profile.unit_id(unit, ortholog_gene_id)
                    selected_units.add(group_id)
                    all_units.add(group_id)
                ortholog_counts.append(len(selected_orthologs))
                unit_counts.append(len(selected_units))
            if not ortholog_counts:
                continue
            rows.append(
                {
                    "taxonomic_scope": scope,
                    "evidence_unit": unit,
                    "gene_count": len(ortholog_counts),
                    "ortholog_count": sum(ortholog_counts),
                    "taxon_count": len(all_taxa),
                    "unit_count": len(all_units),
                    "orthologs_per_gene_min": min(ortholog_counts),
                    "orthologs_per_gene_median": _format_stat(statistics.median(ortholog_counts)),
                    "orthologs_per_gene_mean": _format_stat(statistics.mean(ortholog_counts)),
                    "orthologs_per_gene_max": max(ortholog_counts),
                    "units_per_gene_min": min(unit_counts),
                    "units_per_gene_median": _format_stat(statistics.median(unit_counts)),
                    "units_per_gene_mean": _format_stat(statistics.mean(unit_counts)),
                    "units_per_gene_max": max(unit_counts),
                }
            )
    return rows


def write_taxonomy_summary(path: Path, rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TAXONOMY_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count
