"""Parameterized SQL compilation for the canonical feature-query grammar."""

from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from ..query import QueryExpressionError, parse_query_ast, validate_query_types
from ..serializers import XREF_NAMESPACES, FeatureRow
from .schema import validate_schema

_MULTI_FIELDS = {
    "ko": "kegg_kos",
    "kegg_ko": "kegg_kos",
    "ec": "ec_numbers",
    "ec_number": "ec_numbers",
    "cog": "cog_ids",
    "pfam": "pfam",
    "rfam": "rfam",
    "go": "go_terms",
    "go_terms": "go_terms",
    "db_xref": "qualifier:db_xref",
}

_SCALAR_FIELDS = {
    "record": "r.record_id",
    "record_id": "r.record_id",
    "record_index": "r.record_index",
    "feature_index": "f.feature_index",
    "type": "f.type",
    "locus_tag": "f.locus_tag",
    "gene": "f.gene",
    "product": "f.product",
    "protein_id": "f.protein_id",
    "start": "f.start",
    "end": "f.end",
    "strand_value": "f.strand",
    "length": "f.biological_length",
    "record_length": "r.length",
    "topology": "r.topology",
    "pseudo": "f.pseudo",
    "partial": "f.partial",
    "source": "s.display_path",
    "sample": "s.sample_key",
    "sample_key": "s.sample_key",
}


@dataclass(frozen=True)
class _Operand:
    sql: str | None
    params: tuple[object, ...]
    multi: str | None = None
    casefold: bool = False
    literal: object = None
    is_literal: bool = False


def _strand_sql() -> str:
    return "CASE f.strand WHEN 1 THEN '+' WHEN -1 THEN '-' WHEN 0 THEN '?' ELSE '.' END"


def _operand(node: object) -> _Operand:
    kind = node[0]  # type: ignore[index]
    if kind == "literal":
        return _Operand("?", (node[1],), literal=node[1], is_literal=True)  # type: ignore[index]
    if kind == "field":
        field = str(node[1])  # type: ignore[index]
        if field == "strand":
            return _Operand(_strand_sql(), ())
        if field in _MULTI_FIELDS:
            return _Operand(None, (), multi=_MULTI_FIELDS[field])
        sql = _SCALAR_FIELDS.get(field)
        if sql is None:
            raise QueryExpressionError(f"indexed query does not support field {field!r}")
        return _Operand(sql, ())
    if kind == "qualifier":
        return _Operand(None, (), multi="qualifier:" + str(node[1]))  # type: ignore[index]
    if kind == "casefold":
        inner = _operand(node[1])  # type: ignore[index]
        if inner.multi is not None:
            return _Operand(None, inner.params, multi=inner.multi, casefold=True)
        assert inner.sql is not None
        return _Operand(f"gb_casefold({inner.sql})", inner.params, literal=inner.literal, is_literal=inner.is_literal, casefold=True)
    if kind == "list":
        values = tuple(_operand(item) for item in node[1])  # type: ignore[index]
        if any(item.multi is not None for item in values):
            raise QueryExpressionError("indexed query lists cannot contain multivalue fields")
        return _Operand(
            "(" + ", ".join(item.sql or "?" for item in values) + ")",
            tuple(param for item in values for param in item.params),
        )
    raise QueryExpressionError(f"invalid indexed query operand {kind!r}")


def _multi_sql(namespace: str, expression: str, params: list[object]) -> tuple[str, list[object]]:
    if namespace.startswith("qualifier:"):
        return (
            "EXISTS (SELECT 1 FROM qualifiers q WHERE q.feature_pk = f.feature_pk AND q.key = ? AND "
            + expression
            + ")",
            [namespace[len("qualifier:") :], *params],
        )
    return (
        "EXISTS (SELECT 1 FROM xrefs q WHERE q.feature_pk = f.feature_pk AND q.namespace = ? AND "
        + expression
        + ")",
        [namespace, *params],
    )


def _multi_value_expression(operand: _Operand) -> tuple[str, str, list[object]]:
    if operand.multi is None:
        raise QueryExpressionError("expected a multivalue operand")
    value = "gb_casefold(q.value)" if operand.casefold else "q.value"
    return operand.multi, value, []


def _equal_sql(left: str, right: str) -> str:
    """Use Python-like equality without SQLite numeric/text affinity."""

    return f"gb_equal({left}, {right}) = 1"


def _contains_sql(left: str, right: str) -> str:
    return f"gb_contains({left}, {right}) = 1"


def _list_operands(node: object) -> tuple[_Operand, ...]:
    if node[0] != "list":  # type: ignore[index]
        return ()
    return tuple(_operand(item) for item in node[1])  # type: ignore[index]


def _comparison(operator: str, left_node: object, right_node: object) -> tuple[str, list[object]]:
    left = _operand(left_node)
    right = _operand(right_node)

    if left.multi is not None and right.multi is not None:
        lns = left.multi
        rns = right.multi
        ltable = "qualifiers" if lns.startswith("qualifier:") else "xrefs"
        rtable = "qualifiers" if rns.startswith("qualifier:") else "xrefs"
        lkey = lns.removeprefix("qualifier:")
        rkey = rns.removeprefix("qualifier:")
        lfilter = "l.key = ?" if ltable == "qualifiers" else "l.namespace = ?"
        rfilter = "r.key = ?" if rtable == "qualifiers" else "r.namespace = ?"
        lvalue = "gb_casefold(l.value)" if left.casefold else "l.value"
        rvalue = "gb_casefold(r.value)" if right.casefold else "r.value"
        if operator in {"==", "!=", "in"}:
            relation = _equal_sql(lvalue, rvalue)
        elif operator == "contains":
            relation = _contains_sql(lvalue, rvalue)
        elif operator == "~=":
            relation = f"gb_regexp({rvalue}, {lvalue}) = 1"
        else:
            relation = f"gb_compare('{operator}', {lvalue}, {rvalue}) = 1"
        exists = (
            f"EXISTS (SELECT 1 FROM {ltable} l JOIN {rtable} r ON l.feature_pk = r.feature_pk "
            f"WHERE l.feature_pk = f.feature_pk AND {lfilter} AND {rfilter} "
            + ("AND r.ordinal = 0 AND " if operator in {"contains", "~=", "<", "<=", ">", ">="} else "AND ")
            + f"{relation})"
        )
        params = [lkey, rkey]
        return (f"NOT ({exists})", params) if operator == "!=" else (exists, params)

    if left.multi is not None:
        namespace = left.multi
        value = "gb_casefold(q.value)" if left.casefold else "q.value"
        rhs_operands = _list_operands(right_node)
        if right_node[0] == "list" and not rhs_operands:  # type: ignore[index]
            return ("1" if operator == "!=" else "0"), []
        if operator in {"==", "!=", "in"} and rhs_operands:
            relations = [_equal_sql(value, item.sql or "?") for item in rhs_operands]
            rhs_params = [param for item in rhs_operands for param in item.params]
            relation = "(" + " OR ".join(relations) + ")"
        else:
            item = rhs_operands[0] if rhs_operands else right
            rhs_sql = item.sql or "?"
            rhs_params = list(item.params)
            if operator == "contains":
                relation = _contains_sql(value, rhs_sql)
            elif operator == "~=":
                relation = f"gb_regexp({rhs_sql}, {value}) = 1"
            elif operator in {"==", "in"} or operator == "!=":
                relation = _equal_sql(value, rhs_sql)
            else:
                relation = f"gb_compare('{operator}', {value}, {rhs_sql}) = 1"
        sql, params = _multi_sql(namespace, relation, rhs_params)
        return (f"NOT ({sql})", params) if operator == "!=" else (sql, params)

    if right.multi is not None:
        namespace = right.multi
        value = "gb_casefold(q.value)" if right.casefold else "q.value"
        assert left.sql is not None
        rhs_filter = "q.ordinal = 0 AND "
        if operator in {"==", "in"}:
            relation = _equal_sql(left.sql, value)
        elif operator == "!=":
            relation = _equal_sql(left.sql, value)
            sql, params = _multi_sql(namespace, relation, list(left.params))
            return f"NOT ({sql})", params
        elif operator == "contains":
            relation = _contains_sql(left.sql, value)
        elif operator == "~=":
            relation = f"gb_regexp({value}, {left.sql}) = 1"
        else:
            relation = f"gb_compare('{operator}', {left.sql}, {value}) = 1"
        if operator in {"contains", "~=", "<", "<=", ">", ">="}:
            relation = rhs_filter + relation
        return _multi_sql(namespace, relation, list(left.params))

    assert left.sql is not None and right.sql is not None
    left_sql = left.sql
    right_operands = _list_operands(right_node)
    if right_node[0] == "list" and not right_operands:  # type: ignore[index]
        if operator == "!=":
            return "1", []
        return "0", []
    if right_operands:
        if operator in {"==", "in"}:
            relations = [_equal_sql(left_sql, item.sql or "?") for item in right_operands]
            return f"({left_sql} IS NOT NULL AND ({' OR '.join(relations)}))", list(left.params) + [param for item in right_operands for param in item.params]
        if operator == "!=":
            relations = [_equal_sql(left_sql, item.sql or "?") for item in right_operands]
            return f"({left_sql} IS NULL OR NOT ({' OR '.join(relations)}))", list(left.params) + [param for item in right_operands for param in item.params]
        item = right_operands[0]
        right_sql = item.sql or "?"
        params = list(left.params) + list(item.params)
    else:
        right_sql = right.sql
        params = list(left.params) + list(right.params)
    if operator == "in":
        return f"({left_sql} IS NOT NULL AND {_equal_sql(left_sql, right_sql)})", params
    if operator == "==":
        return _equal_sql(left_sql, right_sql), params
    if operator == "!=":
        return f"NOT ({_equal_sql(left_sql, right_sql)})", params
    if operator == "contains":
        return _contains_sql(left_sql, right_sql), params
    if operator == "~=":
        return f"gb_regexp({right_sql}, {left_sql}) = 1", params
    return f"gb_compare('{operator}', {left_sql}, {right_sql}) = 1", params


def _compile(node: object) -> tuple[str, list[object]]:
    kind = node[0]  # type: ignore[index]
    if kind == "or":
        left, lp = _compile(node[1])  # type: ignore[index]
        right, rp = _compile(node[2])  # type: ignore[index]
        return f"({left} OR {right})", lp + rp
    if kind == "and":
        left, lp = _compile(node[1])  # type: ignore[index]
        right, rp = _compile(node[2])  # type: ignore[index]
        return f"({left} AND {right})", lp + rp
    if kind == "not":
        child, params = _compile(node[1])  # type: ignore[index]
        return f"NOT ({child})", params
    if kind == "truth":
        operand = _operand(node[1])  # type: ignore[index]
        if operand.multi is not None:
            namespace = operand.multi
            expression = "gb_truth(q.value) = 1"
            return _multi_sql(namespace, expression, [])
        assert operand.sql is not None
        return f"gb_truth({operand.sql}) = 1", list(operand.params)
    if kind == "compare":
        return _comparison(str(node[1]), node[2], node[3])  # type: ignore[index]
    raise QueryExpressionError(f"invalid indexed query expression node {kind!r}")


def compile_query_sql(expression: str) -> tuple[str, tuple[object, ...]]:
    """Compile the canonical query AST to SQL with only bound values."""

    tree = validate_query_types(parse_query_ast(expression))
    sql, params = _compile(tree)
    return sql, tuple(params)


def _register_functions(connection: sqlite3.Connection) -> None:
    def regexp(pattern: object, value: object) -> int:
        if not isinstance(pattern, str) or not isinstance(value, str):
            return 0
        return int(re.search(pattern, value) is not None)

    def casefold(value: object) -> str | None:
        return None if value is None else str(value).casefold()

    def truth(value: object) -> int:
        return int(bool(value))

    def equal(left: object, right: object) -> int:
        if left is None or right is None:
            return int(left is None and right is None)
        if isinstance(left, str) != isinstance(right, str):
            return 0
        if isinstance(left, (int, float)) and isinstance(right, (int, float)):
            return int(left == right)
        return int(left == right)

    def contains(value: object, needle: object) -> int:
        return int(isinstance(value, str) and isinstance(needle, str) and needle in value)

    def compare(operator: str, left: object, right: object) -> int:
        if left is None or right is None:
            return 0
        if isinstance(left, bool) or isinstance(right, bool):
            return 0
        if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
            return 0
        return int(
            {
                "<": left < right,
                "<=": left <= right,
                ">": left > right,
                ">=": left >= right,
            }[operator]
        )

    connection.create_function("gb_regexp", 2, regexp)
    connection.create_function("gb_casefold", 1, casefold)
    connection.create_function("gb_truth", 1, truth)
    connection.create_function("gb_equal", 2, equal)
    connection.create_function("gb_contains", 2, contains)
    connection.create_function("gb_compare", 3, compare)


def _row_from_database(
    connection: sqlite3.Connection,
    row: sqlite3.Row,
    *,
    qualifiers_by_feature: dict[int, dict[str, tuple[str, ...]]] | None = None,
    xrefs_by_feature: dict[int, tuple[dict[str, tuple[str, ...]], dict[str, tuple[tuple[str, str], ...]]]] | None = None,
    segments_by_feature: dict[int, tuple[tuple[int, int], ...]] | None = None,
) -> FeatureRow:
    feature_pk = int(row["feature_pk"])
    if qualifiers_by_feature is None:
        qualifiers: dict[str, tuple[str, ...]] = {}
        for key, value in connection.execute(
            "SELECT key, value FROM qualifiers WHERE feature_pk = ? ORDER BY key, ordinal",
            (feature_pk,),
        ):
            qualifiers.setdefault(str(key), ())
            qualifiers[str(key)] = (*qualifiers[str(key)], str(value))
    else:
        qualifiers = qualifiers_by_feature.get(feature_pk, {})
    if xrefs_by_feature is None:
        xrefs: dict[str, tuple[str, ...]] = {key: () for key in XREF_NAMESPACES}
        xref_sources: dict[str, tuple[tuple[str, str], ...]] = {key: () for key in XREF_NAMESPACES}
        for namespace, value, source_field in connection.execute(
            "SELECT namespace, value, source_field FROM xrefs WHERE feature_pk = ? ORDER BY namespace, ordinal",
            (feature_pk,),
        ):
            key = str(namespace)
            xrefs[key] = (*xrefs.get(key, ()), str(value))
            xref_sources[key] = (*xref_sources.get(key, ()), (str(value), str(source_field)))
    else:
        xrefs, xref_sources = xrefs_by_feature.get(
            feature_pk,
            ({key: () for key in XREF_NAMESPACES}, {key: () for key in XREF_NAMESPACES}),
        )
    if segments_by_feature is None:
        segments = tuple(
            (int(start), int(end))
            for start, end in connection.execute(
                "SELECT start, end FROM segments WHERE feature_pk = ? ORDER BY ordinal",
                (feature_pk,),
            )
        )
    else:
        segments = segments_by_feature.get(feature_pk, ())
    strand = row["strand"]
    symbol = "+" if strand == 1 else "-" if strand == -1 else "?" if strand == 0 else "."
    return FeatureRow(
        source=str(row["display_path"]),
        sample_key=str(row["sample_key"]),
        record=str(row["record_id"]),
        record_index=int(row["record_index"]),
        feature_index=int(row["feature_index"]),
        type=str(row["type"]),
        locus_tag=row["locus_tag"],
        gene=row["gene"],
        product=row["product"],
        protein_id=row["protein_id"],
        start=int(row["start"]),
        end=int(row["end"]),
        strand=strand,
        strand_symbol=symbol,
        length=int(row["biological_length"]),
        record_length=int(row["record_length"]),
        topology=row["topology"],
        partial=bool(row["partial"]),
        pseudo=bool(row["pseudo"]),
        segments=segments,
        qualifiers=qualifiers,
        xrefs=xrefs,
        xref_sources=xref_sources,
    )


def _bulk_hydrate(
    connection: sqlite3.Connection,
    rows: list[sqlite3.Row],
) -> tuple[
    dict[int, dict[str, tuple[str, ...]]],
    dict[int, tuple[dict[str, tuple[str, ...]], dict[str, tuple[tuple[str, str], ...]]]],
    dict[int, tuple[tuple[int, int], ...]],
]:
    """Load related projections in three bounded bulk queries."""

    qualifiers: dict[int, dict[str, list[str]]] = {}
    xrefs: dict[int, dict[str, list[str]]] = {}
    xref_sources: dict[int, dict[str, list[tuple[str, str]]]] = {}
    segments: dict[int, list[tuple[int, int]]] = {}
    feature_ids = [int(row["feature_pk"]) for row in rows]
    for offset in range(0, len(feature_ids), 900):
        chunk = feature_ids[offset : offset + 900]
        placeholders = ",".join("?" for _ in chunk)
        for feature_pk, key, value in connection.execute(
            f"SELECT feature_pk, key, value FROM qualifiers WHERE feature_pk IN ({placeholders}) ORDER BY feature_pk, key, ordinal",
            chunk,
        ):
            qualifiers.setdefault(int(feature_pk), {}).setdefault(str(key), []).append(str(value))
        for feature_pk, namespace, _ordinal, value, source_field in connection.execute(
            f"SELECT feature_pk, namespace, ordinal, value, source_field FROM xrefs WHERE feature_pk IN ({placeholders}) ORDER BY feature_pk, namespace, ordinal",
            chunk,
        ):
            feature_key = int(feature_pk)
            namespace_key = str(namespace)
            xrefs.setdefault(feature_key, {}).setdefault(namespace_key, []).append(str(value))
            xref_sources.setdefault(feature_key, {}).setdefault(namespace_key, []).append((str(value), str(source_field)))
        for feature_pk, start, end in connection.execute(
            f"SELECT feature_pk, start, end FROM segments WHERE feature_pk IN ({placeholders}) ORDER BY feature_pk, ordinal",
            chunk,
        ):
            segments.setdefault(int(feature_pk), []).append((int(start), int(end)))
    qualifier_result = {
        feature_pk: {key: tuple(values) for key, values in by_key.items()}
        for feature_pk, by_key in qualifiers.items()
    }
    xref_result = {}
    for feature_pk in feature_ids:
        by_namespace = xrefs.get(feature_pk, {})
        by_source = xref_sources.get(feature_pk, {})
        xref_result[feature_pk] = (
            {key: tuple(by_namespace.get(key, ())) for key in XREF_NAMESPACES},
            {key: tuple(by_source.get(key, ())) for key in XREF_NAMESPACES},
        )
    return qualifier_result, xref_result, {key: tuple(value) for key, value in segments.items()}


def query_index(
    database: str | Path,
    where: str,
    *,
    limit: int | None = None,
) -> tuple[FeatureRow, ...]:
    """Query an index with the same bounded grammar as direct ``gbparse query``."""

    if limit is not None and limit <= 0:
        raise ValueError("limit must be positive")
    path = Path(database)
    if not path.is_file():
        raise FileNotFoundError(f"index does not exist: {path}")
    predicate, params = compile_query_sql(where)
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    _register_functions(connection)
    try:
        validate_schema(connection)
        sql = (
            "SELECT s.display_path, s.sample_key, r.record_id, r.record_index, r.length AS record_length, "
            "r.topology, f.feature_pk, f.feature_index, f.type, f.locus_tag, f.gene, f.product, "
            "f.protein_id, f.start, f.end, f.strand, f.biological_length, f.partial, f.pseudo "
            "FROM features f JOIN records r ON r.record_pk = f.record_pk "
            "JOIN sources s ON s.source_pk = r.source_pk "
            f"WHERE {predicate} ORDER BY s.display_path, r.record_index, f.feature_index"
        )
        if limit is not None:
            sql += " LIMIT ?"
            params = (*params, limit)
        rows = list(connection.execute(sql, params))
        qualifiers, xrefs, segments = _bulk_hydrate(connection, rows)
        return tuple(
            _row_from_database(
                connection,
                row,
                qualifiers_by_feature=qualifiers,
                xrefs_by_feature=xrefs,
                segments_by_feature=segments,
            )
            for row in rows
        )
    finally:
        connection.close()


__all__ = ["compile_query_sql", "query_index"]
