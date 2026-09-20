"""Parameterized SQL compilation for the canonical feature-query grammar."""

from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from ..query import QueryExpressionError, parse_query_ast
from ..serializers import FeatureRow
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


def _comparison(operator: str, left_node: object, right_node: object) -> tuple[str, list[object]]:
    left = _operand(left_node)
    right = _operand(right_node)

    if left.multi is not None or right.multi is not None:
        if left.multi is not None and right.multi is not None:
            # Pairwise existential semantics for two multivalue operands.
            lns, _lvalue, _ = _multi_value_expression(left)
            rns, _rvalue, _ = _multi_value_expression(right)
            params: list[object] = []
            ltable = "qualifiers" if lns.startswith("qualifier:") else "xrefs"
            rtable = "qualifiers" if rns.startswith("qualifier:") else "xrefs"
            lkey = lns.removeprefix("qualifier:")
            rkey = rns.removeprefix("qualifier:")
            lfilter = "l.key = ?" if ltable == "qualifiers" else "l.namespace = ?"
            rfilter = "r.key = ?" if rtable == "qualifiers" else "r.namespace = ?"
            params.extend((lkey, rkey))
            if operator == "contains":
                relation = "instr(CAST(l.value AS TEXT), r.value) > 0"
            elif operator == "~=":
                relation = "gb_regexp(r.value, l.value) = 1"
            elif operator == "!=":
                relation = "l.value = r.value"
            else:
                relation = f"l.value {operator if operator in {'=', '!=', '<', '<=', '>', '>='} else '='} r.value"
            if operator == "in":
                relation = "l.value = r.value"
            exists = (
                f"EXISTS (SELECT 1 FROM {ltable} l JOIN {rtable} r ON l.feature_pk = r.feature_pk "
                f"WHERE l.feature_pk = f.feature_pk AND {lfilter} AND {rfilter} AND {relation})"
            )
            return (f"NOT ({exists})", params) if operator == "!=" else (exists, params)
        if left.multi is not None:
            namespace, value, _ = _multi_value_expression(left)
            if right.multi is not None:
                raise AssertionError("handled above")
            assert right.sql is not None
            rhs_sql = right.sql
            rhs_params = list(right.params)
            if right_node[0] == "list" and operator in {"contains", "~="}:  # type: ignore[index]
                first = _operand(right_node[1][0])  # type: ignore[index]
                assert first.sql is not None
                rhs_sql = first.sql
                rhs_params = list(first.params)
            if right.is_literal and isinstance(right.literal, bool) and operator in {"<", "<=", ">", ">="}:
                raise QueryExpressionError("numeric comparisons do not accept booleans")
            if left.casefold:
                rhs_sql = f"gb_casefold({rhs_sql})"
            if operator == "==":
                relation = f"{value} IN {rhs_sql}" if right_node[0] == "list" else f"{value} = {rhs_sql}"  # type: ignore[index]
            elif operator == "!=":
                relation = f"{value} IN {rhs_sql}" if right_node[0] == "list" else f"{value} = {rhs_sql}"  # type: ignore[index]
                sql, params = _multi_sql(namespace, relation, rhs_params)
                return f"NOT ({sql})", params
            elif operator == "in":
                relation = f"{value} IN {rhs_sql if rhs_sql.startswith('(') else '(' + rhs_sql + ')'}"
            elif operator == "contains":
                relation = f"instr(CAST({value} AS TEXT), {rhs_sql}) > 0"
            elif operator == "~=":
                relation = f"gb_regexp({rhs_sql}, {value}) = 1"
            else:
                relation = f"{value} {operator} {rhs_sql}"
            return _multi_sql(namespace, relation, rhs_params)
        # Scalar left, multivalue right: ``in`` and equality both mean that
        # at least one right-hand value satisfies the relation.
        namespace, value, _ = _multi_value_expression(right)
        assert left.sql is not None
        if operator == "contains":
            relation = f"instr(CAST({left.sql} AS TEXT), {value}) > 0"
        elif operator == "~=":
            relation = f"gb_regexp({value}, {left.sql}) = 1"
        elif operator == "!=":
            relation = f"{left.sql} = {value}"
            sql, params = _multi_sql(namespace, relation, list(left.params))
            return f"NOT ({sql})", params
        else:
            relation = f"{left.sql} {'=' if operator == 'in' else operator} {value}"
        return _multi_sql(namespace, relation, list(left.params))

    assert left.sql is not None and right.sql is not None
    params = list(left.params) + list(right.params)
    if operator in {"<", "<=", ">", ">="}:
        for operand in (left, right):
            if operand.is_literal and isinstance(operand.literal, bool):
                raise QueryExpressionError("numeric comparisons do not accept booleans")
    if right_node[0] == "list":  # type: ignore[index]
        if operator in {"==", "in"}:
            return f"{left.sql} IN {right.sql}", params
        if operator == "!=":
            return f"({left.sql} IS NULL OR {left.sql} NOT IN {right.sql})", params
        if operator in {"contains", "~="}:
            first = _operand(right_node[1][0])  # type: ignore[index]
            assert first.sql is not None
            right_sql = first.sql
            params = list(left.params) + list(first.params)
            if operator == "contains":
                return f"({left.sql} IS NOT NULL AND instr(CAST({left.sql} AS TEXT), {right_sql}) > 0)", params
            return f"gb_regexp({right_sql}, {left.sql}) = 1", params
        raise QueryExpressionError(f"operator {operator!r} does not accept a list")
    right_value = right.literal if right.is_literal else object()
    if operator == "==":
        if right.is_literal and right_value is None:
            return f"{left.sql} IS NULL", list(left.params)
        if not right.is_literal:
            if left.is_literal and left.literal is None:
                return f"{right.sql} IS NULL", list(right.params)
            if left.is_literal:
                return f"{right.sql} = {left.sql}", list(right.params) + list(left.params)
            return f"({left.sql} = {right.sql} OR ({left.sql} IS NULL AND {right.sql} IS NULL))", params
        return f"{left.sql} = {right.sql}", params
    if operator == "!=":
        if right.is_literal and right_value is None:
            return f"{left.sql} IS NOT NULL", list(left.params)
        if not right.is_literal:
            if left.is_literal and left.literal is None:
                return f"{right.sql} IS NOT NULL", list(right.params)
            if left.is_literal:
                return f"({right.sql} IS NULL OR {right.sql} != {left.sql})", list(right.params) + list(left.params)
            return f"(({left.sql} IS NULL AND {right.sql} IS NOT NULL) OR ({left.sql} IS NOT NULL AND ({right.sql} IS NULL OR {left.sql} != {right.sql})))", params
        return f"({left.sql} IS NULL OR {left.sql} != {right.sql})", params
    if operator == "contains":
        return f"({left.sql} IS NOT NULL AND instr(CAST({left.sql} AS TEXT), {right.sql}) > 0)", params
    if operator == "~=":
        return f"gb_regexp({right.sql}, {left.sql}) = 1", params
    if operator == "in":
        return f"{left.sql} = {right.sql}", params
    # SQLite's comparison already returns false for NULL operands, matching
    # the direct evaluator while keeping each bound literal parameter once.
    return f"({left.sql} {operator} {right.sql})", params


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

    tree = parse_query_ast(expression)
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

    connection.create_function("gb_regexp", 2, regexp)
    connection.create_function("gb_casefold", 1, casefold)
    connection.create_function("gb_truth", 1, truth)


def _row_from_database(
    connection: sqlite3.Connection,
    row: sqlite3.Row,
) -> FeatureRow:
    qualifiers: dict[str, tuple[str, ...]] = {}
    for key, value in connection.execute(
        "SELECT key, value FROM qualifiers WHERE feature_pk = ? ORDER BY key, ordinal",
        (row["feature_pk"],),
    ):
        qualifiers.setdefault(str(key), ())
        qualifiers[str(key)] = (*qualifiers[str(key)], str(value))
    xrefs: dict[str, list[str]] = {}
    xref_sources: dict[str, list[tuple[str, str]]] = {}
    for namespace, value, source_field in connection.execute(
        "SELECT namespace, value, source_field FROM xrefs WHERE feature_pk = ? ORDER BY namespace, value, source_field",
        (row["feature_pk"],),
    ):
        xrefs.setdefault(str(namespace), []).append(str(value))
        xref_sources.setdefault(str(namespace), []).append((str(value), str(source_field)))
    segments = tuple(
        (int(start), int(end))
        for start, end in connection.execute(
            "SELECT start, end FROM segments WHERE feature_pk = ? ORDER BY ordinal",
            (row["feature_pk"],),
        )
    )
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
        xrefs={key: tuple(values) for key, values in xrefs.items()},
        xref_sources={key: tuple(values) for key, values in xref_sources.items()},
    )


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
        return tuple(_row_from_database(connection, row) for row in connection.execute(sql, params))
    finally:
        connection.close()


__all__ = ["compile_query_sql", "query_index"]
