"""Feature search and query engine supporting gene, product, KO, COG, EC, Pfam and regex filters."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

from .io import extract_xrefs, read_genbank
from .model import GenBankFeature
from .serializers import FeatureRow


class QueryExpressionError(ValueError):
    """A query is syntactically invalid or exceeds a resource limit."""


@dataclass(frozen=True)
class _Token:
    kind: str
    value: object
    position: int


@dataclass(frozen=True)
class QueryMatch:
    """A matching feature with both its canonical row and source objects."""

    row: FeatureRow
    feature: GenBankFeature
    record: Any


_QUERY_FIELDS = frozenset(
    {
        "record",
        "record_id",
        "record_index",
        "feature_index",
        "type",
        "locus_tag",
        "gene",
        "product",
        "protein_id",
        "start",
        "end",
        "strand",
        "strand_value",
        "length",
        "record_length",
        "topology",
        "pseudo",
        "partial",
        "ko",
        "kegg_ko",
        "ec",
        "ec_number",
        "cog",
        "pfam",
        "rfam",
        "go",
        "go_terms",
        "db_xref",
        "source",
    }
)
_MAX_QUERY_LENGTH = 8_192
_MAX_QUERY_DEPTH = 32
_MAX_REGEX_LENGTH = 1_024
_MAX_REGEX_COUNT = 32


def _decode_query_string(text: str, position: int) -> str:
    """Decode the small quoted-string grammar without evaluating Python."""

    chars: list[str] = []
    index = 0
    escapes = {"n": "\n", "r": "\r", "t": "\t", "\\": "\\", '"': '"', "'": "'"}
    while index < len(text):
        char = text[index]
        if char != "\\":
            chars.append(char)
            index += 1
            continue
        index += 1
        if index >= len(text):
            raise QueryExpressionError(f"unterminated escape at position {position}")
        escaped = text[index]
        if escaped == "u":
            digits = text[index + 1 : index + 5]
            if len(digits) != 4 or not re.fullmatch(r"[0-9A-Fa-f]{4}", digits):
                raise QueryExpressionError(f"invalid unicode escape at position {position}")
            chars.append(chr(int(digits, 16)))
            index += 5
            continue
        chars.append(escapes.get(escaped, escaped))
        index += 1
    return "".join(chars)


def _tokenize_query(expression: str) -> tuple[_Token, ...]:
    if len(expression) > _MAX_QUERY_LENGTH:
        raise QueryExpressionError(
            f"query exceeds the {_MAX_QUERY_LENGTH}-character limit"
        )
    tokens: list[_Token] = []
    index = 0
    two_char_ops = {"==", "!=", "<=", ">=", "~=", "<", ">"}
    while index < len(expression):
        char = expression[index]
        if char.isspace():
            index += 1
            continue
        if char in "()[],":
            kind = {
                "(": "LPAREN",
                ")": "RPAREN",
                "[": "LBRACKET",
                "]": "RBRACKET",
                ",": "COMMA",
            }[char]
            tokens.append(_Token(kind, char, index))
            index += 1
            continue
        if expression[index : index + 2] in two_char_ops:
            tokens.append(_Token("OP", expression[index : index + 2], index))
            index += 2
            continue
        if char in "<>":
            tokens.append(_Token("OP", char, index))
            index += 1
            continue
        if char in "\"'":
            quote = char
            start = index
            index += 1
            raw: list[str] = []
            while index < len(expression):
                current = expression[index]
                if current == quote:
                    break
                if current == "\\" and index + 1 < len(expression):
                    raw.extend((current, expression[index + 1]))
                    index += 2
                else:
                    raw.append(current)
                    index += 1
            if index >= len(expression) or expression[index] != quote:
                raise QueryExpressionError(f"unterminated string at position {start}")
            tokens.append(_Token("STRING", _decode_query_string("".join(raw), start), start))
            index += 1
            continue
        number = re.match(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", expression[index:])
        if number:
            raw_number = number.group(0)
            value: int | float = (
                float(raw_number) if "." in raw_number else int(raw_number)
            )
            tokens.append(_Token("NUMBER", value, index))
            index += len(raw_number)
            continue
        word = re.match(r"[A-Za-z_][A-Za-z0-9_]*", expression[index:])
        if word:
            value = word.group(0)
            tokens.append(_Token("WORD", value, index))
            index += len(value)
            continue
        raise QueryExpressionError(f"unexpected character {char!r} at position {index}")
    tokens.append(_Token("EOF", None, len(expression)))
    return tuple(tokens)


class _QueryParser:
    def __init__(self, expression: str) -> None:
        self.tokens = _tokenize_query(expression)
        self.index = 0
        self.depth = 0
        self.regex_count = 0

    @property
    def current(self) -> _Token:
        return self.tokens[self.index]

    def advance(self) -> _Token:
        token = self.current
        self.index += 1
        return token

    def accept(self, kind: str, value: str | None = None) -> _Token | None:
        if self.current.kind != kind:
            return None
        if value is not None and str(self.current.value).casefold() != value.casefold():
            return None
        return self.advance()

    def expect(self, kind: str, value: str | None = None) -> _Token:
        token = self.accept(kind, value)
        if token is None:
            actual = self.current.value if self.current.kind != "EOF" else "end of expression"
            expected = value or kind
            raise QueryExpressionError(f"expected {expected} near {actual!r}")
        return token

    def parse(self) -> object:
        tree = self.parse_or()
        if self.current.kind != "EOF":
            raise QueryExpressionError(f"unexpected token {self.current.value!r}")
        return tree

    def parse_or(self) -> object:
        node = self.parse_and()
        while self.accept("WORD", "or"):
            node = ("or", node, self.parse_and())
        return node

    def parse_and(self) -> object:
        node = self.parse_not()
        while self.accept("WORD", "and"):
            node = ("and", node, self.parse_not())
        return node

    def parse_not(self) -> object:
        if self.accept("WORD", "not"):
            return ("not", self.parse_not())
        return self.parse_primary()

    def parse_primary(self) -> object:
        if self.accept("LPAREN"):
            self.depth += 1
            if self.depth > _MAX_QUERY_DEPTH:
                raise QueryExpressionError("query nesting exceeds the resource limit")
            node = self.parse_or()
            self.expect("RPAREN")
            self.depth -= 1
            return node
        left = self.parse_operand()
        if self.current.kind == "WORD" and str(self.current.value).casefold() in {
            "contains",
            "in",
        }:
            operator = str(self.advance().value).casefold()
        elif self.current.kind == "OP":
            operator = str(self.advance().value)
        else:
            return ("truth", left)
        right = self.parse_operand(allow_list=True)
        if operator == "~=":
            regex_value = right[1] if isinstance(right, tuple) and right[0] == "literal" else None
            if not isinstance(regex_value, str):
                raise QueryExpressionError("~= requires a string regular expression")
            if len(regex_value) > _MAX_REGEX_LENGTH:
                raise QueryExpressionError("regular expression exceeds the resource limit")
            self.regex_count += 1
            if self.regex_count > _MAX_REGEX_COUNT:
                raise QueryExpressionError("too many regular expressions in query")
            try:
                re.compile(regex_value)
            except re.error as exc:
                raise QueryExpressionError(f"invalid regular expression: {exc}") from exc
        return ("compare", operator, left, right)

    def parse_operand(self, *, allow_list: bool = False) -> object:
        if allow_list and self.accept("LBRACKET"):
            values: list[object] = []
            if not self.accept("RBRACKET"):
                while True:
                    values.append(self.parse_operand())
                    if self.accept("RBRACKET"):
                        break
                    self.expect("COMMA")
            return ("list", tuple(values))
        if allow_list and self.accept("LPAREN"):
            values: list[object] = []
            if not self.accept("RPAREN"):
                while True:
                    values.append(self.parse_operand())
                    if self.accept("RPAREN"):
                        break
                    self.expect("COMMA")
            return ("list", tuple(values))
        token = self.advance()
        if token.kind == "STRING" or token.kind == "NUMBER":
            return ("literal", token.value)
        if token.kind != "WORD":
            raise QueryExpressionError(f"expected a field or literal near {token.value!r}")
        word = str(token.value)
        folded = word.casefold()
        if folded == "true":
            return ("literal", True)
        if folded == "false":
            return ("literal", False)
        if folded in {"null", "none"}:
            return ("literal", None)
        if self.accept("LPAREN"):
            if folded == "qualifier":
                key = self.expect("STRING").value
                self.expect("RPAREN")
                return ("qualifier", str(key))
            if folded == "casefold":
                value = self.parse_operand()
                self.expect("RPAREN")
                return ("casefold", value)
            raise QueryExpressionError(f"unsupported query function {word!r}")
        if word not in _QUERY_FIELDS:
            raise QueryExpressionError(f"unknown query field {word!r}")
        return ("field", word)


def parse_query(expression: str) -> object:
    """Parse a safe declarative query into a private immutable AST."""

    return _QueryParser(expression).parse()


def _evaluate_operand(node: object, row: FeatureRow) -> object:
    kind = node[0]  # type: ignore[index]
    if kind == "literal":
        return node[1]  # type: ignore[index]
    if kind == "field":
        return row.value(node[1])  # type: ignore[index]
    if kind == "qualifier":
        return row.qualifier(node[1])  # type: ignore[index]
    if kind == "casefold":
        value = _evaluate_operand(node[1], row)  # type: ignore[index]
        if isinstance(value, (tuple, list)):
            return tuple(str(item).casefold() for item in value)
        return None if value is None else str(value).casefold()
    if kind == "list":
        return tuple(_evaluate_operand(item, row) for item in node[1])  # type: ignore[index]
    raise QueryExpressionError(f"invalid query node {kind!r}")


def _as_values(value: object) -> tuple[object, ...]:
    if value is None:
        return ()
    if isinstance(value, (tuple, list)):
        return tuple(value)
    return (value,)


def _compare_values(operator: str, left: object, right: object) -> bool:
    left_values = _as_values(left)
    right_values = _as_values(right)
    if operator == "in":
        if not right_values:
            return False
        return any(item == candidate for item in left_values for candidate in right_values)
    if operator == "contains":
        if not right_values:
            return False
        needle = right_values[0]
        if not isinstance(needle, str):
            raise QueryExpressionError("contains requires a string right operand")
        return any(isinstance(item, str) and needle in item for item in left_values)
    if operator == "~=":
        pattern = right_values[0] if right_values else None
        if not isinstance(pattern, str):
            return False
        compiled = re.compile(pattern)
        return any(isinstance(item, str) and compiled.search(item) is not None for item in left_values)
    if operator in {"<", "<=", ">", ">="}:
        if not left_values or not right_values:
            return False
        left_value = left_values[0]
        right_value = right_values[0]
        if isinstance(left_value, bool) or isinstance(right_value, bool):
            raise QueryExpressionError("numeric comparisons do not accept booleans")
        if not isinstance(left_value, (int, float)) or not isinstance(right_value, (int, float)):
            raise QueryExpressionError("numeric comparisons require numeric operands")
        return {
            "<": left_value < right_value,
            "<=": left_value <= right_value,
            ">": left_value > right_value,
            ">=": left_value >= right_value,
        }[operator]
    if operator in {"==", "!="}:
        equal = bool(left_values) and any(item == right for item in left_values for right in right_values)
        if left is None and right is None:
            equal = True
        return equal if operator == "==" else not equal
    raise QueryExpressionError(f"unsupported query operator {operator!r}")


def evaluate_query(tree: object, row: FeatureRow) -> bool:
    kind = tree[0]  # type: ignore[index]
    if kind == "or":
        return evaluate_query(tree[1], row) or evaluate_query(tree[2], row)  # type: ignore[index]
    if kind == "and":
        return evaluate_query(tree[1], row) and evaluate_query(tree[2], row)  # type: ignore[index]
    if kind == "not":
        return not evaluate_query(tree[1], row)  # type: ignore[index]
    if kind == "truth":
        value = _evaluate_operand(tree[1], row)  # type: ignore[index]
        return any(bool(item) for item in _as_values(value))
    if kind == "compare":
        return _compare_values(
            tree[1],  # type: ignore[index]
            _evaluate_operand(tree[2], row),  # type: ignore[index]
            _evaluate_operand(tree[3], row),  # type: ignore[index]
        )
    raise QueryExpressionError(f"invalid query expression node {kind!r}")


def query_features(
    source: str | Path | TextIO,
    where: str,
) -> list[QueryMatch]:
    """Return all features satisfying a safe declarative expression."""

    tree = parse_query(where)
    document = read_genbank(source)
    source_label = document.source_label or str(document.path or "")
    matches: list[QueryMatch] = []
    for record in document.records:
        for feature in record.features:
            row = FeatureRow.from_feature(feature, source=source_label)
            if evaluate_query(tree, row):
                matches.append(QueryMatch(row=row, feature=feature, record=record))
    return matches


def parse_select_fields(select: str | None) -> tuple[str, ...] | None:
    if select is None:
        return None
    fields = tuple(item.strip() for item in select.split(",") if item.strip())
    if not fields:
        raise QueryExpressionError("--select must contain at least one field")
    for field in fields:
        if field.startswith("qualifier(") and field.endswith(")"):
            if len(field) <= len("qualifier()"):
                raise QueryExpressionError("qualifier() requires a key")
            continue
        if field not in _QUERY_FIELDS and field not in {"segments", "qualifiers", "xrefs", "xref_sources"}:
            raise QueryExpressionError(f"unknown selected field {field!r}")
    return fields


def search_features(
    filepath: str | Path,
    gene: str | None = None,
    product: str | None = None,
    ko: str | None = None,
    ec: str | None = None,
    cog: str | None = None,
    pfam: str | None = None,
    ftype: str | None = None,
    gene_regex: str | None = None,
    product_regex: str | None = None,
    format_type: str = 'text',
    output_path: str | Path | None = None,
) -> list[dict[str, Any]]:
    doc = read_genbank(filepath)
    results: list[dict[str, Any]] = []

    gene_pat = re.compile(gene_regex, re.IGNORECASE) if gene_regex else None
    prod_pat = re.compile(product_regex, re.IGNORECASE) if product_regex else None

    for f in doc.all_features:
        if ftype and f.type.casefold() != ftype.casefold():
            continue

        f_gene = f.gene or ''
        f_prod = f.product or ''
        xrefs = extract_xrefs(f)

        if gene and gene.casefold() != f_gene.casefold():
            continue
        if product and product.casefold() not in f_prod.casefold():
            continue
        if gene_pat and not gene_pat.search(f_gene):
            continue
        if prod_pat and not prod_pat.search(f_prod):
            continue
        if ko and not any(ko.casefold() == k.casefold() for k in xrefs['kegg_kos']):
            continue
        if ec and not any(ec in e for e in xrefs['ec_numbers']):
            continue
        if cog and not any(cog.casefold() == c.casefold() for c in xrefs['cog_ids']):
            continue
        if pfam and not any(pfam.casefold() == p.casefold() for p in xrefs['pfam']):
            continue

        res = {
            'record': f.record_id,
            'locus_tag': f.locus_tag or '-',
            'type': f.type,
            'gene': f.gene or '-',
            'product': f.product or '-',
            'start': f.start,
            'end': f.end,
            'strand': f.strand_symbol,
            'length': f.length,
            'kegg_ko': ';'.join(xrefs['kegg_kos']),
            'ec_number': ';'.join(xrefs['ec_numbers']),
            'cog': ';'.join(xrefs['cog_ids']),
            'pfam': ';'.join(xrefs['pfam']),
        }
        results.append(res)

    if format_type == 'json':
        out_json = json.dumps(results, indent=2)
        if output_path:
            Path(output_path).write_text(out_json, encoding='utf-8')
        else:
            print(out_json)
        return results

    if format_type in ('tsv', 'csv'):
        delim = '\t' if format_type == 'tsv' else ','
        fields = ['record', 'locus_tag', 'type', 'gene', 'product', 'start', 'end', 'strand', 'length', 'kegg_ko', 'ec_number', 'cog', 'pfam']
        if output_path:
            with Path(output_path).open('w', newline='', encoding='utf-8') as fh:
                writer = csv.DictWriter(fh, fieldnames=fields, delimiter=delim)
                writer.writeheader()
                writer.writerows(results)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter=delim)
            writer.writeheader()
            writer.writerows(results)
        return results

    # Text mode
    print("=" * 80)
    print(f"  FEATURE SEARCH RESULTS: {filepath}")
    print("=" * 80)
    print(f"  Matches found: {len(results)}")
    print()
    if results:
        print(f"{'Record':16s}  {'Strand':6s}  {'Start':>9s}..{'End':<9s}  {'Locus Tag':16s}  {'Gene':8s}  {'Product'}")
        print("-" * 80)
        for r in results:
            print(f"  {r['record']:16s}  [{r['strand']}]    {r['start']:>9,d}..{r['end']:<9,d}  {r['locus_tag']:16s}  {r['gene']:8s}  {r['product'][:35]}")
        print("=" * 80)
    else:
        print("  (No features matched search criteria)")
        print("=" * 80)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Search GenBank annotations by gene, product, KO, EC, Pfam, or regex.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--gene', help="Exact gene name match (e.g. 'ladA')")
    parser.add_argument('--product', help="Substring product match (e.g. 'monooxygenase')")
    parser.add_argument('--ko', help="KEGG KO match (e.g. 'K20938')")
    parser.add_argument('--ec', help="EC number match (e.g. '1.14.14.1')")
    parser.add_argument('--cog', help="COG match (e.g. 'COG0596')")
    parser.add_argument('--pfam', help="Pfam ID match (e.g. 'PF00067')")
    parser.add_argument('--feature', dest='ftype', help="Feature type (e.g. 'CDS', 'tRNA')")
    parser.add_argument('--gene-regex', help="Regex pattern on /gene")
    parser.add_argument('--product-regex', help="Regex pattern on /product")
    parser.add_argument('--format', choices=['text', 'tsv', 'csv', 'json'], default='text', help="Output format")
    parser.add_argument('--output', help="Output file path")
    args = parser.parse_args()

    search_features(
        args.input,
        gene=args.gene,
        product=args.product,
        ko=args.ko,
        ec=args.ec,
        cog=args.cog,
        pfam=args.pfam,
        ftype=args.ftype,
        gene_regex=args.gene_regex,
        product_regex=args.product_regex,
        format_type=args.format,
        output_path=args.output,
    )


if __name__ == '__main__':
    main()
