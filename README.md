# ensemblGraphQLr

[![R-CMD-check](https://github.com/XxMarce22/ensemblGraphQLr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/XxMarce22/ensemblGraphQLr/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

> Query Ensembl annotation through Ensembl's **GraphQL** API, using a
> `biomaRt::getBM()`-style interface of attributes, filters and values.

`ensemblGraphQLr` is a modern successor to `biomaRt`. It retrieves Ensembl
annotation from the GraphQL service that powers Ensembl's own web interface,
instead of the legacy BioMart and REST engines, and returns results as flattened
tibbles.

```r
library(ensemblGraphQLr)

get_ensembl_data(
  attributes = c("stable_id", "symbol", "name", "metadata.biotype.value"),
  filters    = "hgnc_symbol",
  values     = "BRCA2"
)
#> # A tibble: 1 x 4
#>   stable_id          symbol name                        metadata_biotype_value
#>   <chr>              <chr>  <chr>                       <chr>
#> 1 ENSG00000139618.19 BRCA2  BRCA2 DNA repair associated protein_coding
```

---

## ⚠️ Status

Ensembl's GraphQL service is **beta**. Ensembl's own documentation states that
breaking changes may be deployed without notice and advises against building
anything that requires stability on it. Accordingly this package is marked
*experimental*: treat the API as liable to change, and pin a commit if you need
reproducibility.

This package is **not affiliated with or endorsed by EMBL-EBI**.

## Installation

Not on CRAN; install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("XxMarce22/ensemblGraphQLr")
```

Requires R >= 4.1.0. Dependencies: `ghql`, `jsonlite`, `tibble`, `rlang`.

## Quick start

Three genes in **one** HTTP request:

```r
get_ensembl_data(
  attributes = c("stable_id", "symbol", "slice.region.name"),
  filters    = "hgnc_symbol",
  values     = c("BRCA2", "TP53", "BRCA1")
)
#> # A tibble: 3 x 3
#>   stable_id          symbol slice_region_name
#>   <chr>              <chr>  <chr>
#> 1 ENSG00000139618.19 BRCA2  13
#> 2 ENSG00000141510.21 TP53   17
#> 3 ENSG00000012048.28 BRCA1  17
```

Coordinates, which live under `slice`:

```r
get_ensembl_data(
  attributes = c("stable_id", "slice.location.start",
                 "slice.location.end", "slice.strand.code"),
  filters    = "ensembl_gene_id",
  values     = "ENSG00000139618"
)
#> # A tibble: 1 x 4
#>   stable_id          slice_location_start slice_location_end slice_strand_code
#>   <chr>                             <int>              <int> <chr>
#> 1 ENSG00000139618.19             32315086           32400268 forward
```

A query that matches nothing is a normal empty result, not an error:

```r
get_ensembl_data(
  attributes = c("stable_id", "symbol", "slice.location.start"),
  filters    = "hgnc_symbol",
  values     = "NOT_A_REAL_GENE"
)
#> # A tibble: 0 x 3
#> # i 3 variables: stable_id <chr>, symbol <chr>, slice_location_start <int>
```

Note the column is typed `int`, so empty and non-empty results bind together.

## How it works

1. `attributes` are validated and expanded into a GraphQL selection set. Dotted
   names become nested blocks, so `"metadata.biotype.value"` is queried as
   `metadata { biotype { value } }`.
2. `filters` and `values` become GraphQL **variables**. User values are never
   interpolated into the query document, so no input can alter its shape.
3. `build_graphql_query()` assembles the document, which `ghql` validates
   locally before anything is sent.
4. The request is executed with a per-request timeout and retries with
   exponential back-off for transient failures.
5. The JSON payload is parsed with `jsonlite::fromJSON(flatten = TRUE)` and
   reshaped into a tibble, one column per requested attribute, in the order
   requested.

Inspect the generated query without sending anything:

```r
cat(build_graphql_query(
  attributes = c("stable_id", "symbol"),
  filters    = "hgnc_symbol",
  values     = c("BRCA2", "TP53"),
  genome_id  = "59871324-7803-4234-856e-2a2bd96d7b3c"
))
```

```graphql
query EnsemblQuery(
  $genomeId: String!,
  $value1: String!,
  $value2: String!
) {
  result1: genes(by_symbol: {genome_id: $genomeId, symbol: $value1}) {
    stable_id
    symbol
  }
  result2: genes(by_symbol: {genome_id: $genomeId, symbol: $value2}) {
    stable_id
    symbol
  }
}
```

Each value becomes its own aliased block, so several lookups share one request.

## Scoping: `genome_id`

Every Ensembl core lookup is scoped to a genome, identified by a **UUID** that
covers both assembly *and* release. `ensemblGraphQLr` resolves it for you once
per session — by default to Homo sapiens GRCh38.p14 — and caches the result.

```r
ensembl_genome_id()                    # "59871324-7803-4234-856e-2a2bd96d7b3c"
ensembl_genomes("GCA_000001405.29")    # discover genomes by accession
```

To query another assembly or release, resolve it and pass it in:

```r
mouse <- ensembl_genomes("GCA_000001635.9")$genome_id

get_ensembl_data(
  attributes = c("stable_id", "symbol"),
  filters    = "hgnc_symbol",
  values     = "Trp53",
  genome_id  = mouse
)
```

## Attributes

Use dots for nesting. These paths are verified against the live core API:

| Group | Attribute paths |
|---|---|
| Identity | `stable_id`, `unversioned_stable_id`, `version`, `symbol`, `name`, `alternative_symbols` |
| Biotype | `metadata.biotype.value`, `metadata.biotype.label`, `metadata.biotype.definition` |
| Coordinates | `slice.location.start`, `slice.location.end`, `slice.location.length`, `slice.strand.code`, `slice.strand.value`, `slice.region.name`, `slice.region.length`, `slice.region.code` |
| References | `metadata.name.accession_id`, `metadata.name.url`, `external_references.accession_id`, `external_references.name` |
| Transcripts | `transcripts.stable_id` |

Structural attributes must name a sub-field. Asking for `"transcripts"` alone is
rejected by GraphQL; use `"transcripts.stable_id"`. The package detects this and
tells you so.

Nested **arrays** are returned as list-columns rather than being flattened into
one row per record, because an attribute set such as
`c("stable_id", "transcripts.stable_id")` has no single rectangular shape:

```r
get_ensembl_data("alternative_symbols", "ensembl_gene_id", "ENSG00000139618")
#> # A tibble: 1 x 1
#>   alternative_symbols
#>   <list>
#> 1 <chr [7]>
```

## Filters

The core schema accepts **one identifier field per lookup**, so one filter per
call. `ensembl_filters()` reports what is available:

| Filter | GraphQL field | Argument | Returns |
|---|---|---|---|
| `ensembl_gene_id` | `gene` | `by_id` | one |
| `hgnc_symbol` | `genes` | `by_symbol` | many |
| `ensembl_transcript_id` | `transcript` | `by_id` | one |

Pass a vector of `values` to look up many records at once; they are batched into
one request and the results concatenated. Missing entries are dropped rather
than failing the batch:

```r
get_ensembl_data(c("stable_id", "symbol"), "hgnc_symbol", c("BRCA2", "NOT_A_GENE"))
#> # A tibble: 1 x 2
#>   stable_id          symbol
#>   <chr>              <chr>
#> 1 ENSG00000139618.19 BRCA2
```

`hgnc_symbol` targets the schema's `symbol` field, so it also resolves symbols
in other species — `"Trp53"` against the mouse genome works. The name reflects
the common human case, not a restriction.

## Error handling

Every condition inherits from `ensemblGraphQLr_error`, so all of them can be
caught at once.

| Class | Meaning |
|---|---|
| `ensembl_input_error` | Malformed arguments |
| `ensembl_schema_error` | Unknown filter name |
| `ensembl_genome_error` | Genome UUID could not be resolved |
| `ensembl_connection_error` | Timeout / DNS failure / dropped connection |
| `ensembl_http_error` | Non-success status; `condition$status` |
| `ensembl_graphql_error` | Server rejected the query; `condition$errors` |
| `ensembl_response_error` | Body was not JSON (e.g. a proxy error page) |
| `ensembl_empty_response` | Body was empty |

```r
tryCatch(
  get_ensembl_data("biotype", "hgnc_symbol", "BRCA2"),
  ensembl_graphql_error = function(e) e$errors,
  ensembl_http_error    = function(e) e$status,
  ensemblGraphQLr_error = function(e) conditionMessage(e)
)
```

Transient failures (timeouts, connection resets, HTTP 408/425/429/5xx) are
retried automatically. `*_NOT_FOUND` responses from Ensembl are treated as empty
results, *not* errors.

## Configuration

Set via `options()`:

| Option | Default | Purpose |
|---|---|---|
| `ensemblGraphQLr.endpoint` | `https://www.ensembl.org/api/graphql/core` | Service URL |
| `ensemblGraphQLr.genome_id` | auto-resolved | Genome UUID |
| `ensemblGraphQLr.timeout` | `30` | Seconds per attempt |
| `ensemblGraphQLr.retries` | `3` | Retries for transient failures |
| `ensemblGraphQLr.backoff` | `0.5` | Base seconds of exponential back-off; `0` disables waiting |

Ensembl's GraphQL service is split per data type. This package targets the
**core** service (`/api/graphql/core`); `/api/graphql/variation` and
`/api/graphql/compara` are siblings. The bare `/api/graphql` path is *not* the
GraphQL service. For hand-written queries, `connect_ensembl()` returns a
configured client:

```r
con <- connect_ensembl()
q <- ghql::Query$new()
q$query("g", '{ version { api { major minor } } }')
con$exec(q$queries$g)
```

## Compared with biomaRt

| `biomaRt` | `ensemblGraphQLr` |
|---|---|
| `useEnsembl(biomart = "ensembl")` | `connect_ensembl()` |
| dataset `hsapiens_gene_ensembl` | `genome_id` UUID |
| `getBM(attributes, filters, values, mart)` | `get_ensembl_data(attributes, filters, values)` |
| multiple filters combined with AND | one filter per call |
| `listFilters()` / `listAttributes()` | `ensembl_filters()` |
| legacy BioMart engine | Ensembl GraphQL API |

### Known limitations

- **One filter per call.** The core schema accepts a single identifier field per
  lookup, so `getBM()`-style multi-filter AND semantics are not expressible.
- **Nested values stay nested**, as list-columns, rather than one row per
  transcript or exon. Query transcripts directly via
  `filters = "ensembl_transcript_id"` instead.
- **No expression or variation data** through the core service.
- **`stable_id` is versioned.** Expect `ENSG00000139618.19`. Request
  `unversioned_stable_id`, or strip with `sub("\\..*$", "", x)`, when joining
  against unversioned identifiers.
- **Beta upstream.** See the status note above.

## Development

```bash
# Run the test suite (offline; the live test is skipped by default)
Rscript -e 'testthat::test_local()'

# Include the live tests against the real Ensembl API
NOT_CRAN=true ENSEMBL_GRAPHQL_LIVE=1 Rscript -e 'devtools::test()'

# Prove the suite is hermetic: everything must still pass
ENSEMBL_GRAPHQL_ENDPOINT="http://127.0.0.1:9/closed" Rscript -e 'devtools::test()'

# Full package check
Rscript -e 'devtools::check(args = "--no-manual")'
```

The suite makes no network calls unless `ENSEMBL_GRAPHQL_LIVE` is set.

## License

MIT © Marcelo Bertuol. See [LICENSE.md](LICENSE.md).
