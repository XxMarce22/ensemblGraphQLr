# ensemblGraphQLr 0.0.0.9000

## Breaking changes

* Initial development version; the API is not yet stable.

## Features

* `get_ensembl_data()` queries Ensembl annotation through the Ensembl GraphQL
  API using a `biomaRt::getBM()`-style interface of attributes, filters and
  values, and returns a flattened tibble.
* `build_graphql_query()` and `build_selection_set()` expose the generated
  GraphQL document for inspection and reuse without sending a request.
* `connect_ensembl()` returns a configured `ghql` GraphQL client for
  hand-written queries.
* `ensembl_filters()` lists the supported filters;
  `ensembl_genome_id()` and `ensembl_genomes()` resolve the genome UUID that
  scopes every core lookup.
* Dotted attribute names expand into nested GraphQL selection sets, so
  `"metadata.biotype.value"` is queried as `metadata { biotype { value } }`.
* Several `values` for one filter are batched into a single HTTP request using
  aliased root fields.
* Transport failures are retried with exponential back-off; timeouts, HTTP
  errors, server-side GraphQL errors, non-JSON bodies and empty bodies are
  raised as distinct classed conditions.
* Query results that match nothing return a zero-row tibble with the requested
  columns and sensible types, rather than `NULL` or an error. Lookups that
  Ensembl reports as `*_NOT_FOUND` are treated the same way, and partial
  matches are preserved.
* Filter values are always sent as GraphQL variables, never interpolated into
  the query document.

## Known limitations

* The Ensembl core schema accepts exactly one identifier field per lookup, so
  one filter per call is supported. Multiple filters with AND semantics, as
  accepted by `biomaRt::getBM()`, are not expressible.
* Nested structural attributes such as `transcripts` are returned as
  list-columns rather than being flattened into one row per record.
