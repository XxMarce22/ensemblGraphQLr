#' ensemblGraphQLr: modern Ensembl queries through the GraphQL API
#'
#' `ensemblGraphQLr` is a modern, high-performance successor to `biomaRt`. It
#' retrieves Ensembl annotation through the Ensembl GraphQL API instead of the
#' legacy BioMart and REST engines, while keeping the familiar `getBM()`
#' mental model of attributes and filters.
#'
#' The three functions most users need are:
#'
#' * [get_ensembl_data()] - the `getBM()` equivalent.
#' * [build_graphql_query()] - the generated query, for inspection or reuse.
#' * [connect_ensembl()] - a configured `ghql` client, for hand-written
#'   queries.
#'
#' [ensembl_filters()] lists the supported filters, [ensembl_genomes()] and
#' [ensembl_genome_id()] resolve the genome a query is scoped to, and
#' [empty_ensembl_tibble()] describes the shape of an empty result.
#'
#' @section Configuration options:
#' * `ensemblGraphQLr.endpoint` - overrides the service URL.
#' * `ensemblGraphQLr.genome_id` - overrides the genome UUID.
#' * `ensemblGraphQLr.timeout` - seconds allowed per HTTP attempt.
#' * `ensemblGraphQLr.retries` - additional attempts for transient failures.
#' * `ensemblGraphQLr.backoff` - base seconds of the exponential back-off; set
#'   to `0` to disable waiting.
#'
#' @keywords internal
"_PACKAGE"
