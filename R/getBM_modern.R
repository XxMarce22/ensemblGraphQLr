# ---------------------------------------------------------------------------
# getBM_modern.R -- the public `getBM()`-style interface
#
# This file owns the user-facing workflow: validate the request, build the
# query, execute it, flatten the payload and coerce it into a tibble.
# ---------------------------------------------------------------------------

#' Query Ensembl through the GraphQL API
#'
#' `get_ensembl_data()` is the modern counterpart of `biomaRt::getBM()`. It
#' takes the same mental model - a set of *attributes* to retrieve and a
#' *filter* restricting which records are returned - builds a GraphQL query
#' from it, executes the query against the Ensembl GraphQL service and returns
#' the result as a flattened [tibble][tibble::tibble].
#'
#' @section Request workflow:
#' 1. `attributes` are validated and expanded into a GraphQL selection set
#'    (see [build_selection_set()]). Dotted names become nested selection
#'    blocks, so `"metadata.biotype.value"` is queried as
#'    `metadata { biotype { value } }`.
#' 2. `filters` and `values` are translated into GraphQL variables. Values are
#'    never interpolated into the query document, so no user input can change
#'    the shape of the query.
#' 3. The document is assembled by [build_graphql_query()] and validated
#'    locally by the `graphql` parser inside `ghql`.
#' 4. The request is executed by the `ghql` client, with a per-request timeout
#'    and retries with exponential back-off for transient failures.
#' 5. The JSON response is parsed with `jsonlite::fromJSON(flatten = TRUE)` and
#'    reshaped into a tibble with one column per requested attribute, in the
#'    order requested.
#'
#' @section Scoping to a genome:
#' Every Ensembl core lookup is scoped to a genome, identified by a UUID that
#' covers both assembly and release. `genome_id` overrides it; otherwise the
#' `ensemblGraphQLr.genome_id` option is used, and failing that the genome is
#' resolved once per session by [ensembl_genome_id()] (Homo sapiens GRCh38.p14
#' by default).
#'
#' @section One filter per call:
#' The core schema accepts exactly one identifier field per lookup, so exactly
#' one `filter` may be given. Pass a vector of `values` to look up several
#' genes at once: each value becomes its own aliased query block within a
#' single HTTP request, and the results are concatenated.
#'
#' @section Error handling:
#' Every condition signalled by this package inherits from
#' `ensemblGraphQLr_error`, so all of them can be caught at once. The concrete
#' classes are:
#'
#' * `ensembl_input_error` - malformed arguments, such as `values` without
#'   `filters` or more than one filter.
#' * `ensembl_schema_error` - an unknown filter name.
#' * `ensembl_genome_error` - a genome UUID could not be resolved.
#' * `ensembl_connection_error` - a timeout, DNS failure or dropped
#'   connection, after the retries were exhausted.
#' * `ensembl_http_error` - the service answered with a non-success HTTP
#'   status. The status code is available as `condition$status`.
#' * `ensembl_graphql_error` - the service rejected the query; the server-side
#'   messages are available as `condition$errors`. Ensembl reports validation
#'   failures as HTTP 400 with a GraphQL `errors` body, and those are surfaced
#'   here rather than as an HTTP error.
#' * `ensembl_response_error` - the body was not JSON at all, for example an
#'   HTML error page from an intermediate proxy.
#' * `ensembl_empty_response` - the body was empty.
#'
#' @section Empty results:
#' A query that matches nothing returns a zero-row tibble whose columns match
#' `attributes` exactly, with numeric attributes typed as `integer`, rather
#' than `NULL` or an error. See [empty_ensembl_tibble()].
#'
#' @param attributes A character vector of attributes to retrieve. Nested
#'   attributes are expressed with dots, for example
#'   `"metadata.biotype.value"` or `"slice.location.start"`.
#' @param filters `NULL`, or a single biomaRt-style filter name, such as
#'   `"ensembl_gene_id"`, `"hgnc_symbol"` or `"ensembl_transcript_id"`. See
#'   [ensembl_filters()].
#' @param values `NULL`, or a vector of values to match against `filters`.
#' @param genome_id `NULL` to resolve the genome automatically, or a genome
#'   UUID from [ensembl_genomes()].
#' @param endpoint `NULL` to use the resolved default endpoint, or a URL. See
#'   [ensembl_graphql_endpoint()].
#' @param client `NULL` to use the default (cached) client, or a
#'   `ghql::GraphqlClient` from [connect_ensembl()] to reuse a configured
#'   connection.
#' @param timeout Number of seconds to allow for each HTTP attempt. Defaults to
#'   the `ensemblGraphQLr.timeout` option, or 30.
#' @param retries Number of *additional* attempts after the first failure, used
#'   only for transient failures. Defaults to the `ensemblGraphQLr.retries`
#'   option, or 3.
#' @param root `NULL` to take the GraphQL root field from the filter, or an
#'   override such as `"gene"`.
#' @param verbose Logical; if `TRUE`, report each attempt and retry.
#'
#' @return A [tibble][tibble::tibble] with one row per record and one column
#'   per requested attribute. Columns are named after the attributes, with dots
#'   replaced by underscores (`"metadata.biotype.value"` becomes
#'   `metadata_biotype_value`). Genuinely nested values are returned as
#'   list-columns.
#'
#' @examples
#' \dontrun{
#' # One gene, by stable ID
#' get_ensembl_data(
#'     attributes = c("stable_id", "symbol", "name"),
#'     filters = "ensembl_gene_id",
#'     values = "ENSG00000139618"
#' )
#'
#' # Several genes in a single request, looked up by HGNC symbol
#' get_ensembl_data(
#'     attributes = c("stable_id", "symbol", "metadata.biotype.value"),
#'     filters = "hgnc_symbol",
#'     values = c("BRCA2", "BRCA1", "TP53")
#' )
#'
#' # Coordinates live under `slice`
#' get_ensembl_data(
#'     attributes = c("stable_id", "slice.location.start", "slice.location.end"),
#'     filters = "ensembl_gene_id",
#'     values = "ENSG00000139618"
#' )
#'
#' # A query that matches nothing returns a zero-row tibble
#' nrow(get_ensembl_data(
#'     attributes = c("stable_id", "symbol"),
#'     filters = "hgnc_symbol",
#'     values = "NOT_A_REAL_GENE"
#' ))
#'
#' # Inspect the query without sending it
#' cat(build_graphql_query(
#'     attributes = c("stable_id", "symbol"),
#'     filters = "hgnc_symbol",
#'     values = "BRCA2",
#'     genome_id = "59871324-7803-4234-856e-2a2bd96d7b3c"
#' ))
#' }
#'
#' @seealso [connect_ensembl()], [build_graphql_query()],
#'   [build_selection_set()], [ensembl_filters()], [ensembl_genome_id()],
#'   [empty_ensembl_tibble()]
#' @export
get_ensembl_data <- function(attributes, filters = NULL, values = NULL,
                             genome_id = NULL, endpoint = NULL, client = NULL,
                             timeout = NULL, retries = NULL,
                             root = NULL, verbose = FALSE) {
    # Validated up-front so a bad request fails before any network call.
    attributes <- .check_attributes(attributes)

    client <- client %||% connect_ensembl(endpoint = endpoint)
    genome_id <- .resolve_genome_id(
        genome_id, client = client, timeout = timeout, retries = retries,
        verbose = verbose
    )

    plan <- .plan_query(
        attributes = attributes, filters = filters, values = values,
        genome_id = genome_id, root = root
    )

    if (verbose) {
        rlang::inform(c("i" = "Ensembl GraphQL query", " " = plan$query))
    }

    response <- .ensembl_request(
        query = plan$query,
        variables = plan$variables,
        client = client,
        timeout = timeout,
        retries = retries,
        verbose = verbose
    )

    parsed <- jsonlite::fromJSON(response, flatten = TRUE)

    data <- parsed$data
    if (is.null(data) || !is.list(data)) {
        data <- list()
    }

    # Aliases the service reported as not-found are expected to be missing, so
    # they are not evidence of a schema mismatch.
    not_found <- attr(response, "ensembl_not_found") %||% character()
    absent <- setdiff(plan$root_fields, c(names(data), not_found))
    if (length(absent) > 0L) {
        rlang::warn(
            c(
                sprintf(
                    "The response contains no %s field, so an empty result is returned.",
                    paste(encodeString(absent, quote = "`"), collapse = ", ")
                ),
                "i" = "This usually means the GraphQL schema does not expose this root field; the query itself succeeded."
            ),
            class = "ensembl_missing_field"
        )
    }

    # Each requested value produced its own aliased root field, so the records
    # are gathered back into a single result set.
    rows <- list()
    for (field in plan$root_fields) {
        rows <- c(rows, .to_record_rows(data[[field]]))
    }

    out <- .rows_to_tibble(rows, plan$attributes)
    .warn_missing_attributes(out, plan$attributes)

    out
}

# Warn when the service accepted the query but never populated an attribute.
# Silently returning an all-`NA` column would look like real missing data, so
# the mismatch is surfaced explicitly. Empty results carry no information about
# which attributes exist, so they are exempt.
.warn_missing_attributes <- function(out, attributes) {
    if (nrow(out) == 0L) {
        return(invisible(out))
    }
    columns <- names(out)
    empty <- vapply(out, function(column) {
        !is.list(column) && all(is.na(column))
    }, logical(1))
    if (any(empty)) {
        rlang::warn(
            sprintf(
                "Attribute(s) not returned by the Ensembl GraphQL API: %s.",
                paste(columns[empty], collapse = ", ")
            ),
            class = "ensembl_missing_attributes"
        )
    }
    invisible(out)
}
