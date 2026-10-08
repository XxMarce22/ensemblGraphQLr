# ---------------------------------------------------------------------------
# client.R -- Ensembl GraphQL transport layer
#
# This file owns everything related to *transport*: resolving the service
# endpoint, constructing a `ghql` client, executing requests, retrying
# transient failures and translating transport problems into classed,
# actionable R conditions.
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Known endpoints ------------------------------------------------------

# Ensembl's new architecture exposes a GraphQL service per data type under the
# main host. The core (genes, transcripts, genomes, regions) service lives at
# `/api/graphql/core`; `.../variation` and `.../compara` are siblings.
#
# The endpoint is a POST endpoint expecting a JSON body of the form
# `{"query": "...", "variables": {...}}`, which is exactly what
# `ghql::GraphqlClient$exec()` produces. No authentication is required.
#
# Verified 2026-10-08 against core API v0.2.0-beta (Ensembl release 116).
# Note that the bare `/api/graphql` and `/api/v1/graphql` paths are *not* the
# GraphQL service: they are handled by a legacy-URL resolver that advertises
# `allow: GET` and rejects POST with HTTP 405.
.ensembl_default_endpoint <- "https://www.ensembl.org/api/graphql/core"

# Environment used to cache resolved settings between calls, so that repeated
# `get_ensembl_data()` calls reuse a single client instead of rebuilding one.
.ensembl_cache <- new.env(parent = emptyenv())

#' Resolve the Ensembl GraphQL endpoint
#'
#' Endpoint resolution happens in the following order, first match wins:
#'
#' 1. the `endpoint` argument, when supplied;
#' 2. the `ensemblGraphQLr.endpoint` R option;
#' 3. the `ENSEMBL_GRAPHQL_ENDPOINT` environment variable;
#' 4. the package default, `https://www.ensembl.org/api/graphql/core`.
#'
#' The indirection exists because Ensembl moves its GraphQL entry points
#' between releases, and because the service is split per data type: the core
#' service used by [get_ensembl_data()] lives at `/api/graphql/core`, with
#' `/api/graphql/variation` and `/api/graphql/compara` as siblings. Setting the
#' option (for example in your `.Rprofile`) re-points the package at a new host
#' without reinstalling it:
#'
#' ```r
#' options(ensemblGraphQLr.endpoint = "https://example.org/graphql")
#' ```
#'
#' @section Transport requirements:
#' The client sends a GraphQL request as an HTTP `POST` with a JSON body of the
#' form `{"query": "...", "variables": {...}}`, which is what `ghql` produces.
#' An endpoint must accept that, and must accept *arbitrary* query documents
#' rather than only pre-registered (persisted) queries. Ensembl's core,
#' variation and compara services all satisfy this; the bare `/api/graphql` and
#' `/api/v1/graphql` paths do not, as they advertise `allow: GET` and answer
#' `POST` with HTTP 405.
#'
#' @param endpoint `NULL`, or a single non-empty string giving the full URL of
#'   the GraphQL service, including scheme and path.
#'
#' @return A single string: the endpoint that will be used.
#'
#' @examples
#' ensembl_graphql_endpoint()
#' ensembl_graphql_endpoint("https://example.org/graphql")
#'
#' @seealso [connect_ensembl()]
#' @export
ensembl_graphql_endpoint <- function(endpoint = NULL) {
    if (!is.null(endpoint)) {
        return(.check_endpoint(endpoint))
    }
    from_option <- getOption("ensemblGraphQLr.endpoint")
    if (!is.null(from_option)) {
        return(.check_endpoint(from_option))
    }
    from_env <- Sys.getenv("ENSEMBL_GRAPHQL_ENDPOINT", unset = NA_character_)
    if (!is.na(from_env) && nzchar(from_env)) {
        return(.check_endpoint(from_env))
    }
    .ensembl_default_endpoint
}

.check_endpoint <- function(endpoint) {
    if (!is.character(endpoint) || length(endpoint) != 1L || is.na(endpoint)) {
        .ensembl_abort("`endpoint` must be a single non-missing string.")
    }
    endpoint <- trimws(endpoint)
    if (!nzchar(endpoint)) {
        .ensembl_abort("`endpoint` must not be an empty string.")
    }
    if (!grepl("^https?://", endpoint)) {
        .ensembl_abort(sprintf(
            "`endpoint` must be an absolute `http://` or `https://` URL, not %s.",
            encodeString(endpoint, quote = "\"")
        ))
    }
    sub("/+$", "", endpoint)
}

# ---- Client construction --------------------------------------------------

#' Connect to the Ensembl GraphQL service
#'
#' Creates a `ghql` GraphQL client pointed at the Ensembl GraphQL endpoint. The
#' returned object is an ordinary [`ghql::GraphqlClient`], so it can also be
#' used directly for hand-written queries.
#'
#' The client itself carries no timeout or retry configuration: `ghql` forwards
#' extra arguments of `$exec()` down to `crul`, so timeouts and retries are
#' applied per request by `get_ensembl_data()` (see its `timeout` and `retries`
#' arguments).
#'
#' @param endpoint `NULL` to use the resolved default, or a URL. See
#'   [ensembl_graphql_endpoint()].
#' @param headers A named list of additional HTTP headers. A `User-Agent`
#'   header identifying the package is always sent and can be overridden here.
#' @param cache Logical. If `TRUE` (the default) the client is cached and
#'   reused by later calls that use the same endpoint and headers. Set to
#'   `FALSE` to always build a fresh client.
#'
#' @return A `ghql::GraphqlClient` object.
#'
#' @examples
#' \dontrun{
#' con <- connect_ensembl()
#' con
#' }
#'
#' @seealso [get_ensembl_data()], [ensembl_graphql_endpoint()]
#' @export
connect_ensembl <- function(endpoint = NULL, headers = list(), cache = TRUE) {
    endpoint <- ensembl_graphql_endpoint(endpoint)

    if (!is.list(headers)) {
        .ensembl_abort("`headers` must be a named list of HTTP headers.")
    }
    if (length(headers) > 0L && (is.null(names(headers)) || any(!nzchar(names(headers))))) {
        .ensembl_abort("`headers` must be a *named* list, e.g. `list(Authorization = \"...\")`.")
    }
    headers <- utils::modifyList(list(`User-Agent` = .ensembl_user_agent()), headers)

    # Only an unmodified header set can safely share the cached client.
    cacheable <- isTRUE(cache) && identical(headers, list(`User-Agent` = .ensembl_user_agent()))
    if (cacheable) {
        cached <- .ensembl_cache$client
        if (!is.null(cached) && identical(.ensembl_cache$endpoint, endpoint)) {
            return(cached)
        }
    }

    client <- ghql::GraphqlClient$new(url = endpoint, headers = headers)

    if (cacheable) {
        .ensembl_cache$client <- client
        .ensembl_cache$endpoint <- endpoint
    }
    client
}

.ensembl_user_agent <- function() {
    version <- tryCatch(
        as.character(utils::packageVersion("ensemblGraphQLr")),
        error = function(e) "dev"
    )
    paste0("ensemblGraphQLr/", version)
}

.ensembl_client_endpoint <- function(client) {
    if (is.null(client$url)) "<unknown>" else client$url
}

# ---- Request execution ----------------------------------------------------

# Transient conditions worth retrying. Matched case-insensitively against the
# message raised by `curl`/`crul` or by `ghql`'s HTTP status handler.
.ensembl_transient_messages <- c(
    "timeout was reached", "operation timed out", "resolving timed out",
    "connection timed out", "timed out",
    "could not resolve host", "name or service not known",
    "failed to connect", "connection refused", "connection reset",
    "recv failure", "send failure", "empty reply from server",
    "broken pipe", "network is unreachable", "ssl connect error",
    "unexpected eof", "transfer closed"
)

# HTTP statuses that are worth another attempt.
.ensembl_retry_status <- c(408L, 425L, 429L, 500L, 502L, 503L, 504L)

.check_timeout <- function(timeout, default = 30) {
    if (is.null(timeout)) {
        timeout <- getOption("ensemblGraphQLr.timeout", default)
    }
    if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
        .ensembl_abort("`timeout` must be a single positive number of seconds.")
    }
    as.numeric(timeout)
}

.check_retries <- function(retries, default = 3L) {
    if (is.null(retries)) {
        retries <- getOption("ensemblGraphQLr.retries", default)
    }
    if (!is.numeric(retries) || length(retries) != 1L || is.na(retries) || retries < 0) {
        .ensembl_abort("`retries` must be a single non-negative number.")
    }
    as.integer(retries)
}

# Exponential back-off, capped, and disabled entirely by setting the
# `ensemblGraphQLr.backoff` option to 0 (useful in tests).
.backoff_seconds <- function(attempt) {
    base <- getOption("ensemblGraphQLr.backoff", 0.5)
    if (!is.numeric(base) || length(base) != 1L || is.na(base) || base <= 0) {
        return(0)
    }
    min(base * 2^(attempt - 1L), 15)
}

#' Execute a GraphQL query against the Ensembl service
#'
#' Returns the raw JSON response body as a string. Transport failures are
#' retried with exponential back-off when they are transient; all failures are
#' re-raised as classed conditions (see [get_ensembl_data()] for the taxonomy).
#'
#' @param query A single string containing a syntactically valid GraphQL
#'   document. Syntax is validated locally by `ghql` before the request is
#'   sent.
#' @param variables A named list of GraphQL variable values.
#' @param client A `ghql::GraphqlClient`, or `NULL` to build the default one.
#' @param timeout,retries Request settings, or `NULL` for the resolved
#'   defaults.
#' @param verbose Logical; if `TRUE`, report each attempt and retry.
#'
#' @return A single string: the JSON response body.
#'
#' @noRd
.ensembl_request <- function(query, variables = list(), client = NULL,
                             timeout = NULL, retries = NULL, verbose = FALSE) {
    client <- client %||% connect_ensembl()
    endpoint <- .ensembl_client_endpoint(client)
    timeout <- .check_timeout(timeout)
    retries <- .check_retries(retries)

    attempt <- 0L
    repeat {
        attempt <- attempt + 1L
        if (verbose) {
            rlang::inform(sprintf("Ensembl GraphQL request to %s (attempt %d).", endpoint, attempt))
        }

        result <- tryCatch(
            .ensembl_exec_once(client, query, variables, timeout),
            error = function(e) e
        )

        if (!inherits(result, "error")) {
            return(.ensembl_validate_response(result, endpoint))
        }

        failure <- .classify_transport_error(result, endpoint)

        # A GraphQL validation failure can arrive as an HTTP 400 whose body
        # still carries the server-side `errors` array. Reporting those
        # messages beats reporting the status code, and such failures are
        # deterministic, so they are never retried.
        gql_errors <- .ensembl_graphql_errors_from_message(conditionMessage(result))
        if (length(gql_errors) > 0L) {
            .ensembl_raise_graphql_error(gql_errors, endpoint)
        }

        if (!failure$retryable || attempt > retries) {
            .ensembl_raise_transport_error(failure, result, endpoint, attempt, retries)
        }

        wait <- .backoff_seconds(attempt)
        if (verbose) {
            rlang::inform(sprintf(
                "Transient failure (%s); retrying in %.1fs.", failure$kind, wait
            ))
        }
        if (wait > 0) {
            Sys.sleep(wait)
        }
    }
}

# `ghql` only accepts a `Query` object (not a bare string), and registering the
# query runs it through the `graphql` parser, which gives us free local syntax
# validation. The `timeout` is forwarded through `...` to `crul`'s curl handle.
.ensembl_exec_once <- function(client, query, variables, timeout) {
    if (!is.character(query) || length(query) != 1L || is.na(query)) {
        .ensembl_abort("`query` must be a single string holding a GraphQL document.")
    }
    qry <- ghql::Query$new()
    qry$query("ensemblGraphQLrQuery", query)

    if (length(variables) > 0L) {
        client$exec(qry$queries$ensemblGraphQLrQuery, variables = variables, timeout = timeout)
    } else {
        client$exec(qry$queries$ensemblGraphQLrQuery, timeout = timeout)
    }
}

# ---- Response validation --------------------------------------------------

# A full parse of every payload would double the cost of large responses, so
# the cheap `grepl()` guard is checked first and the payload is only parsed
# when it actually advertises GraphQL-level errors.
.ensembl_validate_response <- function(txt, endpoint) {
    if (is.null(txt) || !nzchar(trimws(txt))) {
        .ensembl_abort(
            sprintf("The Ensembl GraphQL endpoint (%s) returned an empty response body.", endpoint),
            class = "ensembl_empty_response"
        )
    }

    body <- trimws(txt)
    if (!startsWith(body, "{") && !startsWith(body, "[")) {
        .ensembl_abort(
            sprintf(
                "The Ensembl GraphQL endpoint (%s) returned a non-JSON response. First 200 characters:\n%s",
                endpoint, substr(body, 1L, 200L)
            ),
            class = "ensembl_response_error"
        )
    }

    if (grepl("\"errors\"", body, fixed = TRUE)) {
        parsed <- tryCatch(
            jsonlite::fromJSON(body, simplifyVector = FALSE),
            error = function(e) NULL
        )
        errors <- parsed$errors
        messages <- .ensembl_error_messages(errors)
        if (length(messages) > 0L) {
            # A lookup that matches nothing is a normal empty result, not a
            # failure: Ensembl reports it as a `*_NOT_FOUND` error and still
            # returns whatever other aliased lookups did match. Those failed
            # aliases are recorded so callers can tell "no such gene" apart from
            # "the schema has no such field".
            if (.ensembl_not_found(errors)) {
                attr(txt, "ensembl_not_found") <- .ensembl_error_paths(errors)
                return(txt)
            }
            .ensembl_raise_graphql_error(messages, endpoint)
        }
    }

    txt
}

.ensembl_error_codes <- function(errors) {
    vapply(errors, function(e) as.character(e$extensions$code %||% NA_character_), character(1))
}

.ensembl_not_found <- function(errors) {
    if (length(errors) == 0L) {
        return(FALSE)
    }
    codes <- .ensembl_error_codes(errors)
    if (all(!is.na(codes))) {
        return(all(grepl("NOT_FOUND$", codes)))
    }
    # Fall back to the message when the service omits the extension code.
    all(grepl("NOT_FOUND", .ensembl_error_messages(errors), fixed = TRUE))
}

# Every GraphQL error carries the field it failed on, which lines up with the
# alias names used in the query document.
.ensembl_error_paths <- function(errors) {
    paths <- lapply(errors, function(e) e$path)
    unique(unlist(
        lapply(paths, function(p) if (is.null(p)) character() else as.character(unlist(p, use.names = FALSE))),
        use.names = FALSE
    ))
}

# `ghql` raises HTTP failures as `stop(sprintf("[HTTP %s]", code), " ", body)`,
# so the JSON body is recoverable by stripping the status prefix.
.ensembl_graphql_errors_from_message <- function(message) {
    body <- sub("^\\[HTTP [0-9]{3}\\]\\s*", "", message)
    if (identical(body, message) || !startsWith(trimws(body), "{")) {
        return(character())
    }
    parsed <- tryCatch(
        jsonlite::fromJSON(body, simplifyVector = FALSE),
        error = function(e) NULL
    )
    .ensembl_error_messages(parsed$errors)
}

.ensembl_raise_graphql_error <- function(messages, endpoint) {
    hints <- c(
        "Inspect the generated query with `build_graphql_query()` and check the attribute and filter names.",
        "Attribute paths must match the schema, for example \"metadata.biotype.value\" rather than \"biotype\"."
    )

    # GraphQL requires a sub-selection for object-valued fields, which is the
    # one schema rule users hit most often when naming a structural attribute
    # on its own.
    if (any(grepl("selection of subfields", messages, fixed = TRUE))) {
        hints <- append(
            hints,
            "Structural attributes need a sub-field, for example \"transcripts.stable_id\" or \"external_references.accession_id\".",
            after = 1L
        )
    }

    .ensembl_abort(
        c(
            sprintf("The Ensembl GraphQL service (%s) rejected the query.", endpoint),
            stats::setNames(messages, rep("x", length(messages))),
            stats::setNames(hints, rep("i", length(hints)))
        ),
        class = "ensembl_graphql_error",
        errors = messages,
        endpoint = endpoint
    )
}

.ensembl_error_messages <- function(errors) {
    if (is.null(errors) || length(errors) == 0L) {
        return(character())
    }
    vapply(errors, function(e) {
        if (is.character(e)) {
            return(e)
        }
        msg <- e$message %||% "unspecified GraphQL error"
        code <- e$extensions$code
        if (is.null(code) || !nzchar(code)) msg else sprintf("%s [%s]", msg, code)
    }, character(1))
}

# ---- Error classification -------------------------------------------------

.classify_transport_error <- function(error, endpoint) {
    message <- conditionMessage(error)
    lowered <- tolower(message)

    status <- NA_integer_
    status_match <- regmatches(message, regexpr("\\[HTTP [0-9]{3}\\]", message))
    if (length(status_match) == 1L) {
        status <- as.integer(gsub("[^0-9]", "", status_match))
    }

    if (!is.na(status)) {
        return(list(
            kind = sprintf("HTTP %d", status),
            retryable = status %in% .ensembl_retry_status,
            status = status,
            endpoint = endpoint
        ))
    }

    transient <- any(vapply(
        .ensembl_transient_messages,
        function(pattern) grepl(pattern, lowered, fixed = TRUE),
        logical(1)
    ))

    list(
        kind = if (transient) "network" else "client",
        retryable = transient,
        status = NA_integer_,
        endpoint = endpoint
    )
}

.ensembl_raise_transport_error <- function(failure, error, endpoint, attempt, retries) {
    original <- conditionMessage(error)
    attempts <- sprintf("Gave up after %d attempt(s).", attempt)

    if (!is.na(failure$status)) {
        .ensembl_abort(
            c(
                sprintf(
                    "The Ensembl GraphQL service (%s) returned HTTP %d.",
                    endpoint, failure$status
                ),
                "x" = original,
                "i" = if (failure$retryable) attempts else "This status is not retried.",
                "i" = "Check the endpoint with `ensembl_graphql_endpoint()` and retry later."
            ),
            class = "ensembl_http_error",
            status = failure$status,
            endpoint = endpoint,
            attempts = attempt
        )
    }

    .ensembl_abort(
        c(
            sprintf("Could not reach the Ensembl GraphQL service (%s).", endpoint),
            "x" = original,
            "i" = if (failure$retryable) {
                sprintf("%s The connection appears to have timed out or dropped.", attempts)
            } else {
                "The request could not be completed."
            },
            "i" = "Increase `timeout`, raise `retries`, or check network connectivity."
        ),
        class = "ensembl_connection_error",
        endpoint = endpoint,
        attempts = attempt
    )
}

# ---- Genome discovery -----------------------------------------------------

# Fields of `GenomeBySpecificKeywordInput` verified against the core schema.
# `parlance_name`, `taxon_id` and `release_number` are *not* fields of this
# input type, despite appearing in some third-party notes.
.ensembl_genome_keywords <- c("assembly_accession_id", "scientific_name", "tolid")

# Homo sapiens GRCh38.p14. Used only to resolve a genome UUID on first use.
.ensembl_default_assembly <- "GCA_000001405.29"

.genomes_attributes <- c(
    "genome_id", "genome_tag", "release_number",
    "assembly.accession_id", "assembly.name"
)

# `by` is a schema field name, not user data, so it is interpolated into the
# document after being checked against the whitelist above.
.genomes_query <- function(by) {
    sprintf(
        paste0(
            "query EnsemblGenomes($keyword: String!) {\n",
            "  genomes(by_keyword: {%s: $keyword}) {\n",
            "    genome_id\n",
            "    genome_tag\n",
            "    release_number\n",
            "    assembly {\n",
            "      accession_id\n",
            "      name\n",
            "    }\n",
            "  }\n",
            "}"
        ),
        by
    )
}

#' List Ensembl genomes matching a keyword
#'
#' Every Ensembl core GraphQL lookup is scoped to a `genome_id`, which is a
#' UUID identifying an assembly *and* release. This function finds those UUIDs.
#'
#' @param keyword One or more search terms. What they mean depends on `by`.
#' @param by Which field to match, one of `"assembly_accession_id"`,
#'   `"scientific_name"` or `"tolid"`.
#' @param endpoint `NULL` to use the resolved default endpoint, or a URL.
#' @param client `NULL` to use the default client, or a `ghql::GraphqlClient`.
#' @param timeout,retries,verbose Passed through to the transport; see
#'   [get_ensembl_data()].
#'
#' @return A [tibble][tibble::tibble] with the columns `genome_id`,
#'   `genome_tag`, `release_number`, `assembly_accession_id` and
#'   `assembly_name`. Genomes that cannot be found are omitted, so an
#'   unmatched keyword yields a zero-row tibble rather than an error.
#'
#' @examples
#' \dontrun{
#' ensembl_genomes("GCA_000001405.29")
#' ensembl_genomes(c("GCA_000001405.29", "GCA_000001635.9"))
#' }
#'
#' @seealso [ensembl_genome_id()]
#' @export
ensembl_genomes <- function(keyword = .ensembl_default_assembly,
                            by = "assembly_accession_id", endpoint = NULL,
                            client = NULL, timeout = NULL, retries = NULL,
                            verbose = FALSE) {
    if (!is.character(by) || length(by) != 1L || is.na(by) ||
        !by %in% .ensembl_genome_keywords) {
        .ensembl_abort(
            sprintf(
                "`by` must be one of %s.",
                paste(encodeString(.ensembl_genome_keywords, quote = "\""), collapse = ", ")
            ),
            class = "ensembl_input_error"
        )
    }
    if (!is.character(keyword) || length(keyword) == 0L || anyNA(keyword)) {
        .ensembl_abort("`keyword` must be a non-missing character vector.")
    }
    keyword <- trimws(keyword)
    if (any(!nzchar(keyword))) {
        .ensembl_abort("`keyword` must not contain empty or whitespace-only entries.")
    }

    client <- client %||% connect_ensembl(endpoint = endpoint)
    query <- .genomes_query(by)

    results <- lapply(keyword, function(value) {
        response <- tryCatch(
            .ensembl_request(
                query = query,
                variables = list(keyword = value),
                client = client,
                timeout = timeout,
                retries = retries,
                verbose = verbose
            ),
            # A genome that does not exist is a normal empty result, not a
            # failure, so `*_NOT_FOUND` is downgraded to zero rows.
            ensembl_graphql_error = function(e) {
                if (grepl("NOT_FOUND", conditionMessage(e))) NULL else stop(e)
            }
        )
        if (is.null(response)) {
            return(empty_ensembl_tibble(.genomes_attributes))
        }
        parsed <- jsonlite::fromJSON(response, flatten = TRUE)
        .records_to_tibble(parsed$data$genomes, .genomes_attributes)
    })

    out <- do.call(rbind, results)
    tibble::as_tibble(out)
}

#' Resolve the genome UUID for an assembly
#'
#' Returns the `genome_id` UUID that [get_ensembl_data()] needs, resolving it
#' once per session and caching the result.
#'
#' @param assembly_accession_id Assembly accession to look up, defaulting to
#'   Homo sapiens GRCh38.p14 (`"GCA_000001405.29"`).
#' @param endpoint `NULL` to use the resolved default endpoint, or a URL.
#' @param client `NULL` to use the default client, or a `ghql::GraphqlClient`.
#' @param timeout,retries,verbose Passed through to the transport; see
#'   [get_ensembl_data()].
#' @param refresh Logical; if `TRUE`, ignore the cached value and look the
#'   assembly up again.
#'
#' @return A single string: the genome UUID.
#'
#' @examples
#' \dontrun{
#' ensembl_genome_id()
#' get_ensembl_data(
#'     attributes = c("stable_id", "symbol"),
#'     filters = "hgnc_symbol",
#'     values = "BRCA2",
#'     genome_id = ensembl_genome_id()
#' )
#' }
#'
#' @seealso [ensembl_genomes()]
#' @export
ensembl_genome_id <- function(assembly_accession_id = .ensembl_default_assembly,
                              endpoint = NULL, client = NULL, timeout = NULL,
                              retries = NULL, refresh = FALSE, verbose = FALSE) {
    if (!is.character(assembly_accession_id) || length(assembly_accession_id) != 1L ||
        is.na(assembly_accession_id) || !nzchar(trimws(assembly_accession_id))) {
        .ensembl_abort("`assembly_accession_id` must be a single non-empty string.")
    }
    assembly_accession_id <- trimws(assembly_accession_id)

    if (is.null(.ensembl_cache$genome_ids)) {
        .ensembl_cache$genome_ids <- new.env(parent = emptyenv())
    }
    cache <- .ensembl_cache$genome_ids
    if (!isTRUE(refresh) && exists(assembly_accession_id, envir = cache, inherits = FALSE)) {
        return(cache[[assembly_accession_id]])
    }

    genomes <- ensembl_genomes(
        keyword = assembly_accession_id, by = "assembly_accession_id",
        endpoint = endpoint, client = client, timeout = timeout,
        retries = retries, verbose = verbose
    )

    usable <- genomes[!is.na(genomes$genome_id), , drop = FALSE]
    if (nrow(usable) == 0L) {
        .ensembl_abort(
            sprintf(
                "No Ensembl genome found for assembly accession %s.",
                encodeString(assembly_accession_id, quote = "\"")
            ),
            class = "ensembl_genome_error"
        )
    }

    # Several releases of one assembly can be present; the newest wins.
    genome_id <- usable$genome_id[[which.max(usable$release_number)]]
    cache[[assembly_accession_id]] <- genome_id
    genome_id
}

# Resolution order for the genome a call is scoped to: explicit argument, then
# the `ensemblGraphQLr.genome_id` option, then the session cache.
.resolve_genome_id <- function(genome_id = NULL, client = NULL, timeout = NULL,
                               retries = NULL, verbose = FALSE) {
    if (!is.null(genome_id)) {
        if (!is.character(genome_id) || length(genome_id) != 1L || is.na(genome_id) ||
            !nzchar(trimws(genome_id))) {
            .ensembl_abort("`genome_id` must be a single non-empty string.")
        }
        return(trimws(genome_id))
    }

    from_option <- getOption("ensemblGraphQLr.genome_id")
    if (!is.null(from_option)) {
        if (!is.character(from_option) || length(from_option) != 1L || is.na(from_option) ||
            !nzchar(trimws(from_option))) {
            .ensembl_abort("The `ensemblGraphQLr.genome_id` option must be a single non-empty string.")
        }
        return(trimws(from_option))
    }

    ensembl_genome_id(
        client = client, timeout = timeout, retries = retries, verbose = verbose
    )
}

# ---- Condition helpers ----------------------------------------------------

# All conditions signalled by this package inherit from
# `ensemblGraphQLr_error`, so `tryCatch(..., ensemblGraphQLr_error = ...)`
# catches every package-specific failure.
.ensembl_abort <- function(message, class = "ensembl_input_error", ...) {
    rlang::abort(
        message,
        class = c(class, "ensemblGraphQLr_error"),
        ...
    )
}
