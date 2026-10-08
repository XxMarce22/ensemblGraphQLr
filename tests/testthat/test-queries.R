# ---------------------------------------------------------------------------
# Tests for query construction, response shaping and request error handling.
#
# Nothing here touches the network unless ENSEMBL_GRAPHQL_LIVE is set, so the
# suite runs deterministically offline.
# ---------------------------------------------------------------------------

GENOME <- "59871324-7803-4234-856e-2a2bd96d7b3c"

# A fake ghql-compatible client: `.ensembl_request()` only ever reads `$url`
# and calls `$exec()`, so a plain list is enough to exercise the transport.
fake_client <- function(handler, url = "https://fake.test/graphql") {
    list(url = url, exec = handler)
}

# Pins the settings every mocked test depends on: a fixed genome UUID so no
# genome lookup is attempted, and no back-off waits. `.local_envir` is required
# so the options survive until the calling test finishes rather than being
# restored as soon as this helper returns.
local_ensembl_settings <- function(genome_id = GENOME, .local_envir = parent.frame()) {
    withr::local_options(
        ensemblGraphQLr.genome_id = genome_id,
        ensemblGraphQLr.backoff = 0,
        .local_envir = .local_envir
    )
}

# ---- Selection sets -------------------------------------------------------

test_that("build_selection_set maps attribute vectors onto GraphQL fields", {
    expect_equal(
        build_selection_set(c("stable_id", "symbol", "name")),
        "stable_id\nsymbol\nname"
    )
})

test_that("build_selection_set preserves the requested order", {
    expect_equal(build_selection_set(c("symbol", "stable_id")), "symbol\nstable_id")
})

test_that("build_selection_set nests dotted attributes", {
    expect_equal(
        build_selection_set(c("stable_id", "metadata.biotype.value")),
        paste(
            "stable_id",
            "metadata {",
            "  biotype {",
            "    value",
            "  }",
            "}",
            sep = "\n"
        )
    )
})

test_that("build_selection_set collapses repeated prefixes", {
    selection <- build_selection_set(c("slice.location.start", "slice.location.end"))
    expect_equal(regmatches(selection, gregexpr("slice", selection))[[1]], "slice")
    expect_equal(regmatches(selection, gregexpr("location", selection))[[1]], "location")
})

test_that("build_selection_set honours indent", {
    expect_equal(build_selection_set("stable_id", indent = 2), "  stable_id")
})

# ---- Full query generation ------------------------------------------------

test_that("build_graphql_query generates a genome-scoped document", {
    query <- build_graphql_query(
        attributes = c("stable_id", "symbol"),
        filters = "ensembl_gene_id",
        values = "ENSG00000139618",
        genome_id = GENOME
    )
    expect_equal(
        query,
        paste(
            "query EnsemblQuery(",
            "  $genomeId: String!,",
            "  $value1: String!",
            ") {",
            "  gene(by_id: {genome_id: $genomeId, stable_id: $value1}) {",
            "    stable_id",
            "    symbol",
            "  }",
            "}",
            sep = "\n"
        )
    )
})

test_that("build_graphql_query aliases one block per value", {
    query <- build_graphql_query(
        attributes = "symbol",
        filters = "hgnc_symbol",
        values = c("BRCA2", "TP53"),
        genome_id = GENOME
    )
    expect_match(query, "$value1: String!", fixed = TRUE)
    expect_match(query, "$value2: String!", fixed = TRUE)
    expect_match(query, "result1: genes(by_symbol: {genome_id: $genomeId, symbol: $value1})", fixed = TRUE)
    expect_match(query, "result2: genes(by_symbol: {genome_id: $genomeId, symbol: $value2})", fixed = TRUE)
})

test_that("build_graphql_query documents that a genome is required", {
    expect_error(
        build_graphql_query("stable_id", "ensembl_gene_id", "X"),
        class = "ensembl_input_error"
    )
})

test_that("build_graphql_query rejects unsupported requests", {
    expect_error(
        build_graphql_query("stable_id", "nope", "x", genome_id = GENOME),
        class = "ensembl_schema_error"
    )
    expect_error(
        build_graphql_query("stable_id", NULL, NULL, genome_id = GENOME),
        class = "ensembl_input_error"
    )
    expect_error(
        build_graphql_query("stable_id", c("ensembl_gene_id", "hgnc_symbol"), "x", genome_id = GENOME),
        class = "ensembl_input_error"
    )
    expect_error(
        build_graphql_query("stable_id", "hgnc_symbol", "x", genome_id = GENOME, root = "bad-root"),
        class = "ensembl_input_error"
    )
})

test_that("variable values are never interpolated into the document", {
    query <- build_graphql_query(
        attributes = "symbol",
        filters = "hgnc_symbol",
        values = 'BRCA2") { evil } #',
        genome_id = GENOME
    )
    expect_false(grepl("evil", query, fixed = TRUE))
    expect_match(query, "symbol: $value1", fixed = TRUE)
})

# ---- Wire format ----------------------------------------------------------

test_that("variable values are serialised as plain JSON strings", {
    plan <- .plan_query(c("stable_id"), "ensembl_gene_id", "ENSG00000139618", genome_id = GENOME)

    # Mirrors the encoding crul applies for `encode = "json"`.
    body <- jsonlite::toJSON(
        list(query = plan$query, variables = plan$variables),
        auto_unbox = TRUE
    )
    expect_match(as.character(body), sprintf('"genomeId":"%s"', GENOME), fixed = TRUE)
    expect_match(as.character(body), '"value1":"ENSG00000139618"', fixed = TRUE)
    expect_type(plan$variables$value1, "character")
})

test_that("several values produce one variable each", {
    plan <- .plan_query(c("symbol"), "hgnc_symbol", c("BRCA2", "TP53"), genome_id = GENOME)
    expect_equal(names(plan$variables), c("genomeId", "value1", "value2"))
    expect_equal(plan$root_fields, c("result1", "result2"))
})

test_that("a single value keeps the plain root field name", {
    plan <- .plan_query(c("symbol"), "hgnc_symbol", "BRCA2", genome_id = GENOME)
    expect_equal(plan$root_fields, "genes")
})

# ---- Input validation -----------------------------------------------------

test_that("get_ensembl_data validates its inputs", {
    expect_error(get_ensembl_data(character(0)), class = "ensembl_input_error")
    expect_error(get_ensembl_data(NULL), class = "ensembl_input_error")
    expect_error(get_ensembl_data(c("stable_id", NA)), class = "ensembl_input_error")
    expect_error(get_ensembl_data("bad-attribute"), class = "ensembl_input_error")
    expect_error(get_ensembl_data("__schema"), class = "ensembl_input_error")

    # `values` presuppose `filters`, and vice versa.
    local_ensembl_settings()
    expect_error(
        get_ensembl_data("stable_id", values = "BRCA2"),
        class = "ensembl_input_error"
    )
    expect_error(
        get_ensembl_data("stable_id", filters = "hgnc_symbol"),
        class = "ensembl_input_error"
    )

    # The core schema accepts one identifier field per lookup.
    expect_error(
        get_ensembl_data("stable_id", c("ensembl_gene_id", "hgnc_symbol"), c("A", "B")),
        class = "ensembl_input_error"
    )
})

test_that("duplicate attributes are collapsed", {
    expect_equal(.check_attributes(c("stable_id", "stable_id", "symbol")), c("stable_id", "symbol"))
})

test_that("settings are validated", {
    expect_error(.check_timeout(-1), class = "ensembl_input_error")
    expect_error(.check_timeout("soon"), class = "ensembl_input_error")
    expect_error(.check_retries(-1), class = "ensembl_input_error")
    expect_error(ensembl_graphql_endpoint("not-a-url"), class = "ensembl_input_error")
    expect_error(ensembl_graphql_endpoint(c("a", "b")), class = "ensembl_input_error")
})

# ---- Empty results --------------------------------------------------------

test_that("empty_ensembl_tibble matches the requested attributes", {
    out <- empty_ensembl_tibble(c("stable_id", "symbol", "slice.location.start"))
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 0L)
    expect_equal(names(out), c("stable_id", "symbol", "slice_location_start"))
    expect_type(out$stable_id, "character")
    expect_type(out$slice_location_start, "integer")
})

test_that("nested attribute names become underscore column names", {
    out <- empty_ensembl_tibble(c("stable_id", "metadata.biotype.value"))
    expect_equal(names(out), c("stable_id", "metadata_biotype_value"))
})

# ---- Response shaping -----------------------------------------------------

test_that(".to_record_rows normalises every payload shape", {
    expect_equal(.to_record_rows(NULL), list())
    expect_equal(.to_record_rows(list()), list())
    expect_equal(length(.to_record_rows(list(a = 1))), 1L)
    expect_equal(length(.to_record_rows(list(list(a = 1), list(a = 2)))), 2L)

    flat <- jsonlite::fromJSON('{"genes":[{"symbol":"A"},{"symbol":"B"}]}', flatten = TRUE)
    expect_equal(length(.to_record_rows(flat$genes)), 2L)
    expect_equal(length(.to_record_rows(data.frame(symbol = character()))), 0L)
})

test_that("get_ensembl_data returns a flattened tibble for one gene", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(...) {
            '{"data":{"gene":{"stable_id":"ENSG00000139618.19","symbol":"BRCA2","metadata":{"biotype":{"value":"protein_coding"}}}}}'
        }
    )
    out <- get_ensembl_data(
        attributes = c("stable_id", "symbol", "metadata.biotype.value"),
        filters = "ensembl_gene_id",
        values = "ENSG00000139618"
    )
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 1L)
    expect_equal(names(out), c("stable_id", "symbol", "metadata_biotype_value"))
    expect_equal(out$symbol, "BRCA2")
    expect_equal(out$metadata_biotype_value, "protein_coding")
})

test_that("get_ensembl_data concatenates the aliased root fields", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(...) {
            paste0(
                '{"data":{',
                '"result1":[{"stable_id":"A","symbol":"BRCA2"}],',
                '"result2":[{"stable_id":"B","symbol":"TP53"}]',
                '}}'
            )
        }
    )
    out <- get_ensembl_data(c("stable_id", "symbol"), "hgnc_symbol", c("BRCA2", "TP53"))
    expect_equal(nrow(out), 2L)
    expect_equal(out$symbol, c("BRCA2", "TP53"))
})

test_that("get_ensembl_data forwards the planned variables to the transport", {
    local_ensembl_settings()
    seen <- NULL
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(query, variables, ...) {
            seen <<- variables
            '{"data":{"result1":[],"result2":[]}}'
        }
    )
    get_ensembl_data(c("symbol"), "hgnc_symbol", c("BRCA2", "TP53"))
    expect_equal(seen$genomeId, GENOME)
    expect_equal(seen$value1, "BRCA2")
    expect_equal(seen$value2, "TP53")
})

test_that("get_ensembl_data returns an empty tibble when nothing matches", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(...) '{"data":{"genes":[]}}'
    )
    out <- get_ensembl_data(c("stable_id", "symbol"), "hgnc_symbol", "NOT_A_GENE")
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 0L)
    expect_equal(names(out), c("stable_id", "symbol"))
})

test_that("get_ensembl_data warns when an attribute is never returned", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(...) '{"data":{"genes":[{"stable_id":"A"}]}}'
    )
    expect_warning(
        get_ensembl_data(c("stable_id", "not_a_field"), "hgnc_symbol", "X"),
        class = "ensembl_missing_attributes"
    )
})

test_that("get_ensembl_data warns when the root field is absent", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        .ensembl_request = function(...) '{"data":{"somethingElse":{}}}'
    )
    expect_warning(
        get_ensembl_data(c("stable_id"), "hgnc_symbol", "X"),
        class = "ensembl_missing_field"
    )
})

# ---- Not-found handling ---------------------------------------------------

# These go through the real transport (a fake client) so that the
# `*_NOT_FOUND` downgrade is exercised rather than mocked away.
test_that("a not-found lookup yields an empty tibble, not an error", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"data":{"genes":null},"errors":[{"message":"Failed to find gene with ids: symbol=NOPE",',
            '"path":["genes"],"extensions":{"code":"GENE_NOT_FOUND"}}]}'
        )
    })
    out <- get_ensembl_data(
        attributes = c("stable_id", "symbol"),
        filters = "hgnc_symbol", values = "NOPE", client = client
    )
    expect_equal(nrow(out), 0L)
    expect_equal(names(out), c("stable_id", "symbol"))
})

test_that("partial matches are kept and missing aliases are dropped", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"data":{',
            '"result1":[{"stable_id":"A","symbol":"BRCA2"}],',
            '"result2":null',
            '},"errors":[{"message":"Failed to find gene with ids: symbol=NOPE",',
            '"path":["result2"],"extensions":{"code":"GENE_NOT_FOUND"}}]}'
        )
    })
    out <- get_ensembl_data(
        attributes = c("stable_id", "symbol"),
        filters = "hgnc_symbol", values = c("BRCA2", "NOPE"), client = client
    )
    expect_equal(nrow(out), 1L)
    expect_equal(out$symbol, "BRCA2")
})

test_that("a not-found alias does not trigger a missing-field warning", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"data":{"genes":null},"errors":[{"message":"Failed to find gene",',
            '"path":["genes"],"extensions":{"code":"GENE_NOT_FOUND"}}]}'
        )
    })
    expect_no_warning(
        get_ensembl_data(c("symbol"), "hgnc_symbol", "NOPE", client = client)
    )
})

test_that(".ensembl_not_found only accepts NOT_FOUND codes", {
    expect_true(.ensembl_not_found(list(
        list(message = "x", extensions = list(code = "GENE_NOT_FOUND"))
    )))
    expect_false(.ensembl_not_found(list(
        list(message = "x", extensions = list(code = "GRAPHQL_VALIDATION_FAILED"))
    )))
    expect_false(.ensembl_not_found(list()))
})

# ---- Transport failures ---------------------------------------------------

test_that("transient failures are retried and can still succeed", {
    local_ensembl_settings()
    attempts <- 0L
    client <- fake_client(function(...) {
        attempts <<- attempts + 1L
        if (attempts < 2L) {
            stop("Timeout was reached [fake.test]: Operation timed out")
        }
        '{"data":{"genes":[]}}'
    })

    out <- .ensembl_request("query Q { genes { stable_id } }", list(), client = client, retries = 3L)
    expect_equal(attempts, 2L)
    expect_equal(out, '{"data":{"genes":[]}}')
})

test_that("exhausted connection failures raise a classed error", {
    local_ensembl_settings()
    client <- fake_client(function(...) stop("Could not resolve host: fake.test"))

    err <- tryCatch(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = client, retries = 1L),
        error = function(e) e
    )
    expect_s3_class(err, "ensembl_connection_error")
    expect_equal(err$attempts, 2L)
})

test_that("retryable HTTP statuses are retried but client errors are not", {
    local_ensembl_settings()

    server_calls <- 0L
    server <- fake_client(function(...) {
        server_calls <<- server_calls + 1L
        stop("[HTTP 503] Service Unavailable")
    })
    err <- tryCatch(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = server, retries = 2L),
        error = function(e) e
    )
    expect_s3_class(err, "ensembl_http_error")
    expect_equal(err$status, 503L)
    expect_equal(server_calls, 3L)

    client_calls <- 0L
    bad_request <- fake_client(function(...) {
        client_calls <<- client_calls + 1L
        stop("[HTTP 400] Bad Request")
    })
    expect_error(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = bad_request, retries = 2L),
        class = "ensembl_http_error"
    )
    expect_equal(client_calls, 1L)
})

test_that("server-side GraphQL errors are surfaced with their messages", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"errors":[{"message":"Cannot query field \\"biotype\\" on type \\"Gene\\".",',
            '"extensions":{"code":"GRAPHQL_VALIDATION_FAILED"}}]}'
        )
    })
    err <- tryCatch(
        .ensembl_request("query Q { gene { biotype } }", list(), client = client, retries = 0L),
        error = function(e) e
    )
    expect_s3_class(err, "ensembl_graphql_error")
    expect_match(err$errors, "Cannot query field")
    expect_match(err$errors, "GRAPHQL_VALIDATION_FAILED")
})

test_that("validation failures arriving as HTTP 400 report the GraphQL messages", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        stop(paste0(
            '[HTTP 400] {"errors":[{"message":"Field \'parlance_name\' is not defined",',
            '"extensions":{"code":"GRAPHQL_VALIDATION_FAILED"}}]}'
        ))
    })
    err <- tryCatch(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = client, retries = 2L),
        error = function(e) e
    )
    # Reported as a query problem rather than a bare 400, and not retried.
    expect_s3_class(err, "ensembl_graphql_error")
    expect_match(err$errors, "parlance_name")
})

test_that("a structural attribute without a sub-field gets a targeted hint", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"errors":[{"message":"Field \'transcripts\' of type \'[Transcript!]!\'',
            ' must have a selection of subfields."}]}'
        )
    })
    err <- tryCatch(
        .ensembl_request("query Q { gene { transcripts } }", list(), client = client, retries = 0L),
        error = function(e) e
    )
    expect_s3_class(err, "ensembl_graphql_error")
    expect_match(conditionMessage(err), "transcripts.stable_id", fixed = TRUE)
})

test_that("empty and non-JSON bodies are reported distinctly", {
    local_ensembl_settings()

    empty <- fake_client(function(...) "   ")
    expect_error(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = empty),
        class = "ensembl_empty_response"
    )

    proxy <- fake_client(function(...) "<html><title>403 Forbidden</title></html>")
    expect_error(
        .ensembl_request("query Q { genes { stable_id } }", list(), client = proxy),
        class = "ensembl_response_error"
    )
})

# ---- Genome discovery -----------------------------------------------------

test_that("ensembl_genomes validates its keyword field", {
    expect_error(ensembl_genomes("x", by = "parlance_name"), class = "ensembl_input_error")
    expect_error(ensembl_genomes(character(0)), class = "ensembl_input_error")
})

test_that("ensembl_genomes reshapes the payload", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"data":{"genomes":[{"genome_id":"abc","genome_tag":"GRCh38",',
            '"release_number":2,"assembly":{"accession_id":"GCA_000001405.29","name":"GRCh38.p14"}}]}}'
        )
    })
    out <- ensembl_genomes("GCA_000001405.29", client = client)
    expect_equal(nrow(out), 1L)
    expect_equal(out$genome_id, "abc")
    expect_equal(out$assembly_name, "GRCh38.p14")
})

test_that("an unknown assembly yields an empty tibble", {
    local_ensembl_settings()
    client <- fake_client(function(...) {
        paste0(
            '{"data":{"genomes":null},"errors":[{"message":"Failed to find genome with ids: tolid=x",',
            '"path":["genomes"],"extensions":{"code":"GENOME_NOT_FOUND"}}]}'
        )
    })
    out <- ensembl_genomes("x", by = "tolid", client = client)
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 0L)
    expect_equal(
        names(out),
        c("genome_id", "genome_tag", "release_number", "assembly_accession_id", "assembly_name")
    )
})

test_that("ensembl_genome_id caches the resolved UUID", {
    local_ensembl_settings()
    calls <- 0L
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        ensembl_genomes = function(...) {
            calls <<- calls + 1L
            tibble::tibble(
                genome_id = "resolved-uuid",
                genome_tag = "GRCh38",
                release_number = 2,
                assembly_accession_id = "GCA_000001405.29",
                assembly_name = "GRCh38.p14"
            )
        }
    )
    .ensembl_cache$genome_ids <- new.env(parent = emptyenv())

    expect_equal(ensembl_genome_id(), "resolved-uuid")
    expect_equal(ensembl_genome_id(), "resolved-uuid")
    expect_equal(calls, 1L)

    expect_equal(ensembl_genome_id(refresh = TRUE), "resolved-uuid")
    expect_equal(calls, 2L)
})

test_that("ensembl_genome_id prefers the newest release", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        ensembl_genomes = function(...) {
            tibble::tibble(
                genome_id = c("old", "new"),
                genome_tag = c("GRCh38", "GRCh38"),
                release_number = c(1, 3),
                assembly_accession_id = c("a", "a"),
                assembly_name = c("GRCh38.p14", "GRCh38.p14")
            )
        }
    )
    .ensembl_cache$genome_ids <- new.env(parent = emptyenv())
    expect_equal(ensembl_genome_id(), "new")
})

test_that("an unresolvable assembly raises a classed error", {
    local_ensembl_settings()
    local_mocked_bindings(
        .package = "ensemblGraphQLr",
        ensembl_genomes = function(...) empty_ensembl_tibble(.genomes_attributes)
    )
    .ensembl_cache$genome_ids <- new.env(parent = emptyenv())
    expect_error(ensembl_genome_id(), class = "ensembl_genome_error")
})

test_that("the genome can be configured by argument or option", {
    local_ensembl_settings(genome_id = "from-option")
    expect_equal(.resolve_genome_id(), "from-option")
    expect_equal(.resolve_genome_id("from-argument"), "from-argument")
    expect_error(.resolve_genome_id(42), class = "ensembl_input_error")
})

# ---- Client configuration -------------------------------------------------

test_that("the endpoint can be configured by option", {
    withr::local_options(ensemblGraphQLr.endpoint = "https://example.org/graphql")
    expect_equal(ensembl_graphql_endpoint(), "https://example.org/graphql")
    expect_equal(
        ensembl_graphql_endpoint("https://other.org/graphql"),
        "https://other.org/graphql"
    )
})

test_that("the default endpoint targets the Ensembl core GraphQL service", {
    withr::local_options(ensemblGraphQLr.endpoint = NULL)
    withr::local_envvar(ENSEMBL_GRAPHQL_ENDPOINT = NA)
    expect_equal(
        ensembl_graphql_endpoint(),
        "https://www.ensembl.org/api/graphql/core"
    )
})

test_that("connect_ensembl returns a ghql client and caches it", {
    client <- connect_ensembl("https://example.org/graphql")
    expect_s3_class(client, "GraphqlClient")
    expect_equal(client$url, "https://example.org/graphql")
    expect_identical(client, connect_ensembl("https://example.org/graphql"))
    # Custom headers are never shared through the cache.
    other <- connect_ensembl("https://example.org/graphql", headers = list(Authorization = "x"))
    expect_false(identical(client, other))
    expect_error(connect_ensembl(headers = list("no-name")), class = "ensembl_input_error")
})

# ---- Live integration -----------------------------------------------------

test_that("a live query for BRCA2 returns a populated tibble", {
    skip_on_cran()
    skip_if(
        Sys.getenv("ENSEMBL_GRAPHQL_LIVE") == "",
        "Set ENSEMBL_GRAPHQL_LIVE=1 to run live Ensembl queries."
    )

    out <- get_ensembl_data(
        attributes = c("stable_id", "symbol", "metadata.biotype.value"),
        filters = "hgnc_symbol",
        values = "BRCA2",
        timeout = 30
    )

    expect_s3_class(out, "tbl_df")
    expect_true("symbol" %in% names(out))
    expect_true(nrow(out) >= 1L)
    expect_equal(out$symbol[[1L]], "BRCA2")
    expect_equal(out$metadata_biotype_value[[1L]], "protein_coding")

    # A live miss is a zero-row tibble, not an error.
    expect_equal(
        nrow(get_ensembl_data(
            attributes = c("stable_id", "symbol"),
            filters = "hgnc_symbol",
            values = "NOT_A_REAL_GENE_XYZ",
            timeout = 30
        )),
        0L
    )
})
