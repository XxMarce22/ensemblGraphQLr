# ---------------------------------------------------------------------------
# query_builder.R -- translating biomaRt-style requests into GraphQL
#
# This file owns the *schema adapter* between the familiar `biomaRt` vocabulary
# (attributes such as `stable_id`, filters such as `ensembl_gene_id`) and a
# GraphQL document. It also owns the conversion of GraphQL payloads back into
# tibbles, including the empty-result case.
# ---------------------------------------------------------------------------

# ---- Schema adapter -------------------------------------------------------

# The single point of adaptation between biomaRt filter names and the Ensembl
# GraphQL core schema. If Ensembl renames an argument between releases, only
# this table needs to change.
#
# Verified against the Ensembl core API (v0.2.0-beta, release 116) on
# 2026-10-08. Lookups are always scoped to one genome, so the input object
# always carries `genome_id` alongside the identifier field.
.ensembl_filter_registry <- list(
    ensembl_gene_id = list(
        root = "gene",
        argument = "by_id",
        value_field = "stable_id",
        returns = "single"
    ),
    hgnc_symbol = list(
        root = "genes",
        argument = "by_symbol",
        value_field = "symbol",
        returns = "list"
    ),
    ensembl_transcript_id = list(
        root = "transcript",
        argument = "by_id",
        value_field = "stable_id",
        returns = "single"
    )
)

# Attributes whose values are whole numbers. Used to give empty result sets a
# realistic, typed structure instead of an all-character shell. The core schema
# nests coordinates under `slice`, so both the nested and the flattened
# spellings are listed.
.ensembl_attribute_types <- c(
    start = "integer",
    end = "integer",
    strand = "integer",
    length = "integer",
    version = "integer",
    release_number = "integer",
    seq_region_start = "integer",
    seq_region_end = "integer",
    seq_region_strand = "integer",
    transcript_count = "integer",
    gene_count = "integer",
    "slice.location.start" = "integer",
    "slice.location.end" = "integer",
    "slice.location.length" = "integer",
    "slice.strand.code" = "integer"
)

#' Filters supported by the Ensembl core GraphQL API
#'
#' Reports the biomaRt-style filter names that [get_ensembl_data()] knows how to
#' translate, together with the GraphQL root field and argument each one
#' targets and whether the lookup returns a single record or a list.
#'
#' The Ensembl core schema accepts exactly one identifier field per lookup, so
#' exactly one filter may be supplied per call.
#'
#' @return A [tibble][tibble::tibble] with the columns `filter`, `root`,
#'   `argument`, `value_field` and `returns`.
#'
#' @examples
#' ensembl_filters()
#'
#' @seealso [build_graphql_query()]
#' @export
ensembl_filters <- function() {
    tibble::tibble(
        filter = names(.ensembl_filter_registry),
        root = vapply(.ensembl_filter_registry, function(x) x$root, character(1)),
        argument = vapply(.ensembl_filter_registry, function(x) x$argument, character(1)),
        value_field = vapply(.ensembl_filter_registry, function(x) x$value_field, character(1)),
        returns = vapply(.ensembl_filter_registry, function(x) x$returns, character(1))
    )
}

# ---- Input validation -----------------------------------------------------

.graphql_name_pattern <- "^[A-Za-z_][A-Za-z0-9_]*$"
.attribute_path_pattern <- "^[A-Za-z_][A-Za-z0-9_]*(\\.[A-Za-z_][A-Za-z0-9_]*)*$"

.check_attributes <- function(attributes) {
    if (is.null(attributes) || !is.character(attributes) || length(attributes) == 0L) {
        .ensembl_abort("`attributes` must be a character vector with at least one attribute.")
    }
    if (anyNA(attributes)) {
        .ensembl_abort("`attributes` must not contain `NA`.")
    }

    attributes <- trimws(attributes)
    blank <- !nzchar(attributes)
    if (any(blank)) {
        .ensembl_abort("`attributes` must not contain empty or whitespace-only entries.")
    }

    invalid <- !grepl(.attribute_path_pattern, attributes)
    if (any(invalid)) {
        .ensembl_abort(sprintf(
            "Not a valid GraphQL attribute path: %s. Attributes use letters, digits and underscores, optionally nested with dots (e.g. \"transcript.stable_id\").",
            paste(encodeString(attributes[invalid], quote = "\""), collapse = ", ")
        ))
    }

    reserved <- grepl("(^|\\.)__", attributes)
    if (any(reserved)) {
        .ensembl_abort(sprintf(
            "Attributes starting with \"__\" are reserved for GraphQL introspection: %s.",
            paste(encodeString(attributes[reserved], quote = "\""), collapse = ", ")
        ))
    }

    unique(attributes)
}

.check_values <- function(values) {
    if (is.null(values) || length(values) == 0L) {
        .ensembl_abort("`values` must contain at least one value.")
    }
    if (is.list(values)) {
        values <- lapply(values, function(v) {
            if (is.null(v) || length(v) == 0L) {
                .ensembl_abort("Every element of a list `values` must contain at least one value.")
            }
            as.character(v)
        })
        if (any(vapply(values, function(v) anyNA(v), logical(1)))) {
            .ensembl_abort("`values` must not contain `NA`.")
        }
        return(values)
    }
    values <- as.character(values)
    if (anyNA(values)) {
        .ensembl_abort("`values` must not contain `NA`.")
    }
    if (length(values) == 0L) {
        .ensembl_abort("`values` must contain at least one value.")
    }
    values
}

# ---- Filter resolution ----------------------------------------------------

# Turns `filters`/`values` into the GraphQL variable definitions, the aliased
# root-field blocks and the variable values to send. Data always travels as
# GraphQL *variables*, never interpolated into the document, so no user input
# can alter the query shape.
#
# The core schema accepts exactly one identifier field per lookup, so a single
# filter is required; several values are handled by emitting one aliased block
# per value and concatenating the results.
.resolve_filters <- function(filters, values) {
    if (is.null(filters)) {
        if (!is.null(values)) {
            .ensembl_abort("`values` was supplied without `filters`; both must be given together.")
        }
        .ensembl_abort(
            c(
                "`filters` is required.",
                "i" = "The Ensembl core GraphQL schema has no unfiltered gene collection; supply `filters` and `values` to scope a lookup, for example `filters = \"ensembl_gene_id\"`."
            ),
            class = "ensembl_input_error"
        )
    }

    if (!is.character(filters) || length(filters) == 0L || anyNA(filters)) {
        .ensembl_abort("`filters` must be a non-missing character vector.")
    }
    filters <- trimws(filters)
    if (any(!nzchar(filters))) {
        .ensembl_abort("`filters` must not contain empty or whitespace-only entries.")
    }

    unknown <- setdiff(filters, names(.ensembl_filter_registry))
    if (length(unknown) > 0L) {
        .ensembl_abort(
            sprintf(
                "Unknown filter(s): %s.\nSupported filters: %s.",
                paste(encodeString(unknown, quote = "\""), collapse = ", "),
                paste(sort(names(.ensembl_filter_registry)), collapse = ", ")
            ),
            class = "ensembl_schema_error",
            filters = unknown
        )
    }

    if (length(filters) > 1L) {
        .ensembl_abort(
            c(
                sprintf(
                    "Only one filter can be used per call, but %d were supplied.",
                    length(filters)
                ),
                "i" = "The Ensembl core schema accepts a single identifier field per lookup. Call `get_ensembl_data()` once per filter, or pass several values for one filter."
            ),
            class = "ensembl_input_error",
            filters = filters
        )
    }

    if (is.null(values)) {
        .ensembl_abort("`values` must be supplied when `filters` is not `NULL`.")
    }
    if (is.list(values)) {
        if (length(values) != 1L) {
            .ensembl_abort(sprintf(
                "A list `values` must have exactly one element, one per filter, but has %d.",
                length(values)
            ))
        }
        values <- values[[1L]]
    }
    values <- .check_values(values)

    list(
        filter = filters,
        spec = .ensembl_filter_registry[[filters]],
        values = values
    )
}

# ---- Selection sets -------------------------------------------------------

# Attributes are nested paths such as "transcript.stable_id". Expanding them
# into a named list tree first means duplicated prefixes collapse naturally and
# the requested order is preserved, because R lists keep insertion order.
.attribute_tree <- function(attributes) {
    tree <- list()
    for (attribute in attributes) {
        parts <- strsplit(attribute, ".", fixed = TRUE)[[1]]
        tree <- .insert_attribute(tree, parts)
    }
    tree
}

.insert_attribute <- function(node, parts) {
    if (length(parts) == 0L) {
        return(node)
    }
    key <- parts[[1L]]
    child <- node[[key]] %||% list()
    if (length(parts) > 1L) {
        child <- .insert_attribute(child, parts[-1L])
    }
    node[[key]] <- child
    node
}

.render_selection <- function(node, base = 0L, level = 0L) {
    if (length(node) == 0L) {
        return(character())
    }
    padding <- strrep(" ", base + 2L * level)
    lines <- character()
    for (name in names(node)) {
        child <- node[[name]]
        if (length(child) == 0L) {
            lines <- c(lines, paste0(padding, name))
        } else {
            lines <- c(lines, sprintf(
                "%s%s {\n%s\n%s}",
                padding, name,
                paste(.render_selection(child, base, level + 1L), collapse = "\n"),
                padding
            ))
        }
    }
    lines
}

#' Build a GraphQL selection set from attribute names
#'
#' Converts a character vector of attribute names into the selection set of a
#' GraphQL query. Dotted names become nested selection blocks, so
#' `"transcript.stable_id"` is emitted as `transcript { stable_id }`.
#'
#' @param attributes A character vector of attribute names. Names may be nested
#'   using dots, for example `"transcript.stable_id"`.
#' @param indent Number of spaces to indent the block by. Defaults to `0`.
#'
#' @return A single string holding the selection set body, with one field per
#'   line.
#'
#' @examples
#' cat(build_selection_set(c("stable_id", "symbol", "biotype")))
#' cat(build_selection_set(c("stable_id", "transcript.stable_id")))
#'
#' @seealso [build_graphql_query()]
#' @export
build_selection_set <- function(attributes, indent = 0L) {
    attributes <- .check_attributes(attributes)
    if (!is.numeric(indent) || length(indent) != 1L || is.na(indent) || indent < 0) {
        .ensembl_abort("`indent` must be a single non-negative number.")
    }
    paste(.render_selection(.attribute_tree(attributes), as.integer(indent)), collapse = "\n")
}

#' Build a complete GraphQL query for an Ensembl request
#'
#' Assembles the operation, its variable definitions, the aliased root fields
#' with their arguments and the selection set into a single GraphQL document.
#' Filter values are *never* interpolated into the document; they are always
#' passed as GraphQL variables, which makes query-string injection impossible.
#'
#' The core schema accepts a single identifier per lookup, so each requested
#' value becomes its own aliased root field (`result1`, `result2`, ...) sharing
#' one `$genomeId` variable. `get_ensembl_data()` concatenates the per-alias
#' results back into one tibble.
#'
#' @param attributes A character vector of attribute names, possibly nested
#'   with dots. See [build_selection_set()].
#' @param filters A single biomaRt-style filter name, such as
#'   `"ensembl_gene_id"` or `"hgnc_symbol"`. See [ensembl_filters()] for the
#'   supported set.
#' @param values A vector of values to match. Each value produces one aliased
#'   lookup.
#' @param genome_id The UUID of the Ensembl genome to query. See
#'   [ensembl_genome_id()]; it is required because every core lookup is scoped
#'   to a genome.
#' @param root `NULL` to take the root field from the filter, or an override
#'   such as `"gene"`.
#' @param operation_name Name of the GraphQL operation, useful when profiling
#'   requests server-side.
#'
#' @return A single string: the GraphQL document, ready to POST.
#'
#' @examples
#' cat(build_graphql_query(
#'     attributes = c("stable_id", "symbol", "metadata.biotype.value"),
#'     filters = "ensembl_gene_id",
#'     values = "ENSG00000139618",
#'     genome_id = "59871324-7803-4234-856e-2a2bd96d7b3c"
#' ))
#'
#' @seealso [get_ensembl_data()], [build_selection_set()], [ensembl_genome_id()]
#' @export
build_graphql_query <- function(attributes, filters = NULL, values = NULL,
                                genome_id = NULL, root = NULL,
                                operation_name = "EnsemblQuery") {
    .plan_query(
        attributes = attributes, filters = filters, values = values,
        genome_id = genome_id, root = root, operation_name = operation_name
    )$query
}

# The full request plan: the query document plus everything needed to execute
# it (`variables`) and to interpret the response (`root_fields`).
.plan_query <- function(attributes, filters = NULL, values = NULL, genome_id = NULL,
                        root = NULL, operation_name = "EnsemblQuery") {
    attributes <- .check_attributes(attributes)
    resolved <- .resolve_filters(filters, values)

    if (!is.character(genome_id) || length(genome_id) != 1L || is.na(genome_id) ||
        !nzchar(trimws(genome_id))) {
        .ensembl_abort(
            c(
                "`genome_id` must be a single non-empty string.",
                "i" = "Every Ensembl core lookup is scoped to a genome. Resolve one with `ensembl_genome_id()`."
            ),
            class = "ensembl_input_error"
        )
    }
    genome_id <- trimws(genome_id)

    spec <- resolved$spec
    if (!is.null(root) &&
        (!is.character(root) || length(root) != 1L || is.na(root) ||
            !grepl(.graphql_name_pattern, root))) {
        .ensembl_abort(
            c(
                "`root` must be a single valid GraphQL field name, or `NULL`.",
                "i" = "The default comes from the filter; see `ensembl_filters()$root`."
            ),
            class = "ensembl_input_error"
        )
    }
    field <- root %||% spec$root
    values <- resolved$values

    # One aliased block per value. A single value needs no alias at all, which
    # keeps the generated query readable.
    aliases <- if (length(values) == 1L) field else paste0("result", seq_along(values))
    value_variables <- paste0("value", seq_along(values))
    selection <- paste(.render_selection(.attribute_tree(attributes), base = 4L), collapse = "\n")

    calls <- vapply(
        seq_along(values),
        function(i) {
            sprintf(
                "%s(%s: {genome_id: $genomeId, %s: $%s})",
                field, spec$argument, spec$value_field, value_variables[[i]]
            )
        },
        character(1)
    )
    if (length(values) > 1L) {
        calls <- paste0(aliases, ": ", calls)
    }
    blocks <- sprintf("%s {\n%s\n  }", calls, selection)

    definitions <- c(
        "$genomeId: String!",
        sprintf("$%s: String!", value_variables)
    )
    variables <- c(
        list(genomeId = genome_id),
        stats::setNames(as.list(values), value_variables)
    )

    query <- .assemble_query(
        blocks = blocks,
        definitions = definitions,
        operation_name = operation_name
    )

    list(
        query = query,
        variables = variables,
        root_fields = aliases,
        root = field,
        filter = resolved$filter,
        attributes = attributes
    )
}

.assemble_query <- function(blocks, definitions, operation_name) {
    if (!is.character(operation_name) || length(operation_name) != 1L ||
        is.na(operation_name) || !grepl(.graphql_name_pattern, operation_name)) {
        .ensembl_abort("`operation_name` must be a valid GraphQL name.")
    }

    if (length(definitions) == 0L) {
        header <- sprintf("query %s", operation_name)
    } else if (length(definitions) == 1L) {
        header <- sprintf("query %s(%s)", operation_name, definitions)
    } else {
        header <- sprintf(
            "query %s(\n  %s\n)",
            operation_name,
            paste(definitions, collapse = ",\n  ")
        )
    }

    sprintf("%s {\n  %s\n}", header, paste(blocks, collapse = "\n  "))
}

# ---- Result shaping -------------------------------------------------------

# Column names follow `biomaRt` conventions, where nesting is expressed with
# underscores ("transcript.stable_id" -> "transcript_stable_id").
.attribute_column_names <- function(attributes) {
    gsub(".", "_", attributes, fixed = TRUE)
}

.empty_column <- function(type) {
    switch(type,
        integer = integer(0),
        double = numeric(0),
        logical = logical(0),
        character(0)
    )
}

# `[[` on a named vector errors for unknown names, so the type table is probed
# with `[` and falls back to character for anything not listed. The table is
# keyed by attribute path, since nesting is what decides the numeric type.
.attribute_type <- function(attribute) {
    type <- unname(.ensembl_attribute_types[attribute])
    if (length(type) != 1L || is.na(type)) "character" else type
}

#' Build an empty result tibble matching a set of attributes
#'
#' Produces a zero-row tibble whose columns correspond exactly to `attributes`,
#' in the requested order. Columns that are known to be numeric are typed
#' accordingly, so callers can bind the result onto non-empty results without a
#' type clash.
#'
#' @param attributes A character vector of attribute names.
#'
#' @return A zero-row [tibble][tibble::tibble] with one column per attribute.
#'
#' @examples
#' empty_ensembl_tibble(c("stable_id", "symbol", "start"))
#'
#' @seealso [get_ensembl_data()]
#' @export
empty_ensembl_tibble <- function(attributes) {
    attributes <- .check_attributes(attributes)
    columns <- .attribute_column_names(attributes)
    types <- vapply(attributes, .attribute_type, character(1), USE.NAMES = FALSE)
    out <- lapply(types, .empty_column)
    names(out) <- columns
    tibble::as_tibble(out, .name_repair = "minimal")
}

# `jsonlite::flatten = TRUE` turns nested *single* objects into dotted keys
# ("transcript.stable_id") but leaves nested *arrays* as data frames, so both
# spellings have to be probed before walking the path manually.
.pluck_attribute <- function(record, path) {
    if (!is.list(record)) {
        return(NULL)
    }
    if (!is.null(record[[path]])) {
        return(record[[path]])
    }
    flat <- .attribute_column_names(path)
    if (!is.null(record[[flat]])) {
        return(record[[flat]])
    }

    current <- record
    for (part in strsplit(path, ".", fixed = TRUE)[[1]]) {
        if (is.data.frame(current)) {
            current <- current[[part]]
        } else if (is.list(current)) {
            current <- current[[part]]
        } else {
            return(NULL)
        }
        if (is.null(current)) {
            return(NULL)
        }
    }
    current
}

# Collapses a column of per-record values to the most specific rectangular type
# it supports, falling back to a list-column when values are genuinely nested
# (for example one gene with several transcripts).
.simplify_column <- function(values) {
    n <- length(values)
    if (n == 0L) {
        return(character())
    }

    missing <- vapply(values, is.null, logical(1))
    if (all(missing)) {
        return(rep(NA_character_, n))
    }

    if (any(missing)) {
        present_types <- unique(vapply(values[!missing], typeof, character(1)))
        if (length(present_types) == 1L) {
            replacement <- switch(present_types,
                character = NA_character_,
                integer = NA_integer_,
                double = NA_real_,
                logical = NA,
                NULL
            )
            if (!is.null(replacement)) {
                values[missing] <- list(replacement)
            }
        }
    }

    scalars <- vapply(
        values,
        function(v) is.atomic(v) && !is.object(v) && length(v) == 1L,
        logical(1)
    )
    if (all(scalars)) {
        types <- unique(vapply(values, typeof, character(1)))
        if (length(types) == 1L) {
            return(unlist(values, use.names = FALSE))
        }
        if (all(types %in% c("character", "integer", "double", "logical"))) {
            return(vapply(values, as.character, character(1)))
        }
    }

    values
}

#' Normalise a GraphQL response fragment into a list of records
#'
#' A single root field can yield one object (`gene`), an array of objects
#' (`genes`), an empty array, or `NULL`. This flattens all of those shapes into
#' a list of records so that several aliased root fields can be concatenated.
#'
#' @param value The value extracted from the response's `data` member.
#'
#' @return A list of records, possibly empty.
#'
#' @noRd
.to_record_rows <- function(value) {
    if (is.null(value)) {
        return(list())
    }

    if (is.data.frame(value)) {
        if (nrow(value) == 0L) {
            return(list())
        }
        # Attribute paths like "metadata.biotype.value" arrive as dotted column
        # names after flattening, so a row-wise split keeps them reachable.
        return(lapply(seq_len(nrow(value)), function(i) {
            as.list(value[i, , drop = FALSE])
        }))
    }

    if (!is.list(value) || length(value) == 0L) {
        return(list())
    }

    # A named list is a single record; an unnamed list is a collection of them.
    if (is.null(names(value))) {
        return(value[vapply(value, is.list, logical(1))])
    }
    list(value)
}

#' Coerce GraphQL records into a tibble
#'
#' Normalises the many shapes a GraphQL response can take - one object, an
#' array of objects, an empty array, or `NULL` - into a tibble whose columns
#' match the requested attributes exactly, in the requested order.
#'
#' @param records The value extracted from the response's `data` member, or
#'   `NULL`.
#' @param attributes A character vector of requested attribute names.
#'
#' @return A [tibble][tibble::tibble]. Genuinely nested values are preserved as
#'   list-columns; everything else is unlisted to atomic columns.
#'
#' @noRd
.records_to_tibble <- function(records, attributes) {
    .rows_to_tibble(.to_record_rows(records), attributes)
}

#' @noRd
.rows_to_tibble <- function(rows, attributes) {
    if (length(rows) == 0L) {
        return(empty_ensembl_tibble(attributes))
    }

    columns <- .attribute_column_names(attributes)
    out <- lapply(attributes, function(attribute) {
        .simplify_column(lapply(rows, .pluck_attribute, path = attribute))
    })
    names(out) <- columns
    tibble::as_tibble(out, .name_repair = "minimal")
}
