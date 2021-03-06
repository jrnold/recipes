

#' @name selections
#' @aliases selections
#' @aliases selection
#' @title Methods for Select Variables in Step Functions
#' @description When selecting variables or model terms in `step`
#'  functions, `dplyr`-like tools are used. The *selector* functions
#'  can choose variables based on their name, current role, data
#'  type, or any combination of these. The selectors are passed as
#'  any other argument to the step. If the variables are explicitly
#'  stated in the step function, this might be similar to:
#'
#' \preformatted{
#'   recipe( ~ ., data = USArrests) \%>\%
#'     step_pca(Murder, Assault, UrbanPop, Rape, num = 3)
#' }
#'
#'   The first four arguments indicate which variables should be
#'  used in the PCA while the last argument is a specific argument
#'  to [step_pca()].
#'
#' Note that:
#'
#'   \enumerate{
#'   \item The selector arguments should not contain functions
#'    beyond those supported (see below).
#'   \item These arguments are not evaluated until the `prep`
#'    function for the step is executed.

#'   \item The `dplyr`-like syntax allows for negative sings to
#'    exclude variables (e.g. `-Murder`) and the set of selectors will
#'    processed in order.

#'   \item A leading exclusion in these arguments (e.g. `-Murder`)
#'   has the effect of adding all variables to the list except the
#'   excluded variable(s).

#'   }
#'
#' Also, select helpers from the `tidyselect` package can also be used:
#'   [tidyselect::starts_with()], [tidyselect::ends_with()],
#'   [tidyselect::contains()], [tidyselect::matches()],
#'   [tidyselect::num_range()], [tidyselect::everything()], and
#'   [tidyselect::one_of()].
#'   For example:
#'
#' \preformatted{
#'   recipe(Species ~ ., data = iris) \%>\%
#'     step_center(starts_with("Sepal"), -contains("Width"))
#' }
#'
#' would only select `Sepal.Length`
#'
#' **Inline** functions that specify computations, such as
#'  `log(x)`, should not be used in selectors and will produce an
#'  error. A list of allowed selector functions is below.
#'
#' Columns of the design matrix that may not exist when the step
#'  is coded can also be selected. For example, when using
#'  `step_pca`, the number of columns created by feature extraction
#'  may not be known when subsequent steps are defined. In this
#'  case, using `matches("^PC")` will select all of the columns
#'  whose names start with "PC" *once those columns are created*.
#'
#' There are sets of functions that can be used to select
#'  variables based on their role or type: [has_role()] and
#'  [has_type()]. For convenience, there are also functions that are
#'  more specific: [all_numeric()], [all_nominal()],
#'  [all_predictors()], and [all_outcomes()]. These can be used in
#'  conjunction with the previous functions described for selecting
#'  variables using their names:

#'
#' \preformatted{
#'   data(biomass)
#'   recipe(HHV ~ ., data = biomass) \%>\%
#'     step_center(all_numeric(), -all_outcomes())
#' }
#'
#'   This results in all the numeric predictors: carbon, hydrogen,
#'  oxygen, nitrogen, and sulfur.
#'
#'   If a role for a variable has not been defined, it will never be
#'  selected using role-specific selectors.
#'
#' Selectors can be used in [step_interact()] in similar ways but
#'  must be embedded in a model formula (as opposed to a sequence
#'  of selectors). For example, the interaction specification
#'  could be `~ starts_with("Species"):Sepal.Width`. This can be
#'  useful if `Species` was converted to dummy variables
#'  previously using [step_dummy()].
#'
#' The complete list of allowable functions in steps:
#'
#'   \itemize{
#'     \item **By name**: [tidyselect::starts_with()],
#'       [tidyselect::ends_with()], [tidyselect::contains()],
#'       [tidyselect::matches()], [tidyselect::num_range()], and
#'       [tidyselect::everything()]
#'     \item **By role**: [has_role()],
#'       [all_predictors()], and [all_outcomes()]
#'     \item **By type**: [has_type()], [all_numeric()],
#'       and [all_nominal()]
#'   }
NULL

## These are the allowable functions for formulas in the the `terms` arguments
## to the steps or to `recipes.formula`.
name_selectors <- c("starts_with",
                    "ends_with",
                    "contains",
                    "matches",
                    "num_range",
                    "everything",
                    "one_of")

role_selectors <-
  c("has_role", "all_predictors", "all_outcomes")

type_selectors <- c("has_type", "all_numeric", "all_nominal")

selectors <-
  unique(c(name_selectors, role_selectors, type_selectors))

## Get the components of the formula split by +/-. The
## function also returns the sign
f_elements <- function(x) {
  trms_obj <- terms(x)
  ## Their order will change here (minus at the end)
  clls <- attr(trms_obj, "variables")
  ## Any formula element with a minus prefix will not
  ## have an colname in the `factor` attribute of the
  ## terms object. We will check these against the
  ## list of calls
  tmp <- colnames(attr(trms_obj, "factors"))
  kept <- vector(mode = "list", length = length(tmp))
  for (j in seq_along(tmp))
    kept[[j]] <- as.name(tmp[j])
  
  term_signs <- rep("", length(clls) - 1)
  for (i in seq_along(term_signs)) {
    ## Check to see if the elements are in the `factors`
    ## part of `terms` and these will have a + sign
    retained <- any(unlist(lapply(kept,
                                  function(x, y)
                                    any(y == x),
                                  y = clls[[i + 1]])))
    term_signs[i] <- if (retained)
      "+"
    else
      "-"
  }
  list(terms  = clls, signs = term_signs)
}

## This adds the appropriate argument based on whether the call is for
## a variable name, role, or data type.
add_arg <- function(cl) {
  func <- fun_calls(cl)
  if (func %in% name_selectors) {
    cl$vars <- quote(var_vals)
  } else {
    if (func %in% role_selectors) {
      cl$roles <- quote(role_vals)
    } else
      cl$types <- quote(type_vals)
  }
  cl
}

## This flags formulas that are not allowed. When called from `recipe.formula`
## `allowed` is NULL.
element_check <- function(x, allowed = selectors) {
  funs <- fun_calls(x)
  funs <- funs[!(funs %in% c("~", "+", "-"))]
  if (!is.null(allowed)) {
    # when called from a step
    not_good <- funs[!(funs %in% allowed)]
    if (length(not_good) > 0)
      stop(
        "Not all functions are allowed in step function selectors (e.g. ",
        paste0("`", not_good, "`", collapse = ", "),
        "). See ?selections.",
        call. = FALSE
      )
  } else {
    # when called from formula.recipe
    if (length(funs) > 0)
      stop(
        "No in-line functions should be used here; use steps to define ",
        "baking actions", call. = FALSE
      )
  }
  invisible(NULL)
}

has_selector <- function(x, allowed = selectors) {
  res <- rep(NA, length(x) - 1)
  for (i in 2:length(x))
    res[[i - 1]] <- isTRUE(fun_calls(x[[i]]) %in% allowed)
  res
}

#' Select Terms in a Step Function.
#'
#' This function bakes the step function selectors and might be
#'  useful when creating custom steps.
#'
#' @param info A tibble with columns `variable`, `type`, `role`,
#'  and `source` that represent the current state of the data. The
#'  function [summary.recipe()] can be used to get this information
#'  from a recipe.
#' @param terms A list of formulas whose right-hand side contains
#'  quoted expressions. See [rlang::quos()] for examples.
#' @keywords datagen
#' @concept preprocessing
#' @return A character string of column names or an error of there
#'  are no selectors or if no variables are selected.
#' @seealso [recipe()] [summary.recipe()]
#'   [prep.recipe()]
#' @importFrom purrr map_lgl map_if map_chr map
#' @importFrom rlang names2
#' @export
#' @examples
#' library(rlang)
#' data(okc)
#' rec <- recipe(~ ., data = okc)
#' info <- summary(rec)
#' terms_select(info = info, quos(all_predictors()))
terms_select <- function(terms, info) {
  vars <- info$variable
  roles <- info$role
  types <- info$type
  
  if (is_empty(terms)) {
    stop("At least one selector should be used", call. = FALSE)
  }
  
  ## check arguments against whitelist
  lapply(terms, element_check)
  
  # Set current_info so available to helpers
  old_info <- set_current_info(info)
  on.exit(set_current_info(old_info), add = TRUE)
  
  sel <- with_handlers(tidyselect::vars_select(vars, !!! terms),
                       tidyselect_empty = abort_selection
  )
  
  unname(sel)
}

abort_selection <- exiting(function(cnd) {
  abort("No variables or terms were selected.")
})

#' Role Selection
#'
#' `has_role`, `all_predictors`, and `all_outcomes` can be used to
#'  select variables in a formula that have certain roles.
#'  Similarly, `has_type`, `all_numeric`, and `all_nominal` are used
#'  to select columns based on their data type. See [selections()]
#'  for more details. `current_info` is an internal function that is
#'  unlikely to help users while the others have limited utility
#'  outside of step function arguments.
#'
#' @param match A single character string for the query. Exact
#'  matching is used (i.e. regular expressions won't work).
#' @param roles A character string of roles for the current set of
#'  terms.
#' @param types A character string of roles for the current set of
#'  data types.
#' @return Selector functions return an integer vector while
#'   `current_info` returns an environment with vectors `vars`,
#'   `roles`, and `types`.
#' @keywords datagen
#' @examples
#' data(biomass)
#'
#' rec <- recipe(biomass) %>%
#'   add_role(carbon, hydrogen, oxygen, nitrogen, sulfur,
#'            new_role = "predictor") %>%
#'   add_role(HHV, new_role = "outcome") %>%
#'   add_role(sample, new_role = "id variable") %>%
#'   add_role(dataset, new_role = "splitting indicator")
#' recipe_info <- summary(rec)
#' recipe_info
#'
#' has_role("id variable", roles = recipe_info$role)
#' all_outcomes(roles = recipe_info$role)
#' @export

has_role <-
  function(match = "predictor",
           roles = current_info()$roles)
    which(roles %in% match)

#' @export
#' @rdname has_role
#' @inheritParams has_role
all_predictors <- function(roles = current_info()$roles)
  has_role("predictor", roles = roles)

#' @export
#' @rdname has_role
#' @inheritParams has_role
all_outcomes <- function(roles = current_info()$roles)
  has_role("outcome", roles = roles)

#' @export
#' @rdname has_role
#' @inheritParams has_role
has_type <-
  function(match = "numeric",
           types = current_info()$types)
    which(types %in% match)

#' @export
#' @rdname has_role
#' @inheritParams has_role
all_numeric <- function(types = current_info()$types)
  has_type("numeric", types = types)

#' @export
#' @rdname has_role
#' @inheritParams has_role
all_nominal <- function(types = current_info()$types)
  has_type("nominal", types = types)

## functions to get current variable info for selectors modeled after
## dplyr versions

#' @import rlang
cur_info_env <- child_env(env_parent(env))

set_current_info <- function(x) {
  # stopifnot(!is.environment(x))
  old <- cur_info_env
  cur_info_env$vars <- x$variable
  cur_info_env$roles <- x$role
  cur_info_env$types <- x$type
  
  invisible(old)
}

#' @export
#' @rdname has_role
current_info <- function() {
  cur_info_env %||% stop("Variable context not set", call. = FALSE)
}
