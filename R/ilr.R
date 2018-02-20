#' Isometric Log Ratio Transformation of Compositional Data
#'
#' `step_ilr` creates a *specification* of a recipe
#'  step that will apply the isometric log-ratio transformation,
#'  (ILR) compositional data to have a mean of zero.
#'
#' @param recipe A recipe object. The step will be added to the
#'  sequence of operations for this recipe.
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Not used by this step since no new variables are
#'  created.
#' @param trained A logical to indicate if the quantities for
#'  preprocessing have been estimated.
#' @param means A named numeric vector of means. This is
#'  `NULL` until computed by [prep.recipe()].
#' @param na.rm A logical value indicating whether `NA`
#'  values should be removed during computations.
#' @param skip A logical. Should the step be skipped when the
#'  recipe is baked by [bake.recipe()]? While all operations are baked
#'  when [prep.recipe()] is run, some operations may not be able to be
#'  conducted on new data (e.g. processing the outcome variable(s)).
#'  Care should be taken when using `skip = TRUE` as it may affect
#'  the computations for subsequent operations
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any). For the
#'  `tidy` method, a tibble with columns `terms` (the
#'  selectors or variables selected) and `value` (the means).
#'
#' @keywords datagen
#' @concept preprocessing normalization_methods
#' @export
#' @details ilring data means that the average of a variable is
#'  subtracted from the data. `step_ilr` estimates the
#'  variable means from the data used in the `training`
#'  argument of `prep.recipe`. `bake.recipe` then applies
#'  the ilring to new data sets using these means.
#'
#' @examples
#' data(biomass)
#'
#' library('compositions')
#' data(Hydrochem)
#'
#' recipe(~ ., data = Hydrochem) %>%
#'   step_ilr(all_numeric(), -Code, -Site) %>%
#'   prep() %>%
#'   bake(newdata = Hydrochem)
#' @seealso [recipe()] [prep.recipe()]
#'   [bake.recipe()]
step_ilr <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           V = compositions::ilrBase,
           prefix = "ilr_",
           skip = FALSE) {
    add_step(
      recipe,
      step_ilr_new(
        terms = ellipse_check(...),
        trained = trained,
        V = V,
        prefix = prefix,
        role = role,
        skip = skip
      )
    )
  }

## Initializes a new object
step_ilr_new <-
  function(terms = NULL,
           role = NA,
           trained = FALSE,
           V = compositions::ilrBase,
           prefix = "ilr_",
           skip = FALSE) {
    step(
      subclass = "ilr",
      terms = terms,
      role = role,
      trained = trained,
      V = V,
      prefix = prefix,
      skip = skip
    )
  }

prep.step_ilr <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(x$terms, info = info)
  if (is.matrix(x$V) ||
      is.data.frame(x$V) ||
       is.numeric(x$V)) {
    V <- as.matrix(x$V)
    if (nrow(V) != length(col_names)) {
      stop("`V` must be have the same number of rows as selected columns.",
           .call = FALSE)
    }
  } else {
    V <- rlang::as_function(x$V)(training[ , col_names, drop = FALSE])
  }
  rownames(V) <- col_names
  step_ilr_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    V = V,
    prefix = x$prefix,
    skip = x$skip
  )
}

#' @importFrom dplyr bind_cols
bake.step_ilr <- function(object, newdata, ...) {
  col_names <- rownames(object$V)
  res <- compositions::ilr(newdata[, col_names, drop = FALSE],
                           V = object$V)
  colnames(res) <- paste0("ilr_", seq_len(ncol(res)))
  for (i in col_names) {
    newdata[[i]] <- NULL
  }
  newdata <- bind_cols(newdata, as_tibble(res))
  as_tibble(newdata)
}

print.step_ilr <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Isometric log ratio transform for ", sep = "")
    printer(rownames(x$V), x$terms, x$trained, width = width)
    invisible(x)
  }

#' @rdname step_ilr
#' @param x A `step_ilr` object.
tidy.step_ilr <- function(x, ...) {
  simple_terms(x, ...)
}
