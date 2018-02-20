#' Additive Log Ratio Transformation of Compositional Data
#'
#' `step_alr` creates a *specification* of a recipe
#'  step that will normalize numeric data to have a mean of zero.
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
#'  subtracted from the data. `step_alr` estimates the
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
#'   step_alr(all_numeric(), -Code, -Site) %>%
#'   prep() %>%
#'   bake(newdata = Hydrochem)
#' @seealso [recipe()] [prep.recipe()]
#'   [bake.recipe()]
step_alr <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           columns = NULL,
           ivar = NULL,
           prefix = "ilr_",
           skip = FALSE) {
    add_step(
      recipe,
      step_alr_new(
        terms = ellipse_check(...),
        trained = trained,
        columns = columns,
        ivar = ivar,
        prefix = prefix,
        role = role,
        skip = skip
      )
    )
  }

## Initializes a new object
step_alr_new <-
  function(terms = NULL,
           role = NA,
           trained = FALSE,
           columns = NULL,
           ivar = NULL,
           prefix = "alr_",
           skip = FALSE) {
    step(
      subclass = "alr",
      terms = terms,
      role = role,
      trained = trained,
      ivar = ivar,
      columns = columns,
      prefix = prefix,
      skip = skip
    )
  }

prep.step_alr <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(x$terms, info = info)
  ivar <- x$ivar %||% col_names[length(col_names)]
  step_alr_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    columns = col_names,
    ivar = ivar,
    prefix = x$prefix,
    skip = x$skip
  )
}

bake.step_alr <- function(object, newdata, ...) {
  col_names <- object$columns
  res <- compositions::alr(newdata[, col_names, drop = FALSE],
                           ivar = object$ivar)
  new_colnames <- setdiff(colnames, object$ivar)
  newdata[[object$ivar]] <- NULL
  newdata[ , new_colnames] <- res
  as_tibble(newdata)
}

print.step_alr <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Additive log ratio transform for ", sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }


#' @rdname step_alr
#' @param x A `step_alr` object.
tidy.step_alr <- function(x, ...) {
  simple_terms(x, ...)
}
