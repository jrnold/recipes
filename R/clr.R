#' Centered Log Ratio Transformation of Compositional Data
#'
#' `step_clr` creates a *specification* of a recipe
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
#' @details clring data means that the average of a variable is
#'  subtracted from the data. `step_clr` estimates the
#'  variable means from the data used in the `training`
#'  argument of `prep.recipe`. `bake.recipe` then applies
#'  the clring to new data sets using these means.
#'
#' @examples
#' data(biomass)
#'
#' library('compositions')
#' data(Hydrochem)
#'
#' recipe(~ ., data = Hydrochem) %>%
#'   step_clr(all_numeric(), -Code, -Site) %>%
#'   prep() %>%
#'   bake(newdata = Hydrochem)
#' @seealso [recipe()] [prep.recipe()]
#'   [bake.recipe()]
step_clr <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           columns = NULL,
           skip = FALSE) {
    add_step(
      recipe,
      step_clr_new(
        terms = ellipse_check(...),
        trained = trained,
        columns = columns,
        role = role,
        skip = skip
      )
    )
  }

## Initializes a new object
step_clr_new <-
  function(terms = NULL,
           role = NA,
           trained = FALSE,
           columns = NULL,
           skip = FALSE) {
    step(
      subclass = "clr",
      terms = terms,
      role = role,
      trained = trained,
      columns = columns,
      skip = skip
    )
  }

prep.step_clr <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(x$terms, info = info)
  step_clr_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    columns = col_names,
    skip = x$skip
  )
}

bake.step_clr <- function(object, newdata, ...) {
  res <- compositions::clr(newdata[, object$columns, drop = FALSE])
  newdata[, object$columns] <- res
  as_tibble(newdata)
}

print.step_clr <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Centered log ratio transform for ", sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }

#' @rdname step_clr
#' @param x A `step_clr` object.
tidy.step_clr <- function(x, ...) {
  simple_terms(x, ...)
}
