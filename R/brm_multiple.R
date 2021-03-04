#' Run the same \pkg{brms} model on multiple datasets
#' 
#' Run the same \pkg{brms} model on multiple datasets and then combine the
#' results into one fitted model object. This is useful in particular for
#' multiple missing value imputation, where the same model is fitted on multiple
#' imputed data sets. Models can be run in parallel using the \pkg{future}
#' package.
#' 
#' @inheritParams brm
#' @param data A \emph{list} of data.frames each of which will be used to fit a
#'   separate model. Alternatively, a \code{mids} object from the \pkg{mice}
#'   package.
#' @param data2 A \emph{list} of named lists each of which will be used to fit a
#'   separate model. Each of the named lists contains objects representing data
#'   which cannot be passed via argument \code{data} (see \code{\link{brm}} for
#'   examples). The length of the outer list should match the length of the list
#'   passed to the \code{data} argument.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled for every imputed data set. Defaults to \code{FALSE}. If
#'   \code{NULL}, \code{brm_multiple} tries to figure out internally, if recompilation
#'   is necessary, for example because data-dependent priors have changed.
#'   Using the default of no recompilation should be fine in most cases.
#' @param combine Logical; Indicates if the fitted models should be combined
#'   into a single fitted model object via \code{\link{combine_models}}.
#'   Defaults to \code{TRUE}.
#' @param fit An instance of S3 class \code{brmsfit_multiple} derived from a
#'   previous fit; defaults to \code{NA}. If \code{fit} is of class
#'   \code{brmsfit_multiple}, the compiled model associated with the fitted
#'   result is re-used and all arguments modifying the model code or data are
#'   ignored. It is not recommended to use this argument directly, but to call
#'   the \code{\link[brms:update.brmsfit_multiple]{update}} method, instead.
#' @param file Either \code{NULL} or a character string. In the latter case, the
#'   fitted model(s) object is saved via \code{\link{saveRDS}} in a file named
#'   after the string supplied in \code{file}. The \code{.rds} extension is
#'   added automatically. If the file already exists, \code{brm} will load and
#'   return the saved model object instead of refitting the model. 
#'   Unless you specify the \code{file_refit} argument as well, the existing
#'   files won't be overwritten, you have to manually remove the file in order
#'   to refit and save the model under an existing file name. The file name
#'   is stored in the \code{brmsfit} object for later usage.
#' @param file_refit Modifies when the fits stored via the \code{file} parameter
#'   is re-used. For \code{"never"} (default) the fit is always loaded if it
#'   exists and fitting is skipped.  The fits for individual datasets will be
#'   refit if the model, data or algorithm as passed to Stan differ from
#'   what is stored in the file. This also covers changes in priors,
#'   \code{sample_prior}, \code{stanvars}, covariance structure, etc. If you
#'   believe there was a false positive, you can use
#'   \code{\link{brmsfit_needs_refit}} to see why refit is deemed necessary.
#'   The decision is done for each fit, so if you e.g. add more imputed 
#'   datasets, only the new (or modified) datasets will be refit.
#'   Refit will not be triggered for changes in additional parameters of the fit
#'   (e.g., initial values, number of iterations, control arguments, ...). A
#'   known limitation is that a refit will be triggered if within-chain
#'   parallelization is switched on/off.
#'   If set to \code{"on_change"}, 
#'   \code{combine} needs to be \code{FALSE}.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @details The combined model may issue false positive convergence warnings, as
#'   the MCMC chains corresponding to different datasets may not necessarily
#'   overlap, even if each of the original models did converge. To find out
#'   whether each of the original models converged, investigate
#'   \code{fit$rhats}, where \code{fit} denotes the output of
#'   \code{brm_multiple}.
#' 
#' @return If \code{combine = TRUE} a \code{brmsfit_multiple} object, which
#'   inherits from class \code{brmsfit} and behaves essentially the same. If
#'   \code{combine = FALSE} a list of \code{brmsfit} objects.
#'   
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes2)
#' 
#' # fit the model using mice and lm
#' fit_imp1 <- with(lm(bmi ~ age + hyp + chl), data = imp)
#' summary(pool(fit_imp1))
#' 
#' # fit the model using brms
#' fit_imp2 <- brm_multiple(bmi ~ age + hyp + chl, data = imp, chains = 1)
#' summary(fit_imp2)
#' plot(fit_imp2, pars = "^b_")
#' # investigate convergence of the original models
#' fit_imp2$rhats
#' 
#' # use the future package for parallelization
#' library(future)
#' plan(multiprocess)
#' fit_imp3 <- brm_multiple(bmi~age+hyp+chl, data = imp, chains = 1)
#' summary(fit_imp3)
#' }
#' 
#' @export
brm_multiple <- function(formula, data, family = gaussian(), prior = NULL, 
                         data2 = NULL, autocor = NULL, cov_ranef = NULL, 
                         sample_prior = c("no", "yes", "only"), 
                         sparse = NULL, knots = NULL, stanvars = NULL,
                         stan_funs = NULL, recompile = FALSE,
                         combine = TRUE, fit = NA,
                         seed = NA, file = NULL, file_refit = "never",
                         silent = TRUE,  ...) {
  
  combine <- as_one_logical(combine)
  file_refit <- match.arg(file_refit, c("never", "on_change"))
  fits_from_file <- NULL
  if (!is.null(file)) {
    if(combine) {
      # optionally load saved model object
      if (file_refit == "on_change") {
        stop2("Cannot use 'file_refit = \"on_change\"' with `combine = TRUE` yet.")
      }
      fits <- read_brmsfit(file)
      if (!is.null(fits)) {
        return(fits)
      }
    } else {
      fits_from_file <- read_brmsfit_list(file)
      if(file_refit == "never" && !is.null(fits_from_file)) {
        return(fits_from_file)
      }
    }
  }
  
  recompile <- as_one_logical(recompile)
  data_name <- substitute_name(data)
  if (inherits(data, "mids")) {
    require_package("mice", version = "3.0.0")
    data <- lapply(seq_len(data$m), mice::complete, data = data)
  } else if (!is_data_list(data)) {
    stop2("'data' must be a list of data.frames.")
  }
  if (!is.null(data2)) {
    if (!is_data2_list(data2)) {
      stop2("'data2' must be a list of named lists.")
    }
    if (length(data2) != length(data)) {
      stop2("'data2' must have the same length as 'data'.")
    }
  }

  brm_args <- nlist(
    formula, data = data[[1]], family, prior, data2 = data2[[1]], 
    autocor, cov_ranef, sample_prior, sparse, knots, stanvars, 
    stan_funs, seed, silent, ...
  )
  
  # Extracted into a function to lazily compile only once actually needed
  compile_model_when_needed <- function() {
    brm_args$chains <- 0
    message("Compiling the C++ model")
    suppressMessages(do_call(brm, brm_args))
  }  
    
  if (is.brmsfit(fit)) {
    # avoid complications when updating the model
    class(fit) <- setdiff(class(fit), "brmsfit_multiple")
  } 
  dots <- list(...)
  # allow compiling the model without sampling (#671)
  if (isTRUE(dots$chains == 0) || isTRUE(dots$iter == 0)) {
    if(!is.brmsfit(fit)) {
      fit <- compile_model_when_needed()
    }
    class(fit) <- c("brmsfit_multiple", class(fit))
    return(fit)
  }
  
  brm_args_for_empty_fit <- brm_args
  brm_args_for_empty_fit$empty <- TRUE
  any_needed_refit <- FALSE
  if(!is.null(fits_from_file)) {
    if(is.null(recompile) || !recompile) {
      # If recompilation is off, I need to generate the Stan code for 
      # checking cache only once, using the first dataset
      # otherwise some aspects of the code might change - notably the
      # default priors for intercept.
      empty_fit <- suppressMessages(do.call(brm, brm_args_for_empty_fit))
      shared_code <- empty_fit$model
    }
  } 
  
  fits <- futures <- rhats <- vector("list", length(data))
  
  for (i in seq_along(data)) {
    needs_refit <- TRUE
    if(!is.null(fits_from_file) && !is.null(fits_from_file[[i]])) {
      
      # Build an empty fit to reconstruct standata (and potentially code)
      brm_args_for_empty_fit$data <- data[[i]]
      brm_args_for_empty_fit$data2 <- data2[[i]]
      empty_fit <- suppressMessages(do.call(brm, brm_args_for_empty_fit))
      
      if(is.null(recompile) || !recompile) {
        code_for_cache <- shared_code
      } else {
        code_for_cache <- empty_fit$model
      }
      
      needs_refit <- brmsfit_needs_refit(fits_from_file[[i]],
                                         sdata = standata(empty_fit),
                                         scode = code_for_cache,
                                         algorithm = empty_fit$algorithm,
                                         silent = silent,
                                         verbose = TRUE)
    }
      
    if(needs_refit) {
      any_needed_refit <- TRUE
      if(!is.brmsfit(fit)) {
        fit <- compile_model_when_needed()
      }
      futures[[i]] <- future::future(
        update(fit, newdata = data[[i]], data2 = data2[[i]],
               recompile = recompile, silent = silent, ...),
        packages = "brms", seed = TRUE
      )
    } else {
      futures[[i]] <- NA
      fits[[i]] <- fits_from_file[[i]]
    }
  }
  for (i in seq_along(data)) {
    if(!is.logical(futures[[i]])) {
      message("Fitting imputed model ", i)
      fits[[i]] <- future::value(futures[[i]]) 
    }
    
    rhats[[i]] <- data.frame(as.list(rhat(fits[[i]])))
    if (any(rhats[[i]] > 1.1, na.rm = TRUE)) {
      warning2("Imputed model ", i, " did not converge.")
    }
  }
  if (combine) {
    fits <- combine_models(mlist = fits, check_data = FALSE)
    attr(fits$data, "data_name") <- data_name
    fits$rhats <- do_call(rbind, rhats)
    class(fits) <- c("brmsfit_multiple", class(fits))
    if (!is.null(file)) {
      write_brmsfit(fits, file)
    }
  } else {
    # Write the fit only when we actually did refit at least some models
    if (!is.null(file) && any_needed_refit) {
      write_brmsfit_list(fits, file)
    }
  }
  fits
}

#' Combine Models fitted with \pkg{brms}
#' 
#' Combine multiple \code{brmsfit} objects, which fitted the same model.
#' This is usefully for instance when having manually run models in parallel.
#' 
#' @param ... One or more \code{brmsfit} objects.
#' @param mlist Optional list of one or more \code{brmsfit} objects.
#' @param check_data Logical; indicates if the data should be checked
#'   for being the same across models (defaults to \code{TRUE}).
#'   Setting it to \code{FALSE} may be useful for instance
#'   when combining models fitted on multiple imputed data sets.
#'   
#' @details This function just takes the first model and replaces 
#'   its \code{stanfit} object (slot \code{fit}) by the combined 
#'   \code{stanfit} objects of all models.
#'   
#' @return A \code{brmsfit} object.
#' 
#' @export
combine_models <- function(..., mlist = NULL, check_data = TRUE) {
  models <- c(list(...), mlist)
  check_data <- as_one_logical(check_data)
  if (!length(models)) {
    stop2("No models supplied to 'combine_models'.")
  }
  for (i in seq_along(models)) {
    if (!is.brmsfit(models[[i]])) {
      stop2("Model ", i, " is no 'brmsfit' object.")
    }
    models[[i]] <- restructure(models[[i]])
  }
  ref_formula <- formula(models[[1]])
  ref_pars <- parnames(models[[1]])
  ref_mf <- model.frame(models[[1]]) 
  for (i in seq_along(models)[-1]) {
    if (!is_equal(formula(models[[i]]), ref_formula)) {
      stop2("Models 1 and ", i, " have different formulas.")
    }
    if (!is_equal(parnames(models[[i]]), ref_pars)) {
      stop2("Models 1 and ", i, " have different parameters.")
    }
    if (check_data && !is_equal(model.frame(models[[i]]), ref_mf)) {
      stop2(
        "Models 1 and ", i, " have different data. ", 
        "Set 'check_data' to FALSE to turn off checking of the data."
      )
    }
  }
  sflist <- lapply(models, "[[", "fit")
  models[[1]]$fit <- rstan::sflist2stanfit(sflist)
  models[[1]]
}

# validity check for 'data' input of 'brm_multiple'
is_data_list <- function(x) {
  is.list(x) && is.vector(x)
}

# validity check for 'data2' input of 'brm_multiple'
is_data2_list <- function(x) {
  is.list(x) && all(ulapply(x, function(y) is.list(y) && is_named(y)))
}

warn_brmsfit_multiple <- function(x, newdata = NULL) {
  if (is.brmsfit_multiple(x) && is.null(newdata)) {
    warning2(
      "Using only the first imputed data set. Please interpret the results ", 
      "with caution until a more principled approach has been implemented."
    )
  }
  invisible(x)
}
