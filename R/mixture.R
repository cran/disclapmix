#' Generate a mixture
#' 
#' This function can generate a mixture given a list of contributors.
#' 
#' 
#' @param profiles A list with profiles to mix.
#' @return A list, e.g. for use with \code{\link{contributor_pairs}}. See
#' example usage at \code{\link{rank_contributor_pairs}}.
#' @seealso \code{\link{contributor_pairs}}
#' \code{\link{rank_contributor_pairs}} \code{\link{disclapmix-package}}
#' \code{\link{disclapmix}} \code{\link{disclapmixfit}}
#' \code{\link{clusterprob}} \code{\link{predict.disclapmixfit}}
#' \code{\link{print.disclapmixfit}} \code{\link{summary.disclapmixfit}}
#' \code{\link{simulate.disclapmixfit}} %\code{\link{haplotype_diversity}}
#' \code{\link[disclap:disclap-package]{disclap}}
#' @keywords mixture separation deconvolution
#' @export
generate_mixture <- function(profiles) {
  if (!is.list(profiles)) stop("profiles must be a list")
  if (length(profiles) <= 0L) stop("No profiles provided")
  if (any(!unlist(lapply(profiles, is.integer)))) stop("Elements of profiles must be integer vectors")
  if (length(unique(unlist(lapply(profiles, length)))) != 1L) stop("All elements of profiles must have the same length")
  
  n <- length(profiles)
  loci <- length(profiles[[1L]])
  
  mixture <- lapply(1L:loci, function(k) sort(unique(unlist(lapply(1L:n, function(i) profiles[[i]][k])))))
  
  names(mixture) <- names(profiles[[1L]])
  
  return(mixture)
}

# R saves matrices in column major order, hence each column is a haplotype.
# FIXME: Optimise?
.get_cand2_set <- function(mixture, org_mixture, fixed_index, fixed_allele, candidate_pers1, loci) {
  candidate_pers2 <- matrix(0L, nrow = loci, ncol = ncol(candidate_pers1),
    dimnames = list(rownames(candidate_pers1), NULL))
  
  for (j in 1L:ncol(candidate_pers2)) {
    pers1 <- candidate_pers1[, j]
    pers2 <- rep(0L, loci)
    
    for (k in 1L:loci) {
      if (k == fixed_index) {
        pers2[fixed_index] <- fixed_allele
        next
      }
      
      mix_k <- mixture[[k]]
      
      if (length(mix_k) == 2L) {
        if (mix_k[1L] == pers1[k]) {
          pers2[k] <- mix_k[2L]
          next
        }
        
        pers2[k] <- mix_k[1L]
        next
      }
      
      pers2[k] <- mix_k
    }
    
    candidate_pers2[, j] <- pers2
    
    #mix_j <- lapply(seq_along(pers1), function(k) sort(unique(c(pers1[k], pers2[k]))))
    #if (any(unlist(lapply(seq_along(pers1), function(k) (length(mix_j[[k]]) != length(org_mixture[[k]]) || (any(mix_j[[k]] != org_mixture[[k]]))))))) {
    #  print(mixture)
    #  print(pers1)
    #  print(pers2)
    #  stop("Something is wrong with mixture synthesis.")
    #}
  }
  
  return(candidate_pers2)
}

# R saves matrices in column major order, hence each column is a haplotype.


#' Contributor pairs from a 2 person mixture
#' 
#' Get all possible contributor pairs from a 2 person mixture
#' 
#' @param mixture A list of integer vectors. The k'th element in the list is an
#' integer vector with the alleles in the mixture at locus k.
#' @return A \code{contrib_pairs} object that is a unordered list of pairs.
#' Note, that contributor order is disregarded so that each contributor pair is
#' only present once (and not twice as would be the case if taking order into
#' consideration). See example usage at \code{\link{rank_contributor_pairs}}.
#' @seealso \code{\link{rank_contributor_pairs}} \code{\link{generate_mixture}}
#' \code{\link{disclapmix-package}} \code{\link{disclapmix}}
#' \code{\link{disclapmixfit}} \code{\link{clusterprob}}
#' \code{\link{predict.disclapmixfit}} \code{\link{print.disclapmixfit}}
#' \code{\link{summary.disclapmixfit}} \code{\link{simulate.disclapmixfit}}
#' %\code{\link{haplotype_diversity}}
#' \code{\link[disclap:disclap-package]{disclap}}
#' @keywords mixture separation deconvolution
#' @export
contributor_pairs <- function(mixture) {
  if (!is.list(mixture)) stop("mixture must be a list")
  if (any(!unlist(lapply(mixture, is.integer)))) stop("Elements of mixture must be integer vectors")

  mixture <- lapply(mixture, unique)
  mixture <- lapply(mixture, sort)
  
  alleles_count <- unlist(lapply(mixture, length))
  
  if (any(alleles_count < 1L | alleles_count > 2L)) {
    stop("Elements of mixture must have one or two integer elements")
  }
  
  org_mixture <- mixture
      
  # To avoid doing double the number of required calculations
  fixed_index <- which(alleles_count > 1L)[1L]
  fixed_allele <- NA
  
  if (!is.na(fixed_index)) {
    #fixed_allele <- mixture[[fixed_index]][1L]
    #mixture[[fixed_index]] <- mixture[[fixed_index]][2L]
    fixed_allele <- mixture[[fixed_index]][2L]
    mixture[[fixed_index]] <- mixture[[fixed_index]][1L]
  } else {
    warning("All loci have one allele, skipping")
    return(NULL)
  }
  
  number_of_candidate_pairs <- prod(sapply(mixture, length))
  
  # Column major order, each column is a haplotype
  # FIXME: Optimise?
  candidate_pers1 <- t(as.matrix(expand.grid(mixture)))

  if (ncol(candidate_pers1) != number_of_candidate_pairs) {
    stop("Unexpected number of candidates for person 1")  
  }
  
  candidate_pers2 <- .get_cand2_set(mixture = mixture, org_mixture = org_mixture, fixed_index = fixed_index, fixed_allele = fixed_allele, candidate_pers1 = candidate_pers1, loci = nrow(candidate_pers1))
  
  ans <- list(
    mixture = org_mixture, 
    pairs_count = number_of_candidate_pairs, 
    pairs_person1 = candidate_pers1,
    pairs_person2 = candidate_pers2
  )
  
  class(ans) <- c("contrib_pairs", class(ans))
  return(ans)
}



#' Separate a 2 person mixture
#' 
#' Separate a 2 person mixture by ranking the possible contributor pairs.
#' 
#' @param contrib_pairs A \code{contrib_pairs} object obtained from
#' \code{\link{contributor_pairs}}.
#' @param fit A \code{\link{disclapmixfit}} object.
#' @param max_rank Not used. Reserved for future use.
#' @return A \code{ranked_contrib_pairs} object that is basically an order
#' vector and the probabilities for each pair (in the same order as given in
#' \code{contrib_pairs}), found by using \code{fit}. Note, that contributor
#' order is disregarded so that each contributor pair is only present once (and
#' not twice as would be the case if taking order into consideration).
#' @seealso \code{\link{contributor_pairs}} \code{\link{generate_mixture}}
#' \code{\link{disclapmix-package}} \code{\link{disclapmix}}
#' \code{\link{disclapmixfit}} \code{\link{clusterprob}}
#' \code{\link{predict.disclapmixfit}} \code{\link{print.disclapmixfit}}
#' \code{\link{summary.disclapmixfit}} \code{\link{simulate.disclapmixfit}}
#' %\code{\link{haplotype_diversity}}
#' \code{\link[disclap:disclap-package]{disclap}}
#' @keywords mixture separation deconvolution
#' @examples
#' 
#' data(danes)
#' db <- as.matrix(danes[rep(1L:nrow(danes), danes$n), 1L:(ncol(danes) - 1L)])
#' 
#' set.seed(1)
#' true_contribs <- sample(1L:nrow(db), 2L)
#' h1 <- db[true_contribs[1L], ]
#' h2 <- db[true_contribs[2L], ]
#' db_ref <- db[-true_contribs, ]
#' 
#' h1h2 <- c(paste(h1, collapse = ";"), paste(h2, collapse = ";"))
#' tab_db <- table(apply(db, 1, paste, collapse = ";"))
#' tab_db_ref <- table(apply(db_ref, 1, paste, collapse = ";"))
#' tab_db[h1h2]
#' tab_db_ref[h1h2]
#' 
#' rm(db) # To avoid use by accident
#' 
#' mixture <- generate_mixture(list(h1, h2))
#' 
#' possible_contributors <- contributor_pairs(mixture)
#' possible_contributors
#' 
#' fits <- lapply(1L:5L, function(clus) disclapmix(db_ref, clusters = clus))
#' 
#' best_fit_BIC <- fits[[which.min(sapply(fits, function(fit) fit$BIC_marginal))]]
#' best_fit_BIC
#' 
#' ranked_contributors_BIC <- rank_contributor_pairs(possible_contributors, best_fit_BIC)
#' ranked_contributors_BIC
#' 
#' plot(ranked_contributors_BIC, top = 10L, type = "b")
#' 
#' get_rank(ranked_contributors_BIC, h1)
#' 
#' @export
rank_contributor_pairs <- function(contrib_pairs, fit, max_rank = NULL) {
  if (!is(contrib_pairs, "contrib_pairs")) stop("contrib_pairs must be a contrib_pairs object")
  if (!is(fit, "disclapmixfit")) stop("fit must be a disclapmixfit object")
  if (!is.null(max_rank)) stop("max_rank must be NULL (reserved for future use)")
  
  pers1_prob <- predict.disclapmixfit(fit, newdata = t(contrib_pairs$pairs_person1))
  pers2_prob <- predict.disclapmixfit(fit, newdata = t(contrib_pairs$pairs_person2))
  
  pers12_prob <- pers1_prob * pers2_prob
  pers12_prob_order <- order(pers12_prob, decreasing = TRUE)
  
  #pers1_prob_ordered <- pers1_prob[pers12_prob_order]
  #pers2_prob_ordered <- pers2_prob[pers12_prob_order]
  #pers12_prob_ordered <- pers12_prob[pers12_prob_order]
  
  #candidate_pers1_ordered <- candidate_pers1[pers12_prob_order, , drop = FALSE]
  #candidate_pers2_ordered <- candidate_pers2[pers12_prob_order, , drop = FALSE]
  
  # Rearrange (takes more memory, but should be easier to understand and deal with output)
  #ranked_pairs <- vector("list", nrow(candidate_pers1_ordered))
  
  #for (i in 1L:nrow(candidate_pers1)) {
  #  mat <- rbind(candidate_pers1_ordered[i, ], candidate_pers2_ordered[i, ])
  #  colnames(mat) <- names(contrib_pairs$mixture)
  #  rownames(mat) <- c("H1", "H2")
  #  
  #  probs <- c("P(H1)" = pers1_prob_ordered[i], "P(H2)" = pers2_prob_ordered[i], "P(H1)P(H2)" = pers12_prob_ordered[i])
  #  
  #  ranked_pairs[[i]] <- list(contrib_pair = mat, probs = probs)
  #}
  
  ans <- list(
    mixture = contrib_pairs$mixture,
    fit = fit,
    info = list(
      max_rank = max_rank, 
      pairs_count = contrib_pairs$pairs_count,
      pairs_prod_sum = sum(pers12_prob)
    ), 
    pairs_person1 = contrib_pairs$pairs_person1,
    pairs_person2 = contrib_pairs$pairs_person2,
    probs_person1 = pers1_prob,
    probs_person2 = pers2_prob,
    probs_pair = pers12_prob,
    #ranked_pairs = ranked_pairs
    order = pers12_prob_order
  )
  
  class(ans) <- c("ranked_contrib_pairs", class(ans))
      
  return(ans)
}

#' Print contributor pairs
#' 
#' @param x A \code{contrib_pairs} object.
#' @param \dots Ignored
#' 
#' @export
print.contrib_pairs <- function(x, ...) {
  if (!is(x, "contrib_pairs")) stop("x must be a contrib_pairs")
  
  cat("Mixture:\n")
  print(data.frame(do.call(cbind, lapply(x$mixture, paste, collapse = ","))))
  cat("\n")
  
  cat("Number of possible contributor pairs = ", x$pairs_count, "\n", sep = "")

  return(invisible(NULL))
}

#' Print ranked contributor pairs
#' 
#' @param x A \code{ranked_contrib_pairs} object.
#' @param top The top ranked number of pairs to print/plot. \code{NULL} for
#' all.
#' @param hide_non_varying_loci Whether to hide alleles on loci that do not
#' vary.
#' @param \dots Ignored
#' 
#' @export
print.ranked_contrib_pairs <- function(x, top = 5L, hide_non_varying_loci = TRUE, ...) {

  if (!is(x, "ranked_contrib_pairs")) stop("x must be a ranked_contrib_pairs")
  
  if (is.null(top)) {
    top <- x$info$pairs_count
  } else {
    if (!is.integer(top) || length(top) != 1L || top < 1L) {
      stop("top must be one single integer (L postfix), e.g. 10L")
    }
  }
  
  if (length(hide_non_varying_loci) != 1L || !is.logical(hide_non_varying_loci)) stop("hide_non_varying_loci must be TRUE or FALSE")
  
  not_printed <- NA
  top <- min(top, x$info$pairs_count)
    
  if (x$info$pairs_count > top) {
    not_printed <- x$info$pairs_count - top
  }
  
  cat("Mixture:\n")
  print(data.frame(do.call(cbind, lapply(x$mixture, paste, collapse = ","))))
  cat("\n")
  
  cat("Contributor pairs = ", x$info$pairs_count, "\n", sep = "")
  cat("\n")
  
  cat("Sum of all (product of contributor pair haplotypes) = ", x$info$pairs_prod_sum, "\n", sep = "")
  cat("\n")
  
  if (top == 1L) {
    cat("Showing rank 1:\n", sep = "")
  } else {
    cat("Showing rank 1-", top, ":\n", sep = "")
  }
  
  hide_loci <- rep(FALSE, length(x$mixture))
  if (hide_non_varying_loci) {
    hide_loci <- unlist(lapply(x$mixture, length)) == 1L
  }
  
  if (sum(hide_loci) == 0L) {
    hide_loci <- rep(FALSE, length(x$mixture))
    hide_non_varying_loci <- FALSE
  }
  
  for (i in 1L:top) {
    index <- x$order[i]
    cat("\nRank ", i, " [ P(H1)*P(H2) = ", x$probs_pair[index], " ]:\n", sep = "")
    
    h1 <- data.frame(rbind(x$pairs_person1[, index]), Prob = x$probs_person1[index])
    h2 <- data.frame(rbind(x$pairs_person2[, index]), Prob = x$probs_person2[index])
    
    index_df <- data.frame(rbind(H1 = h1, H2 = h2))

    if (hide_non_varying_loci) {
      index_df[, c(hide_loci, FALSE)] <- "."
    }
    
    print(index_df)
  }
  
  if (!is.na(not_printed) && not_printed >= 1L) {
    cat("\n")
    cat(" (", not_printed, " contributor pairs hidden.)\n", sep = "")
  } 
  
  return(invisible(NULL))
}

#' Plot ranked contributor pairs
#' 
#' @param x A \code{ranked_contrib_pairs} object.
#' @param top The top ranked number of pairs to print. \code{NULL} for
#' all.
#' @param \dots Delegated to the generic \code{\link[graphics]{plot}} function.
#' @param xlab Graphical parameter.
#' @param ylab Graphical parameter.
#' @export
plot.ranked_contrib_pairs <- function(x, top = NULL, ..., xlab = "Rank", ylab = "P(H1)P(H2)") {

  if (!is(x, "ranked_contrib_pairs")) stop("x must be a ranked_contrib_pairs")
  
  if (is.null(top)) {
    top <- x$info$pairs_count
  } else {
    if (!is.integer(top) || length(top) != 1L || top < 1L) {
      stop("top must be one single integer (L postfix), e.g. 10L")
    }
  }
  
  top <- min(top, x$info$pairs_count)
  
  xs <- numeric(top)
  ys <- numeric(top)
  
  for (i in 1L:top) {
    index <- x$order[i]
    xs[i] <- i
    ys[i] <- x$probs_pair[index]
  }

  return(plot(xs, ys, xlab = xlab, ylab = ylab, ...))
}

#' Get rank of pair
#' 
#' @param x A \code{ranked_contrib_pairs} object.
#' @param haplotype A haplotype.
#' @export
get_rank <- function(x, haplotype) {
  if (!is(x, "ranked_contrib_pairs")) stop("x must be a ranked_contrib_pairs")
  
  mixture <- x$mixture
  loci <- length(mixture)
  pairs_count <- x$info$pairs_count
  max_rank <- x$info$max_rank

  if (is.null(haplotype) || !is.integer(haplotype) || length(haplotype) != loci) {
    stop("haplotype must be an integer vector of length ", loci)
  }
  
  # Check if haplotype is consistent with mixture
  consistent <- TRUE
  
  for (k in 1L:loci) {
    if (!(haplotype[k] %in% mixture[[k]])) {
      consistent <- FALSE
      break
    }
  }
  
  if (!consistent) {
    stop("The provided haplotype is not consistent with the mixture")
  }
  
  # Check the rank:  
  rank <- NA

  for (i in 1L:ncol(x$pairs_person1)) {
    index <- x$order[i]
    h1 <- x$pairs_person1[, index]
    h2 <- x$pairs_person2[, index]

    if (all(h1 == haplotype) || all(h2 == haplotype)) {
      rank <- i
      break
    }
  }
  
  if (is.na(rank)) {  
    if (is.null(max_rank)) {
      stop("Unexpected condition occured")
    } else {
      warning("Only a subset of all possible contributor pairs was created and the haplotype was not among these. Hence, it has a rank somewhere between ", (max_rank + 1L), " and ", pairs_count, ".")
    }
  }
  
  return(rank)
}

