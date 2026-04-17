# Code to compute efficient DiD estimator with variation in treatment timing, 
# in the simulations

efficient_did_unc_stagg <- function(
    data,
    yname,
    tname,
    gname,
    idname,
    return_weights = FALSE
){
  #-----------------------------------------------------------------------------
  df <- as.data.frame(data)
  Y <- df[[yname]]
  Time <- df[[tname]]
  G <- df[[gname]]
  sample_size = length(unique(df[[idname]]))
  G_cs <- G[Time==min(Time)]
  
  # list of dates from smallest to largest
  tlist <- unique(df[,tname])[order(unique(df[,tname]))]
  glist <- unique(df[,gname])[order(unique(df[,gname]))]
  #-----------------------------------------------------------------------------
  g_treated <- glist[glist>0]
  #-----------------------------------------------------------------------------
  pi2 <- rep(0, length(glist))
  for (i in 1:length(glist)) {
    pi2[i] = mean(G==glist[i])
  }
  #-----------------------------------------------------------------------------
  # Rule: (g',t''): if g'=g, then t_min <= t'' <g'
  # else, g' \in G' :g' <> g and t_min < t'' <g' 
  
  # Collect all (g,t_post) in the data
  gt_gprime_tprime <- expand.grid(g = g_treated, 
                                  tt = tlist, 
                                  g_prime = g_treated, 
                                  t_prime = tlist)
  # Sort gt_gprime_tprime by g and time
  gt_gprime_tprime <- gt_gprime_tprime[order(gt_gprime_tprime$g,
                                             gt_gprime_tprime$tt,
                                             gt_gprime_tprime$g_prime,
                                             gt_gprime_tprime$t_prime),]
  # Keep only observations that satisfy selection rule
  gt_gprime_tprime <- gt_gprime_tprime |>
    filter(tt>=g) |> 
    dplyr::mutate(keeper = ifelse((g == g_prime ) & (t_prime < g_prime),
                           1, 0)) |>
    dplyr::mutate(keeper = ifelse((g != g_prime ) &
                             ((min(t_prime) < t_prime)*(t_prime < g_prime))==1,
                           1, keeper))  |>
    filter(keeper == 1)
  
  #-----------------------------------------------------------------------------
  # Function to compute influence function
  IF_gprime_tprime <- function(g, 
                               tt,
                               g_prime,
                               t_prime, 
                               Y, 
                               G,
                               Time){
    
    # First, compute the average terms we will need for the influence function
    mean_g_t <- mean(Y[G==g & Time==tt]) - 
      mean(Y[G==g & Time==min(Time)])
    
    mean_inf_t_tprime <- mean(Y[G==0 & Time==tt]) - 
      mean(Y[G==0 & Time==t_prime])
    
    mean_gprime_tprime <- mean(Y[G==g_prime & Time==t_prime]) - 
      mean(Y[G==g_prime & Time==min(Time)])
    
    # Now we are ready to start collecting the IFs of each piece
    IF_g <- ((G_cs == g)/pi2[which(glist == g)]) *
      (Y[Time == tt] - Y[Time == min(Time)] - mean_g_t)
    
    IF_inf <-  ((G_cs == 0)/pi2[which(glist == 0)]) *
      (Y[Time == tt] - Y[Time == t_prime] - mean_inf_t_tprime)
    
    IF_g_prime <- ((G_cs == g_prime)/pi2[which(glist == g_prime)]) *
      (Y[Time == t_prime] - Y[Time == min(Time)] - mean_gprime_tprime)
    
    # Combine the IFs
    IFs <- IF_g - (IF_inf + IF_g_prime)
    att <- mean_g_t - (mean_inf_t_tprime + mean_gprime_tprime)
    return(list(att = att, IFs = IFs))
  }
  #-----------------------------------------------------------------------------
  # Compute all IFs
  #-----------------------------------------------------------------------------
  # Check how many pairs we have
  n_pairs <- nrow(gt_gprime_tprime)
  # We will populate this with att_gts
  att_gt_list <- list()
  # place holder in lists
  counter <- 1
  for (k in 1:n_pairs){
    # Get the influence function
    IFs <- IF_gprime_tprime(
      g = gt_gprime_tprime$g[k],
      tt = gt_gprime_tprime$tt[k],
      g_prime = gt_gprime_tprime$g_prime[k],
      t_prime = gt_gprime_tprime$t_prime[k],
      Y = Y,
      G = G,
      Time = Time
    )
    
    att = IFs$att
    inf_function = IFs$IFs
    
    att_gt_list[[counter]] <- 
      list(
        group = gt_gprime_tprime$g[k],
        tt = gt_gprime_tprime$tt[k],
        g_prime = gt_gprime_tprime$g_prime[k],
        t_prime = gt_gprime_tprime$t_prime[k],
        att = att,
        inf_function = inf_function
      )
    
    counter <- counter + 1
  }
  
  
  # Process the att_gt_list
  process_attgt <- function(att_gt_list) {
    size_list <- length(att_gt_list)
    # create vectors to hold the results
    group <- c()
    att <- c()
    tt <- c()
    g_prime <- c()
    t_prime <- c()
    inf_function <- matrix(0, nrow = sample_size, ncol = size_list)
    
    # populate result vectors and matrices
    for (i in 1:size_list) {
      group[i] <- att_gt_list[[i]]$group
      tt[i] <- att_gt_list[[i]]$tt
      g_prime[i] <- att_gt_list[[i]]$g_prime
      t_prime[i] <- att_gt_list[[i]]$t_prime
      inf_function[,i] <- att_gt_list[[i]]$inf_function
      att[i] <- att_gt_list[[i]]$att
    }
    
    list(group=group, 
         tt = tt,
         g_prime = g_prime,
         t_prime = t_prime,
         inf_function = inf_function,
         att=att)
  }
  att_gt_results <- process_attgt(att_gt_list)
  
  # Now compute the efficient IF for each group-time pair
  unique_gt_pairs <- data.frame(g = att_gt_results$group, 
                                tt = att_gt_results$tt) |>
    unique()|>
    dplyr::mutate(estimate = 0,
           std_error = 0)
  
  #Length of unique g-t pairs
  n_gt <- nrow(unique_gt_pairs)
  eff_inf_function <- matrix(0, nrow = sample_size, ncol = n_gt)
  weights_out <- matrix(0, nrow = n_pairs, ncol = 7) 
  for (k in 1:n_gt) {
    
    IF_g_t <- att_gt_results$inf_function[, 
                                          which(att_gt_results$group == unique_gt_pairs$g[k] & 
                                                  att_gt_results$tt == unique_gt_pairs$tt[k])]
    
    theta <- att_gt_results$att[which(att_gt_results$group == unique_gt_pairs$g[k] & 
                                        att_gt_results$tt == unique_gt_pairs$tt[k])]
    if(return_weights){
      # Get the g_prime associated with the g and tt
      g_prime_g_t <- att_gt_results$g_prime[which(att_gt_results$group == unique_gt_pairs$g[k] & 
                                          att_gt_results$tt == unique_gt_pairs$tt[k])]
      
      t_prime_g_t <- att_gt_results$t_prime[which(att_gt_results$group == unique_gt_pairs$g[k] & 
                                                    att_gt_results$tt == unique_gt_pairs$tt[k])]
    }
    #-----------------------------------------------------------------------------
    # Compute Omega and Theta
    #-----------------------------------------------------------------------------
    Omega <- cov(IF_g_t)
    Omega_inv <- solve(Omega)
    w <- colSums(Omega_inv) / sum(Omega_inv)
    unique_gt_pairs$estimate[k] <- sum(w*theta)
    Asy_var <- 1/sum(Omega_inv)
    unique_gt_pairs$std_error[k] <- sqrt(Asy_var/sample_size)
    eff_inf_function[,k] <- IF_g_t %*% w
    
    weights = NULL
    if(return_weights){
      weights <- cbind(
        group = unique_gt_pairs$g[k],
        time = unique_gt_pairs$tt[k],
        e = unique_gt_pairs$tt[k]- unique_gt_pairs$g[k],
        g_prime = g_prime_g_t, 
        t_prime = t_prime_g_t,
        theta = theta,
        w_eff = as.numeric(w))
    }
    # rbind ther weights into weights_out
    if(k==1){
      weights_out = weights
    }
    else{
      weights_out = rbind(weights_out, weights)
    }
  }
  weights_out <- as.data.frame(weights_out)
  
  
  #-----------------------------------------------------------------------------
  return(list(
    estimate = unique_gt_pairs$estimate,
    std.error = unique_gt_pairs$std_error,
    eff_inf_function = eff_inf_function,
    group = unique_gt_pairs$g,
    time = unique_gt_pairs$tt,
    weights = weights_out,
    return_weights = return_weights
  ))
  
}

# Function to compute aggregated event-study

es_efficient <- function(eff_did, 
                         data,
                         yname,
                         tname,
                         gname,
                         idname
){
  
  df <- as.data.frame(data)
  # Reduce to cross-sections
  df <- data[data[,tname]==eff_did$time[1], ]
  
  group <- eff_did$group
  time <- eff_did$time
  att <- eff_did$estimate
  inf_func_attgt <- eff_did$eff_inf_function
  return_weights <- eff_did$return_weights
  weights = NULL
  if(return_weights){
  weights <- as.data.frame(eff_did$weights)
  weights[,"e"] <- weights[,"time"] - weights[,"group"]
  }
  glist<- unique(group)
  # Probability of each group
  pg <- sapply(glist, function(g) mean((df[,gname]==g)))
  
  # length of this is equal to number of groups
  pgg <- pg
  
  # same but length is equal to the number of ATT(g,t)
  pg <- pg[match(group, glist)]
  
  # n x 1 vector of group variable
  G <-  data.frame(df[,gname])
  # Event times
  eseq <- unique(time - group) # compute atts that are specific to each event time
  eseq <- eseq[order(eseq)]

  # ---------------------------------------------------------------------
  # Some aux functions
  wif <- function(keepers, 
                  pg,
                  G,
                  group) {
    # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
    # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here
    
    # effect of estimating weights in the numerator
    if1 <- sapply(keepers, function(k) {
      ( 1*BMisc::TorF(G==group[k]) - pg[k]) /
        sum(pg[keepers])
    })
    # effect of estimating weights in the denominator
    if2 <- base::rowSums( sapply( keepers, function(k) {
      1*BMisc::TorF(G==group[k]) - pg[k]
    })) %*%
      t(pg[keepers]/(sum(pg[keepers])^2))
    
    # return the influence function for the weights
    if1 - if2
  }
  
  get_agg_inf_func <- function(att, 
                               inffunc1, 
                               whichones,
                               weights_agg, 
                               wif=NULL) {
    # enforce weights are in matrix form
    weights_agg <- as.matrix(weights_agg)
    
    # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    thisinffunc <- inffunc1[,whichones]%*%weights_agg
    
    # Incorporate influence function of the weights
    if (!is.null(wif)) {
      thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
    }
    
    # return influence function
    return(thisinffunc)
  }
  #---------------------------------------------------------------------
  # Compute ES
  # compute atts that are specific to each event time
  es_att <- sapply(eseq, function(e) {
    # keep att(g,t) for the right g&t as well as ones that
    # are not trimmed out from balancing the sample
    whiche <- which( (time - group == e))
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })
  
  # Weights for each event time
  weights_event_study = NULL
  if (return_weights) {
  weights_event_study <- sapply(eseq, function(e) {
    
    whiche <- which((weights$e ==e)  )
    which_g <- weights$group[whiche]
    # compute the ES-weights for each group
    pgg_e <- pgg[match(which_g, glist)]/
      sum(pgg[match(unique(which_g), glist)])
    # compute the weights for each group
    es_weights <- pgg_e * weights$w_eff[whiche]
    weights_es <- cbind(
      group = which_g,
      time = weights$time[whiche],
      e = e,
      g_prime = weights$g_prime[whiche], 
      t_prime = weights$t_prime[whiche],
      theta = weights$theta[whiche],
      w_eff = weights$w_eff[whiche],
      w_eff_es = as.numeric(es_weights)
    )
  })
  
  # rbind all list elements in weights_event_study
  weights_event_study <- as.data.frame(do.call(rbind, weights_event_study))
  # Weights for the overall ES
  weights_event_study$w_eff_overall <- weights_event_study$w_eff_es/length(unique(weights_event_study$e))
  }
  
  # compute standard errors for ES
  es_se_inner <- lapply(eseq, function(e) {
    whiche <- which( (time - group == e) )
    pge <- pg[whiche]/(sum(pg[whiche]))
    
    wif_e <- wif(whiche, pg, G, group)
    
    inf_func_e <- as.numeric(get_agg_inf_func(att = att,
                                              inffunc1 = inf_func_attgt,
                                              whichones = whiche,
                                              weights_agg = pge,
                                              wif = wif_e))
    se_e <- sqrt(mean(inf_func_e^2)/length(inf_func_e))
    list(inf_func = inf_func_e, se = se_e)
  })
  
  es_se_e <- unlist(BMisc::getListElement(es_se_inner, "se"))
  es_se_e[es_se_e <= sqrt(.Machine$double.eps)*10] <- NA
  
  es_inf_func_e <- simplify2array(BMisc::getListElement(es_se_inner, "inf_func"))
  
  # get overall average treatment effect
  # by averaging over positive dynamics
  epos <- eseq >= 0
  overall_es <- mean(es_att[epos])
  overall_inf_func <- get_agg_inf_func(att = overall_es[epos],
                                       inffunc1 = as.matrix(es_inf_func_e[,epos]),
                                       whichones = (1:sum(epos)),
                                       weights_agg = (rep(1/sum(epos), sum(epos))),
                                       wif=NULL)
  
  overall_inf_func <- as.numeric(overall_inf_func)
  overall_se <- sqrt(mean(overall_inf_func^2)/length(overall_inf_func))
  if(!is.na(overall_se)){
    if (overall_se <= sqrt(.Machine$double.eps)*10) overall_se <- NA
  }
  
  return(list(es = es_att,
              se = es_se_e,
              inf_func = es_inf_func_e,
              overall_es = overall_es,
              overall_se = overall_se,
              overall_inf_func = overall_inf_func,
              eseq = eseq,
              weights_es = weights_event_study,
              return_weights = return_weights))
  
  
  
}