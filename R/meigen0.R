meigen0		<- function(meig, coords0, coords_z0 = NULL, s_id0 = NULL) {
  if(meig$other$model == "other") {
    stop("meigen0 is not supported for user-specified spatial proximity matrix")
  }

  nz          <- 0
  exception   <- FALSE
  if( !is.null(coords_z0) ){
    if( !is.null(meig$other$z_use[1]) ){
      if( sum(meig$other$z_use)==0 ){
        coords_z0 <- NULL
        exception <- TRUE
      } else {
        coords_z0 <- as.matrix(as.matrix(coords_z0)[,meig$other$z_use])
      }
    } else {
      coords_z0   <- as.matrix(coords_z0)
    }

    if( !is.null(coords_z0) ){
      nz        <- dim(coords_z0)[2]
      if( nz != length(meig$ev_z) ){
        stop("Dimension mismatch: coords_z0 (meigen0) and coords_z (meigen)")
      }
    }
  }
  if( !is.null(meig$ev_z)&is.null(coords_z0) & !exception ){
    stop("coords_z0 is missing")
  }


  if( !is.null(s_id0) ){
      coords0_0 <- coords0
      coords0_x <- tapply(coords0[, 1], s_id0, mean)
      coords0_y <- tapply(coords0[, 2], s_id0, mean)
      coords0 <- as.matrix(cbind(coords0_x, coords0_y))
      s_id2 <- data.frame(s_id = rownames(coords0_x), s_id_num = 1:length(coords0_x))
      s_id_dat <- data.frame(s_id = s_id0, id = 1:length(s_id0))
      s_id_dat2 <- merge(s_id_dat, s_id2, by = "s_id",
                         all.x = TRUE)
      s_id_dat2 <- s_id_dat2[order(s_id_dat2$id), c("s_id", "s_id_num")]
  }
  if(meig$other$fast == 0) {
    if (!is.null(meig$other$s_id)[1]) {
      sfk <- meig$other$sfk
      evk <- meig$ev
      coordk <- meig$other$coordk
      Cmean <- meig$other$Cmeank
    } else {
      sfk <- meig$sf
      evk <- meig$ev
      coordk <- meig$other$coords
      Cmean <- meig$other$Cmean
    }
  } else {
    sfk <- meig$other$sfk
    evk <- meig$other$evk
    coordk <- meig$other$coordk
    Cmean <- meig$other$Cmean
  }
  h <- meig$other$h
  sfk <- sfk %*% diag(1/(evk + 1))
  Dk <- rdist(coordk, coords0)
  if (meig$other$model == "exp") {
    sfk <- t(exp(-Dk/h) - Cmean) %*% sfk
  } else if (meig$other$model == "gau") {
    sfk <- t(exp(-(Dk/h)^2) - Cmean) %*% sfk
  } else if (meig$other$model == "sph") {
    sfk <- t(ifelse(Dk < h, 1 - 1.5 * (Dk/h) + 0.5 * (Dk/h)^3, 0) - Cmean) %*% sfk
  }
  if( !is.null(s_id0) ){#if (!is.null(meig$other$s_id)[1]) {
    sfk <- sfk[s_id_dat2$s_id_num, ]
  }

  sfk_z <- list(NULL)
  if( nz > 0 ){
    for(zz in 1:nz){
      id_z    <- meig$other$id_z[[zz]]
      sf_z    <- as.matrix( aggregate(meig$sf_z[[zz]],by=list(id_z),function(x) head(x,1))[,-1] )
      coordk_z<- meig$other$coordk_z[[zz]]#aggregate( meig$other$coordk_z[[zz]],by=list(id_z),function(x) head(x,1))[,-1]
      h_z     <- meig$other$h_z[[zz]]
      ev_z    <- meig$ev_z[[zz]]
      Cmean_z <- meig$other$Cmean_z[[zz]]

      sf_z    <- as.matrix( sf_z %*% diag(1/(ev_z+1)) )
      Dk_z    <- rdist(coordk_z,coords_z0[,zz])
      if (meig$other$model == "exp") {
        sfk_z[[zz]]<- t(exp(-Dk_z/h_z) - Cmean_z) %*% sf_z
      } else if (meig$other$model == "gau") {
        sfk_z[[zz]]<- t(exp(-(Dk_z/h_z)^2) - Cmean_z) %*% sf_z
      } else if (meig$other$model == "sph") {
        sfk_z[[zz]]<- t(ifelse(Dk_z < h_z, 1 - 1.5 * (Dk_z/h_z) + 0.5 * (Dk_z/h_z)^3, 0)- Cmean_z ) %*% sf_z#
      }
    }
  }

  other  <-list(coords0=coords0)
  return(list(sf = sfk, ev = meig$ev, sf_z = sfk_z, ev_z = meig$ev_z, other= other))
}
