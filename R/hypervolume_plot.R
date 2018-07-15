do_outline_alpha <- function(rp, alpha)
{
  ah = alphahull::ashape(rp,alpha=alpha)
  return(ah)
}


#do_outline_concave <- function(rp, concavity)
#{
#  cm = concaveman::concaveman(rp,concavity=concavity, length_threshold=0)
#  return(cm)
#}

do_outline_ball <- function(rp, radius)
{
  gb = rgeos::gBuffer(sp::SpatialPoints(rp), quadsegs=2, width=radius)
  return(gb)
}

do_outline_raster <- function(pts,res)
{
  pts <- as.matrix(pts)
  
  pr <- padded_range(pts,multiply.interval.amount=0.25)
  
  e <- extent(t(pr))
  
  r <- raster::raster(e, ncol=res, nrow=res)
  
  x <- raster::rasterize(pts, r, rep(1, nrow(pts)), fun=mean,background=NA)
  
  w <- raster::rasterToPolygons(x,dissolve=TRUE)
  
  return(w)	
}

plot.Hypervolume <- function(x, ...)
{
  templist = new("HypervolumeList")
  templist@HVList=list(x)
  plot.HypervolumeList(templist, ...)
}

extendrange <- function(x,factor=0.5)
{
  xmin <- min(x,na.rm=TRUE)
  xmax <- max(x,na.rm=TRUE)
  
  xminf <- xmin - (xmax - xmin)*factor
  xmaxf <- xmax + (xmax - xmin)*factor
  
  result <- c(xminf, xmaxf)
  
  return(result)
}

plot.HypervolumeList <- function(x, 
                                 show.3d=FALSE,plot.3d.axes.id=NULL,
                                 show.axes=TRUE, show.frame=TRUE,
                                 show.random=TRUE, show.density=TRUE,show.data=TRUE,
                                 names=NULL, show.legend=TRUE, limits=NULL, 
                                 show.contour=TRUE, contour.lwd=1.5, 
                                   contour.type='kde', 
                                   contour.alphahull.alpha=0.25,
                                   contour.ball.radius.factor=1, 
                                   contour.kde.level=0.01,
                                   contour.raster.resolution=100,
                                 show.centroid=TRUE, cex.centroid=2,
                                 colors=rainbow(floor(length(x@HVList)*1.5),alpha=0.8), 
                                 point.alpha.min=0.2, point.dark.factor=0.5,
                                 cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
                                 num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
                                 plot.function.additional=NULL,
                                 verbose=FALSE,
                                 ...)
{
  sapply(x@HVList, function(z)
  {
    if (verbose==TRUE)
    {
      cat(sprintf("Showing %d random points of %d for %s\n",min(nrow(z@RandomPoints), num.points.max.random), nrow(z@RandomPoints), z@Name))
    }
    if (show.data && length(z@Data) > 0)
    {
      npd <- ifelse(all(is.nan(z@Data)), 0, nrow(z@Data))
      if (verbose==TRUE)
      {
        cat(sprintf("Showing %d data points of %d for %s\n",min(num.points.max.data, npd), npd, z@Name))
      }
    }    
    
  })
  
  if (!requireNamespace("alphahull", quietly = TRUE)) {
    warning("The package 'alphahull' is needed for contour plotting with contour.type='alphahull'. Please install it to continue.\n\n *** Temporarily setting contour.type='kde'.", call. = FALSE)
    contour.type <- 'kde'
  }
  
  alldims = sapply(x@HVList, function(z) { z@Dimensionality })
  allnames = sapply(x@HVList, function(z) { z@Name })
  stopifnot(all(alldims[1] == alldims))
  
  all <- NULL
  alldata <- NULL
  for (i in 1:length(x@HVList))
  {
    ivals = sample(nrow(x@HVList[[i]]@RandomPoints), min(c(num.points.max.random, nrow(x@HVList[[i]]@RandomPoints))))
    subsampledpoints = data.frame(x@HVList[[i]]@RandomPoints[ivals,,drop=FALSE])
    densityvals = x@HVList[[i]]@ValueAtRandomPoints[ivals]
    
    if (nrow(subsampledpoints) > 0)
    {  
      subsampledpoints = cbind(subsampledpoints, ID=rep(i, nrow(subsampledpoints)), Density=(densityvals-min(densityvals,na.rm=TRUE))/(max(densityvals,na.rm=TRUE)-min(densityvals,na.rm=TRUE)))
      subsampledpoints[is.nan(subsampledpoints[,"Density"]),"Density"] <- 1
      all <- rbind(all, subsampledpoints)
    }
    
    thisdata=x@HVList[[i]]@Data
    alldata <- rbind(alldata, cbind(thisdata, ID=rep(i,nrow(thisdata))))
  }  

  alldata <- as.data.frame(alldata)
  if (num.points.max.data < nrow(alldata) && !is.null(num.points.max.data))
  {
  	alldata <- alldata[sample(nrow(alldata), min(c(num.points.max.data, nrow(alldata)))),]
  }
  
  if (is.null(all))
  {
    warning('No random points to plot.')
    if (is.null(dimnames(x@HVList[[1]]@RandomPoints)[[2]]))
    {
      all <- matrix(0,ncol=2+alldims,nrow=1,dimnames=list(NULL,c(paste("X",1:alldims,sep=""),"ID","Density")))
    }
    else
    {
      all <- matrix(0,ncol=2+alldims,nrow=1,dimnames=list(NULL,c(dimnames(x@HVList[[1]]@RandomPoints)[[2]],"ID","Density")))
    }
    all <- as.data.frame(all)
  }
  
  if (reshuffle==TRUE)
  {
    all <- all[sample(nrow(all),replace=FALSE),,drop=FALSE] # reorder to shuffle colors
    alldata <- alldata[sample(nrow(alldata),replace=FALSE),,drop=FALSE]
  }
  
  no_names_supplied = FALSE
  
  if (is.null(names))
  {
    dn = dimnames(all)[[2]]
    names = dn[1:(ncol(all)-2)]
    
    no_names_supplied = TRUE
  }  
  
  if (!is.null(limits) & !is.list(limits))
  {
    varlimlist = vector('list',ncol(all)-2)
    for (i in 1:length(varlimlist))
    {
      varlimlist[[i]] <- limits
    }
    limits = varlimlist
  }
  
  colorlist <- colors[all$ID]
  alphavals <- (all$Density - quantile(all$Density, 0.025, na.rm=T)) / (quantile(all$Density, 0.975, na.rm=T) - quantile(all$Density,0.025, na.rm=T))
  alphavals[is.nan(alphavals)] <- 0.5 # in case the quantile is un-informative
  alphavals[alphavals < 0] <- 0
  alphavals[alphavals > 1] <- 1
  alphavals <- point.alpha.min + (1 - point.alpha.min)*alphavals
  
  if (show.density==FALSE)
  {
    alphavals <- rep(1, length(colorlist))
  }
  
  for (i in 1:length(colorlist))
  {
    colorlist[i] <- rgb_2_rgba(colorlist[i], alphavals[i])
  }
  
  colorlistdata = colors[alldata$ID]
  for (i in 1:length(colorlistdata))
  {
    colorlistdata[i] <- rgb_2_set_hsv(colorlistdata[i], v=1-point.dark.factor)
  }
  
  if (ncol(all) < 2)
  {
    stop('Plotting only available in n>=2 dimensions.')
  }
  
  if (show.3d==FALSE)
  {
    op = par(no.readonly = T)
    
    par(mfrow=c(ncol(all)-2, ncol(all)-2))
    par(mar=c(0,0,0,0))
    par(oma=c(0.5,0.5,0.5,0.5))
    
    for (i in 1:(ncol(all)-2))
    {
      for (j in 1:(ncol(all)-2))  
      {
        if (j > i)
        {
          # set up axes with right limits
          plot(all[,j], all[,i],type="n",axes=FALSE,xlim=limits[[j]], ylim=limits[[i]],bty='n')
          
          # draw random points
          if(show.random==TRUE)
          {
            points(all[,j], all[,i], col=colorlist,cex=cex.random,pch=16)
          }
          
          # show data
          if (show.data & nrow(alldata) > 0)
          {
            points(alldata[,j], alldata[,i], col=colorlistdata,cex=cex.data,pch=16)
          }
          
          if (show.centroid == TRUE)
          {
            for (whichid in 1:length(unique(all$ID)))
            {
              allss <- subset(all, all$ID==whichid)
              centroid_x <- mean(allss[,j],na.rm=TRUE) 
              centroid_y <- mean(allss[,i],na.rm=TRUE)
              
              # draw point
              points(centroid_x, centroid_y, col=colors[whichid],cex=cex.centroid,pch=16)
              # add a white boundary for clarity
              points(centroid_x, centroid_y, col='white',cex=cex.centroid,pch=1,lwd=1.5)
            }
          }
          
          
          # calculate contours
          if (show.contour==TRUE)
          {
            # draw shaded centers
            for (whichid in 1:length(unique(all$ID)))
            {
              allss <- subset(all, all$ID==whichid) # remove oversampling of some points
              
              if (nrow(allss) > 0)
              {     
                contourx <- allss[,j]
                contoury <- allss[,i]
                
                rp = cbind(contourx, contoury)
                vol_this = x@HVList[[whichid]]@Volume
                density_this = nrow(rp) / vol_this
                dim_this = x@HVList[[whichid]]@Dimensionality
                radius_critical <- density_this^(-1/dim_this) * contour.ball.radius.factor
                
                if (contour.type=='alphahull')
                {
                  poly_outline = do_outline_alpha(rp=rp, alpha=contour.alphahull.alpha)
                  plot(poly_outline,add=TRUE,wpoints=FALSE,wlines='none',lwd=contour.lwd,col=colors[whichid])
                }
                else if (contour.type=='ball')
                {
                  poly_outline <- do_outline_ball(rp=rp, radius=radius_critical)
                  
                  sp::plot(poly_outline, add=TRUE,lwd=contour.lwd,col=colors[whichid])
                }
                else if (contour.type=='kde')
                {
                  m_kde = kde2d(rp[,1], rp[,2], n=50, h=radius_critical)
                  contour(m_kde, add=TRUE, levels=contour.kde.level,drawlabels=FALSE,lwd=contour.lwd,col=colors[whichid])
                }
                else if (contour.type=='raster')
                {
                  poly_raster <- do_outline_raster(as.matrix(rp),res=contour.raster.resolution)
                  sp::plot(poly_raster, add=TRUE, lwd=contour.lwd,col=colors[whichid])
                }
              }
            }
          }
          
          if (!is.null(plot.function.additional))
          {
            plot.function.additional(j,i)
          }
          
          if (show.frame==TRUE)
          {
            box()
          }
        }
        else if (j == i)
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=FALSE)
          text(0.5, 0.5, names[j],cex=cex.names)
        }
        else if (j==1 & i == (ncol(all) - 2))
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=FALSE)
          
          if (show.legend == TRUE)
          {
            legend('topleft',legend=allnames,text.col=colors,bty='n',cex=cex.legend)
          }
        }
        else
        {
          plot(0,0,type="n",axes=FALSE)    
        }
        
        if (j==i+1)
        {
          if (show.axes==TRUE)
          {
            axis(side=1,cex.axis=cex.axis)
            axis(side=2,cex.axis=cex.axis)
          }
        }
      }
    }  
    par(op)
  }
  else
  {
    if (is.null(plot.3d.axes.id))
    {
      plot.3d.axes.id=1:3  
    }
    
    if (no_names_supplied==TRUE)
    {
      axesnames <- names[plot.3d.axes.id]
    }
    else
    {
      axesnames <- names
    }
    
    if(length(plot.3d.axes.id)!=3) { stop('Must specify three axes') }

    if (show.density==TRUE)
    {
      for (i in 1:length(colorlist))
      {
        colorlist[i] <- rgb_2_set_hsv(colorlist[i], s=(alphavals[i]^2)) # should do this with alpha, but a workaround for non-transparent OpenGL implementations...
      }
    }
    rgl::plot3d(all[,plot.3d.axes.id],col=colorlist,expand=1.05, xlab=axesnames[1], ylab=axesnames[2], zlab=axesnames[3], xlim=limits[[1]],ylim=limits[[2]],zlim=limits[[3]],size=cex.random,type='p',box=show.frame,axes=show.axes)
    
    if (show.legend==TRUE)
    {
      for (i in 1:length(allnames))
      {
        rgl::mtext3d(allnames[i],edge='x-+',line=1+i*cex.legend*1.25,color=colors[i],cex=cex.legend)  
      }
    }
    
    if (show.data)
    {
      if (!any(is.nan(as.matrix(alldata[,plot.3d.axes.id]))))
      {
        rgl::points3d(x=alldata[,plot.3d.axes.id[1]], y=alldata[,plot.3d.axes.id[2]], z=alldata[,plot.3d.axes.id[3]], col=colorlistdata,cex=cex.data,pch=16)
      }
    }
    
    if (show.centroid == TRUE)
    {
      for (whichid in 1:length(unique(all$ID)))
      {
        allss <- subset(all, all$ID==whichid)
        centroid_1 <- mean(allss[,plot.3d.axes.id[1]],na.rm=TRUE)
        centroid_2 <- mean(allss[,plot.3d.axes.id[2]],na.rm=TRUE)
        centroid_3 <- mean(allss[,plot.3d.axes.id[3]],na.rm=TRUE)
        
        # draw point
        rgl::points3d(x=centroid_1, y=centroid_2, z=centroid_3, col=colors[whichid],cex=cex.centroid,pch=16)
      }
    }
  }
}  
