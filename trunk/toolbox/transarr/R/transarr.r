# transarr -- a simple class for dealing with data from series of microarrays
# (dealing means mostly plotting at this point).

# $Id$

# transarr is subclassed from data.frame. Names of the factors provide
# the row.names. The column names are the "default" ones (V2, V3, ...),
# they could be taken from the input file, though (FIXME)
# transarr may have two attributes, expspec ("experiment specification")
# and ggroups ("gene groups"). Both are lists that contain ragged arrays.
# The rows of ggroups contain character vectors, each representing a
# list of names of genes forming a group. The rows of expspec contain
# integer vectors, each representing a set of columns belongint to
# an experiment (e.g. a time series).

# File formats:
# The microarray data themselves are basically loaded from a
# tab-separated values file via read.table(). The ggroups attribute
# and the expspec attribute are loaded from tab-separated files in
# which each line contains the name of the gene group (or experiment
# spec), followed by the gene names (or column indices).

# plot method:
# The plot method generates a "expspec x ggroups" series of plots, i.e.
# for each experiment, a series of plot with one plot per gene group
# is generated. Plots for one experiments are grouped on one page, this
# be overridden by providing an mfrow argument.

# Future plans:
# * additional parameters for the plot function for selecting experiments
#   and gene groups
# * an interface for pulling out a data.frame specified by a list of
#   experiments and gene groups
# * interface for computing clusterings
# * ... and for rendering clusterings in parenthesized format
# * facility to generate bootstrap experiment specs and gene groups

# $Log$
# Revision 1.1  2005/03/08 17:12:02  jtk
# Initial revision
#
# Revision 1.3  2003/03/14 01:20:18  kim
# added stuff for Pearson-like distance accroding to Eisen et.al.
#
# Revision 1.2  2003/02/28 00:46:46  kim
# added proper xlim and ylim handling to transarr plot method
#
# Revision 1.1  2003/02/17 21:40:31  kim
# Transformed microarr.r into a R library called transarr. Added array_background
#     to regstruct system
#
# Revision 1.1  2003/02/14 16:31:27  kim
# regstruct project is nearing completion. regstruct_transsys for
#     generation of transsys programs and regstruct_transarray for
#     array simulation are (more or less) properly implemented.
#     Analysis with R has been assembled in microarr.r
#


# library(mva)


# compute parenthesized string representation of cluster tree
# ("New Hamapshire" format like) for a cluster structure as produced
# by the hclust method provided by the mva library.
# This is the default method for now, since I don't know any other
# cluster stuff yet...

clusterstring.default <- function(hc, nodesep = "")
{

  nodestring <- function(s, d)
  {
    paste(s, ":", d, sep = "");
  }

  clist <- list();
  for (i in 1:nrow(hc$merge))
  {
    j <- hc$merge[i, 1];
    if (j < 0)
    {
      c0 <- hc$labels[-j];
      attr(c0, "height") <- 0;
    }
    else
    {
      c0 <- clist[[j]];
    }
    j <- hc$merge[i, 2];
    if (j < 0)
    {
      c1 <- hc$labels[-j];
      attr(c1, "height") <- 0;
    }
    else
    {
      c1 <- clist[[j]];
    }
    n0 <- nodestring(c0, hc$height[i] - attr(c0, "height"));
    n1 <- nodestring(c1, hc$height[i] - attr(c1, "height"));
    cnew <- paste("(", n0, ",", nodesep, n1, ")", sep = "");
    attr(cnew, "height") <- hc$height[i];
    clist[[i]] <- cnew;
  }
  topnode <- clist[[length(clist)]];
  paste(topnode, ":0;", sep = "");
}


# The generic function

clusterstring <- function(clust, ...)
{
  UseMethod("clusterstring");
}


# recursive version of parenthesized cluster string -- not recommended,
# since R seems to float on a shallow stacking sea...

hcluststringrek <- function(hc)
{

  hcstring <- function(hcnode)
  {
    if (is.null(attr(hcnode, "descendants")))
    {
      s <- paste(hcnode, ":", attr(hcnode, "distance"), sep = "");
    }
    else
    {
      c0 <- attr(hcnode, "descendants")[[1]];
      c1 <- attr(hcnode, "descendants")[[2]];
      s <- paste("(", hcstring(c0), ",",  hcstring(c1), "):", attr(hcnode, "distance"), sep = "");
    }
    s;
  }

  leaflist <- list();
  for (i in 1:length(hc$labels))
  {
    leaf <- hc$labels[i];
    class(leaf) <- "clustnode";
    attr(leaf, "height") <- 0.0;
    attr(leaf, "distance") <- 0.0;
    leaflist[[i]] <- leaf;
  }
  clist <- list();
  for (i in 1:nrow(hc$merge))
  {
    j <- hc$merge[i, 1];
    if (j < 0)
    {
      c0 <- leaflist[[-j]];
    }
    else
    {
      c0 <- clist[[j]];
    }
    j <- hc$merge[i, 2];
    if (j < 0)
    {
      c1 <- leaflist[[-j]];
    }
    else
    {
      c1 <- clist[[j]];
    }
    cnew <- paste("bletch", i);
    class(cnew) <- "clustnode";
    attr(cnew, "height") <- hc$height[i];
    attr(c0, "distance") <- attr(cnew, "height") - attr(c0, "height");
    attr(c1, "distance") <- attr(cnew, "height") - attr(c1, "height");
    attr(cnew, "distance") <- 0.0;
    desc <- list();
    desc[[1]] <- c0;
    desc[[2]] <- c1;
    attr(cnew, "descendants") <- desc;
    clist[[i]] <- cnew;
  }
  topnode <- clist[[length(clist)]];
  hcstring(topnode);
}


# Pearson-like distance according to Eisen et.al. "Cluster Analysis and
# Display of Genome-Wide Expression Patterns", PNAS 95: 14863-14868 (1998).

dist.eisen <- function(arr, offset = 0.0)
{
  phi <- function(expvec, offset)
  {
    v <- expvec - offset;
    s <- sum(v * v);
    sqrt(s / length(v));
  }

  if (length(offset) == 1)
  {
    offset <- double(nrow(arr)) + offset;
  }
  if (length(offset) != nrow(arr))
  {
    stop("length of offset does not match number of genes");
  }
  rl <- list();
  for (r in 1:nrow(arr))
  {
    rl[[row.names(arr)[r]]] <- as.double(arr[r, 1:ncol(arr)]);
  }
  ph <- c()
  for (i in 1:nrow(arr))
  {
    ph[i] <- phi(rl[[i]], offset[i]);
  }
  dm <- matrix(nrow = nrow(arr), ncol = nrow(arr));
  for (i in 1:nrow(arr))
  {
    for (j in 1:nrow(arr))
    {
      if ((ph[i] > 0.0) && (ph[j] > 0.0))
      {
	x <- (rl[[i]] - offset[i]) / ph[i];
	y <- (rl[[j]] - offset[j]) / ph[j];
	dm[i,j] <- sum(x * y) / ncol(arr);
      }
      else
      {
	dm[i,j] <- 0.0;
      }
    }
  }
  dm <- as.dist(1.0 - dm, upper = TRUE);
  names(dm) <- row.names(arr);
  dm;
}


cluster.transarr <- function(m, distfunc, clustfunc, expnames = NULL, ggnames = NULL)
{
  expname <- function(en)
  {
    ifelse(is.null(en), "all_experiments", en);
  }

  ggname <- function(gn)
  {
    ifelse(is.null(gn), "all_ggroups", gn);
  }

  if (is.null(expnames))
  {
    expnames <- list(NULL);
  }
  if (is.null(ggnames))
  {
    ggnames <- list(NULL);
  }
  clustermatrix <- list();
  for (en in expnames)
  {
    clustercolumn <- list();
    for (gn in ggnames)
    {
      # print(paste("clustering experiment '", en, "', gene group '", gn, "'", sep = ""));
      me <- transarrExtract(m, en, gn);
      # print(me);
      dmatrix <- distfunc(me);
      cl <- clustfunc(dmatrix);
      class(cl) <- append("cluster.transarr", class(cl));
      attr(cl, "expname") <- expname(en);
      attr(cl, "ggname") <- ggname(gn);
      clustercolumn[[ggname(gn)]] <- cl;
    }
    clustermatrix[[expname(en)]] <- clustercolumn;
  }
  class(clustermatrix) <- c("clustermatrix", class(clustermatrix));
  clustermatrix;
}


clusterstring.clustermatrix <- function(cm, ...)
{
  cstring <- "";
  for (en in names(cm))
  {
    for (gn in names(cm[[en]]))
    {
      cstring <- paste(cstring, clusterstring(cm[[en]][[gn]]), "\n", sep = "");
    }
  }
  cstring;
}


cluster.transarr.allexp <- function(m, distfunc, clustfunc, ggnames = NULL)
{
  cluster.transarr(m, distfunc, clustfunc, names(attr(m, "expspec")), ggnames);
}


plot.cluster.transarr <- function(x, ...)
{
  NextMethod(plot, main = toString(c(attr(x, "expname"), attr(x, "ggname"))), ...);
}


plot.clustermatrix <- function(x, ...)
{
  for (en in names(x))
  {
    for (gn in names(x[[en]]))
    {
      plot(x[[en]][[gn]]);
      if (dev.interactive())
      {
	invisible(readline("Hit return"));
      }
    }
  }
}


# read a gene groups file in which each gene group is described by a
# line of the format
#
#     groupname ( "\t" genename) +

read.genegroups <- function(file = "")
{
  if (is.character(file))
  {
    if (file == "")
    {
      file <- stdin();
    }
    else
    {
      file <- file(file, "r");
      open(file);
      on.exit(close(file));
    }
  }
  else
  {
    if (!isOpen(file))
    {
      open(file, "r");
      on.exit(close(file));
    }
  }
  genegroups <- list();
  groupname <- scan(file, "", n = 1, sep = "\t", quiet = TRUE);
  while (length(groupname) > 0)
  {
    genes <- scan(file, what = "", nlines = 1, sep = "\t", quiet = TRUE);
    genegroups[[groupname]] <- genes;
    groupname <- scan(file, "", n = 1, sep = "\t", quiet = TRUE);
  }
  class(genegroups) <- c("ggroups", class(genegroups));
  genegroups;
}


# read a file of experiment specification, each experiment is described
# by a line with the structure
#
#     experimentname ( "\t" columnname ) +

read.experimentspecs <- function(file = "")
{
  if (is.character(file))
  {
    if (file == "")
    {
      file <- stdin();
    }
    else
    {
      file <- file(file, "r");
      open(file);
      on.exit(close(file));
    }
  }
  else
  {
    if (!isOpen(file))
    {
      open(file, "r");
      on.exit(close(file));
    }
  }
  expspec <- list();
  expspecname <- scan(file, "", n = 1, sep = "\t", quiet = TRUE);
  while (length(expspecname) > 0)
  {
    expspeccols <- scan(file, what = "", nlines = 1, sep = "\t", quiet = TRUE);
    expspec[[expspecname]] <- expspeccols;
    expspecname <- scan(file, "", n = 1, sep = "\t", quiet = TRUE);
  }
  class(expspec) <- c("expspec", class(expspec));
  expspec;
}


# Read a microarray file, optionally accompanied by experimentspecs in
# espfile and gene groups in ggfile. The microarray data itself is
# in a tab-separated file suitable for read.table().
#
# FIXME: If the optional experiment specs and gene groups are also
# read, this function should check that all columns specified by the
# experiments and all genes in the gene groups are indeed present in
# the microarray data.

read.transarr <- function(file = "", espfile = NULL, ggfile = NULL, ...)
{
  t <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, ...);
  for (r in 1:nrow(t))
  {
    row.names(t)[r] <- as.character(t[r,1]);
  }
  t <- t[1:nrow(t), 2:ncol(t)];
  class(t) <- append("transarr", class(t));
  if (is.null(ggfile))
  {
    ggroups <- list();
    ggroups[["all_genes"]] <- row.names(t);
    attr(t, "ggroups") <- ggroups;
  }
  else
  {
    attr(t, "ggroups") <- read.genegroups(ggfile);
  }
  if (is.null(espfile))
  {
    expspec <- list();
    expspec[["all_arrays"]] <- 1:ncol(t);
    attr(t, "expspec") <- expspec;
  }
  else
  {
    attr(t, "expspec") <- read.experimentspecs(espfile);
  }
  t;
}


# Extract a subset of experiments and gene groups from a transarr
# object. By default, all experiments and all gene groups are extracted.
# Thus, without additional arguments, a copy of the transarr object
# is generated.
# Note that the order of the experiments (columns) in the microarray
# data frame will depend on the order with which they are given in the
# experiment specs, and likewise, the row order depends on the order in
# the gene groups. Thus, if such order is informative, specify it by an
# experiment spec or gene group, respectively.
# Extraction should work properly with non-exclusive experiments or
# gene group (e.g. with some genes being members of several gene groups).

transarrExtract <- function(m, expnames = NULL, ggnames = NULL)
{
  # print(paste("expnames:", expnames));
  # print(paste("ggnames:", ggnames));
  if (is.null(expnames))
  {
    expnames <- names(attr(m, "expspec"));
  }
  if (is.null(ggnames))
  {
    ggnames <- names(attr(m, "ggroups"));
  }
  # print(paste("expnames:", expnames));
  # print(paste("ggnames:", ggnames));
  fnames <- c();
  for (gn in ggnames)
  {
    fnames <- union(fnames, attr(m, "ggroups")[[gn]])
  }
  fnames <- intersect(fnames, row.names(m));
  colvec <- c();
  for (en in expnames)
  {
    colvec <- union(colvec, attr(m, "expspec")[[en]]);
  }
  colvec <- intersect(colvec, names(m));
  # print(paste("colvec:", colvec));
  # print(paste("fnames:", fnames));
  mnew <- m[fnames, colvec]
  attr(mnew, "expspec") <- list();
  for (en in expnames)
  {
    attr(mnew, "expspec")[[en]] <- attr(m, "expspec")[[en]];
  }
  attr(mnew, "ggroups") <- list();
  for (gn in ggnames)
  {
    attr(mnew, "ggroups")[[gn]] <- attr(m, "ggroups")[[gn]];
  }
  mnew;
}


# function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL,
# log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
# ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL,
# panel.last = NULL, col = par("col"), bg = NA, pch = par("pch"),
# cex = 1, lty = par("lty"), lab = par("lab"), lwd = par("lwd"),
# asp = NA, ...)

# Plot method for transarr object

plot.transarr <- function(x, mfrow = NULL, type = NULL, xlim = NULL, ylim = NULL, axes = TRUE, frame.plot = axes, ann = par("ann"), ...)
{
  expspec <- attr(x, "expspec");
  ggroups <- attr(x, "ggroups");
  if (is.null(type))
  {
    type <- "l";
  }
  if (is.null(ylim))
  {
    ylim <- range(x);
  }
  if (is.null(mfrow))
  {
    mfrow <- c(length(ggroups), 1);
    mar <- c(2, 4, 2, 2);
  }
  else
  {
    mfrow <- par()$mfrow;
    mar <- par()$mar;
  }
  opar <- par(mfrow = mfrow, mar = mar);
  on.exit(par(opar));
  mfg <- par()$mfg;
  plotnum <- 0;
  plotsperpage <- mfrow[1] * mfrow[2];
  for (expspecname in names(expspec))
  {
    for (ggname in names(ggroups))
    {
      if (is.null(xlim))
      {
	xlim1 <- c(1, ncol(x[expspec[[expspecname]]]));
      }
      else
      {
        xlim1 <- xlim;
      }
      plot.new();
      plot.window(xlim1, ylim, ...);
      xv <- 1:length(expspec[[expspecname]]);
      for (r in ggroups[[ggname]])
      {
	yv <- c();
	i <- 1;
	for (c in expspec[[expspecname]])
	{
	  yv[i] <- x[r, c];
	  i <- i + 1;
	}
	plot.xy(xy.coords(xv, yv), type = type, ...)
      }
      if (axes)
      {
	axis(1, ...);
	axis(2, ...);
      }
      if (frame.plot)
      {
	box(...);
      }
      if (ann)
      {
	title(main = toString(c(expspecname, ggname)), ...)
      }
      plotnum <- plotnum + 1;
      if (dev.interactive())
      {
	if ((plotnum %% plotsperpage) == 0)
	{
	  invisible(readline("Hit return"));
	}
      }
    }
  }
  invisible();
}

