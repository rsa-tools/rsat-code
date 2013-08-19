################################################################
##
## General utilities, to be used by other R files.
##
##  source(file.path(dir.util, "util.R"))

################################################################
# Treatment of verbosity
#
verbosity <- 1
verbose <- function(x, level=0)  {
  if (level <= verbosity) {
    print (paste(date(), " - ", x))
  }
}

################################################################
#
# Export the current device to a set of files with different
# formats.
# Example: 
#     export.plot(file.prefix="test",
#                 export.formats=c("postscript","jpg","png","bmp","pdf"))
export.plot <- function (file.prefix="PlotExport",
	                 export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
			 width=11, # in inches
			 height=8, # in inches
			 horizontal=TRUE,
                         ... ## Additional parameters are passed to the export method
                         ) {

   for (f in export.formats) {
     from.dev <- dev.cur();
     open.plot.device(file.prefix,format=f, width=width,height=height,horizontal=horizontal, ...)
     to.dev <- dev.cur()
     dev.set(which=from.dev)
     dev.copy(which=to.dev)
     dev.set(which=to.dev)
     dev.off()
     dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
   }
}

open.plot.device <- function (file.prefix,
                              format='pdf',
                              width=8,
                              height=8,
                              horizontal=TRUE,
                              ppi=72, ## Points per Inch (screen resolution by default)
                              ...) {  

  format <- tolower(format)  
  file.ext <- c(
                x11 = "x11",
                postscript = "ps",
                pdf = "pdf",
                ps = "ps",
                eps = "eps",
                jpeg="jpg",
                jpg="jpg",
                bmp="bmp",
                png="png")
  file.name <- paste(file.prefix,file.ext[format], sep=".")
  verbose(paste("Exporting plot to file ",file.name, " ", f, " format", sep=""), 2)
  if ((format == "postscript") || (format == "ps")) {
    postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
  } else if (format == "eps") {
    postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
  } else if (format == "pdf") {
    pdf(file.name, paper="special",width=width,height=height, ...)
  } else if ((format == "jpg") || (format == "jpeg")) {
    jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
  } else if (format == "png") {
    png(file.name,width=width*ppi,height=height*ppi, ...)
  } else if (format == "bmp") {
    bitmap(file.name,width=width*ppi,height=height*ppi, ...)
  } else if (format == "x11") {
    x11(width=width,height=height)
  } else {
    print(paste("Error: format ", f, " is not supported", sep=""))
    return()
  }
}

################################################################
#
# Export an objet in dofferent formats
#
export.object <- function (x,
			   file.prefix="export", 
			   export.formats=c("save"), # supported: "save", "print", "table", "dput"
			   ... ## Additional argumens are passed to the exporting function
			   ) {
  file.ext <- c(
	      save="Rdata",
	      table="tab",
	      dput="R",
	      print="txt"
	      )
  if ((is.na(export.formats)) || (is.na(export.formats))) {
    return
  }

  for (f in export.formats) {
    file.name <- paste(file.prefix,file.ext[f],sep=".")
    verbose(paste("Exporting object to file", file.name),2)
    if (f=="save") {
      names(x) 
      ## Tricky way to circumvent a problem : I don't know why, but
      ## when I do nothing with x before saving, the resulting file
      ## does not contain x
      save(x, file=file.name, ...)

    } else if (f == "print") {
      sink(file.name, ...)
      print(x)
      sink()

    } else if (f == "dput") {
      dput(x,file=file.name, ...)

    } else if (f == "table") {

      ## Tricky way to circumvent a bug with header in write.table
      ## (headers are shifted one cell left when the row.names are
      ## exported)
      x.export <- cbind(row.names(x), x)
      names(x.export)[1] <- "row.ID"
      write.table(x.export,file=file.name,quote=F,sep='\t', row.names=F, ...)
      rm(x.export)
      
      ##      x.fixed <- cbind(row.names(x),x)
      ##      names(x.fixed) <- c("ID",names(x))
      ##      write.table(x.fixed,file=file.name,quote=F,sep='\t', row.names=F, ...)

    } else {
      print (paste("ERROR: format", f, "is not supported by export.object()"))
      return()

    }
  }
}



## ##############################################################
## automatic creation of output directories
## Check if the directory exists, warn when creating it
dir.create.check <- function(x ## vector with the dicetories to be created
			     ) {
  for (dir in x) {
    if (!file.exists(dir)) {
      print (paste("Creating directory", dir))
      dir.create(dir)
      if (!file.exists(dir)) {
        stop(paste("Cannot create directory", dir))
      }
    }
  }
}
