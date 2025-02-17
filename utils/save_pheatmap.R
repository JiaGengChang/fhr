### slightly adapted from here to save png or pdf and specify resolution etc: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file

# function detects '.png' or '.pdf' in the declared filename and assigns that file type

### pheatmap save function
save_pheatmap <- function(x, filename, width=12, height=12){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(grepl(".png",filename)){
    png(filename, width=width, height=height, units = "in", res=300)
    grid::grid.newpage()
    gridExtra::grid.arrange(x$gtable, vp=grid::viewport(width=0.9, height=1))
    # grid::grid.draw(x$gtable)
    dev.off()
  }
  else if(grepl(".pdf",filename)){
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    gridExtra::grid.arrange(x$gtable, vp=grid::viewport(width=0.9, height=1))
    # grid::grid.draw(x$gtable)
    dev.off()
  }
  else{
    print("Filename did not contain '.png' or '.pdf'")
  }
}