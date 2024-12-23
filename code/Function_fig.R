###function for figure
########combine two plot function########
library(gtable)
library(grid)
library(jtools)
library(ggplot2)
p_out <- function(p){
  if(p>=0.01){p_value=sprintf("%0.2f",p)}
  if(p>=0.001 & p<0.01){p_value=sprintf("%0.3f",p)}
  if(p<0.001){p_value="<0.001"}
  p_value
}
p_out_fig <- function(p){
  if(p < 0.001){p <- "< 0.001"}
  if(p >= 0.001 & p < 0.01){p <- paste("=",sprintf("%0.3f",p))}
  if(p >= 0.01){p <- paste("=",sprintf("%0.2f",p))}
  p
}

hinvert_title_grob <- function(grob){
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
}
