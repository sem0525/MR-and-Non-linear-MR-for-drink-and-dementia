####function for figure
library(ggplot2)
library(gtable)
###combind the two figure
combind_figs <- function(F1,F2){
  g1 <- ggplotGrob(F1)
  g2 <- ggplotGrob(F2)
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp1 <- c(subset(g1$layout, name == 'panel', se = t:r))
  pp2 <- c(subset(g2$layout, name == 'panel', se = t:r))
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  
  # Overlap panel for second plot on that of the first plot
  g1_2 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp1$t, pp1$l, pp1$b, pp1$l)
  
  ###y title
  ##from g1
  indexy_g1 <- which(g1$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab_g1 <- g1$grobs[[indexy_g1]]        # Extract that grob
  ylab_g1_new <- hinvert_title_grob(ylab_g1)     # Swap margins and fix justifications
  # Put the transformed label on the right side of g1
  g1_3 <- gtable_add_cols(g1_2, g1$widths[g1$layout[indexy_g1, ]$l], pp$r)
  
  g1_4 <- gtable_add_grob(g1_3, ylab_g1, 7, 7, 7, 7
                          ,clip = 'off', name = 'ylab-r')
  ########y-axix
  ##g1
  index_g1 <- which(g1$layout$name == 'axis-l')  # Which grob
  yaxis <- g1_4$grobs[[index_g1]]          # Extract the grob
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1_5 <- gtable_add_cols(g1_4, g1$widths[g1$layout[index_g1, ]$l], pp$r)
  g1_5 <- gtable_add_grob(g1_5, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  
  g1_6 <- g1_5
  ylab_l_index <- which(g1_5$layout$name == 'ylab-l')
  ylab_l_index_g2 <- which(g2$layout$name == 'ylab-l')
  g1_6$grobs[ylab_l_index] <- g2$grobs[ylab_l_index_g2]
  
  
  
  g1_7 <- g1_6
  yaxis_g1_7 <- which(g1_7$layout$name == 'axis-l')
  yaxis_g2 <- which(g2$layout$name == 'axis-l')
  
  g1_7$grobs[yaxis_g1_7] <- g2$grobs[yaxis_g2]
  
  g <- g1_7
  axis_label_index <- which(g$layout$name == "ylab-l")
  axis_label_grob <- g$grobs[[axis_label_index]]
  
  # Reduce the left margin of the y-axis label
  axis_label_grob <- editGrob(
    axis_label_grob, 
    vp = viewport(layout.pos.row = axis_label_grob$vp$layout.pos.row,
                  layout.pos.col = axis_label_grob$vp$layout.pos.col,
                  width = unit(1, "npc") - unit(10, "mm"), # Adjust the width here
                  height = unit(1, "npc"))
  )
  
  # Put the adjusted label back into the gtable
  g$grobs[[axis_label_index]] <- axis_label_grob
  
  # Redraw the plot
  grid.newpage()
  grid.draw(g)
  g
}
