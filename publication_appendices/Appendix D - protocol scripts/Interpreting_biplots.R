# License notice:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################
# cerUB - Appendix on interpreting biplots 
##########################################################


# NOTE: "Initial procedures.R", "Protocol_1.R", and "Protocol_2.R" must be ran once first.

# Protocol 1 example----

# Let's recover protocol 1 override for variable label positions 
arrows_label_adj <- rbind( c(.5, -.5), c(1, .5), c(1.2, 1.2), c(1.2, .4),
                           c(.8, .5),
                           c(0, 0), c(-.2, 1), c(.5, 1.2), c(-.5, .5),
                           c(-.2, .5), c(0, .5), c(0, 0))
row.names(arrows_label_adj) <- c("Fe2O3", "Al2O3", "SiO2", "TiO2",
                                 "MgO",
                                 "Th", "Nb", "Cr", "Ce",
                                 "Ga", "Zn", "Y")

# Protocol 1, representing and testing chemical reference groups
biplot_2d(prot1, 
          groups = factor_list$ChemGroup, 
          group_color = color_list$ChemGroup,
          group_label_cex = 0.6,
          invert_coordinates = c(TRUE, TRUE),
          arrow_label_cex = 0.7,
          test_text = prot1_tests$text(prot1_tests),
          test_cex = 0.8,
          test_fig = c(0, 0.5, 0.65, .99))


# Interpret Protocol 1 in terms of CaO

# Create factor variable containing the classification (5 categories)
CaO_level <- cut(cleanAmphorae$CaO[!isShipwreck], 5)

# Select 5 colours from the 'topo.colors' palette
CaO_level_colors <- topo.colors(nlevels(CaO_level))

# Test the classification
prot1_tests_CaO <- test_groups(prot1$dist_matrix, CaO_level)

# This is for highlighting CaO arrow
arrow_colors <- rep("darkorange", nrow(prot1$loadings))
arrow_colors[row.names(prot1$loadings) == "CaO"] <- "red"

# Protocol 1, grouping by level of CaO content
biplot_2d(prot1, 
          groups = CaO_level,
          group_color = CaO_level_colors,
          group_star_cex = 0,
          group_label_cex = 0,
          show_group_legend = TRUE,
          group_legend_title = "CaO",
          group_legend_title_pos = c(0.5,0.9),
          group_legend_text_cex = 0.8,
          group_legend_fig = c(0.7,0.99,0.68,0.95),
          invert_coordinates = c(TRUE, TRUE),
          arrow_label_cex = 0.8,
          arrow_fig = c(.6,.95,0,.35),
          arrow_label_adj_override = arrows_label_adj,
          arrow_color = arrow_colors,
          test_text = prot1_tests_CaO$text(prot1_tests_CaO),
          test_cex = 0.8,
          test_fig = c(0, 0.5, 0.65, .99),
          output_type = "preview")

#----

# Protocol 2 examples----

# Let us recover protocol 2a override for variable label positions 
arrows_label_adj <- rbind(c(.5,.8),c(.5,1),c(.5,1),c(.5,0),c(.5,1),
                          c(.5,0),c(0,.5))
row.names(arrows_label_adj) <- c("L48","L24","L5","L36","S7",
                                 "S8","S11")

# This will help us select different arrow colours
isDisplayed <- 
  row.names(prot2a_2d$loadings) %in% row.names(
    filter_arrows(prot2a_2d$loadings, min_dist = 0.5))

# protocol 2a, representing and testing fabric groups
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      ylim = c(-.3,.29),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 0.6,
                      arrow_mim_dist = 0.5,
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2a_2d$sub2D,
                      test_text = prot2a_tests$text(prot2a_tests),
                      test_fig = c(0, 0.5, 0.65, .99),
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      output_type = c("preview"))

# Interpret Protocol 2a in terms of INCLUS_ORIENT

# You may want to assure that the true categories are corectly represented:
cleanAmphorae <- order_petro(cleanAmphorae)

levels(cleanAmphorae$INCLUS_ORIENT[!isShipwreck])

# Let us declare this factor separately as an object for clearness
I2 <- cleanAmphorae$INCLUS_ORIENT[!isShipwreck]

# Select colours from the 'topo.colors' palette
I2_colors <- topo.colors(nlevels(I2))

# Test the classification
prot1_tests_I2 <- test_groups(prot2a_2d$dist_matrix, I2)

# This is for highlighting CaO arrow
arrow_colors <- rep("darkorange", nrow(prot2a_2d$loadings))
arrow_colors[row.names(prot2a_2d$loadings) == "I2"] <- "red"
# filter arrows colours, since not all variables are displayed
arrow_colors <- arrow_colors[isDisplayed]

# Protocol 2a, grouping by INCLUS_ORIENT
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      ylim = c(-.3,.29),
                      groups = I2,
                      group_color = I2_colors,
                      group_star_cex = 0,
                      group_label_cex = 0,
                      show_group_legend = TRUE,
                      group_legend_title = "INCLUS_ORIENT",
                      group_legend_title_pos = c(0.5,0.9),
                      group_legend_text_cex = 0.8,
                      group_legend_fig = c(0.6,0.99,0.68,0.95),
                      arrow_mim_dist = .5,
                      arrow_label_cex = 0.8,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      arrow_color = arrow_colors,
                      subtitle = prot2a_2d$sub2D,
                      test_text = prot1_tests_I2$text(prot1_tests_I2),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      output_type = c("preview"))

# Interpret Protocol 2a in terms of COAR_R_CHERT

# Let us declare this factor separately as an object for clearness
L33 <- cleanAmphorae$COAR_R_CHERT[!isShipwreck]

# Select colours from the 'topo.colors' palette
L33_colors <- topo.colors(nlevels(L33))

# Test the classification
prot1_tests_L33 <- test_groups(prot2a_2d$dist_matrix, L33)

# This is for highlighting CaO arrow
arrow_colors <- rep("darkorange", nrow(prot2a_2d$loadings))
arrow_colors[row.names(prot2a_2d$loadings) == "L33"] <- "red"
# filter arrows colours, since not all variables are displayed
arrow_colors <- arrow_colors[isDisplayed]

# Protocol 2a, grouping by COAR_R_CHERT
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      ylim = c(-.3,.29),
                      groups = L33,
                      group_color = L33_colors,
                      group_star_cex = 0,
                      group_label_cex = 0,
                      show_group_legend = TRUE,
                      group_legend_title = "COAR_R_CHERT",
                      group_legend_title_pos = c(0.5,0.9),
                      group_legend_text_cex = 0.8,
                      group_legend_fig = c(0.6,0.99,0.68,0.95),
                      arrow_mim_dist = .5,
                      arrow_label_cex = 0.8,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      arrow_color = arrow_colors,
                      subtitle = prot2a_2d$sub2D,
                      test_text = prot1_tests_L33$text(prot1_tests_L33),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      output_type = c("preview"))

#----