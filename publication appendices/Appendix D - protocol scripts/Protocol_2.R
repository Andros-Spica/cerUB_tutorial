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
# cerUB - Protocol 2
#
# 1. Ordinal petrographic data
# 2. transformed to ranks
# 3. (a) Relative ranking
#    (b) Neighbor interchange
# 4. (a) Principal Coordinates Analysis (PCoA)
#    (b) Non-metric Dimensional Scaling (NMDS)
# 4. PERMANOVA & PERMDISP tests
#
##########################################################


# NOTE: "Initial procedures.R" must be ran once first.


# Apply protocol 2 to identify workshops ----------------------------------

# project in 2d (default)

prot2a_2d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                              "2a", # select protocol 2a (RRD & PCoA)
                              exception_columns = excep_cols,
                              variable_tags = varCode)

prot2b_2d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                              "2b", # select protocol 2a (NI & NMDS)
                              exception_columns = excep_cols,
                              variable_tags = varCode)

# project in 3d (dimensions = 3)

prot2a_3d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                              "2a", # select protocol 2a (RRD & PCoA)
                              exception_columns = excep_cols,
                              variable_tags = varCode,
                              dimensions = 3)

prot2b_3d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                              "2b", # select protocol 2a (NI & NMDS)
                              exception_columns = excep_cols,
                              variable_tags = varCode,
                              dimensions = 3)

# Test the given fabric groups --------------------------------

# Note1: Tests depend only on the distance matrix and the groups, not on the projections
# Note2: This test batch may take several minutes depending on the size of the
# data matrix and the number of groups.
prot2a_tests <- test_groups(prot2a_2d$dist_matrix,
                            factor_list$FabricGroup)
prot2b_tests <- test_groups(prot2b_2d$dist_matrix,
                            factor_list$FabricGroup)

# Biplot 2D ---------------------------------------------------------------

## prot2a
arrows_label_adj <- rbind(c(.5,.8),c(.5,1),c(.5,1),c(.5,0),c(.5,1),c(.5,0),c(0,.5))
row.names(arrows_label_adj) <- c("L48","L24","L5","L36","S7","S8","S11")

# better "preview" (R UI device) version
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      #ylim = c(-.3,.27),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 0.6,
                      arrow_mim_dist = 0,
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2a_2d$sub2D,
                      test_text = prot2a_tests$text(prot2a_tests),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      output_type = "preview")

# better PNG version
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.27),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 2,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2a_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot2a_tests$text(prot2a_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 2,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot2a_Biplot2D",
                      directory = directories$prot2,
                      width = 1000, height = 1000,
                      output_type = "png")

# better EPS version
biplot2d3d::biplot_2d(prot2a_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.27),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 1.5,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2a_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot2a_tests$text(prot2a_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 1.5,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot2a_Biplot2D",
                      directory = directories$prot2,
                      width = 1000, height = 1000,
                      output_type = "eps")

## prot2b

arrows_label_adj <- rbind(c(.5,1),c(.5,0),c(.5,1),c(.5,1),c(.5,0),c(0,.5),c(1,.5))
row.names(arrows_label_adj) <- c("S7","S8","CLAY","L24","L43","L5","L36")

# better "preview" (R UI device) version
biplot2d3d::biplot_2d(prot2b_2d,
                      ordination_method = "NMDS",
                      invert_coordinates = c (FALSE,TRUE),
                      ylim = c(-.33,.25),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 0.6,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2b_2d$sub2D,
                      test_text = prot2b_tests$text(prot2b_tests),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      fitAnalysis_fig = c(.1,.5,.1,.3),
                      output_type = "preview")

# better PNG version
biplot2d3d::biplot_2d(prot2b_2d,
                      ordination_method = "NMDS",
                      invert_coordinates = c (FALSE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.25),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 2,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2b_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot2b_tests$text(prot2b_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 2,
                      fitAnalysis_fig = c(.1,.5,.1,.3),
                      file_name = "Prot2b_Biplot2D",
                      directory = directories$prot2,
                      width = 1000, height = 1000,
                      output_type = "png")

# better EPS version
biplot2d3d::biplot_2d(prot2b_2d,
                      ordination_method = "NMDS",
                      invert_coordinates = c (FALSE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.25),
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 1.5,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot2b_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot2b_tests$text(prot2b_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 1.5,
                      fitAnalysis_fig = c(.1,.5,.1,.3),
                      file_name = "Prot2b_Biplot2D",
                      directory = directories$prot2,
                      width = 1000, height = 1000,
                      output_type = "eps")

# Biplot 3D ---------------------------------------------------------------

## prot2a

biplot2d3d::biplot_3d(prot2a_3d,
                      ordination_method = "PCoA",
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_representation = "stars",
                      star_centroid_radius = 0,
                      star_label_cex = .8,
                      arrow_min_dist = .5,
                      arrow_body_length = .025,
                      subtitle = prot2a_3d$sub3D,
                      test_text = prot2a_tests$text(prot2a_tests),
                      test_cex = 1.25,
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      view_zoom = 0.9)

# save animated GIF or PNG snapshot
biplot2d3d::animation(directory = directories$prot2,
                       file_name = "Prot2a_Biplot3D")


biplot2d3d::biplot_3d(prot2b_3d,
                      ordination_method = "NMDS",
                      point_type = "point",
                      groups = factor_list$FabricGroup,
                      group_color = color_list$FabricGroup,
                      group_representation = "stars",
                      star_centroid_radius = 0,
                      star_label_cex = .8,
                      arrow_min_dist = .5,
                      arrow_body_length = .025,
                      subtitle = prot2b_3d$sub3D,
                      test_text = prot2b_tests$text(prot2b_tests),
                      test_cex = 1.25,
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      view_zoom = 0.9)

# save animated GIF or PNG snapshot
biplot2d3d::animation(directory = directories$prot2,
                       file_name = "Prot2b_Biplot3D")

