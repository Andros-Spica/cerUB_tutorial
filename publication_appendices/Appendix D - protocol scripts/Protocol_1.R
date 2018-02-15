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
# cerUB - Protocol 1
#
# 1. geochemical compositional data (CoDa)
# 2. isometric log-ratio transformation (ilr)
# 3. robust Principal Components Analysis (robPCA), implicitly using Euclidean distance
# 4. PERMANOVA & PERMDISP tests
# 5. (outlier detection)
#
##########################################################


# NOTE: "Initial procedures.R" must be ran once first.



# Apply protocol 1 to confirm workshops' chemical reference groups --------

prot1 <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                          "1", # select protocol 1
                          coda_override = chemVars16,
                          coda_transformation = "ILR")

# Simplify CoDa names for plot clearness
prot1 <- simplify_coda_names(prot1)


# Test the given chemical reference groups --------------------------------

# Note: This test batch may take several minutes depending on the size of the
# data matrix and the number of groups.
prot1_tests <- test_groups(prot1$dist_matrix, factor_list$ChemGroup)


# Biplot 2D ---------------------------------------------------------------

arrows_label_adj <- rbind( c(.5, -.5), c(1, .5), c(1.2, 1.2), c(1.2, .4),
                           c(.8, .5),
                           c(0, 0), c(-.2, 1), c(.5, 1.2), c(-.5, .5),
                           c(-.2, .5), c(0, .5), c(0, 0))
row.names(arrows_label_adj) <- c("Fe2O3", "Al2O3", "SiO2", "TiO2",
                                 "MgO",
                                 "Th", "Nb", "Cr", "Ce",
                                 "Ga", "Zn", "Y")

# better "preview" (R UI device) version
biplot2d3d::biplot_2d(prot1, 
                      groups = factor_list$ChemGroup, 
                      group_color = color_list$ChemGroup,
                      group_ellipse_cex = 2.5,
                      group_label_cex = 0.6,
                      invert_coordinates = c(FALSE, TRUE),
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      test_text = prot1_tests$text(prot1_tests),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      output_type = "preview")

# better PNG version
biplot2d3d::biplot_2d(prot1,
                      ordination_method = "PCA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.1,.1),
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_label_cex = 1.5,
                      arrow_label_cex = 2,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot1$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot1_tests$text(prot1_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 2,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot1_Biplot2D",
                      directory = directories$prot1,
                      width = 1000, height = 1000,
                      output_type = "png")

# better EPS version
biplot2d3d::biplot_2d(prot1,
                      ordination_method = "PCA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.1,.1),
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_label_cex = 1.5,
                      arrow_label_cex = 1.5,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot1$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot1_tests$text(prot1_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 1.5,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot1_Biplot2D",
                      directory = directories$prot1,
                      width = 1000, height = 1000,
                      output_type = c("eps"))

# Biplot 3D ---------------------------------------------------------------

biplot2d3d::biplot_3d(prot1,
                      ordination_method = "PCA",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      point_type = "point",
                      group_representation = "stars",
                      star_centroid_radius = 0,
                      star_label_cex = .8,
                      test_text = prot1_tests$text(prot1_tests),
                      test_cex = 1.25,
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99))

# save animated GIF or PNG snapshot
biplot2d3d::animation(file_name = "Prot1_Biplot3D")


