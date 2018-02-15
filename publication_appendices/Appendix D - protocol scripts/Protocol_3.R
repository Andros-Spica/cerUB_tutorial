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
# cerUB - Protocol 3
#
# 1. geochemical compositional data (CoDa) & Ordinal petrographic data
# 2. centred log-ratio transformation (clr) & transformed to ranks
# 3. Extended Gower Coefficient: Euclidean distance & Relative ranking
# 4. Principal Coordinates Analysis (PCoA)
# 4. PERMANOVA & PERMDISP tests
#
##########################################################


# NOTE: "Initial procedures.R" must be ran once first.


# Apply protocol 3 to identify workshops ----------------------------------

# project in 2d (default)

prot3_2d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                             "3", # select protocol 3
                             exception_columns = excep_cols,
                             variable_tags = varCode,
                             coda_override = chemVars16,
                             coda_transformation = "CLR")

# Simplify CoDa names for plot clearness
prot3_2d <- simplify_coda_names(prot3_2d)


# project in 3d (dimensions = 3)

prot3_3d <- apply_ordination(cleanAmphorae[!isShipwreck,], # no shipwrecks
                             "3", # select protocol 3
                             exception_columns = excep_cols,
                             variable_tags = varCode,
                             coda_override = chemVars16,
                             coda_transformation = "CLR",
                             dimensions = 3)

# Simplify CoDa names for plot clearness
prot3_3d <- simplify_coda_names(prot3_3d)

# Test the given provenance groups --------------------------------

# Note1: Tests depend only on the distance matrix and the groups, not on the projections
# Note2: This test batch may take several minutes depending on the size of the
# data matrix and the number of groups.
prot3_tests <- test_groups(prot3_2d$dist_matrix,
                            factor_list$ChemGroup)

# Biplot 2D ---------------------------------------------------------------

arrows_label_adj <- rbind(c(0,0),c(.5,0),c(.8,.2),c(.5,0),c(.5,1),c(.5,1.1),c(1,.5),
                          c(0,.5))
row.names(arrows_label_adj) <- c("L6","Ga","S11","S7","S8","Nb","SiO2",
                                 "Y")

# better "preview" (R UI device) version
biplot2d3d::biplot_2d(prot3_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (FALSE,FALSE),
                      ylim = c(-.3,.29),
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_label_cex = 0.6,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot3_2d$sub2D,
                      test_text = prot3_tests$text(prot3_tests),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      output_type = "preview")

# better PNG version
biplot2d3d::biplot_2d(prot3_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.1,.1),
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 2,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot3_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot3_tests$text(prot3_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 2,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot3_Biplot2D",
                      directory = directories$prot3,
                      width = 1000, height = 1000,
                      output_type = "png")

# better EPS version
biplot2d3d::biplot_2d(prot3_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.1,.1),
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 1.5,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot3_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot3_tests$text(prot3_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 1.5,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot3_Biplot2D",
                      directory = directories$prot3,
                      width = 1000, height = 1000,
                      output_type = "eps")

# Biplot 3D ---------------------------------------------------------------

biplot2d3d::biplot_3d(prot3_3d,
                      ordination_method = "PCoA",
                      point_type = "point",
                      groups = factor_list$ChemGroup,
                      group_color = color_list$ChemGroup,
                      group_representation = "stars",
                      star_centroid_radius = 0,
                      star_label_cex = .8,
                      arrow_min_dist = .5,
                      arrow_body_length = .025,
                      subtitle = prot3_3d$sub3D,
                      test_text = prot3_tests$text(prot3_tests),
                      test_cex = 1.25,
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      view_zoom = 0.9)

# save animated GIF or PNG snapshot
biplot2d3d::animation(directory = directories$prot3,
                       file_name = "Prot3_Biplot3D")


