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
# cerUB - Protocol 4 (with shipwrecks, no true outlier)
#
# 1. geochemical compositional data (CoDa) & Ordinal petrographic data (Provenance selection)
# 2. centred log-ratio transformation (clr) & transformed to ranks
# 3. Extended Gower Coefficient: Euclidean distance & Relative ranking
# 4. Principal Coordinates Analysis (PCoA)
# 4. PERMANOVA & PERMDISP tests
#
##########################################################


# NOTE: "Initial procedures.R" must be ran once first.


# Apply protocol 4 to identify the provenance of shipwreck samples --------------------

# project in 2d (default)

prot4_Shipwreck_2d <- apply_ordination(cleanAmphorae[!isTrueIND,], # with shipwrecks, not true outliers
                                       "4", # select protocol 4
                                       exception_columns = excep_cols,
                                       variable_tags = varCode,
                                       coda_override = chemVars16,
                                       coda_transformation = "CLR")

# Simplify CoDa names for plot clearness
prot4_Shipwreck_2d <- simplify_coda_names(prot4_Shipwreck_2d)


# project in 3d (dimensions = 3)

prot4_Shipwreck_3d <- apply_ordination(cleanAmphorae[!isTrueIND,], # no shipwrecks
                                       "4", # select protocol 4
                                       exception_columns = excep_cols,
                                       variable_tags = varCode,
                                       coda_override = chemVars16,
                                       coda_transformation = "CLR",
                                       dimensions = 3)

# Simplify CoDa names for plot clearness
prot4_Shipwreck_3d <- simplify_coda_names(prot4_Shipwreck_3d)

# Test the given provenance groups --------------------------------

# Note1: Tests depend only on the distance matrix and the groups, not on the projections
# Note2: This test batch may take several minutes depending on the size of the
# data matrix and the number of groups.
prot4_Shipwreck_tests <- test_groups(prot4_Shipwreck_2d$dist_matrix,
                                     factor_list_Shipwreck$ProvGroup)

# Biplot 2D ---------------------------------------------------------------

arrows_label_adj <- rbind(c(.5,0),c(.5,1),c(.5,0),c(.5,1),c(.5,0),c(.5,1),c(.8,0),
                          c(1,.5),c(.5,0),c(1,.2),c(.5,1),c(.2,.7))
row.names(arrows_label_adj) <- c("S7","S8","S4","CaO","MgO","S11","L48",
                                 "SiO2","Ce","Nb","Th","TiO2")

# better "preview" (R UI device) version
biplot2d3d::biplot_2d(prot4_Shipwreck_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      ylim = c(-.3,.25),
                      point_type = "point",
                      groups = factor_list_Shipwreck$ProvGroup,
                      group_color = color_list_Shipwreck$ProvGroup,
                      group_label_cex = 0.6,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 0.6,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot4_Shipwreck_2d$sub2D,
                      test_text = 
                        prot4_Shipwreck_tests$text(prot4_Shipwreck_tests),
                      test_cex = 0.8,
                      test_fig = c(0, 0.5, 0.65, .99),
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      output_type = "preview")

# better PNG version
biplot2d3d::biplot_2d(prot4_Shipwreck_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.25),
                      point_type = "point",
                      groups = factor_list_Shipwreck$ProvGroup,
                      group_color = color_list_Shipwreck$ProvGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 2,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot4_Shipwreck_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = 
                        prot4_Shipwreck_tests$text(prot4_Shipwreck_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 2,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot4_Shipwreck_Biplot2D",
                      directory = directories$prot4_Shipwreck,
                      width = 1000, height = 1000,
                      output_type = "png")


# better EPS version
biplot2d3d::biplot_2d(prot4_Shipwreck_2d,
                      ordination_method = "PCoA",
                      invert_coordinates = c (TRUE,TRUE),
                      grid_cex = 2.5,
                      ylim = c(-.3,.25),
                      point_type = "point",
                      groups = factor_list_Shipwreck$ProvGroup,
                      group_color = color_list_Shipwreck$ProvGroup,
                      group_label_cex = 1.5,
                      arrow_mim_dist = .5,
                      arrow_label_cex = 1.5,
                      arrow_cex = 0.2,
                      arrow_lwd = 2.5,
                      arrow_fig = c(.6,.95,0,.35),
                      arrow_label_adj_override = arrows_label_adj,
                      subtitle = prot4_Shipwreck_2d$sub2D,
                      subtitle_cex = 2.5,
                      test_text = prot4_Shipwreck_tests$text(prot4_Shipwreck_tests),
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      test_cex = 1.5,
                      fitAnalysis_fig = c(0,.7,.05,.5),
                      file_name = "Prot4_Shipwreck_Biplot2D",
                      directory = directories$prot4_Shipwreck,
                      width = 1000, height = 1000,
                      output_type = "eps")


# Biplot 3D ---------------------------------------------------------------

biplot2d3d::biplot_3d(prot4_Shipwreck_3d,
                      ordination_method = "PCoA",
                      point_type = "point",
                      groups = factor_list_Shipwreck$FabricGroup,
                      group_color = color_list_Shipwreck$FabricGroup,
                      group_representation = "stars",
                      star_centroid_radius = 0,
                      star_label_cex = .8,
                      arrow_min_dist = .5,
                      arrow_body_length = .025,
                      subtitle = prot4_Shipwreck_3d$sub3D,
                      test_text = 
                        prot4_Shipwreck_tests$text(prot4_Shipwreck_tests),
                      test_cex = 1.25,
                      test_spacing_line = 0.9,
                      test_fig =c(0, 0.5, 0.72, .99),
                      view_zoom = 0.9)

# save animated GIF or PNG snapshot
biplot2d3d::animation(directory = directories$prot4_Shipwreck,
                       file_name = "Prot4_Shipwreck_Biplot3D")


