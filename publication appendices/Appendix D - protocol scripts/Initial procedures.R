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
# cerUB - Initial procedures (preparing for protocols)
# This script should be executed before applying a cerUB protocol
##########################################################

library(cerUB)

# Directories -------------------------------------------------------------

directories <- list(
  data = "data",
  transCoDa = "transformed CoDa",

  prot1 = "Protocol_1_geochemical_data",
  prot2 = "Protocol_2_petrographic_data",
  prot3 = "Protocol_3_geochemical_and_petrographic_data",
  prot4 = "Protocol_4_provenance_data",
  prot4_Shipwreck = "Protocols_4_provenance_data_with_shipwrecks"
)

lapply(directories, dir.create, showWarnings = FALSE)


# Read data ---------------------------------------------------------------

data(amphorae)

# alternativelly, using your own data (e.g., CSV)
# dt <- cbind(read.csv(paste(directories$data,
#                            "petrographic_data.csv",sep="/"),
#                      row.names=1), # assuming the first column contains row names
#             read.csv(paste(directories$data,
#                            "geochemical_data.csv",sep="/"),
#                      row.names=1)) # assuming the first column contains row names


# Codify petrographical variables -------------------------------------

varCode <- code_variables(amphorae)


# Clean and format data ---------------------------------------------------

cleanAmphorae <- clean_and_format(amphorae,
                                  completion_variable = c("CHARAC", "complete"),
                                  categorical_columns = 1:112,
                                  numerical_columns = 113:ncol(amphorae),
                                  as_na = c("NULL", "indeterminate", "unfired"),
                                  method = NULL,
                                  # don't use the following variables
                                  columns_to_exclude = c("VOID_VESIC_MEGA",
                                                         "VOID_VUGH_MEGA",
                                                         "VOID_CHAN_MEGA",
                                                         "VOID_PLAN_MEGA",
                                                         "COAR_R_DAC_AND",
                                                         "COAR_R_EVAP",
                                                         "COAR_R_CONGBREC",
                                                         "COAR_R_SERP",
                                                         "COAR_C_SPL",
                                                         "COAR_C_OPX",
                                                         "COAR_C_OL",
                                                         "COAR_C_SIL",
                                                         "COAR_C_ST",
                                                         "COAR_C_ZRN",
                                                         "COAR_C_PY",
                                                         "FINE_C_OPX",
                                                         "FINE_C_ZRN"),
                                  # don't use the following observations
                                  rows_to_exclude = c("PV4033", # PV4-IND4
                                                      "PV4017", # PV4-CAMP
                                                      # PV4-ITT
                                                      "PV4021", "PV4023",
                                                      "PV4024", "PV4025",
                                                      "PV4035", "PV4037",
                                                      # PV4-NAP
                                                      "PV4022", "PV4026",
                                                      "PV4027", "PV4028",
                                                      "PV4029", "PV4030",
                                                      "PV4036")
)


# Subsetting criteria -----------------------------------------------------

# Build vector indicating wheter each observation is from a shipwreck
isShipwreck <-
  cleanAmphorae$Site_Name=="Cap del Vol" |
  cleanAmphorae$Site_Name=="Ullastres I" |
  cleanAmphorae$Site_Name=="Port-Vendres 4"


# Build vectors indicating provenance group and whether observations are
# true outliers (IND). Also, reformat FabricGroup and ChemReferenceGroup,
# so true outliers are single out separatelly and not as a extra group.
ProvenanceGroup <- c()
isTrueIND <- c()

# coerce the original group variables (factors) into character vectors
# so we can use stringr package to operate on them.
cleanAmphorae$FabricGroup <- as.character(cleanAmphorae$FabricGroup)
cleanAmphorae$ChemReferenceGroup <- as.character(cleanAmphorae$ChemReferenceGroup)

for (i in 1:nrow(cleanAmphorae)){

  groupChem <-
    stringr::str_split(cleanAmphorae$ChemReferenceGroup[i], "-")[[1]]
  groupFabric <-
    stringr::str_split(cleanAmphorae$FabricGroup[i], "-")[[1]]
  group <- ""
  isATrueInd <- FALSE

  if (groupChem[2] == "IND" || groupFabric[2] == "IND") {

    group <- cleanAmphorae$ChemReferenceGroup[i]

    if (!isShipwreck[i]) isATrueInd <- TRUE

    index <- 1
    for (j in 1:length(ProvenanceGroup)){

      if (ProvenanceGroup[j] == paste(group, index, sep = ""))
        index <- index + 1

    }

    group <- paste(group, index, sep = "")
    cleanAmphorae$ChemReferenceGroup[i] <- group
    cleanAmphorae$FabricGroup[i] <- group
  }
  else {
    if (groupChem[1] == "ULL" || groupChem[1] == "PV4" || groupChem[1] == "CDV"){
      group <- cleanAmphorae$ChemReferenceGroup[i]
    }
    else if (groupChem[1] == groupFabric[1]){
      group <- groupChem[1]
    }
  }
  ProvenanceGroup <- c(ProvenanceGroup, group[1])
  isTrueIND <- c(isTrueIND, isATrueInd)
}


# Build lists of named group factors for easiness of reference ------------

# list aiming to define workshops productions, so no shipwrecks
factor_list <-
  list(
    Site = factor(cleanAmphorae$Site_Name[!isShipwreck]),
    FabricGroup = factor(cleanAmphorae$FabricGroup[!isShipwreck]),
    ChemGroup = factor(cleanAmphorae$ChemReferenceGroup[!isShipwreck]),
    ProvGroup = factor(ProvenanceGroup[!isShipwreck])
  )

# list aiming to assign shipwreck observations to workshop productions,
# so with shipwreck samples but no true outliers
factor_list_Shipwreck <-
  list(
    Site = factor(cleanAmphorae$Site_Name[!isTrueIND]),
    FabricGroup = factor(cleanAmphorae$FabricGroup[!isTrueIND]),
    ChemGroup = factor(cleanAmphorae$ChemReferenceGroup[!isTrueIND]),
    ProvGroup = factor(ProvenanceGroup[!isTrueIND])
  )


# Build lists of named point types vectors for easiness of reference -------------

# point type full vectors
labels_code <- as.character(row.names(cleanAmphorae))
labels_cross <- rep("+", nrow(cleanAmphorae))
labels_x <- rep(4, nrow(cleanAmphorae)) # using pch code
labels_point <- rep(20, nrow(cleanAmphorae)) # using pch code

# list aiming to define workshops productions, so no shipwrecks
labels_list <- list(
  Code = labels_code[!isShipwreck],
  Cross = labels_cross[!isShipwreck],
  X = labels_x[!isShipwreck],
  Point = labels_point[!isShipwreck]
)

# list aiming to assign shipwreck observations to workshop productions,
# so with shipwreck samples but no true outliers
labels_list_Shipwreck <- list(
  Code = labels_code[!isTrueIND],
  Cross = labels_cross[!isTrueIND],
  X = labels_x[!isTrueIND],
  Point = labels_point[!isTrueIND]
)


# Build lists of named group color vectors for easiness of reference ------------

color_list <- list()
for (i in 1:length(factor_list)){
  cv <- rainbow(nlevels(factor_list[[i]]), v=.8)
  color_list[[i]] <- cv
  names(color_list)[i] = names(factor_list)[i]
}

color_list_Shipwreck <- list()
for (i in 1:length(factor_list_Shipwreck)){
  cv <- rainbow(nlevels(factor_list_Shipwreck[[i]]), v=.8)
  color_list_Shipwreck[[i]] <- cv
  names(color_list_Shipwreck)[i] = names(factor_list_Shipwreck)[i]
}


# Enunciate exception columns ---------------------------------------------

# enunciates which ordinal variables have "none" as a exceptional value
# when calculating the distance between values.

excep_cols <- c("INCLUS_DISTRIB","INCLUS_ORIENT","COAR_ROUNDNESS",
                "COAR_FORM","COAR_SPACING","COAR_SORTING","FINE_FORM")


# Check order of petrographic variables ------------------------------------

# no need to save it, because apply_protocol will do it internally
str(order_petro(cleanAmphorae))


# Select (and save transformed geochemical data) ----------------------------

chemVars16 <- c("Fe2O3","Al2O3","TiO2","MgO","CaO","SiO2",
                "Th","Nb","Zr","Y","Ce","Ga","V","Zn","Ni","Cr")

# Save transformed CoDa to file (optional)
# no need to save it in the environment, because apply_protocol will transform
# the data internally and save the results in ordination_object$logratio_data,
# when applicable
# NOTE: In the output table, CoDa columns will be ordered as:
# (1) CoDa variables not transformed,
# (2) Raw version of the selected CoDa variables,
# (3) Transformed version of the selected CoDa variables.

write(transform_coda(cleanAmphorae,
                     coda_variables = chemVars16,
                     method = c("CLR")),
      file = paste(directories$transCoDa, "transAmphorae_clr.csv", sep = "/"))

