# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
#                    Set configurations for modelling pipeline
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

gc(reset=TRUE)

# ------------------- 

results.directory <- "results/"
data.folder <- file.path("data", "Input")

# -------------------

project.name <- "Dictyota_dichotoma"
models.to.use <- c("BRT") 
# -------------------

rasters.names.s <- c("light.at.bottom.Benthic.Min.Var.Mean","Present.Benthic.Min.Depth.Salinity.Lt.min","Present.Benthic.Min.Depth.Temperature.Lt.max","Present.Benthic.Min.Depth.Temperature.Lt.min") 

variable.monotonic.response <- c(+1,+1,-1,+1) #Hardcode relationship between occs and environmental variables 


# -------------------------------

bathymetry.file <- file.path("data", "Input", "Dependencies", "Rasters","BO2BathymetryDepthMin.tif") # DepthMin for Shallow sites

# ------------------------------------------------------------------------------------

n.cores <- 2 #number of cores 

# ------------------------------------------------------------------------------------

region.buffer <- c(12,12,12,12) #"Code Jorge"
depth.buffer <- 90 # max depth for predictions is 140m (50m + 90m)

# ---------------------
# Species

min.depth = 0   #minimum depth at which species occurs
max.depth = 50  # maximum depth at which species occurs

relocate.occ.species.depth <- TRUE
relocate.occ.distance <- 25  #Max number of km that occurence point may be located inland, if lower it is relocated to closest pixel

# ---------------------
# Cross-validation
cv.type <- "blocks.latitudinal"  #latitudinal blocks to do cross validation 
cv.k <- 7  # use 7 blocks
remove.cv.edges <- TRUE # (remove the 2 blocks at the edges for testing: 5 testing blocks remaining)


sre.threhold <- 0.005 # 0.025
simplify.model <- FALSE

# ------------------------------------------------------------------------------------

brt.learning.complex.span <- c(0.01,0.005,0.001)
brt.max.tree.depth <- 1:4 # 1: the number of predictors
