
##Creating enmtools.species objects
setwd("C:/project_sigs/Piaya/modelos/piaya_mexicana/G/Set_1/present_all/")###working directory with environmental variables
env.files <- list.files(".",pattern = "*.asc$",full.names = T)###read the environmental variables within folder.
env<- stack(env.files)##create the stack for environmental variables
env <- setMinMax(env)
plot(env[[1]])##plor the first environmetal variable in the stack.

##specie1: P.c. mexicana
mexicana <- enmtools.species()
mexicana

mexicana.path <- paste(system.file(package="ENMTools"), "/mexicana.csv", sep='')
mexicana <- enmtools.species(species.name = "mexicana", 
                              presence.points = read.csv(mexicana.path))
mexicana$range <- background.raster.buffer(mexicana$presence.points, 50000, mask = env)
mexicana$background.points <- background.points.buffer(points = mexicana$presence.points,
                                                        radius = 20000, n = 1000, mask = env[[1]])
mexicana <- check.species(mexicana)
interactive.plot.enmtools.species(mexicana)


##specie2: P.c. thermophila
thermophila <- enmtools.species()
thermophila

thermophila.path <- paste(system.file(package="ENMTools"), "/thermophila.csv", sep='')
thermophila <- enmtools.species(species.name = "thermophila", 
                             presence.points = read.csv(thermophila.path))
thermophila$range <- background.raster.buffer(thermophila$presence.points, 50000, mask = env)
thermophila$background.points <- background.points.buffer(points = thermophila$presence.points,
                                                       radius = 20000, n = 1000, mask = env[[1]])
thermophila <- check.species(thermophila)
interactive.plot.enmtools.species(thermophila)


##Rangebreak tests
directorio = setwd("C:/project_sigs/Piaya/modelos/")

rbl.glm <- rangebreak.linear(mexicana, thermophila, env, type = "glm", nreps = 1000,
                             nback = 1000, rep.dir = directorio)

rbb.bc <- rangebreak.blob(mexicana, thermophila, env, type = "bc", nreps = 1000,
                          nback = 1000)
rbl.glm
rbb.bc

rbl.glm$replicate.models

##Fin/end