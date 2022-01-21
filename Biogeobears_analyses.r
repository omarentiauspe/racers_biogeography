##########################################################################################
# R script used to run BioGeoBEARS on calibrated trees.
# Wrote for analyses in Grazziotin et al. (in press) - "Phylogenetic position of Amphisbaena ridleyi".
# Felipe G. Grazziotin - fgrazziotin@gmail.com
##########################################################################################
# This is a modular general script combining original functions and commands with several other, derived from blogs, forums, books and regular papers.
# Feel free to change, use and identify mistakes.
# Here I've just modified the script "introductory example script for the BioGeoBEARS" by Nick Matzke (I've kept Nick's comments on some functions)
# Please, check http://phylo.wikidot.com/biogeobears for other informations about the method
##########################################################################################

##########################################################################################
# Install and Load libraries
##########################################################################################

##############################
# install packages
##############################
# install.packages("ape")
# install.packages("rexpokit")
# install.packages("cladoRcpp")
# library(devtools)
# devtools::install_github(repo="nmatzke/BioGeoBEARS", dependencies=FALSE)
# install.packages("snow")

##############################
# load packages
##############################
library(ape)
library(GenSA)  # GenSA is better than optimx (although somewhat slower)
library(FD)     # for FD::maxent() (make sure this is up-to-date)
library(snow)   # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
# library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

library(strap)

##########################################################################################

##########################################################################################
# Edit phylogeny before running 
##########################################################################################

#read tree
t1<-read.tree("TREE1_clean.tre")

# remove multiple individuals ()

out_to_drop<-c("Sordellina_punctata","Taeniophallus_affinis","Drepanoides_anomalus","Pseudoboa_nigra","Siphlophis_pulcher","Thamnodynastes_hypoconia","Tomodon_dorsatus","Helicops_angulatus","Hydrops_triangularis","Apostolepis_flavotorquata","Elapomorphus_quinquelineatus","Cubophis_cantherigerus","Hypsirhynchus_parvifrons","Arrhyton_taeniatum","Erythrolamprus_miliaris","Xenodon_merremii","Lygophis_elegantissimus","Pseudalsophis_biserialis","Pseudalsophis_elegans_CTMZ7428","Psomophis_joberti","Psomophis_obtusus","Atractus_trihedrurus","Dipsas_catesbyi","Imantodes_cenchoa","Synophis_bicolor","Carphophis_amoenus","Contia_tenuis","Heterodon_platirhinos","Diadophis_punctatus","Farancia_erytrogramma","Sticophanes_ningshaanensis","Thermophis_baileyi","Thermophis_zhaoermii","Pseudoxenodon_bambusicola","Pseudoxenodon_karlschmidti","Coluber_constrictor","Lampropeltis_getula","Calamaria_pavimentata","Natrix_natrix","Rhabdophis_subminiatus","Bungarus_fasciatus","Micrurus_surinamensis","Calliophis_maculiceps","Atractaspis_micropholis","Rhamphiophis_oxyrhynchus","Homalopsis_buccata","Azemiops_feae","Bothriechis_schlegelii","Causus_lichtensteinii","Aplopeltura_boa","Pareas_carinatus","Achalinus_rufescens","Xenodermus_javanicus","Acrochordus_granulatus","Boa_constrictor","Eryx_conicus")
in_to_drop<-c("Philodryas_psammophidea-MeloSampaio","Philodryas_viridissima_1-MeloSampaio","Philodryas_laticeps-MeloSampaio")

pt1<-drop.tip(t1,c(out_to_drop,in_to_drop))

# changing tip lables 

pt1$tip.label<-gsub("Philodryas_agassizii_BRA_RS_LavrasdoSul_CTMZ18499","Philodryas_agassizii_-_PA",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_agassizii_Doc","Philodryas_agassizii_-_CE",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_livida_BRA_GO_Mineiros_CTMZ04172","Philodryas_livida",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_psammophidea_ARG_Catamarca_unknown_CTMZ00274","Philodryas_psammophidea",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_patagoniensis_BRA_RS_CachoeiradoSul_CTMZ16411","Philodryas_patagoniensis_-_CE",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_varia_unknown_CTMZ00202","Philodryas_varia",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_aestiva_BRA_RS_DoutorRicardo_CTMZ18480","Philodryas_aestiva",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_olfersii_1-MeloSampaio","Philodryas_olfersii_-_AF",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_olfersii_BRA_SP_Iporanga_CTMZ16513","Philodryas_olfersii_-_CE",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_arnaldoi_BRA_RS_Osorio_CTMZ16404","Philodryas_arnaldoi",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_mattogrossensis-MeloSampaio","Philodryas_erlandi",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_mattogrossensis_BRA_MG_unknown_CTMZ10870","Philodryas_mattogrossensis",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_nattereri_BRA_TO_unknown_CTMZ00916","Philodryas_nattereri",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_baroni_ARG_Tucuman_Tucuman_CTMZ16395","Philodryas_baroni",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_trilineata_ARG_unknown_CTMZ04935","Philodryas_trilineata",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_chamissonis_CHI_unknown_CTMZ16588","Philodryas_chamissonis",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_laticeps_SB0512","Chlorosoma_laticeps",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_viridissima_BRA_PA_Melgaco_CTMZ14603","Chlorosoma_viridissima",pt1$tip.label)
pt1$tip.label<-gsub("Colubroidea_sp_BrJ-2018a_1-MeloSampaio","Chlorosoma_dunupyana",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_argentea_BRA_AC_unknown_CTMZ16538","Xenoxybelis_argentea",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_georgeboulengeri_BRA_AM_Coari_CTMZ13754","Xenoxybelis_georgeboulengeri",pt1$tip.label)
pt1$tip.label<-gsub("Tropidodryas_serra_BRA_SP_unknown_CTMZ00717","Tropidodryas_serra",pt1$tip.label)
pt1$tip.label<-gsub("Tropidodryas_striaticeps_BRA_SP_EmbuGuacu_CTMZ00185","Tropidodryas_striaticeps",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_simonsii_PERU_Arequipa_Yura_CTMZ16424","Incaspis_simonsii",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_tachymenoides_PERU_Ica_Ica_CTMZ16423","Incaspis_tachymenoides",pt1$tip.label)
pt1$tip.label<-gsub("Philodryas_patagoniensis_1-MeloSampaio","Philodryas_patagoniensis_-_PA",pt1$tip.label)


#ladderize tree
pt1<-ladderize(pt1,right=F)
#write read tree to fix node numbers
write.tree(pt1,file="Philodryadini_BEAST-LAD.tre")
pt1<-read.tree("Philodryadini_BEAST-LAD.tre")
#plot tree to check topology
# plot(pt1)


##########################################################################################

##########################################################################################
# BioGeoBEARS Setup
##########################################################################################

##############################
# set Extension data directory
##############################

extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
# extdata_dir
# list.files(extdata_dir)

##############################
# set tree file
##############################

# set tree file name
trfn=("Philodryadini_BEAST-LAD.tre")
# load tree file
tr = read.tree(trfn)

##############################
# set area file
##############################

# load area file name
geogfn = "Philodryadini_geodata_2.txt"

# set tipranges
ranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

# set tip ranges
tipranges=order_tipranges_by_tr(ranges, tr)

# set max range size
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))

# check number of sates
numstates_from_numareas(numareas=length(tipranges@df), maxareas=max_range_size, include_null_range=TRUE)

##########################################################################################

##########################################################################################
# DEC and DEC+J Analyses
##########################################################################################

#############################
# Run DEC
#############################
# 
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
# Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$min_branchlength = 0.000001
# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = TRUE

# Uncomment files you wish to use in time-stratified analyses:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
# # BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
# BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50	# returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE 	# shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE # "GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 2
BioGeoBEARS_run_object$force_sparse = FALSE	# Sparse matrix exponentiation

# loads the dispersal multiplier matrix etc. from the text files into the model object
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)

# The stratified tree is described in this table:
# BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE  # get ancestral states from optim run


# check list of settings
# BioGeoBEARS_run_object
# BioGeoBEARS_run_object$BioGeoBEARS_model_object
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

##############################

# For a slow analysis, run once, then set runslow=FALSE to just load the saved result.
runslow = TRUE
resfn = "Philo_DEC_M0.Rdata"
if (runslow) {
	res = bears_optim_run(BioGeoBEARS_run_object)
	save(res, file=resfn)
	resDEC = res
} else {
	# Loads to "res"
	load(resfn)
	resDEC = res
}

##############################

#############################
# Run DEC+J
#############################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
# BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
# BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE       #"GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)

# The stratified tree is described in this table:
# BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart


# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#############################

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
resfn = "Philo_DEC+J_M0.Rdata"
runslow = TRUE
if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resDECj = res
} else {
    # Loads to "res"
    load(resfn)
    resDECj = res
}

#############################

##########################################################################################
# PDF plots
##########################################################################################
pdffn = "Philo_DEC_vs_DEC+J_M0.pdf"
pdf(pdffn, width=12, height=20)

#############################
# Plot ancestral states - DEC
#############################
analysis_titletxt ="BioGeoBEARS DEC on Philo M0"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#############################
# Plot ancestral states - DECJ
#############################
analysis_titletxt ="BioGeoBEARS DEC+J on Philo M0"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF

#############################

##########################################################################################
# DIVALIKE and DIVALIKE+J Analyses
##########################################################################################

#############################
# Run DIVALIKE
#############################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE       # "GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#############################

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
resfn = "Philo_DIVALIKE_M0.Rdata"
runslow = TRUE
if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resDIVALIKE = res
} else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
}

##############################

##############################
# Run DIVALIKE+J
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE       # "GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

##############################
# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
resfn = "Philo_DIVALIKE+J_M0.Rdata"
runslow = TRUE
if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resDIVALIKEj = res
} else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
}

##############################

##########################################################################################
# PDF plots
##########################################################################################

pdffn = "Philo_DIVALIKE_vs_DIVALIKE+J_M0.pdf"
pdf(pdffn, width=12, height=20)

##############################
# Plot ancestral states - DIVALIKE
##############################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Amphisbaenia M0"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

##############################
# Plot ancestral states - DIVALIKE+J
##############################
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Philo M0"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

##############################

##########################################################################################
# BAYAREALIKE AND BAYAREALIKE+J analyses
##########################################################################################

##############################
# Run BAYAREALIKE
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

#  Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE       # "GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# loads the dispersal multiplier matrix etc. from the text files into the model object
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

##############################

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
resfn = "Philo_BAYAREALIKE_M0.Rdata"
runslow = TRUE
if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resBAYAREALIKE = res
} else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
}

##############################

##############################
# Run BAYAREALIKE+J
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
# BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
# #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE       # "GenSA" or if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

##############################

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
resfn = "Philo_BAYAREALIKE+J_M0.Rdata"
runslow = TRUE
if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resBAYAREALIKEj = res
} else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
}

##############################

##########################################################################################
# PDF plots
##########################################################################################

pdffn = "Philo_BAYAREALIKE_vs_BAYAREALIKE+J_M0.pdf"
pdf(pdffn, width=12, height=20)

##############################
# Plot ancestral states - BAYAREALIKE
##############################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Philodryadini M0"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

##############################
# Plot ancestral states - BAYAREALIKE+J
##############################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Philodryadini M0"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

##############################

##########################################################################################

#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
# 
#########################################################################

##########################################################################################

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

# Load results
load("Philo_DEC_M0.Rdata")
resDEC<-res
load("Philo_DEC+J_M0.Rdata")
resDECj<-res
load("Philo_DIVALIKE_M0.Rdata")
resDIVALIKE<-res
load("Philo_DIVALIKE+J_M0.Rdata")
resDIVALIKEj<-res
load("Philo_BAYAREALIKE_M0.Rdata")
resBAYAREALIKE<-res
load("Philo_BAYAREALIKE+J_M0.Rdata")
resBAYAREALIKEj<-res


#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################
# Loads to "restable"
save(restable, file="restable_v1.Rdata")
load(file="restable_v1.Rdata")

# Loads to "teststable"
save(teststable, file="teststable_v1.Rdata")
load(file="teststable_v1.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")

#######################################################

##########################################################################################

#########################################################################
# 
# PLOT EDITED TREE WITH MOST PROBABLE STATES
# I don't known the origen of the code below... but I've modified it to my needs and taste
# 
#########################################################################

##########################################################################################

#####################
# Prepare infos from results and tree
#####################

# read least table with AICc results
cond_restable_AICc_rellike<-read.csv("restable_AICc_rellike_formatted.txt", sep="\t")

# select best model
best_model<-rownames(cond_restable_AICc_rellike)[cond_restable_AICc_rellike[,1]==max(cond_restable_AICc_rellike[,1])]


if(length(best_model)>1){
	best_model<-best_model[1]
}

# load best model
load(file=paste("Philo_",best_model,"_M0.Rdata",sep=""))

# load tree file
tr = read.tree(res$inputs$trfn)

# set tip and node numbers
tips<-seq(Ntip(tr))
nodes<-seq(Ntip(tr)+1,length.out=Nnode(tr))

# set ranges
ranges = getranges_from_LagrangePHYLIP(lgdata_fn=res$inputs$geogfn)

# set ordered tip ranges
tipranges=order_tipranges_by_tr(ranges, tr)

# set max range size
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))

# set area names based on tiprange
areanames=getareas_from_tipranges_object(tipranges)

# set max range size
states_list=areas_list_to_states_list_new(areanames)

# set area names based on max range
names_areas<-unlist(areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range = TRUE, split_ABC = FALSE))

# get state probabilities for tips and nodes
relprobs_matrix<-res$ML_marginal_prob_each_state_at_branch_top_AT_node

# set probabilities for nodes
relprobs_matrix_for_internal_states <-relprobs_matrix[nodes,]


#####################
# Gather 1st, 2nd and 3rd more probable states for nodes
#####################

#First most probable state
whichmax<-apply(relprobs_matrix[nodes,],1,function(x) which(x==max(x)))
first_prob<-apply(relprobs_matrix[nodes,],1,function(x) max(x))
first_area<-names_areas[whichmax]#Names for the most probable states
# write.csv(first_max, "first_max_states.csv")

#Second most probable state
whichsecondmax<-apply(relprobs_matrix[nodes,],1,function(x) which(x==sort(x, TRUE)[2]))
second_prob<-apply(relprobs_matrix[nodes,],1,function(x) sort(x, TRUE)[2])
second_area<-names_areas[whichsecondmax]
# write.csv(second_max, "second_max_states.csv")

#Third most probable states
whichthirdmax<-apply(relprobs_matrix[nodes,],1,function(x) which(x==sort(x, TRUE)[3]))
third_prob<-apply(relprobs_matrix[nodes,],1,function(x) sort(x, TRUE)[3])
third_area<-names_areas[whichthirdmax]
# write.csv(third_max, "third_max.csv")

#####################
# Combine and remove states with probability lower than 0.05
#####################

temp_data_states<-cbind(s1=first_area,s2=second_area,s3=third_area)

temp_data_probs<-cbind(p1=first_prob,p2=second_prob,p3=third_prob)

for(i in 1:3){
    for(j in seq(temp_data_probs[,i])){
        if(temp_data_probs[j,i]<0.05){
            temp_data_probs[j,i]<-0
            temp_data_states[j,i]<-"XXXX"
        }
    }
}

#####################
# Fill prop info for the rest of the states combined
#####################

fourth_prob<-apply(temp_data_probs,1,function(x) 1-sum(x))

fourth_prob[fourth_prob<0.0000001]<-0

temp_data_probs<-cbind(temp_data_probs,p4=fourth_prob)

temp_data_states<-cbind(temp_data_states,s4=rep("AllOthers",Nnode(tr)))


#####################
# Combine all info for states, probs and nodes in a matrix
#####################

# To identify how many states you have to prepare a legend containing all states 
statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=TRUE, split_ABC=FALSE)

MLstates =  get_ML_states_from_relprobs(relprobs_matrix, names_areas, returnwhat="states", if_ties="takefirst")

unique_states<-unique(MLstates)

all_states<-unlist(statenames[statenames%in%unique_states])

all_states<-c(all_states,"AllOthers")

tot_matrix<-relprobs_matrix[,which(names_areas%in%all_states)]

names_col<-c(names_areas[names_areas%in%all_states],"AllOthers")

for(i in seq(tot_matrix[1,])){
    for(j in seq(tot_matrix[,i])){
        if(tot_matrix[j,i]<0.05){
            tot_matrix[j,i]<-0
        }
    }
}

fourth_toAdd<-c(rep(0,Ntip(tr)),fourth_prob)

tot_matrix<-cbind(tot_matrix,fourth_toAdd)

colnames(tot_matrix)<-names_col

#####################
# set colors for states

colors_matrix = get_colors_for_numareas(length(areanames))

colnames(colors_matrix)<-areanames

write.table(colors_matrix, file="state_color-toSet.txt", quote=FALSE, sep="\t")

# generated colors matrix

# 		AT	CA	CE	CC	PP	NC	SA	AM
# red	255	255	128	0	0	0	128	255
# green	0	191	255	255	255	64	0	0
# blue	0	0	0	64	255	255	255	191

# manually changed colors matrix
#     NC  ME  AT  PA  AL  AM  OD  AF  FN
# red 0   0   255 0   255 0   255 170 170
# green   0   255 0   170 0   255 170 255 0
# blue    255 170 170 255 0   0   0   0   255


colors_matrix2<-read.csv("state_colors.csv",header=T)
colnames(colors_matrix2)<-seq(colors_matrix2[1,])
rownames(colors_matrix2)<-c("red","green","blue")


states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=TRUE)

selected_states<-states_list_0based_index[which(statenames%in%all_states)]

colors_list_for_states2 <- mix_colors_for_states(colors_matrix2, selected_states)

colors_list_for_states2 <- c(colors_list_for_states2,"#FFFFFF")


cols_byNode = rangestxt_to_colors(all_states, colors_list_for_states2, MLstates)


#####################
# Plot tree
#####################

tree_b<-ladderize(tr)
tree_b$root.time<-max(node.depth.edgelength(tr)) 
lim<-tree_b$root.time

# pdffn = "Amphis_DECj_8areas_pretty.pdf"
pdffn = "Philo_BioGeoBEARS_bestResult.pdf"

pdf(pdffn, width=6, height=8)
par(oma=c(2,0,0,0))
par(mai=c(0,0,0,0))
par(mar=c(0,0,0,0))

#plot(tree_b,  no.margin=F,cex=1.2,adj=0, label.offset=2.5) , y.lim= 62
geoscalePhylo(tree=tree_b, units=c("Period", "Epoch", "Age"), boxes="Age", quat.rm=F, cex.age=0.6, cex.ts=0.6, x.lim=c(lim,-15), show.tip.label=TRUE, type= "p", node.pos=1,vers="ICS2013", tick.scale="myr", lwd=1, width=1, cex.tip = 0.6)
nodelabels(node=nodes, pie=tot_matrix[nodes,],cex=0.6, cex.lab=0.3, piecol=colors_list_for_states2)
tiplabels(text=MLstates[tips], tip=tips, bg=cols_byNode[tips], cex=0.4)
#axisPhylo2(side=1, tick=TRUE,las=1, par(cex=2), roundlabels=FALSE, minage=0, pos=0.1, padj=0.8)

dev.off()

#######################################################

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


