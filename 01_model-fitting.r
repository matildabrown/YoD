##############################################################################
#  Model fitting for 'More than 75% of undescribed plant species
#    are likely threatened'
#  MJM Brown, SP Bachman, E Nic Lughadha
##############################################################################

library(tidyverse)
library(brms)
library(tidybayes)

library(phylolm)
library(phytools)
library(caper)


df <- read.csv("df_fit_inclNA.csv")
df_fit <- read.csv("df_fit.csv")

# model fitting

if(!dir.exists("models")) dir.create("models")

fit_and_loo <- function(formula, name, iter.brm=2000) {
  model <- brms::brm(formula, data=df_fit,
                       family = bernoulli(link = "logit"),cores=4,
                       prior = priors,
                       chains=4, iter = iter.brm)
  model <- add_criterion(model, "loo")
  saveRDS(model, file = paste0("models/", name,".rds"))
  return(model)
}

# formulae

model_formulae=c(
  "year" = threatened ~ year,
  "year + climate" =  threatened ~ year*clim_bin,
  "year + lifeform" = threatened ~ year*standard_lifeform,
  "year + lifeform + climate" = threatened ~ year*clim_bin + year*standard_lifeform
)

priors <- c(
  prior(normal(0, 10), class="b"),
  prior(normal(0, 0.05), class="b", coef="year"),
  prior(normal(-35, 4), class="Intercept")
)

# library(doParallel)
# registerDoParallel(cl <- makeCluster(4))
# results_list <- foreach(i = 1:4, .packages = c("brms")) %dopar% {
#
#   name <- gsub( " \\+ ","-", names(model_formulae)[i])
#   formula <- unname(model_formulae[i])[[1]]
#
#   fit_and_loo(formula=formula, name=name, iter.brm=8000)
# }
# stopCluster(cl)


 # if r crashes and you have to load them in again:
results_list <- list()
for (i in 1:4){
  results_list[[i]] <- readRDS(paste0("models/", gsub( " \\+ ","-", names(model_formulae)[i]),".rds"))
  names(results_list)[i] <-  gsub( " \\+ ","-", names(model_formulae)[i])
  }

# check summaries, identify non-converging models
# increased iterations if non-convergence present
summary(results_list[[1]])
summary(results_list[[2]])
summary(results_list[[3]])
summary(results_list[[4]])

#choose best model based on PSIS-LOO

LOOlist <- list()
for (i in 1:4){
  LOOlist[[length(LOOlist)+1]] <- results_list[[i]]$criteria$loo
  names(LOOlist)[length(LOOlist)] <-  gsub( " \\+ ","-", names(model_formulae)[i])
}

loo_compare(LOOlist)
#chose best model as simplest because coeffs are similar across all
best_model <- results_list[[1]]
summary(best_model)

## dose dependency

model_formulae_ENCR=c(
  "year_EN" = threatenedEN ~ year,
  "year_CR" =  threatenedCR ~ year
)
mods_VUENCR <- list()
mods_VUENCR[[1]] <- best_model

# fit the models # comment out when fitted and saved
# for (i in 1:2){
#   name <- gsub( " \\+ ","-", names(model_formulae_ENCR)[i])
#   formula <- unname(model_formulae_ENCR[i])[[1]]
#   mods_VUENCR[[i+1]] <- fit_and_loo(formula=formula, name=name)
# }

#read them back in if already fit
for (i in 1:2){
  mods_VUENCR[[i+1]] <- readRDS(paste0("models/", gsub( " \\+ ","-", names(model_formulae_ENCR)[i]),".rds"))
  names(mods_VUENCR)[i+1] <-  gsub( " \\+ ","-", names(model_formulae_ENCR)[i])
}

endnames <- c(
  "year",
  "extrap",
  "VUepred",
  "VUlower",
  "VUupper",
  "ENepred",
  "ENlower",
  "ENupper",
  ".row",
  "CRepred",
  "CRlower",
  "CRupper",
  ".width",
  ".point",
  ".interval"
)

totals <- df_fit %>% group_by(year) %>%
  summarise(n=n())

obs <- df_fit %>% group_by(year, redlistCategory) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = "redlistCategory", values_from="n", values_fill = 0)

obs[,2:8] <- obs[,2:8]/totals$n

obs <- obs %>%
  pivot_longer(cols=2:8, names_to="Category", values_to = "Percent") %>%
  filter(year<2021) %>%
  mutate(Category=factor(Category, levels=c("Least Concern",
                                                   "Near Threatened",
                                                   "Data Deficient",
                                                   "Vulnerable",
                                                   "Endangered",
                                                   "Critically Endangered",
                                                   "Extinct in the Wild",
                                                   "Extinct")))



preds_plot <- df_fit %>%
  modelr::data_grid(year=seq(from=1753, to=2050)) %>%
  mutate(extrap=year>2020) %>%
  add_epred_draws(best_model) %>%
  median_qi(.epred) %>%
  dplyr::select(-c(.row, .width, .point,.interval)) %>%
  `colnames<-`(endnames[1:5]) %>%
  add_epred_draws(mods_VUENCR[[2]]) %>%
  median_qi(.epred) %>%
  dplyr::select(-c(.row, .width, .point,.interval)) %>%
  `colnames<-`(endnames[1:8]) %>%
  add_epred_draws(mods_VUENCR[[3]]) %>%
  median_qi(.epred) %>%
  `colnames<-`(endnames[1:15])

pal_cat <- c("Least Concern"="#8ccba7",
             "Near Threatened" = "#c8e7af",
             "Data Deficient" = "gray70",
             "Vulnerable" =  "#fee3ba",
             "Endangered" = "#f9b6a1",
             "Critically Endangered" = "#eb9793",
             "Extinct in the Wild" = "#a991a1",
             "Extinct" = "#7f7f7f")


  preds_plot %>%
    ggplot(aes(x=year))+
    # geom_hline(aes(yintercept=0.5), colour="gray60")+
    # geom_hline(aes(yintercept=0.75), colour="gray60")+
    # geom_hline(aes(yintercept=0.90), colour="gray60")+
    # geom_hline(aes(yintercept=0.95), colour="gray60")+
     geom_bar(data=obs,
             aes(y=Percent, fill = Category),
             position = position_stack(),
             stat = "identity", width = 1, alpha=1)+
    geom_line(aes(y=VUepred), colour="#FEC776", linewidth=0.6)+
    geom_ribbon(aes(ymin=VUlower, ymax=VUupper), fill="#FEC776",alpha=0.6)+
    geom_line(aes(y=ENepred), colour="#F46D43", linewidth=0.6)+
    geom_ribbon(aes(ymin=ENlower, ymax=ENupper), fill="#F46D43",alpha=0.6)+
    geom_line(aes(y=CRepred), colour="#D73027", linewidth=0.6)+
    geom_ribbon(aes(ymin=CRlower, ymax=CRupper), fill="#D73027",alpha=0.6)+
  theme_ggdist()+
   geom_vline(aes(xintercept=2020), linetype="dashed")+
    scale_fill_manual(values=pal_cat)+
    coord_cartesian(expand=FALSE)+
    ylab("Predicted probability threatened")+
    xlab("Year of description")

  ggsave("Fig1_yod.pdf", height=11, width =19, units="cm")

### PHYLOGENETIC ANALYSIS ####
## phyloglm

alltrees_red <- read.tree("ANGIOtrees100.tre")  # from F. Forest, unpubl.
goodtips <- gsub(" ","_", df_fit$taxon_name)
setdiff(alltrees_red[[1]]$tip.label, goodtips)
alltrees_red <- drop.tip(alltrees_red, setdiff(alltrees_red[[1]]$tip.label, goodtips))


alltrees_red[[1]]$edge.length<-
  alltrees_red[[1]]$edge.length/max(nodeHeights(alltrees_red[[1]])[,2])


df_phylo <- df_fit
rownames(df_phylo) <- goodtips
df_phylo <- df_phylo[rownames(df_phylo) %in% alltrees_red[[1]]$tip.label,]

saveRDS(alltrees_red, "redlist_phylo.rds")

df_phylo <- df_phylo[alltrees_red[[1]]$tip.label,]

#only modelled on first of 100 trees to ensure convergence
phylo_fit <- phyloglm(formula=threatened~year, data=df_phylo,
                      phy=alltrees_red[[1]], btol=65, log.alpha.bound = 2)

# Warning messages:
#   1: In phyloglm(formula = threatened ~ year, data = df_phylo, phy = alltrees_red[[1]],  :
#      the boundary of the linear predictor has been reached during the optimization procedure.
#      You can increase this bound by increasing 'btol'.
#  2: In phyloglm(formula = threatened ~ year, data = df_phylo, phy = alltrees_red[[1]],  :
#      phyloglm failed to converge.

## Increasing btol gave following: Error in transf.branch.lengths(par[1:dk], par[dk + 1]) :
#  edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.

#check anyway
summary(phylo_fit)


### EXENDED DATA FIG 1 ####
#collated and partitioned from Extended Data Table 1
model_coeffs <- read.csv("model_coefficients.csv")

ggplot(model_coeffs, aes(y=interaction(lifeform, climate), x=est, group=interaction(lifeform, climate, model), colour=factor(model)))+
  geom_point(position = position_dodge(width=-0.5))+
  facet_grid(~type, scales="free")+
  geom_errorbar(aes(xmin=lower, xmax=upper),
                position = position_dodge(width=-0.5),
                width=0.5)+
  theme_ggdist()+
  scale_y_discrete(limits=rev)+
  ylab("Partition")+
  xlab("Coefficient estimate")

ggsave("ED_fig1_coeffs.pdf",height=11, width =19, units="cm")

### EXENDED DATA FIG 2 ####
# tropical species more likely to be assessed?
df %>%
  filter(year<2021, !is.na(clim_bin)) %>%
  mutate(assessed=!is.na(redlistCategory)) %>%
  group_by(year, assessed, clim_bin) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=year,y=n, group=assessed, fill=assessed))+
  geom_bar(position = position_stack(),
           stat = "identity", width = 1)+
 facet_grid(~clim_bin)+
  coord_cartesian(expand=FALSE, ylim=c(0,2500))+
  theme_ggdist()+
  ylab("Number of species")+
  xlab("Year of description")+
  scale_fill_manual(values=c("gray80","gray60","gray10"))

ggsave("ED_fig2A_clim-ass-yod.pdf", height=8, width =12, units="cm")

# by proportion?

df %>%
  filter(year<2021, !is.na(clim_bin)) %>%
  mutate(assessed=!is.na(redlistCategory)) %>%
  group_by(year, assessed, clim_bin, threatened) %>%
  summarise(n=n()) %>%
  group_by(year, clim_bin) %>%
  mutate(pc.ass =  100 *n/sum(n)) %>%
  ungroup() %>%
  filter(assessed==TRUE) %>%
  ggplot(aes(x=year,y=pc.ass, group=threatened, fill=threatened))+
  geom_bar(position = position_stack(),
           stat = "identity", width = 1)+
  facet_grid(~clim_bin)+
  coord_cartesian(expand=FALSE)+
  theme_ggdist()+
  ylab("Percent of species")+
  xlab("Year of description")+
  scale_fill_manual(values=c("gray80","gray60"))

ggsave("ED_fig2B_clim-ass-yod.pdf", height=8, width =12, units="cm")

## bias plots: are we assessing mostly new species?
### EXENDED DATA FIG 3 ####
df %>%
  filter(year<2021) %>%
  mutate(assessed=!is.na(redlistCategory)) %>%
  group_by(year, assessed) %>%
  summarise(n=n()) %>%
  group_by(year) %>%
  mutate(pc.ass =  100 *n/sum(n)) %>%
  ungroup() %>%
  filter(assessed==TRUE) %>%
  ggplot(aes(x=year,y=pc.ass))+
  geom_bar(position = position_stack(),
           stat = "identity", width = 1, fill="gray80")+
  geom_smooth(aes(x=year, y=pc.ass), colour="black", linewidth=0.3, method="gam")+
  coord_cartesian(expand=FALSE)+
  theme_ggdist()+
  ylab("Percent assessed")+
  xlab("Year of description")

ggsave("ED_fig3_ass-yod.pdf", height=8, width =12, units="cm")

### EXENDED DATA FIG 4 ####
# endemism

df %>%
  filter(year<2021) %>%
  group_by(year, endemic) %>%
  summarise(n=n()) %>%
  group_by(year) %>%
  mutate(pc.end =  100 *n/sum(n)) %>%
  ungroup() %>%
  filter(endemic==TRUE) %>%
  ggplot(aes(x=year,y=pc.end))+
  geom_bar(position = position_stack(),
           stat = "identity", width = 1, fill="gray80")+
  geom_smooth(aes(x=year, y=pc.end), colour="black", linewidth=0.3, method="gam")+
  geom_hline(yintercept=56, colour="red", linetype="dashed", linewidth=0.4)+
  coord_cartesian(expand=FALSE, ylim=c(0,100))+
  theme_ggdist()+
  ylab("Percent endemic to a single WGSRPD Area")+
  xlab("Year of description")

ggsave("ED_fig4_endem-yod.pdf", height=8, width =12, units="cm")


