# ANALYSIS: ABUNDANCE ####
## Age ratios ####
ratio <- function(num, denom){
  sum(num, na.rm=T) / sum(denom, na.rm=T)
}
ratioSE <- function(num, denom, reps=9999){
  n <- length(num)
  if(n != length(denom)) stop("num and denom must have equal length")
  rats <- sample(1:n, n*reps, replace=TRUE) %>%
    matrix(nrow=reps) %>%
    apply(1, function(i) ratio(num[i], denom[i]))
  sd(rats)
}
ageRatios <- ctdpbound$data$observations %>%
  group_by(sequenceID, scientificName, year, lifeStage) %>%
  summarise(count = sum(count)) %>% 
  pivot_wider(names_from = lifeStage, values_from = count) %>%
  rename(unknown = 'NA') %>%
  mutate(adult = replace_na(adult, 0),
         subadult = replace_na(subadult, 0),
         juvenile = replace_na(juvenile, 0),
         unknown = replace_na(unknown, 0),
         nonjuv = sum(adult, subadult),
         aged = sum(nonjuv, juvenile),
         total = sum(aged + unknown)) %>%
  group_by(scientificName, year) %>%
  summarise(pJuv = ratio(juvenile, aged),
            pJuvSE = ratioSE(juvenile, aged),
            jpa = ratio(juvenile, adult),
            jpaSE = ratio(juvenile, adult),
            pAged = ratio(aged, total))
View(ageRatios)

ratioPlot <- ggplot(ageRatios, aes(year, jpa, col=scientificName)) +
  geom_point() +
  geom_line() +
  #  geom_errorbar(aes(x=year, ymin=jpa-1.96*jpaSE, ymax=jpa+1.96*jpaSE), width=0.2) + 
  labs(y="Juveniles per adult") +
  theme_classic()

agedPlot <- ggplot(ageRatios, aes(year, pAged, col=scientificName)) +
  geom_point() +
  geom_line() +
  labs(y="Proportion aged") +
  theme_classic()

plot_grid(ratioPlot, agedPlot, nrow=2)


## Trap rate data ####
# Function extracts traprate data for all species from one annual slice
trd <- function(dp){
  lapply(spp, function(sp) get_traprate_data(dp, sp)) %>%
    bind_rows() %>%
    mutate(year = lubridate::year(dp$data$deployments$start[1]))
}
# Trap rate data for all species and years
trdat <- lapply(ctdpslices, trd) %>%
  bind_rows() %>%
  pivot_wider(names_from=scientificName, values_from=n)
View(trdat)
# crude trap rates, not including zeros
trdat2 <- ctdpbound$data$observations %>%
  left_join(select(ctdpbound$data$deployments, deploymentID, locationName), by="deploymentID") %>%
  group_by(scientificName, locationName, year) %>%
  summarise(sumcount = sum(count),
            nobs = n()) %>%
  pivot_wider(names_from = scientificName, values_from = c(sumcount, nobs)) %>% 
  replace(is.na(.), 0)
View(trdat2)

## Compare PAJ v Agouti trap rates ####
check <- deployments %>%
  left_join(trdat2, by=c("year", "locationName"))
View(check)
names(check)
plot(check$roe_deer, check$'sumcount_Capreolus capreolus')
plot(check$red_deer, check$'sumcount_Cervus elaphus')
plot(check$mouflon, check$'sumcount_Ovis ammon')
plot(check$wild_boar, check$'sumcount_Sus scrofa')
plot(check$roe_deer, check$'nobs_Capreolus capreolus')
plot(check$red_deer, check$'nobs_Cervus elaphus')
plot(check$mouflon, check$'nobs_Ovis ammon')
plot(check$wild_boar, check$'nobs_Sus scrofa')
lines(c(0,2000), c(0,2000))
# PAJ deployments table counts observations, not summed counts for trap rate


## Group count corrections ####
dat <- read.csv("./Hoge Veluwe data/DetectionZoneGroupCounts.csv") %>%
  mutate(countCor = ifelse(count>=countDetZone, count, countDetZone),
         sp = as.factor(scientificName),
         nm = as.numeric(as.factor(scientificName)))
View(dat)
n <- nrow(dat)
lim <- range(dat[, c("count", "countDetZone")], na.rm=TRUE)
spp <- levels(dat$sp)

calc_p <- function(num, den, reps=9999){
  est <- function(nm, dn) 100 * sum(nm, na.rm = TRUE) / sum(dn)
  rep <- function(){
    i <- sample(1:n, n, replace=TRUE)
    est(num[i], den[i])
  }
  n <- length(num)
  rnd <- replicate(reps, rep())
  data.frame(est=est(num, den), se=sd(rnd))
}

dzp <- dat %>%
  group_by(scientificName, habitatType) %>%
  reframe(calc_p(countDetZone, count))

ggplot(dat, aes(x = count, y = countDetZone)) +
  geom_point() +
  geom_jitter(height=0.5, width=0.5) +
  geom_line(data=data.frame(lim=lim), aes(x = lim, y = lim))


ggplot(dzp, aes(x = scientificName, y = est, fill = habitatType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = est - 1.96*se, ymax = est + 1.96*se), 
                width = 0.2, position = position_dodge(0.9)) +
  coord_cartesian(ylim=c(50, 100)) +
  labs(x = "Species", y = "% in detection zone", fill = "Habitat type") +
  scale_fill_manual(values = c("closed" = "forestgreen", "open" = "khaki")) +
  theme_classic()

mod <- glm(cbind(countDetZone, countCor-countDetZone) ~ sp + habitatType + countCor, 
           family=quasibinomial, data = dat)
summary(mod)
nd <- expand.grid(sp = spp[2:3], habitatType=unique(dat$habitatType), countCor=100)
prdn <- predict(mod, type="resp", newdata = nd, se.fit=TRUE)
nd$p <- prdn$fit
nd$se <- prdn$se.fit
dzp
nd

sq <- seq(lim[1], lim[2], len=100)
nd <- data.frame(countCor=sq, sp="Cervus elaphus", habitatType="open")
nd$est <- predict(mod, type="resp", newdata = nd)
ggplot(data=dat, aes(x = countCor, y = countDetZone / countCor)) +
  geom_point() +
  geom_jitter(height=0.02) +
  geom_line(data=nd, aes(x = countCor, y = est)) +
  theme_classic()


## Density analysis ####

# Function gets an REM estimate table for a given year and species
yr <- 2020
sp <- "mouflon"
get_estimates <- function(yr, sp, reps=999){
  # speed, radius, angle tables
  sdat <- spd_dat %>%
    subset(species==sp) %>%
    group_by(stratumID) %>%
    summarise(estimate = mean(est),
              se = se_from_ses(est, se)) %>%
    mutate(unit = "km/day",
           param = "speed")
  rdat <- rad_dat %>%
    ungroup() %>%
    subset(species==sp) %>%
    mutate(unit = "km",
           param = "radius") %>%
    select(stratumID, param, estimate=est, se, unit)
  adat <- ang_dat %>%
    group_by(stratumID) %>%
    summarise(estimate = mean(est),
              se = se_from_ses(est, se)) %>%
    mutate(unit = "radian", param="angle")
  if(yr<2016){
    amod <- actmods_14_15
    sdat <- subset(sdat, !grepl("O", stratumID))
    rdat <- subset(rdat, !grepl("O", stratumID))
    adat <- subset(adat, !grepl("O", stratumID))
    str <- subset(strata, !grepl("O", stratumID))
  } else{
    amod <- actmods_16_20
    str <- strata
  }
  
  # activity table
  pdat <- data.frame(stratumID = names(amod[[sp]]$est$act_stratum),
                     param = "activity",
                     estimate = amod[[sp]]$est$act_stratum,
                     se = amod[[sp]]$se$act_stratum,
                     unit = "none")
  
  # traprate table
  trd <- deployments %>%
    filter(year == yr) %>%
    select(deploymentID, stratumID, effort, n=all_of(sp))
  # stratum-specific traprate function
  trdfunc <- function(trd){
    trd %>%
      group_by(stratumID) %>%
      summarise(estimate = sum(n) / sum(effort))
  }
  # bootstrapped stratum-specific trprate function
  trdboot <- function(trd){
    trd %>%
      group_by(stratumID) %>%
      slice_sample(prop=1, replace=TRUE) %>%
      trdfunc() %>%
      select(estimate) %>%
      pull()
  }
  trdboots <- replicate(reps, trdboot(trd))
  tdat <- trdfunc(trd) %>%
    mutate(se = apply(trdboots, 1, sd),
           param = "trap_rate",
           unit = "n/day")
  
  # combine tables with added density and abundance
  prm_tab <- bind_rows(tdat, pdat, sdat, rdat, adat)
  res <- rem_strat(prm_tab, str)
  areas <- c(sum(str$area), str$area) / 100
  res %>%
    filter(param == "density") %>%
    mutate(estimate = estimate * areas,
           se = se * areas,
           param = "abundance",
           unit = "n") %>%
    bind_rows(res) %>%
    arrange(stratumID)
}

reptab <- expand.grid(yr=2014:2020, sp=unique(spd_dat$species))
fitone <- function(i)
  get_estimates(yr=reptab[i,"yr"], sp=reptab[i,"sp"])
rem_list <- sapply(1:nrow(reptab), fitone, simplify=FALSE)
rowcount <- unlist(lapply(rem_list, nrow))
addcols <- reptab[rep(1:nrow(reptab), rowcount), ]
rem_estimates <- rem_list %>%
  bind_rows() %>%
  bind_cols(addcols)

count_data <- read.csv("./Hoge Veluwe data/Ungulate Counts Hoge Veluwe National Park 2012-2023.csv") %>%
  filter(Year>=2014 & Year<=2020) %>%
  rename(year = Year,
         red_deer = Red_deer,
         roe_deer = Roe_deer,
         wild_boar = Wild_boar,
         mouflon = Mouflon_excl_lambs)
sp="roe_deer"
for(sp in unique(rem_estimates$sp)){
  counts <- count_data %>%
    rename(count = all_of(sp)) %>%
    select(year, count)
  dat <- rem_estimates %>%
    rename(species = sp,
           year = yr) %>%
    filter(param=="abundance" & stratumID=="stratified" & species==sp) %>%
    select(year, estimate, se) %>%
    left_join(counts, by = "year") %>%
    pivot_longer(c(2,4)) %>%
    mutate(se = ifelse(name=="count", NA, se))
  
  plt <- ggplot(dat, aes(x=year, y=value, col=name)) +
    geom_point(size=2) + 
    geom_line() +
    geom_errorbar(aes(ymin=value-1.96*se, ymax=value+1.96*se),
                  width=0.2) +
    ylim(0, 510) +
    theme_bw() +
    ggtitle(sp) +
    theme(legend.title = element_blank())
  print(plt)
}
