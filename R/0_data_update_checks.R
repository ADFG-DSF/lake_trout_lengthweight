## This script is intended as another layer of algorithmic check whenever the
## lake trout database is updated, just to make sure the data is structured
## as expected, and that any additions make sense.
##
## Checks are performed on the raw (not R-edited) data, as a check that filtration
## functions etc can be used without changes.

library(tidyverse)

# old & new morphometry datasets
morph_old <- read_csv("flat_data/lake_morphometry3.csv", skip=1) %>%
  filter(!is.na(LakeName))
morph_new <- read_csv("flat_data/lake_morphometry_25_12_23.csv", skip=1) %>%
  filter(!is.na(LakeName))

# old & new length/weight/age datasets
lw_old <- read_csv("flat_data/length_weight4.csv", skip=1) %>%
  filter(!is.na(LakeName))
lw_new <- read_csv("flat_data/length_weight_25_12_23.csv", skip=1) %>%
  filter(!is.na(LakeName))


## ------------------ Morphometry ------------------ ##

# make sure columns are the same
dim(morph_old)
dim(morph_new)

colnames(morph_old)
colnames(morph_new)

# any new columns?
colnames(morph_new)[!(colnames(morph_new) %in% colnames(morph_old))]

# any old columns went away?
colnames(morph_old)[!(colnames(morph_old) %in% colnames(morph_new))]

# same order? doesn't matter too much but makes checks easier
all(colnames(morph_old) == colnames(morph_new))

# if so, make sure the columns are the same class
all(sapply(morph_old, class) == sapply(morph_new, class))

# this is slightly more robust if the column order changed
all(sapply(morph_old, class)[order(names(sapply(morph_old, class)))] ==
      sapply(morph_new, class)[order(names(sapply(morph_new, class)))])

# any new lakes?
morph_new$LakeName[!(morph_new$LakeName %in% morph_old$LakeName)]
# [1] "Chikuminuk Lake" "Nayakuk Lake"    "Tikchik Lake"

# any old lakes went away?
morph_old$LakeName[!(morph_old$LakeName %in% morph_new$LakeName)]
# [1] "Chelle Lake"                              "South Middle Fork Lake of Goodnews River"

# any changes to values?
old_a2a <- morph_old[order(morph_old$LakeName), ] %>%
  filter(LakeName %in% intersect(morph_old$LakeName, morph_new$LakeName))
new_a2a <- morph_new[order(morph_new$LakeName), ]%>%
  filter(LakeName %in% intersect(morph_old$LakeName, morph_new$LakeName))
for(j in names(old_a2a)) {
  check1 <- xor(is.na(old_a2a[j]), is.na(new_a2a[j]))
  check2 <- ifelse(is.na(old_a2a[j] != new_a2a[j]), FALSE, old_a2a[j] != new_a2a[j])
  if(any(check1 | check2)) {
    cat("\n","\n")
    # cat(j)
    print(data.frame(name=old_a2a$LakeName,
                     old=old_a2a[j],
                   new=new_a2a[j])[check1 | check2,])
  }

}





##### look for possibly duplicated lake names (of lakes we care about)
morph_new %>%
  mutate(Model = `Include in "Alaskanizing" Modeling Exercise` %in% c("yes","Yes"))%>%
  mutate(Management = `Potentially Include in Lake Trout Management Plan` != "No") %>%
  mutate(Sample_data = LakeName %in% lw_new$LakeName) %>%
  select(LakeName,  Model, Management, Sample_data) %>%
  filter(Model | Management | Sample_data) %>%
  filter(Sample_data) %>%
  data.frame



### Actually, there are several duplicate-ish named lakes
### Here is a quick function to map them

# identifying which lakes share the same first word in the name string...
firsts <- sapply(strsplit(morph_new$LakeName, split=" "), "[", 1)
seconds <- sapply(strsplit(morph_new$LakeName, split=" "), "[", 2)
modifiers <- c("Big","Lake","Little","Lower","North","South","Upper")
fullfirsts <- ifelse(firsts %in% modifiers, paste(firsts, seconds), firsts)

sums <- sapply(fullfirsts, \(x) sum(fullfirsts==x))
morethan1 <- unique(names(sums)[sums>1])

#### - this section is commented out but results are printed

# # will use this package for interactive mapping
# library(leaflet)
#
# # not the most robust function as it pulls things from the global environment
# check_duplicate <- function(x) {
#   print(morph_new$LakeName[fullfirsts==morethan1[x]])
#   morph_new %>%
#     filter(fullfirsts==morethan1[x]) %>%
#     leaflet() %>%
#     # addTiles() %>%   # default from open street map
#     addProviderTiles("Esri.WorldImagery") %>%  # or satellite imagery!
#     addMarkers(lng=~Longitude_WGS84,
#                lat=~Latitude_WGS84,
#                label=~LakeName)
# }
#
# # how many lakes to check?
# length(morethan1)  # 17
#
# check_duplicate(1)
# # [1] "Beaver Lake (near Chisana)"   "Beaver Lake (near Susitna L)"
# # two different lakes on map
#
# check_duplicate(2)
# # [1] "Big Lake (Bob Johnson Lake)" "Big Lake (Healy)"
# # ONLY ONE LAKE HAS COORDS - Big Lake (Healy)
#
# check_duplicate(3)
# # [1] "Caribou Lake (Cantwell)"    "Caribou Lake (Lake Louise)"
# # two different lakes on map
#
# check_duplicate(4)
# # [1] "Ernie Lake" "Ernie Lake"
# # ONLY ONE LAKE HAS COORDS
#
# check_duplicate(5)
# # [1] "Fish Lake (Chandler)" "Fish Lake (Mt Hayes)"
# # two different lakes on map
#
# check_duplicate(6)
# # [1] "Ghost Lake"                "Ghost Lake (Meadows Road)"
# # two different lakes on map
#
# check_duplicate(7)
# # [1] "Island Lake (Glennallen)"
# # [2] "Island Lake (Philip Smith Mts)"
# # [3] "Island Lake (Richardson Highway, S of Donnelly Dome)"
# # three different lakes on map
#
# check_duplicate(8)
# # [1] "Kisaralik Lake"    "Kisaralik Lake II"
# # two different lakes on map
#
# check_duplicate(9)
# # [1] "Long Lake (McCarthy Rd)"             "Long Lake (Nabesna Road near Slana)"
# # two different lakes on map
#
# check_duplicate(10)
# # [1] "Lost Bob Lake"                                "Lost Haller Lake"
# # [3] "Lost Lake  (Chisholm Lake) (near Birch Lake)"
# # three dots, but Lost Bob and Lost Haller are very close to one another
# # and without lakes on the base layer
# # but a bajillion lakes on the satellite imagery
#
# check_duplicate(11)
# # [1] "Octopus Lake"              "Octopus Lake (Denali Hwy)"
# # two different lakes on map
#
# check_duplicate(12)
# # [1] "Rock Lake" "Rock Lake"
# # ONLY ONE HAS COORDS
#
# check_duplicate(13)
# # [1] "Round Lake (Chandler)" "Round Tangle Lake"
# # two different lakes on map
#
# check_duplicate(14)
# # [1] "Sevenmile Lake"                           "Sevenmile Lake (Denali Hwy)"
# # [3] "Sevenmile Lake (SE of Anderson, no fish)"
# # TWO LAKES ON MAP, ONE (Sevenmile Lake (SE of Anderson, no fish)) HAS NO COORDS
#
# check_duplicate(15)
# # [1] "Summit Lake (Parks Hwy)"                      "Summit Lake (Richardson Highway near Paxson)"
# # two different lakes on map
#
# check_duplicate(16)
# # [1] "Swan Lake (by Lake Louise)" "Swan Lake (Cooper Landing)"
# # two different lakes on map
#
# check_duplicate(17)
# # [1] "Unnamed (F-26)"                   "Unnamed (Philip Smith Mountains)"
# # two different lakes on map

# check differences between those with and without coordinates
morph_new[fullfirsts %in% c("Ernie", "Rock", "Sevenmile", "Big Lake"),] %>%
  filter(is.na(Longitude_WGS84)) %>%
  data.frame

morph_new[fullfirsts %in% c("Ernie", "Rock", "Sevenmile", "Big Lake"),] %>%
  filter(!is.na(Longitude_WGS84)) %>%
  data.frame




# making sure that the components of (inherited) current data script
# are robust
morphometry0 <- morph_new %>%
  # filter(!is.na(LakeName)) %>%
  # filter(!(is.na(Latitude_WGS84) & is.na(`Elevation (m)`) & is.na(`Temp (C)`) & is.na(SurfaceArea_h))) %>%
  mutate(use_fish = `Include in "Alaskanizing" Modeling Exercise` %in% c("Yes","yes")) %>%
  mutate(make_estimates = `Potentially Include in Lake Trout Management Plan` != "No") %>%
  arrange(LakeName)
morphometry1 <- morph_new %>%
  filter(!is.na(LakeName)) %>%
  filter(!(is.na(Latitude_WGS84) & is.na(`Elevation (m)`) & is.na(`Temp (C)`) & is.na(SurfaceArea_h))) %>%
  mutate(use_fish = `Include in "Alaskanizing" Modeling Exercise` %in% c("Yes","yes")) %>%
  mutate(make_estimates = `Potentially Include in Lake Trout Management Plan` != "No") %>%
  arrange(LakeName)

with(morph_new, table(!is.na(LakeName)))
with(morph_new, table(!(is.na(Latitude_WGS84) & is.na(`Elevation (m)`) & is.na(`Temp (C)`) & is.na(SurfaceArea_h))))

with(morphometry0, table(`Include in "Alaskanizing" Modeling Exercise`, use_fish, useNA="ifany"))
with(morphometry0, table(`Potentially Include in Lake Trout Management Plan`, make_estimates, useNA = 'ifany'))
with(morphometry1, table(`Include in "Alaskanizing" Modeling Exercise`, use_fish, useNA="ifany"))
with(morphometry1, table(`Potentially Include in Lake Trout Management Plan`, make_estimates, useNA = 'ifany'))




## ------------------ Length - Weight - Age ------------------ ##
# make sure columns are the same
dim(lw_old)
dim(lw_new)

colnames(lw_old)
colnames(lw_new)

# any new columns?
colnames(lw_new)[!(colnames(lw_new) %in% colnames(lw_old))]

# any old columns went away?
colnames(lw_old)[!(colnames(lw_old) %in% colnames(lw_new))]

# same order? doesn't matter too much but makes checks easier
all(colnames(lw_old) == colnames(lw_new))

# if so, make sure the columns are the same class
all(sapply(lw_old, class) == sapply(lw_new, class))

# # this is a bit more robust if the column order has changed
all(sapply(lw_old, class)[order(names(sapply(lw_old, class)))] ==
  sapply(lw_new, class)[order(names(sapply(lw_new, class)))])

# any new lakes?
unique(lw_new$LakeName[!(lw_new$LakeName %in% lw_old$LakeName)])
# might have to deal with "Other"
# [1] "Other"                          "Chikuminuk Lake"                "Nayakuk Lake"
# [4] "Tikchik Lake"                   "Kagati Lake"                    "Island Lake (Philip Smith Mts)"

# any old lakes went away?
unique(lw_old$LakeName[!(lw_old$LakeName %in% lw_new$LakeName)])

# any lake names that don't match with Morphometry?
unique(lw_old$LakeName[!(lw_old$LakeName %in% morph_old$LakeName)]) # baseline
unique(lw_new$LakeName[!(lw_new$LakeName %in% morph_new$LakeName)]) # what it is now
## NEW: Chelle, Other

# ISSUE: 7 lakes don't match Morphometry
# SOLUTION: same as last round, keep replacement step in data import
# Issue: one new one (Chelle)
# Solution: make it "Chelle Lake (Chaleau Lake or Hot Dog Lake)" to match

# ISSUE: lake "Other"  sum(lw_new$LakeName=="Other") [1] 3
# SOLUTION: take these samples out on import



# new projects?
proj_old <- with(lw_old, paste(LakeName, ProjectTitle, Year))
proj_new <- with(lw_new, paste(LakeName, ProjectTitle, Year))
unique(proj_new[!(proj_new %in% proj_old)])

# projects that went away?
unique(proj_old[!(proj_old %in% proj_new)])

# n samples per lake?
# length, weight, age, length & weight, length & age

alllakes <- sort(unique(c(lw_old$LakeName, lw_new$LakeName)))
# alllakes <- sort(unique(c(morph_old$LakeName, morph_new$LakeName)))

lw_old$lake2 <- factor(lw_old$LakeName, levels=alllakes)
lw_new$lake2 <- factor(lw_new$LakeName, levels=alllakes)
nsamples <- data.frame(nL_old = with(lw_old, tapply(!is.na(ForkLength_mm), lake2, sum)),
           nL_new = with(lw_new, tapply(!is.na(ForkLength_mm), lake2, sum)),
           nW_old = with(lw_old, tapply(!is.na(Weight_g), lake2, sum)),
           nW_new = with(lw_new, tapply(!is.na(Weight_g), lake2, sum)),
           nA_old = with(lw_old, tapply(!is.na(Age), lake2, sum)),
           nA_new = with(lw_new, tapply(!is.na(Age), lake2, sum)),
           nLW_old = with(lw_old, tapply(!is.na(ForkLength_mm) & !is.na(Weight_g), lake2, sum)),
           nLW_new = with(lw_new, tapply(!is.na(ForkLength_mm) & !is.na(Weight_g), lake2, sum)),
           nLA_old = with(lw_old, tapply(!is.na(ForkLength_mm) & !is.na(Age), lake2, sum)),
           nLA_new = with(lw_new, tapply(!is.na(ForkLength_mm) & !is.na(Age), lake2, sum)))
nsamples



nsamples %>%
  select(nL_old, nL_new) %>%
  filter(nL_old != nL_new | xor(is.na(nL_old), is.na(nL_new)))
nsamples %>%
  select(nW_old, nW_new) %>%
  filter(nW_old != nW_new | xor(is.na(nW_old), is.na(nW_new)))
nsamples %>%
  select(nA_old, nA_new) %>%
  filter(nA_old != nA_new | xor(is.na(nA_old), is.na(nA_new)))
nsamples %>%
  select(nLW_old, nLW_new) %>%
  filter(nLW_old != nLW_new | xor(is.na(nLW_old), is.na(nLW_new)))
nsamples %>%
  select(nLA_old, nLA_new) %>%
  filter(nLA_old != nLA_new | xor(is.na(nLA_old), is.na(nLA_new)))



## check: does the length-weight dataframe have the number of samples we expect
## from the morphometry dataframe?
nsamples %>%
  mutate(LakeName = rownames(nsamples)) %>%
  mutate(LakeName = ifelse(LakeName == "Donnelly Lake", "Donnelly Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Four Mile Lake", "Fourmile Lake (Taylor Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Lost Lake", "Lost Lake  (Chisholm Lake) (near Birch Lake)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "North Twin Lake", "North Twin Lake (Meadows Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Paul's Pond", "Pauls Pond (Coal Mine Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Rapids Lake", "Rapids Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Summit Lake (Paxson)", "Summit Lake (Richardson Highway near Paxson)", LakeName)) %>%
  left_join(select(morph_new, c("LakeName", "Count of ForkLength_mm", "Count of Weight_g"))) %>%
  select(LakeName, nL_new, `Count of ForkLength_mm`) %>%
  filter(nL_new != `Count of ForkLength_mm` | xor(is.na(nL_new), is.na(`Count of ForkLength_mm`)))

nsamples %>%
  mutate(LakeName = rownames(nsamples)) %>%
  mutate(LakeName = ifelse(LakeName == "Donnelly Lake", "Donnelly Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Four Mile Lake", "Fourmile Lake (Taylor Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Lost Lake", "Lost Lake  (Chisholm Lake) (near Birch Lake)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "North Twin Lake", "North Twin Lake (Meadows Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Paul's Pond", "Pauls Pond (Coal Mine Road)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Rapids Lake", "Rapids Lake (Richardson Highway)", LakeName)) %>%
  mutate(LakeName = ifelse(LakeName == "Summit Lake (Paxson)", "Summit Lake (Richardson Highway near Paxson)", LakeName)) %>%
  left_join(select(morph_new, c("LakeName", "Count of ForkLength_mm", "Count of Weight_g"))) %>%
  select(LakeName, nW_new, `Count of Weight_g`) %>%
  mutate(nW_new= ifelse(nW_new==0, NA, nW_new)) %>%
  filter(nW_new != `Count of Weight_g` | xor(is.na(nW_new), is.na(`Count of Weight_g`)))
