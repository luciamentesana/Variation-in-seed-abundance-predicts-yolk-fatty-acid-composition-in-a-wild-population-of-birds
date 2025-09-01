Variation-in-seed-abundance-predicts-yolk-fatty-acid-composition-in-a-wild-population-of-birds

In this study, we investigated inter-annual variation in yolk fatty acid composition of free-living great tits (Parus major) in relation to fluctuations in beech (Fagus sylvatica) fructification—their preferred food—across two years differing strongly in seed abundance. Here we place the raw data as well as the R codes used.



Description of the data and file structure

We carried out this study in a nest-box population of great tits (Parus major) in the Dellinger Buchet, in Southern Germany.  We collected data on ambient temperature (°C) every hour from12 i-buttons (DS9093A+ Thermochron iButton) placed at different locations within the forest.

  Beech tree fructification

We obtained information on the annual percentage of beech (Fagus sp.) mast fructification from the Bavarian State Minister for Food, Agriculture and Forestry (Bayerisches Staatsministerium für Ernährung, Landwirtschaft und Forsten 2018), and extracted the values using the R package ‘metaDigitise’ (Pick et al. 2019).

  Seed collection

We obtained tree seeds from the Bavarian Office for Forest Seeding and Planting, including from Norway spruce (Picea abies) and beech (Fagus sylvatica), which together account for 56.39% of the total tree coverage in our study population, as well as Scots pine (Pinus sylvestris).

  Egg collection

During the 2015 and 2016 breeding seasons, we checked each nest box every other day from the start of the breeding season. We collected the fourth egg from first and second clutches.



Files and variables

1) File: Mentesana_etal_Fatty_acid_yolk_content_2025.csv

Description: Database containing all information on great tit egg yolks collected.

  Variables

ID: Unique row number. 

Nest: Nest ID. 

1st egg laid: date (date-month-year) when the first egg of the clutch was laid by the female. 

Year: year when the egg was collected. This database includes only two years: 2015 and 2016.

Date collection: date (date-month-year) when the fourth egg was collected. 

Date collection num: date when the fourth egg was collected transformed into a numeric value so it can be included into the statistical models as a covariate. 

Egg number collected: the position of the egg within the laying sequence. We mostly collected the fourth egg; when this was not possible (e.g., due to harsh environmental conditions), we collected the fifth egg.

Time collection: time of the day when the egg was collected. 

Temperature - time of collection: temperature (in Celsius) when the egg was collected. 

Accumulated rain until collection: total rainfall (mm) from 00:00 h until the egg was collected.

Clutch number: whether the egg was collected from the 1st or 2nd clutch.

Female: unique ID for the mother. 

Prop Lauric (SFA), Prop Myristic acid (SFA), Prop 15_0 (SFA), Prop Palmitic acid (SFA), Prop 17_0 (SFA), Prop Stearic (SFA): proportion ('Prop') of individual saturated fatty acids (SFA). 

Prop Total SFA: sum of the proportions of all saturated fatty acids.

Prop 14_1 (MUFA), Prop Oleic acid (MUFA), Prop 16:1n-9 (MUFA), Prop Palmitoleic acid (MUFA), Prop cis-vaccenic acid (MUFA), Prop Gondoic acid (MUFA), Prop 20:1n-7 (MUFA): proportion ('Prop') of individual monounsaturated fatty acids (MUFA).

Prop Total MUFA: sum of the proportions of all monounsaturated fatty acids.

Prop alfa-linolenic acid (n-3 PUFA), Prop EPA (n-3 PUFA), Prop DPA (n-3 PUFA), Prop DHA (n-3 PUFA): proportion ('Prop') of individual omega-3 polyunsaturated fatty acids (n-3 PUFA).

Prop Total n-3 PUFA: sum of the proportions of all omega-3 polyunsaturated fatty acids.

Prop 16:02 (n-6 PUFA), Prop gamma-linolenic (n-6 PUFA), Prop Linoleic acid (n-6 PUFA), Prop pinoleic acid (n-6 PUFA), Prop Eicosadienoic acid (n-6 PUFA), Prop dihomo-gamma-linolenic acid (n-6 PUFA), Prop eicosatrienoic acid (n-6 PUFA), Prop Arachidonic acid (n-6 PUFA), Prop 22:4n-6 (n-6 PUFA), Prop 22:5n-6 (n-6 PUFA): proportion ('Prop') of individual omega-6 polyunsaturated fatty acids (n-6 PUFA).

Prop Total n-6 PUFA: sum of the proportions of all omega-6 polyunsaturated fatty acids.

Prop 18:2 (n-9 PUFA): proportion ('Prop') of individual omega-9 polyunsaturated fatty acids (n-9 PUFA).

Prop Total n-9 PUFA: sum of the proportions of all omega-9 polyunsaturated fatty acids.

Prop Ratio n6/n3: ratio between the proportion of all omega-6 and all omega-3 polyunsaturated fatty acids.

Temp 3 days before: mean temperature (in Celsius) of the three days prior to egg collection. 

Mean accumulated rain 3 days previous: mean rainfall (mm) of the three days prior to egg collection. 


  2) File: Mentesana_etal_Fatty_acids_yolk_content_2025.R

Description: R code to test for differences in egg fatty acid composition between 2015 and 2016. 


3) File: Mentesana_etal_Fatty_acid_Seed_content_2025.csv

Description: Database containing all information on seeds collected. Seeds belong to three tree species. 

Variables

Tree_spp: scientific name of the tree species for which the seeds were collected. 

Myristic_SFA, Palmitic_SFA, 17.0_SFA, Stearic_SFA, 20.0_SFA, 22.0_SFA: proportion of individual saturated fatty acids (SFA).

Total_SFA: sum of the proportions of all saturated fatty acids.

16:1n-9_MUFA, 16:1n-7 Palmitoleic_MUFA, Oleic_MUFA, cis-vaccenic_MUFA, Eicosenoic acid_MUFA, 22:1n-9_MUFA: proportion of individual monounsaturated fatty acids (MUFA). 

Total_MUFA: sum of the proportions of all monounsaturated fatty acids.

9, 12-Linoleic (18:2n-6)_n6_PUFA, 5, 9, 12-Pinolenic (18:3n-6)_n6_PUFA, 11,14-eicosadienoic acid (20:2n-6)_n6_PUFA, 7, 11, 14-eicosatrienoic acid (20:3n-6)_n6_PUFA: proportion of individual omega-6 polyunsaturated fatty acids (n-6 PUFA).  

Total_n6_PUFA: sum of the proportions of all omega-6 polyunsaturated fatty acids.

A-linolenic_n3_PUFA: proportion of individual omega-3 polyunsaturated fatty acids (n-3 PUFA). 

Total_n3_PUFA: sum of the proportions of all omega-3 polyunsaturated fatty acids.

7, 11-eicosadienoic acid (20:2n-9)_n6_PUFA, 5, 9-octadecadienoic acid (18:2n-9)_n6_PUFA: proportion of individual omega-9 polyunsaturated fatty acids (n-9 PUFA). 

Total_n9_PUFA: sum of the proportions of all omega-9 polyunsaturated fatty acids.

Total Percent: sum of the proportion of all fatty acids (i.e., SFA, MUFA, n-3, n-6 and n-9).


4) File: Mentesana_etal_Fatty_acids_Seed_content_2025.R

  Description: R code to test for differences in fatty acid composition across three seed species. 



5) File: Mentesana_etal_Beech_mass_quantification_2025.R

  Description: R code to obtain data on beech mass fructification across years in Bavaria, Germany.


6) File: Mentesana_etal_Beech_mass_fructification_2025.png

Description: Figure representing the annual percentage of beech (Fagus sp.) mast fructification from the Bavarian State Minister for Food, Agriculture and Forestry (Bayerisches Staatsministerium für Ernährung, Landwirtschaft und Forsten 2018). 



Code/software

In this submission, we included three R codes titled 'Mentesana_etal_Fatty_acids_yolk_content_2025', 'Mentesana_etal_Fatty_acids_Seed_content_2025', and 'Mentesana_etal_Beech_mass_quantification_2025'.

The software used was: R Core Team. 2013. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.


