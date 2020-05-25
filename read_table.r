#Library initiation 
library("readxl")
library("xlsx")

#Reading the main excel files - All data, Reference Panel and Hi Set
my_data = read_excel("Sup_tables_FYI-BRCA1_GIM_APR_2020_Paulo.xlsx", sheet = "All data", col_names = TRUE)
#Appending important informations to the main table
#Adds two columns to the end of the table with numbers of "no functional impact" and "functional impact", respectively, for each missense variant 
my_data$"no_functional_impact" <- rowSums(my_data[,11:141] == 0, na.rm= TRUE)
my_data$"functional_impact" <- rowSums(my_data[,11:141] == 1, na.rm= TRUE)
#Adds one column to the end of the table after sum of the columns "no functional impact" and "functional impact"
my_data$"number of tests" <- my_data$"no_functional_impact"  + my_data$"functional_impact"
#HI SET
#Adds two columns to the end of the table with numbers of "no functional impact" and "functional impact", respectively, for each missense variant 
my_data$"no_functional_impact_Hi_set" <- rowSums(my_data[,c(11,22,28,29,54,55,61,79,81,82,92,93,103,111,116,121,128,132,134,135)] == 0, na.rm= TRUE)
my_data$"functional_impact_Hi_set" <- rowSums(my_data[,c(11,22,28,29,54,55,61,79,81,82,92,93,103,111,116,121,128,132,134,135)] == 1, na.rm= TRUE)
#Adds one column to the end of the table after sum of the columns "no functional impact" and "functional impact"
my_data$"number_of_tests_Hi_set" <- my_data$"no_functional_impact_Hi_set"  + my_data$"functional_impact_Hi_set"
#Load Reference Panel based sheet and creat a data.frame with number of variants tested "no functional impact" and "functional impact" for each track
reference_panel = read_excel("C:/Users/ufes/Desktop/BRCA1_Codes/Sup_tables_FYI-BRCA1_GIM_APR_2020_Paulo.xlsx", sheet = "Reference Panel", col_names = TRUE)
variants_no_functional_impact_reference_panel <- data.frame(colSums(reference_panel[,11:141] == 0, na.rm= TRUE))
variants_functional_impact_reference_panel <- data.frame(colSums(reference_panel[,11:141] == 1, na.rm= TRUE))
df_variants_reference_panel <- merge(variants_no_functional_impact_reference_panel, variants_functional_impact_reference_panel, by.x=0,by.y=0)
colnames(df_variants_reference_panel)<- c("track","Assay_no functional impact", "Assay_functional impact")
df_variants_reference_panel$"total" <- rowSums(df_variants_reference_panel[,2:3])

#START OF NEW TABS

#STable 5 - NUMBERS ABOUT CLASSES OF VARIANTS BEING TESTED

#Missense variants - Calculation for ALL FUNCTIONAL TRACKS and HI-SET of number of times variants were tested
variants_tested<-matrix(c(11009,11009,sum(my_data$"number of tests" == 0),sum(my_data$"number_of_tests_Hi_set" == 0),
                         sum(my_data$"number of tests" >= 1),sum(my_data$"number_of_tests_Hi_set" >= 1),
                         sum(my_data$"number of tests" >= 2),sum(my_data$"number_of_tests_Hi_set" >= 2),
                         sum(my_data$"number of tests" >= 10),sum(my_data$"number_of_tests_Hi_set" >= 10),
                         sum(my_data$"number of tests" >= 15),sum(my_data$"number_of_tests_Hi_set" >= 15),
                         sum(my_data$"number of tests" >= 20),sum(my_data$"number_of_tests_Hi_set" >= 20),
                         sum(my_data$"number of tests" >= 25),sum(my_data$"number_of_tests_Hi_set" >= 25),
                         sum(my_data$"number of tests" >= 30),sum(my_data$"number_of_tests_Hi_set" >= 30),
                         sum(my_data$"number of tests" >= 35),sum(my_data$"number_of_tests_Hi_set" >= 35),
                         sum(my_data$"number of tests" >= 40),sum(my_data$"number_of_tests_Hi_set" >= 40),
                         sum(my_data$"number of tests" >= 50),sum(my_data$"number_of_tests_Hi_set" >= 50),
                         sum(my_data$"number of tests" >= 60),sum(my_data$"number_of_tests_Hi_set" >= 60)),
                         ncol = 2, byrow=TRUE)
rownames(variants_tested)<-c("Possible unique (BRCA1)","Variants not yet tested","Variants tested more than 1 time","Variants tested more than two times","Variants tested more than 10 times",
                            "Variants tested more than 15 times","Variants tested more than 20 times","Variants tested more than 25 times","Variants tested more than 30 times","Variants tested more than 35 times","Variants tested more than 40 times",
                            "Variants tested more than 50 times","Variants tested more than 60 times")
colnames(variants_tested)<-c("ALL FUNCTIONAL TRACKS", "HI-SET")

#VUS only (excluding reference variants) - Calculation for ALL FUNCTIONAL TRACKS and HI-SET of number of times variants were tested
VUS_tested<-matrix(c(sum((my_data$"T7" == "NA") & (my_data$"number of tests" == 1)),0,
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 1)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 1)),  
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 2)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 2)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 10)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 10)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 15)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 15)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 20)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 20)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 25)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 25)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 30)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 30)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 50)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 50)),
                     sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 60)),sum((my_data$"T7" == "NA") & (my_data$"number_of_tests_Hi_set" >= 60))),
                     ncol = 2, byrow=TRUE)
rownames(VUS_tested) <- c("VUS Tested equal one time","VUS Tested more than one time","VUS tested more than two times",
                          "VUS tested more than 10 times", "VUS tested more 15 times", "VUS tested more than 20 times","VUS tested more than 25 times",
                          "VUS tested more than 30 times","VUS tested more than 50 times","VUS tested more than 60 times") 
colnames(VUS_tested)<-c("ALL FUNCTIONAL TRACKS", "HI-SET")

#Documented variants (BRCAExchange)
BRCAExchange <- matrix(c(sum(my_data$"T9" == 1),sum(my_data$"T9" == 1),
                         sum((my_data$"T9" == 1) & (my_data$"number of tests" >= 1)),sum((my_data$"T9" == 1) & (my_data$"number_of_tests_Hi_set" >= 1))),
                         ncol=2, byrow = TRUE)
rownames(BRCAExchange)<- c("Total documented missense variants","Total documented missense variants tested")
colnames(BRCAExchange)<-c("ALL FUNCTIONAL TRACKS", "HI-SET")                    

#STable 6 - FUNCTIONAL TRACK THROUGHPUT
#Number of variants tested by track
Number_of_variants_tested_track <- data.frame(colSums(my_data[,11:141] <= 1, na.rm= TRUE))
colnames(Number_of_variants_tested_track) <-("Total of variants tested")
Tracks <- matrix(c(131, sum(Number_of_variants_tested_track[,1] < 2), sum(Number_of_variants_tested_track[,1] >= 5),
                   sum(Number_of_variants_tested_track[,1] >= 10),sum(Number_of_variants_tested_track[,1] >= 20),
                   sum(Number_of_variants_tested_track[,1] >= 30),sum(Number_of_variants_tested_track[,1] >= 40),
                   sum(Number_of_variants_tested_track[,1] >= 50),sum(Number_of_variants_tested_track[,1] >= 60),
                   sum(Number_of_variants_tested_track[,1] >= 70),sum(Number_of_variants_tested_track[,1] >= 80),
                   sum(Number_of_variants_tested_track[,1] >= 90),sum(Number_of_variants_tested_track[,1] >= 100),
                   sum(Number_of_variants_tested_track[,1] >= 200),sum(Number_of_variants_tested_track[,1] >= 300),
                   sum(Number_of_variants_tested_track[,1] >= 400),sum(Number_of_variants_tested_track[,1] >= 500),
                   sum(Number_of_variants_tested_track[,1] >= 1000),sum(Number_of_variants_tested_track[,1] >= 1500)),
                   byrow = TRUE)
rownames(Tracks) <- c("Total number of functional tracks","Tracks testing 1 variant only","Tracks testing more than 5 variants",
                      "Tracks testing more than 10 variants", "Tracks testing more than 20 variants","Tracks testing more than 30 variants",
                      "Tracks testing more than 40 variants","Tracks testing >= 50 variants","Tracks testing >= 60 variants",
                      "Tracks testing >= 70 variants","Tracks testing >= 80 variants","Tracks testing >= 90 variants",
                      "Tracks testing >= 100 variants","Tracks testing >= 200 variants","Tracks testing >= 300 variants",
                      "Tracks testing >= 400 variants","Tracks testing >= 500 variants","Tracks testing >= 1000 variants",
                      "Tracks testing >= 1500 variants")
colnames(Tracks) <- "Functional track throughput"
                      
#Tracks with number of variants in reference panel non-pathogenic; pathogenic  X:Y
Ratios <- matrix(c(sum((df_variants_reference_panel[,2] >= 3) & (df_variants_reference_panel[,3] >= 3)),
                   sum((df_variants_reference_panel[,2] >= 4) & (df_variants_reference_panel[,3] >= 4)),
                   sum((df_variants_reference_panel[,2] >= 5) & (df_variants_reference_panel[,3] >= 5)),
                   sum((df_variants_reference_panel[,2] >= 10) & (df_variants_reference_panel[,3] >= 10))),
                   byrow = TRUE)
rownames(Ratios) <- c("Tracks with # variants in reference panel [non-pathogenic; pathogenic] equal or greater [3:3]",
                      "Tracks with # variants in reference panel [non-pathogenic; pathogenic] equal or greater [4:4]",
                      "Tracks with # variants in reference panel [non-pathogenic; pathogenic] equal or greater [5:5]",
                      "Tracks with # variants in reference panel [non-pathogenic; pathogenic] equal or greater [10:10]")
colnames(Ratios) <- c("Total")
#Tracks meeting the two criteria ( number of variants tested more than 10 and number of variants in reference panel equal or greater non-pathogenic; pathogenic x:y)
Criteria <- matrix(c(sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 3) & (df_variants_reference_panel[,3] >= 3)),
                     sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 4) & (df_variants_reference_panel[,3] >= 4)),
                     sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 5) & (df_variants_reference_panel[,3] >= 5)),
                     sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 10) & (df_variants_reference_panel[,3] >= 10))),
                     byrow=TRUE)
rownames(Criteria) <- c("Tracks meeting the two criteria ( #variants tested equal or greater 10 and # variants in reference panel equal or greater [non-pathogenic; pathogenic] [3:3])",
                        "Tracks meeting the two criteria ( #variants tested equal or greater 10 and # variants in reference panel equal or greater[non-pathogenic; pathogenic] [4:4])",
                        "Tracks meeting the two criteria ( #variants tested equal or greater 10 and # variants in reference panel equal or greater[non-pathogenic; pathogenic] [5:5])",
                        "Tracks meeting the two criteria ( #variants tested equal or greater 10 and # variants in reference panel equal or greater[non-pathogenic; pathogenic] [10:10])")
colnames(Criteria) <- c("Total")

#STable 8 - SPECIFICITY AND SENSIBILITY
#Calculating and creating data matrix
Espec_Sensi <- matrix(c(colSums((reference_panel[,11:141] == 0) & (reference_panel$"T7" > 3)),
                        colSums((reference_panel[,11:141] == 1) & (reference_panel$"T7" > 3)),
                        colSums((reference_panel[,11:141] == 0) & (reference_panel$"T7" < 3)),
                        colSums((reference_panel[,11:141] == 1) & (reference_panel$"T7" < 3))),
                        ncol = 4)
rownames(Espec_Sensi) <- c("T10",	"T11",	"T12",	"T13",	"T14",	"T15",	"T16",	"T17",	"T18",	
                           "T19",	"T20",	"T21",	"T22",	"T23",	"T24",	"T25",	"T26",	"T27",	
                           "T28",	"T29",	"T30",	"T31",	"T32",	"T33",	"T34",	"T35",	"T36",	
                           "T37",	"T38",	"T39",	"T40",	"T41",	"T42",	"T43",	"T44",	"T45",	
                           "T46",	"T47",	"T48",	"T49",	"T50",	"T51",	"T52",	"T53",	"T54",	
                           "T55",	"T56",	"T57",	"T58",	"T59",	"T60",	"T61",	"T62",	"T63",	
                           "T64",	"T65",	"T66",	"T67",	"T68",	"T69",	"T70",	"T71",	"T72",	
                           "T73",	"T74",	"T75",	"T76",	"T77",	"T78",	"T79",	"T80",	"T81",
                           "T82",	"T83",	"T84",	"T85",	"T86",	"T87",	"T88",	"T89",	"T90",	
                           "T91",	"T92",	"T93",	"T94",	"T95",	"T96",	"T97",	"T98",	"T99",	
                           "T100",	"T101",	"T102",	"T103",	"T104",	"T105",	"T106",	"T107",	"T108",	
                           "T109",	"T110",	"T111",	"T112",	"T113",	"T114",	"T115",	"T116",	"T117",	
                           "T118",	"T119",	"T120",	"T121",	"T122",	"T123",	"T124",	"T125",	"T126",	
                           "T127",	"T128",	"T129",	"T130",	"T131",	"T132",	"T133",	"T134",	"T135",	
                           "T136",	"T137",	"T138",	"T139",	"T140")
colnames(Espec_Sensi) <- c("classified_false_no_impact", "classified_true_impact","classified_true_no_impact", "classified_false_impact")

#Output excel file update
write.xlsx(variants_tested,file="output.xlsx", sheetName="STable5_variants_tested")
write.xlsx(VUS_tested,file="output.xlsx", sheetName="STable5_variants_hi_set", append=TRUE)
write.xlsx(BRCAExchange,file="output.xlsx", sheetName="STable5_Documented_Variants", append=TRUE)
write.xlsx(variants_tested,file="output.xlsx", sheetName="STable6_Variants_by_Track", append=TRUE)
write.xlsx(Ratios,file="output.xlsx", sheetName="STable6_#_variants_classified_by_track", append=TRUE)
write.xlsx(Criteria,file="output.xlsx", sheetName="STable6_#_tracks_meeting_criteria", append=TRUE)
write.xlsx(Espec_Sensi,file="output.xlsx", sheetName="STable8_sensibility_specificity", append=TRUE)

