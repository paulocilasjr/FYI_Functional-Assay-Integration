#Reading the main excel file - STable 1
library("readxl")
my_data = read_excel("C:/Users/ufes/Desktop/BRCA1_Codes/Sup_tables_FYI-BRCA1_GIM_APR_2020_Paulo_V2.xlsx", sheet = "STable 1", col_names = TRUE)

#Adds two columns to the end of the table with numbers of "no functional impact" and "functional impact", respectively, for each missense variant 
my_data$"no_functional_impact" <- rowSums(my_data[,11:141] == 0, na.rm= TRUE)
my_data$"functional_impact" <- rowSums(my_data[,11:141] == 1, na.rm= TRUE)
#Adds one column to the end of the table after sum of the columns "no functional impact" and "functional impact"
my_data$"number of tests" <- my_data$"no_functional_impact"  + my_data$"functional_impact"

#STable 5 

#ALL FUNCTIONAL TRACKS
#Missesense variants
# Total of missense variants - Possible unique (BRCA1)
total_of_missense_variants <- 11009
#Not yet tested
Not_yet_tested <- sum(my_data$"number of tests" == 0)
#tested >= 1 time
tested_more_than_one_time <- sum(my_data$"number of tests" >= 1)
#tested >= two times
tested_more_than_two_time <- sum(my_data$"number of tests" >= 2)
#tested >= 10 times
tested_more_than_ten_time <- sum(my_data$"number of tests" >= 10)
#tested >= 15 times
tested_more_than_15_times <- sum(my_data$"number of tests" >= 15)
#tested >= 20 times
tested_more_than_20_times <- sum(my_data$"number of tests" >= 20)
#tested >= 25 times
tested_more_than_25_times <- sum(my_data$"number of tests" >= 25)
#tested >= 30 times
tested_more_than_30_times <- sum(my_data$"number of tests" >= 30)
#tested >= 35 times
tested_more_than_35_times <- sum(my_data$"number of tests" >= 35)
#tested >= 40 times
tested_more_than_40_times <- sum(my_data$"number of tests" >= 40)
#tested >= 50 times
tested_more_than_50_times <- sum(my_data$"number of tests" >= 50)
#tested >= 60 times
tested_more_than_60_times <- sum(my_data$"number of tests" >= 60)

#HI SET
#Adds two columns to the end of the table with numbers of "no functional impact" and "functional impact", respectively, for each missense variant 
my_data$"no_functional_impact_Hi_set" <- rowSums(my_data[,c(11,22,28,29,54,55,61,79,81,82,92,93,103,111,116,121,128,132,134,135)] == 0, na.rm= TRUE)
my_data$"functional_impact_Hi_set" <- rowSums(my_data[,c(11,22,28,29,54,55,61,79,81,82,92,93,103,111,116,121,128,132,134,135)] == 1, na.rm= TRUE)
#Adds one column to the end of the table after sum of the columns "no functional impact" and "functional impact"
my_data$"number_of_tests_Hi_set" <- my_data$"no_functional_impact_Hi_set"  + my_data$"functional_impact_Hi_set"

#Not yet tested
Not_yet_tested_hi_set <- sum(my_data$"number of tests_Hi_set" == 0)
#tested >= 1 time
tested_more_than_one_time_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 1)
#tested >= two times
tested_more_than_two_time_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 2)
#tested >= 10 times
tested_more_than_ten_time_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 10)
#tested >= 15 times
tested_more_than_15_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 15)
#tested >= 20 times
tested_more_than_20_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 20)
#tested >= 25 times
tested_more_than_25_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 25)
#tested >= 30 times
tested_more_than_30_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 30)
#tested >= 35 times
tested_more_than_35_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 35)
#tested >= 40 times
tested_more_than_40_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 40)
#tested >= 50 times
tested_more_than_50_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 50)
#tested >= 60 times
tested_more_than_60_times_hi_set <- sum(my_data$"number_of_tests_Hi_set" >= 60)

#VUS only (excluding reference variants)

#Tested = one time
VUS_tested_one_time <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" == 1))
#Tested >= one time
VUS_tested_more_than_one_time <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 1))
#tested >= two times
VUS_tested_more_than_two_time <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 2))
#tested >= 10 times
VUS_tested_more_than_ten_time <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 10))
#tested >= 15 times
VUS_tested_more_than_15_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 15))
#tested >= 20 times
VUS_tested_more_than_20_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 20))
#tested >= 25 times
VUS_tested_more_than_25_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 25))
#tested >= 30 times
VUS_tested_more_than_30_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 30))
#tested >= 50 times
VUS_tested_more_than_50_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 50))
#tested >= 60 times
VUS_tested_more_than_60_times <- sum((my_data$"T7" == "NA") & (my_data$"number of tests" >= 60))

#Documented variants (BRCAExchange)
Total_documented_missense_variants <- sum(my_data$"T9" == 1)
Total_documented_missense_variants_tested <- sum((my_data$"T9" == 1) & (my_data$"number of tests" >= 1))

#STable 6 - FUNCTIONAL TRACK THROUGHPUT

#Total number of functional tracks
Total_functional_tracks <- 131

#First, the calculation of variants tested per track
Number_of_variants_tested_track <- data.frame(colSums(my_data[,11:141] <= 1, na.rm= TRUE))

#Tracks testing 1 variant only
tracks_testing_1_variant_only <- sum(Number_of_variants_tested_track[,1] == 1)
#Tracks testing >= 5 variants
tracks_testing_more_than_5_variants <- sum(Number_of_variants_tested_track[,1] >= 5)
#Tracks testing >= 10 variants
tracks_testing_more_than_10_variants <- sum(Number_of_variants_tested_track[,1] >= 10)
#Tracks testing >= 20 variants
tracks_testing_more_than_20_variants <- sum(Number_of_variants_tested_track[,1] >= 20)
#Tracks testing >= 30 variants
tracks_testing_more_than_30_variants <- sum(Number_of_variants_tested_track[,1] >= 30)
#Tracks testing >= 40 variants
tracks_testing_more_than_40_variants <- sum(Number_of_variants_tested_track[,1] >= 40)
#Tracks testing >= 50 variants
tracks_testing_more_than_50_variants <- sum(Number_of_variants_tested_track[,1] >= 50)
#Tracks testing >= 60 variants
tracks_testing_more_than_60_variants <- sum(Number_of_variants_tested_track[,1] >= 60)
#Tracks testing >= 70 variants
tracks_testing_more_than_70_variants <- sum(Number_of_variants_tested_track[,1] >= 70)
#Tracks testing >= 80 variants
tracks_testing_more_than_80_variants <- sum(Number_of_variants_tested_track[,1] >= 80)
#Tracks testing >= 90 variants
tracks_testing_more_than_90_variants <- sum(Number_of_variants_tested_track[,1] >= 90)
#Tracks testing >= 100 variants
tracks_testing_more_than_100_variants <- sum(Number_of_variants_tested_track[,1] >= 100)
#Tracks testing >= 200 variants
tracks_testing_more_than_200_variants <- sum(Number_of_variants_tested_track[,1] >= 200)
#Tracks testing >= 300 variants
tracks_testing_more_than_300_variants <- sum(Number_of_variants_tested_track[,1] >= 300)
#Tracks testing >= 400 variants
tracks_testing_more_than_400_variants <- sum(Number_of_variants_tested_track[,1] >= 400)
#Tracks testing >= 500 variants
tracks_testing_more_than_500_variants <- sum(Number_of_variants_tested_track[,1] >= 500)
#Tracks testing >= 1000 variants
tracks_testing_more_than_1000_variants <- sum(Number_of_variants_tested_track[,1] >= 1000)
#Tracks testing >= 1500 variants
tracks_testing_more_than_1500_variants <- sum(Number_of_variants_tested_track[,1] >= 1500)

#Load Reference Panel based sheet and creat a data.frame with number of variants tested "no functional impact" and "functional impact" for each track
reference_panel = read_excel("C:/Users/ufes/Desktop/BRCA1_Codes/Sup_tables_FYI-BRCA1_GIM_APR_2020_Paulo_V2.xlsx", sheet = "Reference Panel", col_names = TRUE)
variants_no_functional_impact_reference_panel <- data.frame(colSums(reference_panel[,11:141] == 0, na.rm= TRUE))
variants_functional_impact_reference_panel <- data.frame(colSums(reference_panel[,11:141] == 1, na.rm= TRUE))
df_variants_reference_panel <- merge(variants_no_functional_impact_reference_panel, variants_functional_impact_reference_panel, by.x=0,by.y=0)
df_variants_reference_panel$"total" <- rowSums(df_variants_reference_panel[,2:3])

#Tracks with number of variants in reference panel non-pathogenic; pathogenic  X:Y
Ratio_3_3 <- sum((df_variants_reference_panel[,2] >= 3) & (df_variants_reference_panel[,3] >= 3))
Ratio_4_4 <- sum((df_variants_reference_panel[,2] >= 4) & (df_variants_reference_panel[,3] >= 4))
Ratio_5_5 <- sum((df_variants_reference_panel[,2] >= 5) & (df_variants_reference_panel[,3] >= 5))
Ratio_10_10 <- sum((df_variants_reference_panel[,2] >= 10) & (df_variants_reference_panel[,3] >= 10))

#Tracks meeting the two criteria ( number of variants tested more than 10 and number of variants in reference panel equal or greater non-pathogenic; pathogenic x:y)
criterias_Ratio_3_3 <- sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 3) & (df_variants_reference_panel[,3] >= 3))
criterias_Ratio_4_4 <- sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 4) & (df_variants_reference_panel[,3] >= 4))
criterias_Ratio_5_5 <- sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 5) & (df_variants_reference_panel[,3] >= 5))
criterias_Ratio_10_10 <- sum((df_variants_reference_panel[,4] >= 10) & (df_variants_reference_panel[,2] >= 10) & (df_variants_reference_panel[,3] >= 10)) 