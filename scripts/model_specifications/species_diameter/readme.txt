When making decisions about the species + size model there were choices in: 
1. Threshold of removing group: between 12 and 15 
2. Including groups under the threshold as misc small and misc large (problematic because have to remove as they cannot be given a FG) 

The csvs, df_add3.5_t12, and df_add3_5_t15, have no misc and thresholds 12 and 15 respectively. The auc and aic showed no difference between different thresholds but the threshold 12 model had a lot more uncertainty, therefore the threshold 15 was used. 

Misc grouped were included in the modelling but not subsequent as they cannot be given an FG, as they may be important for determining alpha for other groups. Therefore, sp_size_df_t15_updated is used. 

Also the updated csv uses R version 4.6.0 (most up to date at this time), other dfs use 4.5.1. 