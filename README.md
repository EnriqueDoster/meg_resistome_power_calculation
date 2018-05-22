# meg_resistome_power_calculation
Power calculation for shotgun metagenomic studies of AMR

## Current workflow:
1. Use AMR analytic data from NCBA2 dataset (N=60).  
2. Summarize count distribution for each AMR mechanism in the 4 groups of samples (Trt-Arrival,Ctrl-Arrival, Trt-Exit, Ctrl-Exit)
3. Use the count distributions in each sample group to simulate datasets with identical study design, increasing sample numbers, and different scenarios (ie. varying standard deviations, effect sizes)
   - Data can be simulated using; normal distribution, negative binomial distribution (2 different parameterizations),and zero-inflated binomial distribution. 
4. Use the Zero Inflated Gaussian Model to test differential abundance of AMR mechanisms and classes between treatment groups. 
5. Plot the proportion of times that an AMR feature abundance was statistically different between groups. 


### Changes that need to be made:
1. Only use Exit samples
2. Change ZIG model
3. Add theoretic effect for how treatment with a macrolide drug could change abundance for different mechanisms. 
4. Clean up code and make into a usable function or R package. 
