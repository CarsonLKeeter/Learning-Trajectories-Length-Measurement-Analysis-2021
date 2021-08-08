Title: Modeling relative sophistication of problem-solving strategies in early mathematics: a novel hurdle ordinal logit approach
Date: April 2021
Authors: Carson Keeter, MS        -   University of Wyoming  (Dept. of Statistics) 
 	 Pavel Chernyavskiy, PhD  -   University of Virginia (School of Medicine)
 	 Traci Kutaka, PhD        -   University of Virginia   (School of Education and Human Development)
 	 Douglas Clements, PhD    -   University of Denver   (Marsico Institute)
 	 Julie Sarama, PhD        -   University of Denver   (Marsico Institute)
   
This work was supported by Institute for Educational Sciences grant R305A200100 "Evolution of kindergartners' arithmetical problem-solving strategies" to Kutaka, Sarama, and Clements.

Contact: Carson Keeter @ keeterc1@gmail.com or Pavel Chernyavskiy @ pchern@virginia.edu with questions regarding code, analysis, etc. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------

Abstract: To date, interventions have typically focused on increasing correctness as the primary metric of efficacy across education research. 
However, correctness alone fails to comprehensively capture the set of competencies displayed by students. Here, we propose a shift in focus to 
the relative sophistication of problem-solving strategies ordered according to research-based developmental guidelines. We describe a novel 
modeling approach that treats ordinal strategies as the outcome of interest and explicitly accounts for the differential probability of 
detecting a strategy for each experimental condition, item, student, and classroom. Our model is estimated using the efficient No-U-Turn 
Hamiltonian Monte Carlo in Stan. We pilot our analysis on data collected during a kindergarten early geometry intervention situated within an 
urban school district in a Mountain West US state. We investigate the intervention efficacy of two one-on-one learning approaches relative to a 
comparison group and describe how to interpret item, student, and classroom random effects.

-----------------------------------------------------------------------------------------------------------------------------------------------------------
Code: Descriptions and proper order of code are seen below. 

(1) Data_plots.R - Visual summaries of data, mostly the frequency of post-assessment sophistication (y = 0, 1, 2, 3, H) partitioned by the various covariates (experimental condition, public/private school, etc). Constructed with 'ggplot2' and 'RColorBrewer'. Data* is stored as `fig.data` since some factor levels are changed for aesthetic reasons (not to be used in "Model_construction.R"). 

(2) Model_construction.R - Construction of the Ordered Cumulative Logit Hurdle model in Stan (done with 'rstan'), data* import/manipulation/format, prior probability and parameter specification, implementation of Stan model, and construction and use of custom summary function, `stan_sum()`. Additionally, one can find additional analysis regarding random effects and sophistication thresholds. Likelihood and traceplots are extracted to be used in "Model_diagnostics.R". 

(3) Model_diagnostics.R - Rather short, this code produces traceplots from `stan_sum()` and Pareto K/PSIS diagnostic plots from the likelihood extracted in "Model_construction.R". 

(4) Model_based_plots.R - Based on "Model_construction.R", this code constructs plots based on the samples and output taken from the Stan model. Sophistication thresholds are constructed for each experimental condition from the respective intercept and slope estimates. For these interval plots, the 'tidybayes' package is used to construct the visual analogue for these thresholds. Random intercepts/effects for students, classrooms, and items are also extracted** and plotted with 'ggplot2'. Lastly, the hypothetical data generation process is visually shown with 'ggplot2'. This plot is not data nor model based but included here for comparison to the sophistication threshold interval plot. 

(!) Code Order: Data_plots -->  Model_construction.R --> Model_diagnostics.R --> Model_based_plots.R

*There are instances where extra columns are imported that are blank or redundant. However, these columns may not always be present. Ensure that column elimination is necessary ("Data_plots.R": line  4; "Model_construction.R": line 203).  

**The 'tidyverse' package (loaded in "Model_construction.R" as well as 'tidybayes') masks the `extract()` function in 'rstan'. The error "Error in UseMethod("extract_"): no applicable method for 'extract_' applied to an object of class "stanfit"" can be avoided by directly calling `rstan::extract()` (seen in "Model_construction.R": lines 369, 391-393; "Model_based_plots": lines 120, 190, 233)

-----------------------------------------------------------------------------------------------------------------------------------------------------------
Data: "LT_length_data.csv"
- Note: data will be made public according to Institute of Education Sciences data-sharing guidelines following an embargo.

Variable Descriptions

- SID: Student identification

- Private: Private school indicator (1 = private school, 0 = public school)

- Class: Classroom identification

- Sex: Sex indicator (M/F) 

- Condition: Experimental condition (Learning Trajectories, Reverse-order, Business-as-Usual; LT, REV, BAU)

- Item: Questions developed to assess length measurement learning 

- Soph_pre*: Pre-assessment sophistication score 

- Correct_pre*: Pre-assessment correctness indicator for item (0 = incorrect, 1 = correct)

- Soph_post: Post-assessment sophistication score (y)

- Correct_post*: Post-assessment correctness indicator for item (0 = incorrect, 1 = correct)

- Theta: Pre-assessment sophistication z-score 

* Not used in any of the code above 

-----------------------------------------------------------------------------------------------------------------------------------------------------------
