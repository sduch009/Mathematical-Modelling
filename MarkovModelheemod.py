https://arxiv.org/pdf/1702.03252.pdf
Markov Model in R with heemod
Building and analysing a model: an example
Description of the question:
Imaginary disease called shame, assumptions:
asymptomatic at onset,
patients are at risk of being ashamed (marks entry into symptomatic phase)
symptomatic patients are at risk of dying

Compared Strategies:
Base strategy (base): Do nothing, this is the natural evolution of the disease
Medical treatment (med): Patients with asymptomatic disease are treated with ashaminib,
a highly potent shame inhibitor, until progression to symptomatic state in order to lower the risk of being ashamed.
Surgical treatment (surg): patients with asymptomatic disease undergo shamectomy, a surgical procedure that lowers the
risk of being ashamed. The procedure needs to be performed only once.

States (3):
Asymptomatic state (pre): Before the symptomatic state, when treatment can still be provided.
Symptomatic state (symp): Symptomatic disease, after being ashamed.  With degraded health, high hospital costs and
increased probability of dying.
Death (death): Death by natural causes of because of shame.

Model Parameters:
In Markov models values may depend on 2 distinct measurements of time:
time elapsed since the start of the model (called model time), and time spent in a given States
(called state time).  Both situations can co-exist in a same model.
In heemod, time-dependencyu is specified with 2 variables: model_time and state_time.  The pckg-reserved
names return sequential values starting from 1, corresponding to time spent in the model for model_time and
time spent in a given state for state_time.  They can be used in any user-defined expression or function.

In our case the probability of all-cause death depends on age.  Because age increases with time spent since the
beginning of the model, the all-cause death probability is model time dependent.  On the other hand the probability
of dying of shame depends on the time elapsed after being ashamed (i.e.  time spent in the symp state):  this probability
is state time dependent.The probability of being ashamed after surgery, the cost of surgery, and the hospital costs in
the symptomatic state also depend on the time spent in their respective state.In this model we will use a cycle duration
of 1 year:  we must take care that all transitions probabilities, values attached to states, and discount rates are calculated
on this time-frame.

We can now create the global parameters with define_parameters():

R> par_mod <- define_parameters(
R+  age_base = 20,
R+  age_cycle = model_time + age_base)

#The age of individuals for a given cycle age_cycle is the age at the beginning of the model (age_base), plus the time the model
#has run.

R> par_mod <- modify(
R+  par_mod,
R+
R+  sex_indiv = "MLE", #MLE => male in the WHO database
R+  p_death_all = get_who_mr(
R+    age      = age_cycle,
R+    sex      = sex_indiv,
R+    country  = "GBR",
R+    local    = TRUE))

#The death probability p_death_all, as a function of age and sex, is fetched from the WHO database with get_who_mr(), here for
a British population.

R> par_mod <- modify(
R+  par_mod,
R+
R+  p_death_disease = compute_surv(
R+    fit_death_disease,
R+    time     = state_time,
R+    km_limit = 5))

#The probability of dying of shame when the disease is symptomatic p_death_disease is extracted with the get_probs_from_surv() function
# from fit_death_disease,  a model fitted with the flexsurv pckg.  Because this probability depense on time spent with the disease state, the
#state_time model variable is used to specify time.  Here we use non-parametric Kaplan-Meier estimates for the first 5 years in stead of
#model-fitted values with km_limit = 5
#The parametric survival model fit_death_disease used to compute p_death_disease is fitted with the following code:

R> fit_death_disease <- flexsurv::flexsurvreg(
R+   survival::Surv(time, status) ~ 1,
R+   dist = "weibull",
R+   data = tab_surv)
#Where tab_surv is a data-fram containing survival data

R> dput(tab_surv)

structure(list(time = c(0.4, 8.7, 7, 5.1, 9.2, 1, 0.5, 3.3, 1.8,
3, 6.7, 3.7, 1.1, 5.9, 5.1, 10, 10, 10, 10, 10, 10, 10, 10, 10,
10), status = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)), .Names = c("time",
"status"), row.names = c(NA, -25L), class = "data.frame")

R> par_mod <- modify(
R+   par_mod,
R+
R+   p_death_symp = combine_probs(
R+     p_death_all,
R+     p_death_disease))

#The death probability in the symptomatic state p_death_symp is the probability to die either from old age (p_death_all)
#or from the disease (p_death_symp).  Assuming those probabilities are independent, we use the combine_probs() to combine them
#with the formula P(A∪B) = 1−(1−P(A))×(1−P(B)).

R> par_mod <- modify(
R+   par_mod,
R+
R+   p_disease_base = 0.25,
R+   med_effect     = 0.5,
R+   p_disease_med  = p_disease_base * med_effect)

#The probability of disease under medical treatment p_disease is the base probability of disease p_disease_base times the
#protecting effect of the treatment, med_effect.

R> par_mod <- modify(
R+   par_mod,
R+
R+   shape = 1.5,  #We will see later why we need
R+   scale = 5,    #to define these 2 parameters here.
R+   p_disease_surg = define_survival(
R+       distribution = "weibull",
R+       shape        = shape,
R+       scale        = scale) %>%
R+     compute_surv(time = state_time))

#The probability of disease after surgery is extracter with get_probs_from_surv() from a parametric Weibull survival model
#defined with define_survival().  For reason explained in Section 2.10 parameters scale and shape are not written in the
#define_survival() call, but defined separately.

R> par_mod <- modify(
R+   par_mod,
R+
R+   cost_surg     = 20000,
R+   cost_surg_cycle = ifelse(state_time ==1, cost_surg, 0))

#Because surgery is only performed once at the beginning of the pre state, the time-dependant variable state_time
#was  used to limit surgery costs to the first cycle in the pre state.

R> par_mod <- modify(
R+   par_mod,
R+
R+   cost_hospit_start = 11000,
R+   cost_hospit_end   = 9000,
R+   n_years           = 9,
R+   cost_hospit_cycle = ifelse(
R+     state_time < n_years,
R+     cost_hospit_start,
R+     cost_hospit_end))

#After n-years in the symptomatic state the symptoms become milder and hospital costs decrease (from cost_hospit_start to
#cost_hospit_end). We used the time-dependant variable state_time to condition the hospital costs cost_hospit on n_years

R> par_mod <- modify(
R+   par_mod,
R+






#rate for a qaly_disease, the QALY for one year in the symptomatic state.

2.5.  transitions
#We define a transition matrix for the base strategy with the define_transition() function.
#We can reference parameters defined in the previous section.

R> mat_base <- define_transition(
R+   state_names = c("pre", "symp", "death")
R+
R+   C,       p_disease_base,  p_death_all,
R+   p_cured, C,               p_death_symp,
R+   0,       0,               1)

R> mat_base

#A transition matrix, 3 state.
#
#       pre           symp            death
#pre    C             p_disease_base  p_death_all
#symp   p_cured       C               p_death_symp
#death                                1

#p_disease_base is the probability of being ashamed in the base strategy, p_death_all the all
#cause probability of death (not caused by disease) and p_death_symp the death probability
#in the symptomatic state (greater than p_death_all).  p_cured is the unlikely probability to
#revert to the asymptomatic unashamed state. p_cured, p_death_symp and p_death_all do not depend
#on the strategy.  The value of these parameters will be defined later.  C is an alias for the
#probability complement, 1 minus the sum of probabilities in a given row.  Death was modelled as an
#absorbing health state (ie the probability of transitioning from death to other health states was
#set to zero).  The resulting transition diagram is presented in figure 1.  Similarly, transitions can
#be defined for the 2 other strategies.  In our case only the name of the probabilites change:
#p_disease_base becomes p_disease med pr p_disease surg (those parameters will also be defined later)

R> mat_med <- defined_transition(
R+   state_names = c("pre", "symp", "death"),
R+
R+   C,        p_disease_med,     p_death_all,
R+   p_cured,  C,                 p_death_symp,
R+   0,        0,                 1)

2.6 State values
#Next we define the values associated with states using the define_state function.  An arbitraty
#number of values can be attached to a state, here we define: cost_treat the treatment cost (drug cost
# for the med strategy or surgery costs for the surg strategy, there is no treatment in the base strategy),
# cost_hospit the hospitalization costs, cost_total the total cost, and qaly the health_related quality adjusted life-years( QALY),
# where 1 stands for 1 year in perfect health and 0 stands for death.  In the following code we define the state pre:

R> state_pre <- define_state(
R+   cost_treat  = dispatch_strategy(
R+     base = 0, #no treatment => no treatment cost
R+     med  = cost_med,
R+     surg = cost_surg_cycle),
R+   cost_hospit = 0, #good health => no hospital expenses
R+   cost_total = discount(cost_treat + cost_hospit, r = dr),
R+   qaly       = 1)

#To dispatch the cost of treatment according to the strategy we used in the dispatch_strategy()
#function, with arguments named as strategies (we will define the strategy names in section 2.8).
#Another approach would have been to define 3 versions of state_pre, one per strategy, and in the
#next section use the corresponding version in each dinstinct define_strategy() call.

#The total cost is discounted with the discount() function at a givent rate (dr).  The QALY
#attached to 1 year in this state are set to 1, corresponding to 1 year in perfect health.  The variables dr,
# cost_med and cost_surg_cycle will be defined later.
#The 2 other states are defined similarly:

R> state_symp <- define_state(
R+   cost_treat  = 0,
R+   cost_hospit = cost_hospit_cycle,
R+   cost_total  = discount(cost_treat  + cost_hospit, r = dr),
R+   qaly        = qaly_disease)

R> state_death <- define_state(
R+   cost_treat  = 0,
R+   cost_hospit = 0,
R+   cost_total  = 0,
R+   qaly        = 0)

#patients have a degraded quality of life and are hospitalized during the symptomatic disease
#phase, we need to define specific QALYs (qaly_disease) and hospital costs (cost_hospit_cycle)
#for this state.  These variables will be defined later.  Finally dead patients have QALYs at 0,
#and they do not cost anything to the healthcare system.

2.7 strategies
#All the information (states and transitions) is now available to define the strategies.
#For this purpose we use the define_strategy() function.  Only the transition objects differ
#between strat_base, strat_med and strat_surg.

R> strat_base <- define_strategy(
R+   transition  = mat_base,
R+
R+   pre    = state_pre,
R+   symp   = state_symp,
R+   death  = state_death)

R> strat_med <- define_strategy(
R+   transition = mat_med,
R+
R+   pre    = state_pre,
R+   symp   = state_symp,
R+   death  = state_death)

R> strat_surg <- define_strategy(
R+   transition = mat_surg,
R+
R+   pre     = state_pre,
R+   symp    = state_symp,
R+   death   = state_death)

2.8 Running the model
#The model can then be run with run_model():
R> res_mod <- run_model(
R+   parameters = par_mod,
R+
R+   base = strat_base,
R+   med  = strat_med,
R+   surg = strat_surg,
R+
R+   cycles = 10,
R+
R+   cost   = cost_total,
R+   effect = qaly,
R+
R+   method = "life-table")

#base: detected use of'state_time', expanding states: pre, symp.
#Fetching mortality data from package cached data.Using cached data from year 2015.
#Fetching mortality data from package cached data.Using cached data from year 2015.
#med: detected use of'state_time', expanding states: pre, symp.
#surg: detected use of'state_time', expanding states: pre, symp.

#Strategy names are defined at that point by using the argument names provided by the user.
#We define a cost_total and qaly as the respective cost and effectiveness result.  The model is
#run for 10 cycles (ie 10 years), and state membership counts are corrected using the life-table method.
#By default the starting population is made of 1,000 patients in the first state, and no patient in the other states.

2.9 Results interpretation
#How do the strategies compare to each other with regard to their relative cost and effectiveness?
#Result is given by calculating the total expected cost and effectiveness of all strategies, and then computing the ICER.
#ICER defined as [(Cb - Ca) / (Eb - Ea)]
#Where C = total expected cost of a strategy and E its total expected effect (ie sum of life-years of the population).

#The strategies can be presented on a cost-effectiveness plane, where we see that both med and surg are more effective than base, but more costly.

R> summary(res_mod, threshold = c(1000, 5000, 15000))

#3 strategies run for 10 cycles.
#initial state counts:
#pre = 1000L
#symp = 0L
#death = 0L
#Values:
#     cost_treat cost_hospit cost_total     qaly
#base          0    54214446   42615142 5792.258
#med    27619456    37181168   52246211 7224.085
#surg   10074777    47429243   46220058 6553.701

#Net monetary benefit difference:
#      1000     5000     15000
#1 8199.241 2471.932     0.000
#2 5355.769 2674.230  7816.725
#3    0.000    0.000 11846.341

#Efficiency frontier:
#base -> surg -> med

#Differences:
     Cost Diff. Effect Diff.     ICER Ref.
surg   3604.915    0.7614427 4734.322 base
med    6026.153    0.6703846 8989.098 surg

#From the printed model output presented above we see in the ICER column of the Differences
#section that surg is more cost-effective than base if one is willing to pay 4,734 more per QALY gained.
#Furthermore med is more cost-effective than surg if one is willing to pay 8,989 more per QALY gained.
#A net monetary benefit analysis is run by specifying threshold ICER values in the summary() function with the threshhold argument.
#We see in the Net monetary benefit section that at a threshold ICER of 1,000 the strategy with the highest net monetary benefit is base, at 5,000 surg
#and med at 15,000.
##Figure 3 (generated by plot(res_mod, type = "counts", panel = "by_state") and plot(res_mod, type = "values", panel = "by_value")) gives us
#more information about what happens in our model: the effect of surgery seems to wear down with time comapred to the medical treatment.  surgery
# delays the outcome, reporting degraded health status and hospitacl costs further in time.  After a few years the hospital costs
#of the surgery strategy reach similar levels to the base strategy.  Nevertheless these increased hospital costs do not outweigh the important
#treatment costs associated with the medical therapy.

2.10 Uncertainty analysis
#What is the uncertainty of these results?  What strategy is probably the most cost-effective?
#Uncertainty of the results originate from uncertainty regarding the true value of the input parameters (eg treatment effect, hospital costs
#quality of life with teh disease, survival probabilities).  The effect of this uncertainty can be assessed by varying the parameter Values
#and computing the model results with these new inputs.  While multiple methods exist to study uncertainty, deterministic and probabilistic
#sensitivity analysis (DSA and PSA) are the most widely used.

#In a DSA, parameter values are changed one by one, usually to a low and high value (eg the lower and upper bounds of the parameter confidence interval)
#Model results are plotted on a tornado plot to display how a change in the value of one parameter impacts the model results.  A DSA gives
# a good sense of the relative impact of each parameter on the uncertainty of the model outcomes, but does not account for the total uncertainty over all the parameters,
# for skewed or complex parameter distribution, nor for correlations between the errors of different parameter estimeates.

#We define the DSA with define_dsa() by specifying a lowe and upper bound for each parameter of interest:

R> def_dsa <- define_dsa(
R+   age_base,           15,     30,
R+   p_disease_base,     0.2,    0.3,
R+   p_cured,            0.005,  0.02,
R+   med_effect,         0.3,    0.7,
R+   shape,              1.4,    1.6,
R+   scale,              4,      6,
R+   cost_med,           4000,   6000,
R+   cost_surg,          8000,   12000,
R+   cost_hospit_start,  5000,   15000,
R+   dr,                 0,      0.1,
R+  qaly_disease,        0.3,    0.7,
R+  n_years,             8,      10)

R> res_dsa <- run_dsa(res_mod, dsa = def_dsa)
#Running DSA on strategy 'base'...
#Running DSA on strategy 'med'...
#Running DSA on strategy 'surg'...

#Only parameters (eg state values, transition probabilities) defined with define_parameters() can be modified in a DSA (or a PSA)
#Accordingly, many state values, transition probabilities and the shape and scale parameters used in our example were defined as parameters,
#this allowing them to be varied in sensitivity analyses.  Once defined, the analysis can be run using run_dsa().

#Figure 4 (generated by plot(res_dsa, result = "cost", strategy = "med")) shows the impact of varying each parameter individually on total cost for the med strategy.  The results demonstrate
#that the discount rate and hospital costs have a greater impact than other parameters.  Unsurprisingly  parameters used only in the
#strategy and parameters unrelated to costs have no effect on the cost of the med strategy.

#PSA
#In PSA, the model is re-run for a given number of simulations with each parameter being replaced with a value re-sampled from a user-defined
#probability distribution.  These results are then aggregated, allowing us to obtain the probability distribution of model outputs.
#We define the parameter distribution with define_psa(), and optionally their correlation structure with define_correlation().

R> def_psa <- define_psa(
R+   age_base          ~ normal(mean = 20, sd = 5),
R+   p_disease_base    ~ binomial(prob = 0.25, size = 500),
R+   p_cured           ~ binomial(prob = 0.001, size = 500),
R+   med_efect         ~ lognormal(mean = 0.5,  sd = 0.1),
R+   shape             ~ normal(mean = 1.5,  sd = 0.2),
R+   scale             ~ normal(mean = 5,  sd = 1),
R+   cost_med          ~ gamma(mean = 5000,  sd = 1000),
R+   cost_surg         ~ gamma(mean = 200000, sd = 3000),
R+   cost_hospit_start ~ gamma(mean = 11000,  sd = 2000),
R+   dr                ~ binomial(prob = 0.05, sd = 100),
R+   qaly_disease      ~ normal(mean = 0.5,  sd = 0.1),
R+   n_years           ~ poisson(mean = 9),
R+
R+   correlation = define_correlation(
R+      shape,    scale,         -0.5,
R+      age_base, p_disease_base, 0.3,))

R> res_psa <- run_psa(res_mod, psa = def_psa, N = 1000)
#Resampling strategy 'base'...
#Resampling strategy 'med'...
#Resampling strategy 'surg'...

#We then run the PSA with
run_psa()
#, here for 1,000 re-samplings.  The results can be plotted as uncertainty clouds on the cost-effectiveness
#plane (Figure 5) (a similar plot can be generated with
plot(res_psa, type = "ce"))

#The probability of a strategy being cost effective can be plotted for various willingness to pay values on a cost-effectiveness acceptability curve.
#This plot can be generated by:
plot(res_psa, type = "ac")

#It is possible to compute the expected value of perfect information (EVPI) depending on willingness to pay.  This is a quantification of
# potentially choosing the wrong strategy, and thus conversely the price one is ready to pay to reduce the risk of incorrect decisions by
#obtaining more information (conducting more studies).
#In this example, the EVPI generated by:
plot(res_psa, type = "evpi")
#peaks between 1,000 and 10,000, where the uncertainty is high.  It also increases for higher willingness to pay, because even though the Uncertainty
# is not as high, the costs of a wrong decision become higher.
#The EVPI can indicate whether conducting more research is cost-effective.  But it does notinform on the value of getting more information on particular parameters
#(Briggset al.2006).The expected value of perfect information for parameters (EVPPI) is very similar to the EVPI,but returns values by parameters (Adeset  al.2004).
#Unfortunately its computation is nottrivial.  PSA results can be exported to compute EVPPI with the Sheffield Accelerated Valueof InformationSAVIsoftware (Stronget al.2014) with export_savi()

#The individual contribution of parameter uncertainty on the overall uncertainty is illustrated by
plot(res_psa, type = "cov") #We could also perform the same analysis on the differencebetween strategies with the optiondiff = TRUE.
#We can see taht, depending on the   strategy,  different parameters generate the uncertainty on costs and effect.  In all case sdr,cost_hospit_startandqaly_disease explain  a
#high  part  of  variability  for  all  strategies.   Unsurprisingly  the effect  ofscale(the  scale  of  the  post-surgery  Weibull  survival  function),med_effectandcost_medare
#limited to the surg or med strategies.

#n addition, average model values can be computed on the results and presented in a summarysimilar to therun_model()output.  Because of non-linearities in Markov models, averagesover the PSA
#output distribution are more accurate than point estimates (Briggset al.2006).In our case the ICERs changed from 4,734 to 6,453 and 8,989 to 7,059 for the surg and med strategies respectively.

2.11 Heterogeneity analysis --> How does the cost-effectiveness of strategies vary depending on the characteristics of the population?
#If population haracteristics are available, model results can be computed on the different sub-populations to study the heterogeneity of the resulting model outputs (Briggset  al.2006).
#Furthermore, average population-level results can be computed from these distributions.The model we ran in Section 2.9 computed results for a cohort of males aged 20.  To assess how population
#characteristics affect model results we can run a heterogeneity analysis.  We use the update() function to run the model on a table containing population data.


#in this example, we use a table with population characteristics, here named tab_pop, with an optional column .weights giving the relative population weight of each strata
R> head(tab_pop)
   age_base sex_indiv   .weights
 1       10       MLE 0.04242889
 2       10      FMLE 0.86696571
 3       15       MLE 0.69960873
 4       15      FMLE 0.51253057
 5       20       MLE 0.91723545
 6       20      FMLE 0.09685623

R> pop_mod <- update(res_mod, newdata = tab_pop)

Updating strategy'base'...
Updating strategy'med'...
Updating strategy'surg'...

#The summary of the updated model gives the distribution of the values of interest in the population, and the average model values over the entire population.
#Here the average ICERs in the population are 5,052 and 9,150 for the surg and med strategies respectively, quite similar to the values computed in section 2.9
#(4,734 and 8.989).  We can also plot the distribution of model results by generating:
plot(pop_mod, result = "effect", bins = 15)

2.12 Budget impact analysis -->What would be the total cost of a strategy for the health system?
#So far we mostly worked on model results at the scale of the infividual (eg cost per person).
#If we want to implement a strategy at a health system level we also need to know the total cost ofer a given time horizon, in order to assess whether the strategy is sustainable.
#This is called a budget impact analysis (BIA).  The main differences with the classic model are (1) the patient counts at the model start should reflect the population statistics,
#and (2) additional patients may enter the model every year (new disease cases).

#We use the
init and inflow arguments of run_model() to implement BIA, here for the med strategy.
#The inflow of new patients is defined with define_inflow().  Inflow counts can depend on model time (state time dependency is meaningless in this context).

R> res_bia <- run_model(
R+   parameters = par_mod,
R+
R+   med = strat_med,
R+
R+   cycles = 10,
R+
R+   cost = cost_total,
R+   effect = qaly,
R+
R+   method = "life-table",
R+
R+   init =   c(
R+     pre = 25000,
R+     symp = 5000,
R+     death = 0),
R+   inflow = define_inflow(
R+     pre = 8000,
R+     symp = 0,
R+     death = 0))

med: detected use of 'state_time', expanding state: pre, symp.

#At the start of the model there are 25,000 patients with asymptomatic shame and 5,000 with a symptomatic form of the disease in the population.
#Every year 8,000 additional cases of shame are added to the model, starting the disease in the asymptomatic state.

R> summary(res_bia)

1 strategy run for 10 cycles.

Initial state counts:

pre = 25000
symp = 5000
death = 0

Counting method: 'life_table'.

Values:
    cost_treat cost_hospit cost_total     qaly
med 1942366603  2621681008 3531335211 508495.2

The total cost of strategy med over a 10-year time horizon will be 3.5 billions.



3. Other Features and Extentions
#This section introduces features and extensions that were not presented in the previous example

3.1 Survival analysis:  The heemod package provides a nuimber of ways to estimate transition probabilities from survival distributions.
#Survival distributions can come from at least 3 different sources:
•User-definded parametric distributions created using thedefine_survival()function.
•Fitted parametric distributions withflexsurv::flexsurvreg()(Jackson 2016).
•Fitted Kaplan-Meiers withsurvival::survfit()(Therneau 2015; Therneau and Gramb-sch 2000).

#Once defined each of these types of distributions can be combined and modified using a standard set of operations.
#Treatments effects can be applied to any survival distribution:
•Hazard ratio:apply_hr().
•Odds ratio:apply_or().
•Acceleration factor:apply_af().

#The transition or survival probabilities are computed with compute_serv().  Time (usually model_time or state_time) needs to be passed to the function
#as a time argument.
#All these operations can be chained with the %>% piping operator, eg:
R> fit_cov %>%
R+   apply_hr(hr = 2) %>%
R+   join(
R+     fitcov_poor,
R+     at = 3) %>%
R+   pool(
R+     fitcov_medium,
R+     weights = c(0.25, 0.75)) %>%
R+   add_hazards(
R+     fit_w) %>%
R+   compute_surv(time = 1:5)

3.2 Convenience functions
#For reproducibility and ease of use we implemented convenience functions to perform some of the most common calculations needed in health economic evaluation studies
#(eg converting incidence rates, odds ratios, or relative risks to transition priobabilities with
rate_to_prob()
or_to_prob()
rr_to_prob()
#Probabilites and discount rates can be rescaled to fit different time frames (generally the duration of a cycle) with
rescale_prob() and rescale_discount_rate()

3.3 Cluster computing
#PSA and heterogeneity analyses can become time-consuming since they consist in iteratively re-running the model with new parameter inputs.  Because this workload is
#"embarassingly parallel", ie there is no dependency or need for communication between the parallel tasks, it can easily be run on a cluster relying on the parallel (R Core Team 2016) package.
#This is done by calling the
use_cluster()
#function.  This function can either take as an arguments
# 1. A number: a local cluster with the given number of cores will be created.
# 2. A cluster object defined with the
makeCluster()
#function from parallel: the user-defined will be used (eg to use more complex clusters with non-local hosts).

3.5 Extensions to other types of models
#Even though the main focus of heemod is to compute Markov models, other methods that model state changes can be included in the package:
#only the transition argument of
define_strategy()
#and the associated evaluation methods need to be extended.
#For example, partitioned survival models (williams et al 2016b) were added to the package recently.  These models can be computed byt passing an object defined by
define_part_surv()
# to
transition
#In theory most modelling methods that return state counts over time could be integrated into heemod, eg
dynamic models for infectious diseases (Snedecor 2012) https://www.ispor.org/news/articles/Nov-Dec12/understanding-use-dynamic-models.asp
Anderson RM, May RM. Infectious diseases of humans: dynamics and control. Reprinted ed. Oxford etc.: Oxford University Press; 2002.
Elbasha EH, Dasbach EJ. Impact of vaccinating boys and men against HPV in the United States. Vaccine 2010;28:6858-67.
Halloran ME, Cochi SL, Lieu TA, et al. Theoretical epidemiologic and morbidity effects of routine varicella immunization of preschool children in the United States. Am J Epidemiol 1994;140:81-104.
Hethcote HW. The mathematics of infectious diseases. SIAM Rev 2000;42:599-653.
Kermack WO, McKendrick AG. Contributions to the mathematical theory of epidemics--I. 1927. Bull Math Biol 1991;53:33-55.
Kim SY, Goldie SJ. Cost-effectiveness analyses of vaccination programmes : a focused review of modelling approaches. Pharmacoeconomics 2008;26:191-215.
Snedecor SJ, Strutton DR, Ciuryla V, Schwartz EJ, Botteman MF. Transmission-dynamic model to capture the indirect effects of infant vaccination with Prevnar (7-valent pneumococcal conjugate vaccine (PCV7)) in older populations. Vaccine 2009;27:4694-703.
Dynamic Transmission Modeling Working Group. ISPOR 11 A.D. July 11. Available from: http://www.ispor.org/workpaper/modeling_methods/Dynamic-transmission-modeling.asp

4. MATHEMATICAL IMPLEMENTATION
#In this section, we detail the mathematcial implementation of most of the features of heemod.

4.1 Parameter correlation in psa
#Correlation of parameters in PSA was implementated with the following steps:
#   1. A correlation structure is defined with
define_correlation()
#See section 2.10
#   2. Values are sampled from a multi-normal distribution having the required correlation structure with mvnfast (fasiolo 2016)
#   3. The sample values are the mapped to the target distributions on a quantile by quantile basis.

#This approach described in Briggset al.(2006) is an approximation that allows to define cor-relations between arbitrary distributions.
# The final Pearson correlation coefficients betweenthe target distributions may differ slightly from the ones initially defined by the user.
#Thatissue is mostly true if the target distributions are too dissimilar, e.g.  a gamma and a binomialdistribution

4.2 Time-dependency Implementation
#A thorought description of time-dependency in Markov models is given by Hawkinset  al.(2005), with solutions for the computation of both non-homogeneous and semi-Markov mod-els.
#These methods were implemented in theheemodpackage.  Markov models withmodeltimedependency are usually called non-homogeneous Markov models, and models withstatetimedependency are called semi-Markov models

4.3 Implementing budget impact analysis
#...
