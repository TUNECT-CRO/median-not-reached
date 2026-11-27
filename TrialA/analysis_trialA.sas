 /* ----------------------------------------------------------
 Sponsor: TUNECT
 Project Name: Median Not Reached. Supplementary Materials. Trial A.
 
 Program Name: analysis_A
 Description : SAS implementation of the methods presented in the article for survival analysis without competing risks. 

 Programmer : Arsenii Zavalishin
 
 SAS Version : 9.4 (PC) ru 
 Operating System: Windows 10 x64 
----------------------------------------------------------- */

proc delete data = work._all_; run; 

/**********************************************************************************/
/*Please, upload Your ADTTE dataset*/
/*For outputs, please insert your actual Path in the output link*/
/**********************************************************************************/

/*---------------------------------------------------------------------* 
 * Descriptive Summary - Numbers at Risk                               * 
 *                                                                     * 
 * Purpose: Calculates the number of subjects at risk (AVAL >= time)   * 
 *          at specified landmark time points (in months) by treatment * 
 *          arm (TRT) and data cut (DCOID).                            * 
 * Note:    Time points for Interim (0, 6, 12, 15m) and Final (0, 12,  * 
 *          24, 47m) are adjusted based on available follow-up.        * 
 *---------------------------------------------------------------------*/ 

%macro time_points(data=, cutof=, months=, i=);
	/* Count subjects still at risk (time-to-event >= landmark month) */

	proc freq data = &data noprint;
		where DCOID = "&cutof" and aval >= &months;
		table TRT / out = freq_&i;
	run;
	/*Outputs:
		- OneWayFreqs, One-Way Frequencies
	*/

	/* Add data cut and time point identifiers, and rename count to n_risk */
	data freq_&i;
		set freq_&i;
		DCOID = "&cutof";
		time = &months;
		rename COUNT = n_risk;
		drop percent;
	run;

	/* Initialize output dataset on first run (i=11); append on subsequent runs */
	%if &i = 11 %then %do;
		data tA1_at_risk;
			retain DCOID TRT time n_risk;
			set freq_&i;
		run;
	%end;
	%else %do;
		data tA1_at_risk;
			retain DCOID TRT time n_risk;
			set tA1_at_risk freq_&i;
		run;
	%end;

	/* Clean up temporary dataset */
	proc delete data = freq_&i; run; 

%mend;


/* Execute macro for Interim analysis time points */
%time_points(data = adtte, cutof=Interim, months=0,  i=11);
%time_points(data = adtte, cutof=Interim, months=6,  i=12);
%time_points(data = adtte, cutof=Interim, months=12, i=13);
%time_points(data = adtte, cutof=Interim, months=15, i=14);

/* Execute macro for Final analysis time points */
%time_points(data = adtte, cutof=Final, months=0,  i=21);
%time_points(data = adtte, cutof=Final, months=12, i=22);
%time_points(data = adtte, cutof=Final, months=24, i=23);
%time_points(data = adtte, cutof=Final, months=36, i=24);
%time_points(data = adtte, cutof=Final, months=42, i=25);
%time_points(data = adtte, cutof=Final, months=47, i=26);

data tA1_at_risk;
	set tA1_at_risk;
	label
		DCOID = "Analysis Type"
		TRT = "Treatment"
		time = "Timepoint"
		n_risk = "Number at Risk"
	;
run;
/* Sort final at-risk summary table by data cut and treatment arm */
proc sort data = tA1_at_risk; 
	by descending DCOID TRT; 
run;


/*---------------------------------------------------------------------* 
 * Descriptive Summary of Events and Censoring                         * 
 *                                                                     * 
 * Purpose: Summarizes total sample size (N), number and percentage    * 
 *          of events and censored observations, and counts by         * 
 *          specific censoring/event descriptions (EVNTDESC) across    * 
 *          each data cut (DCOID) and treatment arm (TRT).             * 
 *---------------------------------------------------------------------*/ 

/* Sort input ADaM dataset by data cut and treatment group for BY-group processing */
proc sort data = adtte; 
	by DCOID TRT; 
run;

/* Calculate counts and percentages of overall events (CNSR=0) vs. censored (CNSR=1) */
proc freq data = adtte noprint;
	by DCOID TRT;
	table CNSR / out = freq_main_0;
run;

/* Calculate counts for specific event and censoring reasons (e.g., administrative vs. dropout) */
proc freq data = adtte noprint;
	by DCOID TRT;
	table EVNTDESC / out = freq_adm_1;
run;
/*Outputs for both:
	- OneWayFreqs, One-Way Frequencies
*/


/* Format counts as "n (xx.x%)" and assign descriptive labels for transposition */
data freq_main_1;
	length EV $100;
	set freq_main_0;
	res = strip(
			    strip(put(count, best.)) || " (" ||
				strip(put(round(percent, 0.1), best.)) || "%)"
		  );
	if cnsr = 0 then EV = "EVENT"; 
	else EV = "CENSORED"; 
	N = 140;
run;

/* Transpose overall event and censored summary strings into columns */
proc transpose data = freq_main_1 out = freq_main_2(drop = _NAME_); 
	by DCOID TRT N;
	id EV;
	var res;
run;

/* Transpose detailed censoring breakdown counts into columns, dropping total event count */
proc transpose data = freq_adm_1 out = freq_adm_2(drop = _NAME_ _label_ EVENT); 
	by DCOID TRT;
	id EVNTDESC;
	var COUNT;
run;

/* Merge overall summaries with detailed censoring counts into a single descriptive table */
data freq_desc;
	retain DCOID TRT N EVENT CENSORED;
	merge freq_main_2 freq_adm_2;
	by DCOID TRT;
run;

/* Sort final descriptive summary table by data cut and treatment arm */
proc sort data = freq_desc out = tA1_descriptive; 
	by descending DCOID TRT; 
run;

data tA1_descriptive;
	set tA1_descriptive;
	label
		DCOID = "Analysis Type"
		TRT = "Treatment"
		N = "Number of Subjects"
		Event = "Number of Events (%)"
		Censored = "Number of Censored (%)"
		Administrative_censoring = "Number of Administrative Censoring"
		Dropout = "Number of Dropout Censoring"
	;
run;


/*---------------------------------------------------------------------* 
 * Kaplan-Meier Median Time-to-Event & 95% CI                          * 
 *                                                                     * 
 * Purpose: Computes median time-to-event and 95% Brookmeyer-Crowley   * 
 *          confidence limits using a log-log transformation by data   * 
 *          cut (DCOID) and treatment arm (TRT). Formats missing       * 
 *          estimates as "NR" (Not Reached) and missing confidence     * 
 *          limits as "NE" (Not Estimable) in separate variables.      * 
 *---------------------------------------------------------------------*/ 

/* Suppress standard output and capture quartile estimates (including 50th percentile) */
proc sort data = adtte; by DCOID TRT; run;

ods select none;
ods output Quartiles = _raw_medians;
proc lifetest data=adtte conftype=loglog alpha=0.05;
    by DCOID TRT;
    id USUBJID;
    time AVAL*CNSR(1); 
run;
ods select all;
/*Outputs:
                  - OneWayFreqs,  One-Way Frequencies
                  - ProductLimitEstimates,  Product-Limit Estimates
                  - Quartiles,  Quartiles of the Survival Distribution
                  - Means,  Mean
                  - SurvivalPlot,  Survival Curve
                  - CensoredSummary,  Censored Summary
*/ 

/* Extract median (Percent = 50) and format point estimate and 95% CI bounds separately */
data final_summary(keep=DCOID TRT Median Median_LCL Median_UCL);
	retain DCOID TRT Median Median_LCL Median_UCL;
    set _raw_medians;
    where Percent = 50;
    
    length Median Median_LCL Median_UCL $10;
    
    /* Format median point estimate; assign "NR" if KM curve did not reach 0.50 */
    if Estimate = . then do;
        Median = "NR"; 
    end;
    else do;
        Median = put(Estimate, 8.1);
    end;
    
    /* Format lower 95% confidence limit; assign "NE" if not estimable */
    if LowerLimit = . then do;
        Median_LCL = "NE";
    end;
    else do;
        Median_LCL = put(LowerLimit, 8.1);
    end;
    
    /* Format upper 95% confidence limit; assign "NE" if not estimable */
    if UpperLimit = . then do;
        Median_UCL = "NE";
    end;
    else do;
        Median_UCL = put(UpperLimit, 8.1);
    end;

	/* Remove leading spaces from formatted character values */
	median = strip(median);
	median_LCL = strip(median_lcl);
	median_UCL = strip(median_ucl);
run;

/* Sort final median summary table by data cut and treatment arm */
proc sort data=final_summary out = tA2_median;
    by descending DCOID TRT; 
run;

data tA2_median;
	set tA2_median;
	label
		DCOID = "Analysis Type"
		TRT = "Treatment"
		Median = "Median"
		Median_LCL = "95% Lower Limit"
		Median_UCL = "95% Upper Limit"
	;
run;


/*---------------------------------------------------------------------* 
 * Reverse Kaplan-Meier Estimate of Follow-Up                          * 
 *                                                                     * 
 * Purpose: Calculates median, 95% CI, Q1, and Q3 follow-up times by   * 
 *          data cut (DCOID) and treatment arm (TRT) using the reverse * 
 *          Kaplan-Meier method.                                       * 
 * Logic:   Inverts the primary censoring indicator (CNSR) so that     * 
 *          primary censored observations become events (REV_CNSR=0)   * 
 *          and primary clinical events become censored (REV_CNSR=1).  * 
 *---------------------------------------------------------------------*/ 

/* Invert the censoring indicator for reverse Kaplan-Meier analysis */
data adtte_rev_prep;
    set adtte; 
    
    if CNSR = 0 then REV_CNSR = 1;      /* Primary event becomes censored */
    else if CNSR = 1 then REV_CNSR = 0; /* Primary censoring becomes follow-up event */
run;

/* Sort dataset by data cut and treatment group */
proc sort data=adtte_rev_prep;
    by DCOID TRT;
run;

/* Run reverse KM analysis and capture quartile estimates (25th, 50th, 75th) */

proc lifetest data=adtte_rev_prep plots=none;
    by DCOID;              
    strata TRT;           
    time AVAL * REV_CNSR(1);
    
    ods output Quartiles=rev_km_quartiles;
run;
/*Outputs:
                  - ProductLimitEstimates,  Product-Limit Estimates
                  - Quartiles,  Quartiles of the Survival Distribution
                  - Means,  Mean
                  - CensoredSummary,  Censored Summary
                  - HomStats,  Rank Statistics
                  - LogrankHomCov,  Log-Rank Covariance
                  - WilcoxonHomCov,  Wilcoxon Covariance
                  - HomTests,  Homogeneity Tests
*/

/* Extract median follow-up (50th percentile) and its 95% confidence limits */
data followup_50;
    set rev_km_quartiles;
    where Percent = 50;

    FU_median = estimate;
    FU_median_disp = round(estimate, 0.01);

    FU_med_LCL = LowerLimit;
    FU_med_LCL_disp = round(LowerLimit, 0.01);

    FU_med_UCL = UpperLimit;
    FU_med_UCL_disp = round(UpperLimit, 0.01);
    
    keep DCOID TRT FU_median FU_med_LCL FU_med_UCL FU_median_disp FU_med_LCL_disp FU_med_UCL_disp;
run;

/* Extract Q1 follow-up (25th percentile) */
data followup_25;
    set rev_km_quartiles;
    where Percent = 25;

    FU_q1 = estimate;
    FU_q1_disp = round(estimate, 0.01);
    
    keep DCOID TRT FU_q1 FU_q1_disp;
run;

/* Extract Q3 follow-up (75th percentile) */
data followup_75;
    set rev_km_quartiles;
    where Percent = 75;

    FU_q3 = estimate;
    FU_q3_disp = round(estimate, 0.01);
    
    keep DCOID TRT FU_q3 FU_q3_disp;
run;

/* Combine Q1, Median (with 95% CI), and Q3 follow-up estimates into one table */
data followup_all;
    retain DCOID TRT FU_median FU_med_LCL FU_med_UCL FU_q1 FU_q3 
           FU_median_disp FU_med_LCL_disp FU_med_UCL_disp FU_q1_disp FU_q3_disp;
    merge followup_25 followup_50 followup_75;
    by DCOID TRT;
run; 

/* Sort final follow-up summary table by data cut and treatment arm */
proc sort data=followup_all out = tA3_followup;
    by descending DCOID TRT;
run;

data tA3_followup;
	set tA3_followup;
	keep 
		DCOID
		TRT 
		FU_median_disp 
		FU_med_LCL_disp 
		FU_med_UCL_disp
		FU_q1_disp 
		FU_q3_disp 
	;
	label
		DCOID = "Analysis Type"
		TRT = "Treatment"
		FU_median_disp = "Follow-Up Median"
		FU_med_LCL_disp = "Follow-Up Median 95% CI Lower Limit"
		FU_med_UCL_disp = "Follow-Up Median 95% CI Upper Limit"
		FU_q1_disp = "Follow-Up Q1"
		FU_q3_disp = "Follow-Up Q3"
	;
run;



/*---------------------------------------------------------------------* 
 * Landmark Event-Free Survival (EFS) at Fixed Time Points             * 
 *                                                                     * 
 * Purpose: Estimates Kaplan-Meier survival probabilities S(t) and 95% * 
 *          log-log confidence intervals at 12, 24, and 36 months for  * 
 *          the Final data cut by treatment arm (TRT). Also computes   * 
 *          the absolute survival difference (Treatment - Control) and * 
 *          its Wald 95% confidence interval using Greenwood variance. * 
 *---------------------------------------------------------------------*/ 

/* Suppress standard output and extract landmark survival estimates at 12, 24, and 36 months */
ods select none;
ods output ProductLimitEstimates = _raw_efs;
proc lifetest data=adtte alpha=0.05 timelist=(12 24 36) outsurv=survdata reduceout conftype=loglog;
	where DCOID = "Final";
    strata TRT;
    time AVAL*CNSR(1);
run;
/*Outputs:
                  - ProductLimitEstimates,  Product-Limit Estimates
                  - Quartiles,  Quartiles of the Survival Distribution
                  - Means,  Mean
                  - CensoredSummary,  Censored Summary
                  - HomStats,  Rank Statistics
                  - LogrankHomCov,  Log-Rank Covariance
                  - WilcoxonHomCov,  Wilcoxon Covariance
                  - HomTests,  Homogeneity Tests
                  - SurvivalPlot,  Survival Curves
*/

ods select all;

/* Format landmark survival rates and 95% CIs per treatment arm */
data tA4_fixedtime;
    set survdata;
	t = TIMELIST;
	S = round(SURVIVAL, 0.001);
	LCL = round(SDF_LCL, 0.001);
	UCL = round(SDF_UCL, 0.001);
    keep TRT t S LCL UCL;
	label
		TRT = "Treatment"
		t = "Timepoint"
		s = "Event-Free Probability"
		LCL = "95% CI Lower Limit"
		UCL = "95% CI Upper Limit"
	;
run;

/* Sort raw product-limit estimates by time point for transposition */
proc sort data=_raw_efs;
    by Timelist;
run;

/* Transpose survival estimates S(t) to wide format (Treatment vs. Control) */
proc transpose data=_raw_efs out=_trans_diff_surv(drop = _name_ rename = (Treatment = Survival_Treatment Control = Survival_Control));
    by Timelist;
    id TRT;
  	var Survival;
run;

/* Transpose Greenwood standard errors to wide format (Treatment vs. Control) */
proc transpose data=_raw_efs out=_trans_diff_err(drop = _name_ _label_ rename = (Treatment = StdErr_Treatment Control = StdErr_Control));
    by Timelist;
    id TRT;
  	var StdErr;
run;

/* Merge survival probabilities and standard errors by landmark time point */
data _trans_diff;
	merge _trans_diff_surv _trans_diff_err;
	by Timelist;
run;

/* Calculate survival difference (Treatment - Control) and 95% Wald confidence limits */
data tA4_fixedtime_diff;
	retain Parameter t diff lcl ucl;
    set _trans_diff;
    by Timelist;
    
    if Survival_Treatment ne . and Survival_Control ne . then do;
        diff = Survival_Treatment - Survival_Control;
        se_diff = sqrt((StdErr_Treatment**2) + (StdErr_Control**2));
    end;
    else do;
        Diff_Est = "NE";
        Diff_LCL = "NE";
        Diff_UCL = "NE";
    end;

	/* Compute 95% Wald CI using standard normal quantile (z = 1.96) */
	z = probit(1 - 0.05/2); 
    lcl = round(diff - z * se_diff, 0.001);
    ucl = round(diff + z * se_diff, 0.001);
    diff = round(diff, 0.001);
	se_diff = round(se_diff, 0.001);

	t = timelist;

	Parameter = "Difference of Probabilities (Treatment - Control)";
	keep Parameter t diff se_diff lcl ucl;

	label
		Parameter = "Parameter"
		t = "Timepoint"
		diff = "Difference"
		se_diff = "Standard Error"
		LCL = "95% CI Lower Limit"
		UCL = "95% CI Upper Limit"
		se = "Standart Error"
	;
run;




/*---------------------------------------------------------------------* 
 * Restricted Mean Survival Time (RMST) Analysis                       * 
 *                                                                     * 
 * Purpose: Computes RMST by treatment arm (TRT), the RMST difference  * 
 *          (Treatment - Control), and 95% confidence limits across    * 
 *          each data cut (DCOID) up to a common restriction horizon   * 
 *          (tau).                                                     * 
 *---------------------------------------------------------------------*/ 

/* Capture group-level RMST estimates and pairwise RMST differences */
ods output RMST = _raw_rmst RMSTDiff = _raw_rmst_diff;
proc lifetest data=adtte rmst alpha=0.05;
	by DCOID;
    time AVAL*CNSR(1);    
	strata TRT / diff=control('Control')  ADJUST=SIMULATE;
run;
/*Outputs:
                  - ProductLimitEstimates,  Product-Limit Estimates
                  - Quartiles,  Quartiles of the Survival Distribution
                  - Means,  Mean
                  - CensoredSummary,  Censored Summary
                  - HomStats,  Rank Statistics
                  - LogrankHomCov,  Log-Rank Covariance
                  - WilcoxonHomCov,  Wilcoxon Covariance
                  - HomTests,  Homogeneity Tests
                  - SurvivalPlot,  Survival Curves
*/


/* Transpose group-level RMST estimates into wide format (Control vs. Treatment) */
proc transpose data = _raw_rmst out = rmst_1;
	by DCOID tau;
	id TRT;
	var estimate;
run; 

/* Rename and round group RMST estimates to 2 decimal places */
data rmst_2;
	set rmst_1;
	rename Control = RMST_Control Treatment = RMST_Treatment;
	RMST_Control = round(Control, 0.001);
	RMST_Treatment = round(Treatment, 0.001);
	keep DCOID Tau Control Treatment;
run;

/* Merge RMST differences with group estimates and format 95% Wald confidence limits */
data tA5_rmst;
	retain DCOID tau RMST_Control RMST_Treatment diff lcl ucl;
	merge _raw_rmst_diff rmst_2;
	by DCOID;
	RMST_Control = round(RMST_Control, 0.001);
	RMST_Treatment = round(RMST_Treatment, 0.001);
	diff = round(estimate, 0.001);
	z = probit(1 - 0.05/2); 
    lcl = round(estimate - z * StdErr, 0.001);
    ucl = round(estimate + z * StdErr, 0.001);
	tau = round(tau, 0.001);
	format tau 8.3;
	
	keep DCOID tau RMST_Control RMST_Treatment diff lcl ucl;
	label
		DCOID = "Analysis Type"
		tau = "Tau"
		RMST_Control = "RMST, Control"
		RMST_Treatment = "RMST, Treatment"
		diff = "RMST Difference, Treatment - Control"
		lcl = "RMST Difference, 95% CI Lower Limit"
		ucl = "RMST Difference, 95% CI Upper Limit"
	;
run;

proc sort data = tA5_rmst; by descending DCOID; run;


/*---------------------------------------------------------------------* 
 * Parametric Mixture Cure Model (Logistic / Weibull AFT)              * 
 *                                                                     * 
 * Purpose: Fits a two-component parametric mixture cure model to the  * 
 *          Final data cut:                                            * 
 *          1) Incidence: Logistic regression for long-term event-free * 
 *             fraction (cure probability, P_cure).                    * 
 *          2) Latency: Weibull AFT distribution for event time among  * 
 *             susceptible (uncured) subjects.                         * 
 *          Outputs parameter estimates, delta-method 95% CIs, and a   * 
 *          grid of fitted marginal survival curves for KM overlay.    * 
 *---------------------------------------------------------------------*/ 

/* Prepare Final data cut with explicit numeric treatment coding (1=Trt, 0=Ctrl) */
data adtte_final;
    set adtte;
    where lowcase(DCOID) = "final";
    
    if upcase(TRT) = "TREATMENT" then TRT_NUM = 1;
    else if upcase(TRT) = "CONTROL" then TRT_NUM = 0;
run;

/* Capture optimization outputs, parameter estimates, and delta-method quantities */
ods output ParameterEstimates  = nlmixed_raw_output
           AdditionalEstimates = model_estimates;
proc nlmixed data=adtte_final tech=quanew gconv=1e-8 df=1e6;
    /* Define starting values aligned with R optimization setup */
    parms b0=0 b1=0 intercept=1 beta1=0 scale=1;
    
    /* Constrain Weibull scale parameter to be strictly positive */
    bounds scale > 0;
    
    /* 1. INCIDENCE COMPONENT: Logistic model for long-term event-free fraction */
    linear_c = b0 + b1 * TRT_NUM;
    P_cure = 1 / (1 + exp(-linear_c)); 
    
    /* 2. LATENCY COMPONENT: Weibull AFT distribution among susceptible subjects */
    mu = intercept + beta1 * TRT_NUM;
    z = (log(AVAL) - mu) / scale;
    
    f_weibull = (1 / (scale * AVAL)) * exp(z - exp(z));
    S_weibull = exp(-exp(z));
    
    /* 3. MIXTURE LOG-LIKELIHOOD FORMULATION */
    if CNSR = 1 then 
        lik = P_cure + (1 - P_cure) * S_weibull;
    else if CNSR = 0 then 
        lik = (1 - P_cure) * f_weibull;
        
    /* Guard against log(0) numerical underflow during line search */
    if lik > 1e-12 then loglik = log(lik);
    else loglik = -1e8;
    
    model AVAL ~ general(loglik);
    
    /* 4. ESTIMATE TARGET METRICS AND 95% CIs VIA DELTA METHOD */
	estimate 'Long-term event-free fraction, Control'                            -b0                      df=1e6;
	estimate 'Long-term event-free fraction, Treatment'                           b0 + b1                 df=1e6;
	estimate 'Odds ratio for long-term event-free status (Treatment vs Control)' -b1                      df=1e6;
	estimate 'Latency log time ratio (Treatment vs Control)'                      beta1                   df=1e6;
	estimate 'Weibull shape'                                                     -log(scale)              df=1e6;

run;
/*Outputs:
	- Specifications
	- Dimensions
	- Parameters
	- IterHistory
	- ConvergenceStatus
	- FitStatistics
	- ParameterEstimates
	- AdditionalEstimates
*/

ods output close;


/* Format summary table of estimated quantities and 95% confidence limits */
data tA6_curemodel;
  length label $200.;
  set model_estimates;

  est = 0; lcl = 0; ucl = 0;

  if label = "Long-term event-free fraction, Control" then do;   
     est = logistic(-Estimate);
     lcl = logistic(-Upper);  
	 ucl = logistic(-Lower);
  end;

  if label = "Long-term event-free fraction, Treatment"  then do;  
     est = logistic(Estimate);
     lcl = logistic(Lower);  
	 ucl = logistic(Upper);
  end;

   if label = "Odds ratio for long-term event-free status (Treatment vs Control)" then do;                  
     est = exp(-Estimate);
     lcl = exp(-Upper);          
     ucl = exp(-Lower);
  end;

  if label = "Latency log time ratio (Treatment vs Control)" then do;                      
     est = Estimate;
     lcl = Lower;          
	 ucl = Upper;
  end;

  if label = "Weibull shape" then do;
	 est = exp(Estimate);
     lcl = exp(Lower);          
	 ucl = exp(Upper);
  end;

  est = round(est, 0.001);
  lcl = round(lcl, 0.001);
  ucl = round(ucl, 0.001);

  keep
  		Label 
		est
		lcl 
		ucl 
	;

  label 
  		Label = "Parameter"
		est = "Estimate"
		lcl = "95% CI, Lower Limit"
		ucl = "95% CI, Upper Limit"
	;
run;



/**********************************************************************************************************
                                                 OUTPUTS
                                     REPLACE PATH WITH YOUR ACTUAL PATH
**********************************************************************************************************/

proc export data=tA1_at_risk
  outfile='PATH/tA1_at_risk.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA1_descriptive
  outfile='PATH/tA1_descriptive.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA2_median
  outfile='PATH/tA2_median.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA3_followup
  outfile='PATH/tA3_followup.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA4_fixedtime
  outfile='PATH/tA4_fixedtime.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA4_fixedtime_diff
  outfile='PATH/tA4_fixedtime_diff.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA5_rmst
  outfile='PATH/tA5_rmst.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tA6_curemodel
  outfile='PATH/tA6_curemodel.csv'
  label 
  dbms=csv
  replace;
run;
