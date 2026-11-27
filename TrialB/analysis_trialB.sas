 /* ----------------------------------------------------------
 Sponsor: TUNECT
 Project Name: Median Not Reached. Supplementary Materials. Trial B.
 
 Program Name: analysis_trialB
 Description : SAS implementation of the methods presented in the article for survival analysis with competing risks. 

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
 * Derive Event Indicator (EVNT) for Competing Risks Analysis          * 
 *                                                                     * 
 * Note: Under CDISC ADaM standards, the CNSR variable uses 1 to       * 
 * indicate censoring and 0 for events. This step creates a custom     * 
 * numeric variable (EVNT) coded as required for the subsequent        * 
 * cause-specific Cox regression models (0 = Censored, 1 = Progression,* 
 * 2 = Death without progression).                                     * 
 *---------------------------------------------------------------------*/ 
data adtte_event;
	set adtte;
	/* Default to censored (0), then assign event codes based on description */
	evnt = 0;
	if EVNTDESC = "Progression" then EVNT = 1;
	if EVNTDESC = "Death without progression" then EVNT = 2;
run;

/*---------------------------------------------------------------------*
	Frequency of Event Status by TRT
 *---------------------------------------------------------------------*/

ods output CrossTabFreqs = work.freqs
           ChiSq         = work.chisq;
proc freq data=adtte_event;
    tables EVNT*TRT / chisq nopercent norow nocum;
run;
/*Outputs:
	- CrossTabFreqs, Frequency Table of TRT by evnt
	- ChiSq, Statistics for Table of TRT by evnt 
*/

/* 2. Prepare the Frequency and Column Percentage dataset */
data clean_freqs;
    set freqs;
    /* Keep only the actual cross-tabulation cells (exclude row/col totals) */
    where _TYPE_ = '11'; 
    
    /* Combine Frequency and Column Percent into a single formatted string */
    length stat_val $25;
    stat_val = strip(strip(put(Frequency, comma6.)) || ' (' || strip(put(ColPercent, 5.1)) || '%)');
run;

/* 3. Transpose so Treatment groups become columns */
proc sort data=clean_freqs;
    by TRT  EVNT;
run;

proc transpose data=clean_freqs out=trans_freqs(drop=_NAME_);
    by TRT;
    id EVNT;
    var stat_val;
run;

data tB1_counts;
	retain TRT _1 _2 _0;
	set trans_freqs;
	label
		trt = "Treatment"
		_0 = "Administrative censoring"
		_1 = "Progression"
		_2 = "Death without progression"
		;
run;


/*---------------------------------------------------------------------*
	Probability of Progression by Milestones (Aalen-Johansen)
 *---------------------------------------------------------------------*/

/* Compute the Cumulative Incidence Function (CIF) using Aalen-Johansen */
ods output CIF = cif_estimates;
proc lifetest data=adtte_event method=km plots=none timelist=6 12 24 36;
    strata TRT;
    time AVAL * EVNT(0) / failcode=1; /* failcode=1 targets Progression */
run;
/*Outputs:
	- FailureSummary, Failures Summary
	- CIF, Cumulative Incidence Function Estimates
	- GrayTest, Gray's Test
*/

ods output  ProductLimitEstimates = km_estimates;
proc lifetest data=adtte_event method=km plots=none outsurv=km_out timelist=6 12 24 36;
    strata TRT;
    time AVAL * EVNT(0, 2); /* failcode=1 targets Progression */
run;


proc sort data=cif_estimates;
    by TRT Timelist;
run;

proc sort data=km_estimates;
    by TRT Timelist;
run;

data tB2_cif_vs_km;
	retain trt timelist cif_c failure_c;
	length cif_c Failure_c $64;
	merge cif_estimates(keep = TRT Timelist CIF CIF_LCL CIF_UCL)  km_estimates(keep = TRT Timelist Survival Failure StdErr);
	by TRT Timelist;
	cif_c = strip(
				strip(put(round(cif,     0.001), best.)) || " [" ||
				strip(put(round(cif_lcl, 0.001), best.)) || ", " ||
				strip(put(round(cif_ucl, 0.001), best.)) || "]"
				);
	/*1-KM 95% CI interval	*/
	z = probit(1 - 0.05/2);
	if Survival = 1 then do;
            Failure_LCL = 0;
            Failure_UCL = 0;
        end;
    else if Survival = 0 then do;
            Failure_LCL = 1;
            Failure_UCL = 1;
        end;
        else do;
            /* Log-log standard error and bounds for Survival */
            se_loglog = StdErr / (Survival * abs(log(Survival)));
            s_lcl = Survival ** exp( z * se_loglog);
            s_ucl = Survival ** exp(-z * se_loglog);
            
            /* Invert Survival bounds to get Failure bounds */
            Failure_LCL = 1 - s_ucl;
            Failure_UCL = 1 - s_lcl;
    end;

	Failure_c = strip(
				strip(put(round(Failure,     0.001), best.)) || " [" ||
				strip(put(round(Failure_lcl, 0.001), best.)) || ", " ||
				strip(put(round(Failure_ucl, 0.001), best.)) || "]"
				);

	keep trt timelist cif_c Failure_c;
	label 
		TRT = "Treatment"
		Timelist = "Timepoints"
		CIF_c = "CIF [95%]"
		Failure_c = "1 - KM [95%]"
		;
run;

proc sort data = cif_estimates; by timelist; run; 

proc transpose data=cif_estimates out=cif_estimates_trans(drop=_NAME_);
    by timelist;
    id TRT;
    var CIF;
run;

data tB2_cif_diff;
	retain Parameter timelist Treatment Control c_diff_TminusC;
	set cif_estimates_trans;
	c_diff_TminusC = Treatment - Control;
	Parameter = "Cumulative Incidence Function";
	keep Parameter timelist Treatment Control c_diff_TminusC;
	label 
		Parameter = "Parameter"
		timelist = "Timepoints"
		Treatment = "Treatment"
		Control = "Control"
		c_diff_TminusC = "Difference, Treatment - Control"
		;
run;

/*---------------------------------------------------------------------*
	Cause-Specific Cox Proportional Hazards Models
 *---------------------------------------------------------------------*/

ods output ParameterEstimates = cox_1;
proc phreg data=adtte_event;
    /* Define reference group for treatment arm */
    class TRT (ref="Control") / param=ref; 
    
    /* For progression, only EVNT=1 is the event. 
       EVNT=0 and EVNT=2 are treated as censored (indicated by parenthetical values) */
    model AVAL * EVNT(0, 2) = TRT / rl;
run;


ods output ParameterEstimates = cox_2;
proc phreg data=adtte_event;
    /* Define reference group for treatment arm */
    class TRT (ref="Control") / param=ref; 
    
    /* For progression, only EVNT=2 is the event. 
       EVNT=0 and EVNT=1 are treated as censored (indicated by parenthetical values) */
    model AVAL * EVNT(0, 1) = TRT / rl;
run;
/*Outputs for both:
	- ModelInfo, Model Information
	- NObs, Number of Observations
	- ClassLevelInfo, Class Level Information
	- CensoredSummary, Summary of Event and Censored Observations
	- ConvergenceStatus, Convergence Status
	- FitStatistics, Model Fit Statistics
	- GlobalTests, Test of Global Null Hypothesis
	- ModelANOVA, Type 3 Tests
	- ParameterEstimates, Maximum Likelihood Estimates of Model Parameters
*/


data tB3_cs_cox;
	retain Event estimate StdErr HazardRatio HRLowerCL HRUpperCL ProbChiSq;
	length event $200. ;
	set cox_1 (in=a) cox_2 (in=b);
	if a then Event = "Cause-specific Cox: progression (deaths censored)";
	if b then Event = "Cause-specific Cox: death without progression";

	estimate    = round(estimate,    0.001);
	StdErr      = round(StdErr,      0.001);
	HazardRatio = round(HazardRatio, 0.001);
	HRLowerCL   = round(HRLowerCL,   0.001);
	HRUpperCL   = round(HRUpperCL,   0.001);
	ProbChiSq   = round(ProbChiSq,   0.0001);

	format estimate StdErr 8.3;

	label event = "Model"
		  estimate = "Coefficient"
		  StdErr = "Standard Error"
		  HazardRatio = "Hazard Ratio"
		  HRLowerCL = "Hazard Ratio, 95% CI Lower Limit"
		  HRUpperCL = "Hazard Ratio, 95% CI Upper Limit"
		  ProbChiSq = "p-value"
		  ;
	keep Event estimate StdErr HazardRatio HRLowerCL HRUpperCL ProbChiSq;
run;


/*---------------------------------------------------------------------*
	Fine-Gray Subdistribution Hazard Model for Progression
 *---------------------------------------------------------------------*/

ods output ParameterEstimates = fine_gray;
proc phreg data=adtte_event;
    /* Define reference group for treatment arm */
    class TRT (ref="Control") / param=ref; 
    
    /* eventcode=1 specifies Progression as the event of interest.
       SAS automatically identifies any other non-censored value (EVNT=2) 
       as the competing risk. */
    model AVAL * EVNT(0) = TRT / eventcode=1 rl;
run;


/*Outputs:
	- ModelInfo, Model Information
	- NObs, Number of Observations
	- ClassLevelInfo, Class Level Information
	- FailureSummary, Summary of Failure Outcomes
	- ConvergenceStatus, Convergence Status
	- FitStatistics, Model Fit Statistics
	- GlobalTests, Test of Global Null Hypothesis
	- ModelANOVA, Type 3 Tests
	- ParameterEstimates, Maximum Likelihood Estimates of Model Parameters
*/

data tB4_finegray;
	retain model  estimate StdErr HazardRatio HRLowerCL HRUpperCL ProbChiSq;
	length model  $200. ;
	set fine_gray;

	estimate    = round(estimate,    0.001);
	StdErr      = round(StdErr,      0.001);
	HazardRatio = round(HazardRatio, 0.001);
	HRLowerCL   = round(HRLowerCL,   0.001);
	HRUpperCL   = round(HRUpperCL,   0.001);
	ProbChiSq   = round(ProbChiSq,   0.0001);
	format estimate StdErr 8.3;

	model = "Fine-Gray subdistribution model: progression";

	label model = "Model"
		  estimate = "Coefficient"
		  StdErr = "Standard Error"
		  HazardRatio = "Hazard Ratio"
		  HRLowerCL = "Hazard Ratio, 95% CI Lower Limit"
		  HRUpperCL = "Hazard Ratio, 95% CI Upper Limit"
		  ProbChiSq = "p-value"
		  ;
	keep model  estimate StdErr HazardRatio HRLowerCL HRUpperCL ProbChiSq;
run;



/**********************************************************************************************************
                                                 OUTPUTS
                                     REPLACE PATH WITH YOUR ACTUAL PATH
**********************************************************************************************************/

proc export data=tB1_counts
  outfile='PATH/tB1_counts.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tB2_cif_vs_km
  outfile='PATH/tB2_cif_vs_km.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tB2_cif_diff
  outfile='PATH/tB2_cif_diff.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tB3_cs_cox
  outfile='PATH/tB3_cs_cox.csv'
  label 
  dbms=csv
  replace;
run;

proc export data=tB4_finegray
  outfile='PATH/tB4_finegray.csv'
  label 
  dbms=csv
  replace;
run;

