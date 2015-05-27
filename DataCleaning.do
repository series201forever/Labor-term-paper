/* Summary Stats for SIPP Data */

cd "/home/mallick/Desktop/NBER Data/SIPP93"
use "SIPP93.dta", clear

/* Dropping variables to stay within the 2,000 variable limit that Stata imposes
so that we can keep track of our other variables */
drop jobid101-jobid236
drop sc16721-sc16929
drop hh_inc01-hh_oth36
drop clswk101-ind1_01
drop lgtkey01-lgtoth36
drop wavflg01-wavflg36
drop wksper01-wksper36
drop ent_sp01-ent_sp36

destring pp_pnum, replace

////////////////////////////////////////////////////////////////////////////////
// 1. PRELIMINARY CLEANING //
////////////////////////////////////////////////////////////////////////////////
/* Dropping everyone who is not related to the reference person in the HH. So we
are basically left with nuclear and extended families */

/* Dropping "secondary individuals" which are defined as someone who is not
the HH reference person and is not related to anyone in the household */
drop if famtyp01==1 & famtyp02==1 & famtyp03==1 & famtyp04==1 & famtyp05==1 & ///
	famtyp06==1 & famtyp07==1 & famtyp08==1 & famtyp09==1 & famtyp10==1 & ///
	famtyp11==1 & famtyp12==1 & famtyp13==1 & famtyp14==1 & famtyp15==1 & ///
	famtyp16==1 & famtyp17==1 & famtyp18==1 & famtyp19==1 & famtyp20==1 & ///
	famtyp21==1 & famtyp22==1 & famtyp23==1 & famtyp24==1 & famtyp25==1 & ///
	famtyp16==1 & famtyp27==1 & famtyp28==1 & famtyp29==1 & famtyp30==1 & ///
	famtyp31==1 & famtyp32==1 & famtyp33==1 & famtyp34==1 & famtyp35==1 & ///
	famtyp36==1 
	
/* Dropping "unrelated subfamilies" which are nuclear families not related to
the HH reference person */
drop if famtyp01==2 & famtyp02==2 & famtyp03==2 & famtyp04==2 & famtyp05==2 & ///
	famtyp06==2 & famtyp07==2 & famtyp08==2 & famtyp09==2 & famtyp10==2 & ///
	famtyp11==2 & famtyp12==2 & famtyp13==2 & famtyp14==2 & famtyp15==2 & ///
	famtyp16==2 & famtyp17==2 & famtyp18==2 & famtyp19==2 & famtyp20==2 & ///
	famtyp21==2 & famtyp22==2 & famtyp23==2 & famtyp24==2 & famtyp25==2 & ///
	famtyp16==2 & famtyp27==2 & famtyp28==2 & famtyp29==2 & famtyp30==2 & ///
	famtyp31==2 & famtyp32==2 & famtyp33==2 & famtyp34==2 & famtyp35==2 & ///
	famtyp36==2 
	
/* Dropping "primary individuals" which are a reference person that lives alone
or only with unrelated people */
drop if famtyp01==4 & famtyp02==4 & famtyp03==4 & famtyp04==4 & famtyp05==4 & ///
	famtyp06==4 & famtyp07==4 & famtyp08==4 & famtyp09==4 & famtyp10==4 & ///
	famtyp11==4 & famtyp12==4 & famtyp13==4 & famtyp14==4 & famtyp15==4 & ///
	famtyp16==4 & famtyp17==4 & famtyp18==4 & famtyp19==4 & famtyp20==4 & ///
	famtyp21==4 & famtyp22==4 & famtyp23==4 & famtyp24==4 & famtyp25==4 & ///
	famtyp16==4 & famtyp27==4 & famtyp28==4 & famtyp29==4 & famtyp30==4 & ///
	famtyp31==4 & famtyp32==4 & famtyp33==4 & famtyp34==4 & famtyp35==4 & ///
	famtyp36==4 

////////////////////////////////////////////////////////////////////////////////
/* Cleaning marital statuses */
////////////////////////////////////////////////////////////////////////////////
/* Dropping those that were married, spouse absent at some point */
drop if ms_01==2 | ms_02==2 | ms_03==2 | ms_04==2 | ms_05==2 | ///
	ms_06==2 | ms_07==2 | ms_08==2 | ms_09==2 | ms_10==2 | ///
	ms_11==2 | ms_12==2 | ms_13==2 | ms_14==2 | ms_15==2 | ///
	ms_16==2 | ms_17==2 | ms_18==2 | ms_19==2 | ms_20==2 | ///
	ms_21==2 | ms_22==2 | ms_23==2 | ms_24==2 | ms_25==2 | ///
	ms_16==2 | ms_27==2 | ms_28==2 | ms_29==2 | ms_30==2 | ///
	ms_31==2 | ms_32==2 | ms_33==2 | ms_34==2 | ms_35==2 | ///
	ms_36==2 
	
/* Dropping those that were widowed at some point */
drop if ms_01==3 | ms_02==3 | ms_03==3 | ms_04==3 | ms_05==3 | ///
	ms_06==3 | ms_07==3 | ms_08==3 | ms_09==3 | ms_10==3 | ///
	ms_11==3 | ms_12==3 | ms_13==3 | ms_14==3 | ms_15==3 | ///
	ms_16==3 | ms_17==3 | ms_18==3 | ms_19==3 | ms_20==3 | ///
	ms_21==3 | ms_22==3 | ms_23==3 | ms_24==3 | ms_25==3 | ///
	ms_16==3 | ms_27==3 | ms_28==3 | ms_29==3 | ms_30==3 | ///
	ms_31==3 | ms_32==3 | ms_33==3 | ms_34==3 | ms_35==3 | ///
	ms_36==3 
	
/* Dropping those that were divorced at some point */
drop if ms_01==4 | ms_02==4 | ms_03==4 | ms_04==4 | ms_05==4 | ///
	ms_06==4 | ms_07==4 | ms_08==4 | ms_09==4 | ms_10==4 | ///
	ms_11==4 | ms_12==4 | ms_13==4 | ms_14==4 | ms_15==4 | ///
	ms_16==4 | ms_17==4 | ms_18==4 | ms_19==4 | ms_20==4 | ///
	ms_21==4 | ms_22==4 | ms_23==4 | ms_24==4 | ms_25==4 | ///
	ms_16==4 | ms_27==4 | ms_28==4 | ms_29==4 | ms_30==4 | ///
	ms_31==4 | ms_32==4 | ms_33==4 | ms_34==4 | ms_35==4 | ///
	ms_36==4 
	
/* Dropping those that were separated at some point */
drop if ms_01==5 | ms_02==5 | ms_03==5 | ms_04==5 | ms_05==5 | ///
	ms_06==5 | ms_07==5 | ms_08==5 | ms_09==5 | ms_10==5 | ///
	ms_11==5 | ms_12==5 | ms_13==5 | ms_14==5 | ms_15==5 | ///
	ms_16==5 | ms_17==5 | ms_18==5 | ms_19==5 | ms_20==5 | ///
	ms_21==5 | ms_22==5 | ms_23==5 | ms_24==5 | ms_25==5 | ///
	ms_16==5 | ms_27==5 | ms_28==5 | ms_29==5 | ms_30==5 | ///
	ms_31==5 | ms_32==5 | ms_33==5 | ms_34==5 | ms_35==5 | ///
	ms_36==5 
	
////////////////////////////////////////////////////////////////////////////////
/* Culling HH members */		
////////////////////////////////////////////////////////////////////////////////
/* Removing those that are "other relative of HH" for all periods (i.e. not 
reference, spouse, or kid */
drop if rrp_01==5 & rrp_02==5 & rrp_03==5 & rrp_04==5 & rrp_05==5 & ///
	rrp_06==5 & rrp_07==5 & rrp_08==5 & rrp_09==5 & rrp_10==5 & ///
	rrp_11==5 & rrp_12==5 & rrp_13==5 & rrp_14==5 & rrp_15==5 & ///
	rrp_16==5 & rrp_17==5 & rrp_18==5 & rrp_19==5 & rrp_20==5 & ///
	rrp_21==5 & rrp_22==5 & rrp_23==5 & rrp_24==5 & rrp_25==5 & ///
	rrp_16==5 & rrp_27==5 & rrp_28==5 & rrp_29==5 & rrp_30==5 & ///
	rrp_31==5 & rrp_32==5 & rrp_33==5 & rrp_34==5 & rrp_35==5 & ///
	rrp_36==5 

/* Removing those that are "non-relative" of the HH for all periods */
drop if rrp_01==6 & rrp_02==6 & rrp_03==6 & rrp_04==6 & rrp_05==6 & ///
	rrp_06==6 & rrp_07==6 & rrp_08==6 & rrp_09==6 & rrp_10==6 & ///
	rrp_11==6 & rrp_12==6 & rrp_13==6 & rrp_14==6 & rrp_15==6 & ///
	rrp_16==6 & rrp_17==6 & rrp_18==6 & rrp_19==6 & rrp_20==6 & ///
	rrp_21==6 & rrp_22==6 & rrp_23==6 & rrp_24==6 & rrp_25==6 & ///
	rrp_16==6 & rrp_27==6 & rrp_28==6 & rrp_29==6 & rrp_30==6 & ///
	rrp_31==6 & rrp_32==6 & rrp_33==6 & rrp_34==6 & rrp_35==6 & ///
	rrp_36==6 

drop if rrp_01==7 & rrp_02==7 & rrp_03==7 & rrp_04==7 & rrp_05==7 & ///
	rrp_06==7 & rrp_07==7 & rrp_08==7 & rrp_09==7 & rrp_10==7 & ///
	rrp_11==7 & rrp_12==7 & rrp_13==7 & rrp_14==7 & rrp_15==7 & ///
	rrp_16==7 & rrp_17==7 & rrp_18==7 & rrp_19==7 & rrp_20==7 & ///
	rrp_21==7 & rrp_22==7 & rrp_23==7 & rrp_24==7 & rrp_25==7 & ///
	rrp_16==7 & rrp_27==7 & rrp_28==7 & rrp_29==7 & rrp_30==7 & ///
	rrp_31==7 & rrp_32==7 & rrp_33==7 & rrp_34==7 & rrp_35==7 & ///
	rrp_36==7 

/* Removing HH types that were not at least once a "married couple HH" */
drop if lgthht01!=1 & lgthht02!=1 & lgthht03!=1 & lgthht04!=1 & lgthht05!=1 & ///
	lgthht06!=1 & lgthht07!=1 & lgthht08!=1 & lgthht09!=1 & lgthht10!=1 & ///
	lgthht11!=1 & lgthht12!=1 & lgthht13!=1 & lgthht14!=1 & lgthht15!=1 & ///
	lgthht16!=1 & lgthht17!=1 & lgthht18!=1 & lgthht19!=1 & lgthht20!=1 & ///
	lgthht21!=1 & lgthht22!=1 & lgthht23!=1 & lgthht24!=1 & lgthht25!=1 & ///
	lgthht16!=1 & lgthht27!=1 & lgthht28!=1 & lgthht29!=1 & lgthht30!=1 & ///
	lgthht31!=1 & lgthht32!=1 & lgthht33!=1 & lgthht34!=1 & lgthht35!=1 & ///
	lgthht36!=1 
 
////////////////////////////////////////////////////////////////////////////////
/* 2. HOUSEHOLD OBSERVATION CREATION */
////////////////////////////////////////////////////////////////////////////////
/* Generating husband and wife earning variables */
forvalues i=1/36 {
gen hus_ear`i'=0
label variable hus_ear`i' "Husband Earnings"
}

forvalues i=1/36 {
gen wife_ear`i'=0
label variable wife_ear`i' "Wife Earnings"

}

/* Generating husband and wife employment statuses */
forvalues i=1/36 {
gen hus_emp`i'=0
label variable hus_emp`i' "Husband Employment Status"

}

forvalues i=1/36 {
gen wife_emp`i'=0
label variable wife_emp`i' "Wife Employment Status"

}

/* Generating number of kids each period */
forvalues i=1/36 {
gen nkids`i'=0
label variable nkids`i' "Number of Kids"
}

/* Generating age of each kid each period */
forvalues i=1/13 {
	forvalues j=1/36{
	gen kid`i'_`j'=-100
	label variable kid`i'_`j' "Kid's Age"
	}
}

/* Generating race and ethnicity for husband and wife */
gen hus_race=-100
gen hus_eth=-100
gen wife_race=-100
gen wife_eth=-100

/* Now that the variables are generate, we need to fill them in with the proper
values */

/* First we record the spouse's wage and employment status for each period. 
I fix a person and search over all people and each period to find their spouse. 
When I find their spouse, I update the husband/wife's earnings, race, ethnicity,
and employment status. */ 
sort pp_id pp_pnum
forvalues i = 1/40336 {
disp `i'
local start = `i' - 17
local end = `i' + 17
	forvalues j = `start'/`end'{
		forvalues k = 1/36{
			if `k' < 10{
				if pp_id[`i']==pp_id[`j'] & pp_pnum[`j']==pnsp_0`k'[`i'] & sex[`j']==1{
					replace hus_ear`k'=pp_ear0`k'[`j'] in `i'
					replace hus_emp`k'=1 if esr_0`k'[`j'] < 4 in `i'
					replace hus_race=race[`j'] in `i'
					replace hus_eth=ethnicty[`j'] in `i'
				}
				else if pp_id[`i']==pp_id[`j'] & pp_pnum[`j']==pnsp_0`k'[`i'] & sex[`j']==2{
					replace wife_ear`k'=pp_ear0`k'[`j'] in `i'
					replace wife_emp`k'=1 if esr_0`k'[`j'] < 4 in `i'
					replace wife_race=race[`j'] in `i'
					replace wife_eth=ethnicty[`j'] in `i'
				}
				else if pp_id[`i']==pp_id[`j'] & pp_pnum[`i']==pnpt_0`k'[`j']{
					replace nkids`k'=nkids`k'[`i'] + 1 in `i'
					forvalues l=1/36{
						local apple=nkids`k'[`i']
						if `l' < 10 {
								replace kid`apple'_`l'= age_0`l'[`i'] in `i'
						}
						else if `l' >= 10{
								replace kid`apple'_`l'= age_`l'[`i'] in `i'
						}
					}					
				}
			}
			else {
				if pp_id[`i']==pp_id[`j'] & pp_pnum[`j']==pnsp_`k'[`i'] & sex[`j']==1{
					replace hus_ear`k'=pp_ear`k'[`j'] in `i'
					replace hus_emp`k'=1 if esr_`k'[`j'] < 4 in `i'
					replace hus_race=race[`j'] in `i'
					replace hus_eth=ethnicty[`j'] in `i'
				}
				else if pp_id[`i']==pp_id[`j'] & pp_pnum[`j']==pnsp_`k'[`i'] & sex[`j']==2{
					replace wife_ear`k'=pp_ear`k'[`j'] in `i'
					replace wife_emp`k'=1 if esr_`k'[`j'] < 4 in `i'
					replace wife_race=race[`j'] in `i'
					replace wife_eth=ethnicty[`j'] in `i'
				}
				else if pp_id[`i']==pp_id[`j'] & pp_pnum[`i']==pnpt_`k'[`j']{
					replace nkids`k'=nkids`k'[`i'] + 1 in `i'
					forvalues l=1/36{
						local apple=nkids`k'[`i']
						if `l' < 10 {
								replace kid`apple'_`l'= age_0`l'[`i'] in `i'
						}
						else if `l' >= 10{
								replace kid`apple'_`l'= age_`l'[`i'] in `i'
						}
					}
				}
			}
		}
	}
}

save "/home/mallick/Desktop/NBER Data/SIPP93/clean.dta", replace

/* Collapsing to a single HH observation */
use "/home/mallick/Desktop/NBER Data/SIPP93/clean.dta", clear
sort pp_id pp_pnum
bysort pp_id: keep if _n==1
keep if pp_pnum==101  // Keeping only HH reference person

forvalues i = 1/10717{
	forvalues k = 1/36{
		if `k' < 10{
			if sex[`i']==1{
				replace hus_ear`k'=pp_ear0`k'[`i'] in `i'
				replace hus_emp`k'=1 if esr_0`k'[`i'] < 4 in `i'
				replace hus_race=race[`i'] in `i'
				replace hus_eth=ethnicty[`i'] in `i'
			}
			else if sex[`i']==2{
				replace wife_ear`k'=pp_ear0`k'[`i'] in `i'
				replace wife_emp`k'=1 if esr_0`k'[`i'] < 4 in `i'
				replace wife_race=race[`i'] in `i'
				replace wife_eth=ethnicty[`i'] in `i'
			}
		}
		else {
			if sex[`i']==1{
				replace hus_ear`k'=pp_ear`k'[`i'] in `i'
				replace hus_emp`k'=1 if esr_`k'[`i'] < 4 in `i'
				replace hus_race=race[`i'] in `i'
				replace hus_eth=ethnicty[`i'] in `i'
			}
			else if sex[`i']==2{
				replace wife_ear`k'=pp_ear`k'[`i'] in `i'
				replace wife_emp`k'=1 if esr_`k'[`i'] < 4 in `i'
				replace wife_race=race[`i'] in `i'
				replace wife_eth=ethnicty[`i'] in `i'
			}
		}
	}
}	
save "/home/mallick/Desktop/NBER Data/SIPP93/clean.dta", replace

////////////////////////////////////////////////////////////////////////////////
/* 3. MORE CLEANING (SIM TO DEY AND FLINN 2008)*/
////////////////////////////////////////////////////////////////////////////////
/* Dropping couples in the armed forces */
drop if in_af_1==1 | in_af_2==1 | in_af_3==1 | in_af_4==1 | in_af_5==1 | ///
				in_af_6==1 | in_af_7==1 | in_af_8==1 | in_af_9==1
				
/* Dropping couples younger than retirement (i.e. less than 62) */
drop if age_36 >= 62
drop if age_01 >= 62
drop if age_36 == 0

/* Dropping those enrolled in school */
drop if att_sch1!=3 | att_sch2!=3 | att_sch3!=3 | att_sch4!=3 | att_sch5!=3 | ///
				att_sch6!=3 | att_sch7!=3 | att_sch8!=3 | att_sch9!=3

save "/home/mallick/Desktop/NBER Data/SIPP93/dfclean.dta", replace



/* 

by pp_id: egen numberofkids=total(child)
label variable numberofkids "Number of Kids in HH"
drop if numberofkids > 0  // From stats, this is less than 1.2% of sample

/* Dropping people that have negative income */
drop if pp_inc01 < 0

/* Employment rate of each parent by number of kids in age band */
bysort numberofkids relation: tabstat employed, stat(mean sd count)

/* Wage of mothers and fathers by number of kids in age band */
drop if employed==0
bysort numberofkids relation: tabstat pp_inc01, stat(mean sd count)

/* Correlations between each parent's income and number of kids in age band */
bysort relation: correlate pp_inc01 numberofkids

/* Correlations between employment rates and number of kids in age band */
bysort relation: correlate employed numberofkids
*/
