// Moment Conditions //
///////////////////////
// This file gets the appropriate moment conditions for the labor term paper.

cd "/home/mallick/Desktop/NBER Data/SIPP93"
use dfclean.dta, replace

// Generating wage and age bins
gen taxbracket = 0
replace taxbracket = 1 if hus_ear1+wife_ear1 >=20000
replace taxbracket = 2 if hus_ear1+wife_ear1 >=11000 & hus_ear1+wife_ear1 < 20000
replace taxbracket = 3 if hus_ear1+wife_ear1 >=5000 & hus_ear1+wife_ear1 < 11000
replace taxbracket = 4 if hus_ear1+wife_ear1 >=3000 & hus_ear1+wife_ear1 < 5000
replace taxbracket = 5 if hus_ear1+wife_ear1 >=0 & hus_ear1+wife_ear1 < 3000

gen agebracket = 0
replace agebracket = 1 if age_36 >=18 & age_36 < 30
replace agebracket = 2 if age_36 >=30 & age_36 < 40
replace agebracket = 3 if age_36 >=40 & age_36 < 50
replace agebracket = 4 if age_36 >=50 & age_36 < 62

// (1) Childbirth rate conditional on parents' age/wage.
gen births = nkids36 - nkids1
bysort agebracket taxbracket: tab births

// (2) Quitting rate (of male and female,respectively) conditional on 
//     childbirth (say within three months), conditional on parents'age/wage 
//     (if possible. i.e. if we have enough sample size)
// Generating a quit variable that indicates if someone in the HH quit within
// three months of having a child.
gen quit = 0
gen aheadone = 0
gen aheadtwo = 0
gen aheadthree = 0
gen hadkid = 0

forvalues i = 1/5551{
	forvalues j = 1/36{
		replace aheadone = `j' + 1
		replace aheadtwo = `j' + 2
		replace aheadthree = `j' + 3
		replace hadkid = nkids`aheadone' - nkids`j' if `aheadone' <= 36
		}
		}

		if aheadthree <= 36 & aheadtwo <= 36 & aheadone <= 36{
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadone'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadone'[`i']) in `i'
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadtwo'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadtwo'[`i']) in `i'
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadthree'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadthree'[`i']) in `i'
		}
		else if aheadthree > 36 & aheadtwo <= 36 & aheadone <= 36{
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadone'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadone'[`i']) in `i'
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadtwo'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadtwo'[`i']) in `i'
		}
		else if aheadthree > 36 & aheadtwo > 36 & aheadone <= 36{
			replace quit = 1 if nkids`aheadone'[`i'] - nkids`j'[`i'] > 0 & (hus_emp`j'[`i'] > hus_emp`aheadone'[`i'] | wife_emp`j'[`i'] > wife_emp`aheadone'[`i']) in `i'
		}
	}
}

// (3) Mean leaving duration conditional on (2), and on parents'age/wage (again, if possible).
// (4) Wage difference before and after the leave (male and female)
// (5) Employment rate (of male and female respectively) conditional on parents' age/youngest child's age (assume parents' age are the same)
// (6) Mean and standard deviation of wages (both male and female) conditional on parents' age
// (7) Transition probability from ee to uu
// (8) Transition probability from eu to uu, and ee to ue without childbirth 
// (9) Transition probability from ue to uu, and ee to eu without childbirth
// (10) Transition probability from uu to ue, and eu to ee
// (11) Transition probability from uu to eu, and ue to ee
