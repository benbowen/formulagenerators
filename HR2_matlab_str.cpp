#include "mex.h"
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

string tmp_str = "hello";

#define VERSION "20050617"	/* String ! */
#define TRUE 	1
#define FALSE 	0
#define MAXLEN  181          /* max. length of input string */
#define _CRT_SECURE_NO_DEPRECATE 1

typedef struct 	{
		const char *sym;	/* symbol */
		const double mass;	/* accurate mass */
		const float val;	/* to calculate unsaturations */
		const int key;		/* used for decoding cmd line */
		int min,		/* atom count min */
		    max,		/* atom count max */
		    cnt,		/* atom count actual */
		    save;		/* atom count old  - for loop exiting*/
		} Element;


/* --- atomic masses as published by IUPAC, 2002-10-02 ---------- */

Element el[]=
/* array of the elements used here:
   Symbol, exact mass, dbe, keycode, number, default-min, default-max
*/
{
{ "C",  12.000000000,   +2.0, 'C', 0, 41, 0,0 },
{ "13C", 13.0033548378, +2.0, '1', 0, 0, 0 ,0 },
{ "H",   1.0078250321,  -1.0, 'H', 0, 72, 0 ,0},
{ "D",   2.0141017780,  -1.0, 'D', 0, 0, 0 ,0},
{ "N",  14.0030740052,  +1.0, 'N', 0, 34, 0,0 },		//org +1 = valence = 3: now +3 for valence = 5
{ "15N", 15.0001088984,   +1.0, 'M', 0, 0, 0,0 },
{ "O",  15.9949146221,   0.0, 'O', 0, 30, 0 ,0},
{ "F",  18.99840320,    -1.0, 'F', 0, 0, 0 ,0},
{ "Na", 22.98976967,    -1.0, 'A', 0, 0, 0 ,0},
{ "Si", 27.9769265327,  +2.0, 'I', 0, 0, 0 ,0},	
{ "P",  30.97376151,    +3.0, 'P', 0, 0, 0 ,0},		//org +1 valence = 3: now +3 for valence = 5
{ "S",  31.97207069,    +4.0, 'S', 0, 0, 0 ,0},		//org 0 = valence = 2; now +4 for valence = 6
{ "Cl", 34.96885271,    -1.0, 'L', 0, 0, 0 ,0},
{ "Br", 78.9183376,     -1.0, 'B', 0, 0, 0 ,0},
};

const double electron = 0.000549;	/* mass of the electron in amu */


/* --- global variables --- */

double  charge,		/* charge on the molecule */
        tol;		/* mass tolerance in mmu */
char    comment[MAXLEN]="";	/* some text ;-) */
int     single;		/* flag to indicate if we calculate only once and exit */
int     nr_el;		/* number of elements in array (above) */


/* --- some variables needed for reading the cmd line --- */

char   *optarg;		/* global: pointer to argument of current option */

//lear: because input_pr starts from zero position, I changed optind to 0
//attention!!!
int optind = 0;	/* global: index of which argument is next. Is used
                as a global variable for collection of further
                arguments (= not options) via argv pointers. */


/* --- function prototypes ------------------- */

int     input(char *text, double *zahl);
int     readfile(char *whatfile);
double  calc_mass(void);
float   calc_rdb(void);
//char*     do_calculations(double mass, double tolerance);
string*     do_calculations(double mass, double tolerance);
int     clean (char *buf);
int     getopt(int argc, char *argv[], char *optionS);
//you have to compile with C++ or define yourself this bool type (C99 compiler definition)
bool calc_element_ratios(bool element_probability);

// the funcitons is used to convert int to string

string to_string( int val )
{
stringstream stream;
stream << val;
return stream.str();
}

/***************************************************************************
* INPUT:	reads a dataset in "dialog mode".			   *
* Input: 	Pointer to comment text, pointer to mass.		   *
* Returns:	Result of sscanf(), 0 if prog is to be finished. 	   *
* Note:		This is a fairly primitive way to do it, but it works ;-)  *
****************************************************************************/
int input(char *txt, double *zahl)
{
char buf[MAXLEN];				/* input line */

*zahl = 0.0;				    /* reset */

mexPrintf("Comment             : ");	/* display prompt */
fgets(buf, MAXLEN-1, stdin);	/* read line */
buf[MAXLEN] = 0x0;				/* terminate string */
clean (buf);				    /* remove linefeed */
strcpy(txt, buf);			    /* copy text over */

mexPrintf("Mass (ENTER to quit): ");	/* display prompt */
fgets(buf, MAXLEN-1, stdin);       	/* read line */
buf[MAXLEN] = 0x0;			    /* terminate string */
if (!clean (buf))			    /* only a CR ? --> quit */
	return 0;
sscanf(buf,"%lf", zahl);		/* scan string */

return 1;
}


/***************************************************************************
* READFILE:	reads dataset from file.				   *
* Input: 	Pointer to comment text, pointer to mass.		   *
* Returns:	0 if OK, 1 if error.					   *
****************************************************************************/
int readfile(char *whatfile)
{
double mz;		/* measured mass */
char buf[MAXLEN];		/* input line */
FILE *infile;

infile = fopen(whatfile, "r");
if (NULL == infile)
	{
	mexPrintf ("Error: Cannot open.");
	//mexPrintf (stderr, "Error: Cannot open %s.", whatfile);
	return 1;
	}

while (fgets(buf, MAXLEN-1, infile))
	{
	buf[MAXLEN] = 0x0;			/* terminate string */
	if (*buf == ';')		/* comment line */
		continue;
	if (!clean (buf))		/* only a CR ? --> quit */
		return 0;
	sscanf(buf,"%s %lf", comment, &mz);	/* scan string */
	do_calculations(mz, tol);
	mz = 0.0;				/* reset */
	}
return 0;
}


/************************************************************************
* CALC_MASS:	Calculates mass of an ion from its composition.	  	*
* Input: 	nothing (uses global variables) 	      		*
* Returns. 	mass of the ion.	  				*
* Note:		Takes care of charge and electron mass!   		*
* 		(Positive charge means removal of electrons).	 	*
*************************************************************************/
double calc_mass(void)
{
int i;
double sum = 0.0;

for (i=0; i < nr_el; i++)
	sum += el[i].mass * el[i].cnt;

return (sum - (charge * electron));
}


/************************************************************************
* CALC_RDB:	Calculates rings & double bond equivalents.    		*
* Input: 	nothing (uses global variables)			   	*
* Returns. 	RDB.				       			*
*************************************************************************/
float calc_rdb(void)
{
int i;
float sum = 2.0;

for (i=0; i < nr_el; i++)
	sum += el[i].val * el[i].cnt;

return (sum/2.0);
}
/************************************************************************
* Calculates element ratios , CH2 (more than 8 electrons needed is not handled)  		
* Calculations element probabilities if element_probability = true 
* Input: 	nothing (uses global variables)			   	
* Returns. true/false.				       			
*************************************************************************/
bool calc_element_ratios(bool element_probability)
{
bool CHNOPS_ok;	
float HC_ratio;
float NC_ratio;
float OC_ratio;
float PC_ratio;
float SC_ratio;

float C_count = (float)el[0].cnt;
float H_count = (float)el[2].cnt;
float N_count = (float)el[4].cnt;
float O_count = (float)el[6].cnt;
float P_count = (float)el[10].cnt;
float S_count = (float)el[11].cnt;


		/* ELEMENT RATIOS allowed
			MIN		MAX (99.99%)
		H/C	0.07	6.00
		N/C	0.00	4.00
		O/C	0.00	3.00
		P/C	0.00	2.00
		S/C	0.00	6.00
		*/	

	// set CHNOPS_ok = true and assume all ratios are ok
	CHNOPS_ok = true;	
	
	
	if (C_count && H_count >0)					// C and H  must have one count anyway (remove for non-organics//
	{	
		HC_ratio = H_count/C_count;
		if (element_probability)
		{
			if ((HC_ratio <  0.2) || (HC_ratio >  3.0)) // this is the H/C probability check ;
			CHNOPS_ok = false;
		}
		else if (HC_ratio >  6.0) // this is the normal H/C ratio check - type cast from int to float is important
			CHNOPS_ok = false;
	}

	if (N_count >0)	// if positive number of nitrogens then thes N/C ratio else just calc normal
	{
		NC_ratio = N_count/C_count;
		if (element_probability)
		{
			if (NC_ratio >  2.0) // this is the N/C probability check ;
			CHNOPS_ok = false;
		}
		else if (NC_ratio >  4.0)
			CHNOPS_ok = false;
	}	
	
	if (O_count >0)	// if positive number of O then thes O/C ratio else just calc normal
	{	
		OC_ratio = O_count/C_count;
		if (element_probability)
		{
			if (OC_ratio >  1.2) // this is the O/C  probability check ;
			CHNOPS_ok = false;		
		}
		else if (OC_ratio >  3.0)
				CHNOPS_ok = false;
	}	


	if (P_count >0)	// if positive number of P then thes P/C ratio else just calc normal
	{	
		PC_ratio = 	P_count/C_count;
		if (element_probability)
		{
			if (PC_ratio >  0.32) // this is the P/C  probability check ;
			CHNOPS_ok = false;	
		
		}
		else if (PC_ratio >  6.0)
			CHNOPS_ok = false;
	}	

	if (S_count >0)	// if positive number of S then thes S/C ratio else just calc normal
	{	
		SC_ratio = 	S_count/C_count;
		if (element_probability)
		{
			if (SC_ratio >  0.65) // this is the S/C  probability check ;
			CHNOPS_ok = false;	
		}
		else if (SC_ratio >  2.0)
			CHNOPS_ok = false;
	}	

//-----------------------------------------------------------------------------	
		
	// check for multiple element ratios together with probability check 
	//if N<10, O<20, P<4, S<3 then true
	if (element_probability && (N_count > 10) && (O_count > 20) && (P_count > 4) && (S_count > 1))
		CHNOPS_ok = false;	
	
	// NOP check for multiple element ratios together with probability check
	// NOP all > 3 and (N<11, O <22, P<6 then true)
	if (element_probability && (N_count > 3) && (O_count > 3) && (P_count > 3))
		{
		if (element_probability && (N_count > 11) && (O_count > 22) && (P_count > 6))
			CHNOPS_ok = false;	
		}
	
	// OPS check for multiple element ratios together with probability check
	// O<14, P<3, S<3 then true
	if (element_probability && (O_count > 14) && (P_count > 3) && (S_count > 3))
		CHNOPS_ok = false;	

	// PSN check for multiple element ratios together with probability check
	// P<3, S<3, N<4 then true
	if (element_probability && (P_count > 3) && (S_count > 3) && (N_count >4))
		CHNOPS_ok = false;	

	
	// NOS check for multiple element ratios together with probability check
	// NOS all > 6 and (N<19 O<14 S<8 then true)
	if (element_probability && (N_count >6) && (O_count >6) && (S_count >6))
	{
		if (element_probability && (N_count >19) && (O_count >14) && (S_count >8))
			CHNOPS_ok = false;	
	}	


	// function return value;
	if (CHNOPS_ok == true)
		return true;
	else 
		return false;
}

/************************************************************************
* DO_CALCULATIONS: Does the actual calculation loop.			*
* Input: 	   measured mass (in amu), tolerance (in mmu)	    	*
* Returns. 	   number of hits.	       				*
*************************************************************************/
//char* do_calculations (double measured_mass, double tolerance)
string* do_calculations (double measured_mass, double tolerance)
{
 int j;
  //string tmp_string;
  string output_string = "";

  string* output_str_arr = new string [1000000]; 
  int arr_idx = 0;

  /* char* part_str = "hello";
 part_str = "wossshshh";
 int part_int = 5;
 tmp_string = to_string(part_int);
 tmp_string = tmp_string + "hello";
 mexPrintf("the number is %s \n",tmp_string.c_str());
 tmp_string = string(part_str);
  */
time_t start, finish;
double elapsed_time;

double mass;			/* calc'd mass */
double limit_lo, limit_hi;	/* mass limits */
float rdb, lewis;			/* Rings & double bonds */
long i;			
long long hit;		/* counts the hits, with long declaration, overflow after 25h with all formulas < 2000 Da
							    long = FFFFFFFFh = 4,294,967,295d*/
long long counter;
bool elementcheck;
bool set_break; 

time( &start );		// start time
mexPrintf("\n");		/* linefeed */

/* calculate limits */

limit_lo = measured_mass - (tolerance / 1000.0);
limit_hi = measured_mass + (tolerance / 1000.0);

//if (strlen(comment))	/* print only if there is some text to print */
//	mexPrintf ("Text      \t%s\n", comment);

//mexPrintf ("Composition\t");
//for (i=0; i < nr_el; i++)
//	if (el[i].max > 0)
//		mexPrintf("%s:%d-%d ", el[i].sym, el[i].min, el[i].max);
//mexPrintf ("\n");

//mexPrintf ("Tol (mmu)\t%.1f\n",tolerance);
//mexPrintf ("Measured\t%.4lf\n", measured_mass);
//mexPrintf ("Charge  \t%+.1lf\n", charge);

hit = 0;			/* Reset counter */
counter = 0;
set_break = false;	/* set breaker for element counts to false */

/* Now let's run the big big loop ... I'd like to do that
   recursively but did not yet figure out how ;-) 
   TK Adds: the loop is just fine.
*/

/* now comes the "COOL trick" for calculating all formulae:
sorting the high mass elements to the outer loops, the small weights (H)
to the inner loops;

This will reduce the computational time by factor ~10-60-1000
OLD HR: Cangrelor at 1ppm  4465 formulas found in   5866 seconds.
NEW HR2: Cangrelor at 1ppm 4465 formulas found in     96 seconds.
NEW2 HR2: Cangrelor at 1ppm 4465 formulas found in     60 seconds.
NEW3 HR2: Cangrelor at 1ppm 4465 formulas found in     59 seconds.
HR2 Fast: Cangrelor at 1ppm 4465 formulas found in     41 seconds by evaluating 2,003,436,894 formulae.
hr2 -c "Cangrelor" -m  774.948 -t 0.77 -C 1-64 -H 1-112 -N 0-30 -O 0-80 -P 0-12 -S 0-9 -F 0-10 -L 0-10

Another additional trick is to end the 2nd.. 3rd.. 4th.. xth innermost loop
to prevent loops which are just higher and higher in mass.
*/

el[13].cnt = el[13].min - 1;  el[13].save = el[13].cnt; 
while (el[13].cnt++ < el[13].max) /* "Br"*/ { 

el[12].cnt = el[12].min - 1;  el[12].save = el[12].cnt; 
while (el[12].cnt++ < el[12].max) /*"Cl"*/ { 
	 
el[11].cnt = el[11].min - 1;  el[11].save = el[11].cnt; 
while (el[11].cnt++ < el[11].max) /*"S"*/ { 
	 
el[10].cnt = el[10].min - 1;  el[10].save = el[10].cnt; 
while (el[10].cnt++ < el[10].max) /*"P"*/ { 
	 
el[9].cnt = el[9].min - 1;  el[9].save = el[9].cnt; 
while (el[9].cnt++ < el[9].max) /*"Si"*/ { 

el[8].cnt = el[8].min - 1;  el[8].save = el[8].cnt; 
while (el[8].cnt++ < el[8].max) /*"Na"*/{ 

el[7].cnt = el[7].min - 1;  el[7].save = el[7].cnt; 
while (el[7].cnt++ < el[7].max) /*"F"*/ { 
 
el[6].cnt = el[6].min - 1;  el[6].save = el[6].cnt; 
while (el[6].cnt++ < el[6].max) /*"O"*/ { 
	 
el[5].cnt = el[5].min - 1;  el[5].save = el[5].cnt; 
while (el[5].cnt++ < el[5].max) /*"15N"*/{ 

el[4].cnt = el[4].min - 1; el[4].save = el[4].cnt; 
while (el[4].cnt++ < el[4].max) /*"N"*/{ 
	 
el[1].cnt = el[1].min - 1; el[1].save = el[1].cnt; 
while (el[1].cnt++ < el[1].max) /*"13C"*/ { 

el[0].cnt = el[0].min - 1; el[0].save = el[0].cnt; 
while (el[0].cnt++ < el[0].max) /* "C"*/ { 

el[3].cnt = el[3].min - 1; 	el[3].save = el[3].cnt; 
while (el[3].cnt++ < el[3].max) /*"D"*/{ 

el[2].cnt = el[2].min - 1; el[2].save = el[2].cnt; 
while (el[2].cnt++ < el[2].max) /*"H"*/{ 

	mass = calc_mass();
	counter++;

	//just for debug purposes
	//if (mass > limit_hi)  
	//mexPrintf("mass: %f\tC: %d  H: %d  N: %d O: %d P: %d S: %d Cl: %d Br: %d\n",mass,el[0].cnt,el[2].cnt,el[4].cnt,el[6].cnt,el[10].cnt,el[11].cnt,el[12].cnt,el[13].cnt);
    	
	/* if we exceed the upper limit, we can stop the calculation
       for this particular element (JHa 20050227). <-- comment TK that will only bust the innermost while loop, which is "H"*/

	// break H loop 
	if (mass > limit_hi)  break;

    //************************************************************************************************************/	
	//Calculus loop with print out
	//************************************************************************************************************/	

	if ((mass >= limit_lo) && (mass <= limit_hi)) /* within limits? */
	{	
		// element check will be performed always, if variable bool element_probability is true also probabilities will be calculated
		// not an elegant implementation, but fast.
		 elementcheck = calc_element_ratios(true);
		 if (elementcheck)
	{ 
	rdb = calc_rdb();	/* get RDB */
	lewis = (float)(fmod(rdb, 1)); /*calc reminder*/
	if ((rdb >= 0) && (lewis != 0.5) && (lewis !=-0.5))/* less than -0.5 RDB does not make sense */
		{													/* NO(!) CH3F10NS2 exists , RDB =  -4.0   M= 282.9547*/
		output_string = "";
		hit ++;
		for (i = 0; i < nr_el; i++)			/* print composition */
		    if (el[i].cnt > 0)				/* but only if useful */
			 
			 //! comment the following two printlines and add a single ; in a free line for the counting version
		      output_string = output_string + string(el[i].sym) + to_string(el[i].cnt);
		      //mexPrintf("the true formula is %s and the idx is %d \n", output_string.c_str(),arr_idx);
		      output_str_arr[arr_idx] = output_string;
		      //mexPrintf("after this step, the formula is %s\n", output_str_arr[arr_idx].c_str());
		      arr_idx = arr_idx + 1;
		      
		      //mexPrintf("%s%d", el[i].sym, el[i].cnt);	// print formula
			  //mexPrintf("\t\t%.1f\t%.4lf\t%+.1lf mmu \n", rdb, mass, 1000.0 * (measured_mass - mass));
			  
		}	/* end of 'rdb' loop */
	}	// end of elementcheck loop
	
	}	/* end of 'limit' loop */
	//************************************************************************************************************/
		

	/*
	TK: if the current mass is larger than the limit the loop can be exited.
	Each element must point to the element which is in use and before.
	This is a static implementation which can be enhanced with a pointer chain to the lower element.
	Actually now its only allowed for CHNSOP-Fl-Cl-Br-Si !!! Brute-force <> elegance :-)
	*/
		} /*"H"*/
		
		} /*"D"*/
		
		if ((mass >= limit_lo) && (el[2].save == el[2].cnt-1)) break;
		} /* "C"*/
		
		} /*"13C"*/

		if ((mass >= limit_lo) && (el[0].save == el[0].cnt-1)) break;
		} /*"N"*/
		
		} /*"15N"*/

        if ((mass >= limit_lo) && (el[4].save == el[4].cnt-1)) break;
		} /*"O"*/
		
	    if ((mass >= limit_lo) && (el[6].save == el[6].cnt-1)) break;
		} /*"F"*/
		
		} /*"Na"*/
		
	    if ((mass >= limit_lo) && (el[7].save == el[7].cnt-1)) break;
		}  /*"Si"*/
		
		if ((mass >= limit_lo) && (el[9].save == el[9].cnt-1)) break;
		} /*"P"*/
		
		if ((mass >= limit_lo) && (el[10].save == el[10].cnt-1)) break;
		} /*"S"*/
		
		if ((mass >= limit_lo) && (el[11].save == el[11].cnt-1)) break;
		} /*"Cl"*/
		
		if ((mass >= limit_lo) && (el[12].save == el[12].cnt-1)) break;
		} /*"Br" ends*/
/* close that giant loop thing started above */

	

time(&finish);		// stop timer
elapsed_time = difftime(finish , start);	// calulate time difference

//if (!hit)
	//mexPrintf("No matching combination found in %6.0f seconds.\n", elapsed_time );
//else
	//mexPrintf("\n%llu formulas found in %6.0f seconds by evaluating %llu formulae.\n",hit,elapsed_time,counter);
	//mexPrintf("RDBs are overloaded to maximum valence values (N=3,P=5,S=6).\n");
	
//return hit;
 //mexPrintf("%s\n", output_string.c_str());
/*
for (j=0;j<3;j++){
	mexPrintf("At the end of do_cal, the formula is: %s \n",output_str_arr[j].c_str());
}
*/
output_str_arr[arr_idx] = "";
return output_str_arr;
}


/************************************************************************
* CLEAN:	"cleans" a buffer obtained by fgets() 			*
* Input: 	Pointer to text buffer					*
* Returns:	strlen of buffer.   					*
*************************************************************************/
int clean (char *buf)
{
int i;

for(i = 0;  i < strlen(buf);  i++)		/* search for CR/LF */
	{
	if(buf[i] == '\n' || buf[i] == '\r')
		{
		buf[i] = 0;			/* stop at CR or LF */
		break;
		}
	}
return (strlen(buf));
}



/***************************************************************************
* GETOPT: Command line parser, system V style.
*
*  This routine is widely (and wildly) adapted from code that was
*  made available by Borland International Inc.
*
*  Standard option syntax is:
*
*    option ::= SW [optLetter]* [argLetter space* argument]
*
*  where
*    - SW is '-'
*    - there is no space before any optLetter or argLetter.
*    - opt/arg letters are alphabetic, not punctuation characters.
*    - optLetters, if present, must be matched in optionS.
*    - argLetters, if present, are found in optionS followed by ':'.
*    - argument is any white-space delimited string.  Note that it
*      can include the SW character.
*    - upper and lower case letters are distinct.
*
*  There may be multiple option clusters on a command line, each
*  beginning with a SW, but all must appear before any non-option
*  arguments (arguments not introduced by SW).  Opt/arg letters may
*  be repeated: it is up to the caller to decide if that is an error.
*
*  The character SW appearing alone as the last argument is an error.
*  The lead-in sequence SWSW ("--") causes itself and all the rest
*  of the line to be ignored (allowing non-options which begin
*  with the switch char).
*
*  The string *optionS allows valid opt/arg letters to be recognized.
*  argLetters are followed with ':'.  Getopt () returns the value of
*  the option character found, or EOF if no more options are in the
*  command line. If option is an argLetter then the global optarg is
*  set to point to the argument string (having skipped any white-space).
*
*  The global optind is initially 1 and is always left as the index
*  of the next argument of argv[] which getopt has not taken.  Note
*  that if "--" or "//" are used then optind is stepped to the next
*  argument before getopt() returns EOF.
*
*  If an error occurs, that is an SW char precedes an unknown letter,
*  then getopt() will return a '~' character and normally prints an
*  error message via perror().  If the global variable opterr is set
*  to false (zero) before calling getopt() then the error message is
*  not printed.
*
*  For example, if
*
*    *optionS == "A:F:PuU:wXZ:"
*
*  then 'P', 'u', 'w', and 'X' are option letters and 'A', 'F',
*  'U', 'Z' are followed by arguments. A valid command line may be:
*
*    aCommand  -uPFPi -X -A L someFile
*
*  where:
*    - 'u' and 'P' will be returned as isolated option letters.
*    - 'F' will return with "Pi" as its argument string.
*    - 'X' is an isolated option.
*    - 'A' will return with "L" as its argument.
*    - "someFile" is not an option, and terminates getOpt.  The
*      caller may collect remaining arguments using argv pointers.
***************************************************************************/
int getopt(int argc, char *argv[], const char *optionS)
{
static char *letP	= NULL;		/* remember next option char's location */
static char SW		= '-';		/* switch character */

int opterr = 1;				/* allow error message	*/
unsigned char ch;
char *optP;
//tmp modification by lear: argc -> argc - 1
//mexPrintf("argc is %d \n",argc);
//mexPrintf("optind is %d \n",optind);
 if ( argc > optind)
	{
	if (letP == NULL)
		{
		if ((letP = argv[optind]) == NULL || *(letP++) != SW)
			goto gopEOF;
		if (*letP == SW)
			{
			optind++;
			goto gopEOF;
			}
		}
	if (0 == (ch = *(letP++)))
		{
		optind++;
		goto gopEOF;
		}
	if (':' == ch  ||  (optP = strchr(optionS, ch)) == NULL)
		goto gopError;
	if (':' == *(++optP))
		{
		optind++;
		if (0 == *letP)
			{
			if (argc <= optind)
				goto  gopError;
			letP = argv[optind++];
			}
		optarg = letP;
		letP = NULL;
	}
	else
	{
	if (0 == *letP)
		{
		optind++;
		letP = NULL;
		}
	optarg = NULL;
	}
	return ch;
}

gopEOF:
	optarg = letP = NULL;
	return EOF;

gopError:
	optarg = NULL;
	errno  = EINVAL;
	if (opterr)
		perror ("Command line option");
	return ('~');
}

/*
List of elements sorted according to mass
_______________________
INo#	El	Mass
2		H	1.007825032
3		D	2.014101778
0		C	12
1		13C	13.00335484
4		N	14.00307401
5		15N	15.0001089
6		O	15.99491462
7		F	18.9984032
8		Na	22.98976967
9		Si	27.97692653
10		P	30.97376151
11		S	31.97207069
12		Cl	34.96885271
13		Br	78.9183376
------------------------
*/

/* --- main_function --- */

string* main_function (int argc, char **argv)
{
  string* zero_pr = 0;
  string* output_str_arr;
  //int maxFM = 100;
  //output_str_arr = new char*[maxFM];
//mexPrintf("I am in the loop \n");
//flush();
double mz;	/* mass */
char buf[MAXLEN];
int i, tmp;
int j;

char *id =
"hr version %s. Copyright (C) by Joerg Hau 2001...2005 & Tobias Kind 2006 :-).\n";

char *disclaimer =
"\nThis program is free software; you can redistribute it and/or modify it under\n"
"the terms of version 2 of the GNU General Public License as published by the\n"
"Free Software Foundation.\n\n"
"This program is distributed in the hope that it will be useful, but WITHOUT ANY\n"
"WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
"PARTICULAR PURPOSE. See the GNU General Public License for details.\n";

char *msg =
"Calculates possible elemental compositions for a given mass.\n\n"
"usage: hr [options] file\n\nValid command line options are:\n"
"-h      This Help screen.\n"
"-v      Display version information.\n"
"-t tol  Set tolerance to 'tol' mmu (default 5).\n"
"-m mz   Set mass to 'mz'.\n"
"-c txt  Set comment to 'txt' (only useful together with '-m').\n"
"-p      Positive ions; electron mass is removed from the formula.\n"
"-n      Negative ions; electron mass is added to the formula.\n"
"-X a-b  For element X, use atom range a to b. List of valid atoms:\n\n"
"           X    key   mass (6 decimals shown)\n"
"        -------------------------------------\n";

/* initialise variables */

single = FALSE;			/* run continuously */
charge = 0.0;	       	 	/* default charge is neutral */
tol = 5.0;			/* default tolerance in mmu */
nr_el = sizeof(el)/sizeof(el[0]);	/* calculate array size */


/* decode and read the command line */
//mexPrintf("I am before the loop");
 while ((tmp = getopt(argc, argv, "hvpnt:m:c:C:H:N:M:O:D:1:S:F:L:B:P:I:A:")) != EOF)
   {
	switch (tmp)
		{
		case 'h':     	  		/* help me */
			mexPrintf (id, VERSION);
            mexPrintf (msg);
			for (i=0; i < nr_el; i++)
				mexPrintf ("        %4s     -%c %15.6lf\n",
				el[i].sym, el[i].key, el[i].mass);
			mexPrintf(disclaimer, "\n");
			return zero_pr;
		case 'v':     	  		/* version */
			mexPrintf (id, VERSION);
			return zero_pr;
		case 'p':    			/* positive charge */
			charge = +1.0;
			continue;
		case 'n':			/* negative carge */
			charge = -1.0;
			continue;
		case 't':       		/* tolerance */
			strcpy(buf, optarg);
			sscanf(buf, "%lf", &tol);
			continue;
		case 'm':			/* single mass */
			strcpy(buf, optarg);
		        sscanf(buf, "%lf", &mz);
		        single = TRUE;
	        	continue;
		case 'c':			/* comment for single mass */
	   		strcpy(comment, optarg);
		        continue;
		case 'C':      		/* C12 */
 		case 'H':      		/* 1H */
		case 'N':      		/* 14N */
		case 'M':      		/* 15N */
		case 'O':      		/* 16O */
 		case 'D':      		/* 2H */
 		case '1':      		/* 13C */
		case 'A':      		/* Na ('N' is taken for Nitrogen!) */
		case 'S':      		/* 32S */
		case 'F':      		/* 19F */
		case 'L':      		/* 35Cl ('C' is taken!) */
		case 'B':      		/* 79Br */
		case 'P':      		/* 31P */
		case 'I':      		/* 28Si ('S' is taken!) */
		       	i = 0;
			/* compare keys until found */
			while ((i < nr_el) && (el[i].key != tmp))
				i++;
			strcpy(buf, optarg);
			sscanf(buf, "%d-%d", &el[i].min, &el[i].max);	/* copy over */
			if (el[i].min > el[i].max)			/* swap them */
				{
				tmp = el[i].min;
				el[i].min = el[i].max;
				el[i].max = tmp;
				}

			/* mexPrintf ("\n %c = %c ... %s (%d-%d)", \
			tmp, el[i].key, el[i].sym, el[i].min, el[i].max); */
			//mexPrintf("I am here and the value of optind is %d \n", optind);
			//mexPrintf("I am here and single is %d \n", single);
			continue;
		case '~':    	  	/* invalid arg */
		default:
			mexPrintf ("'%s -h' for help.\n", argv[0]);
			return zero_pr;
		}
   }

  //if (argv[optind] != NULL)	 /* remaining parameter on cmd line? */
	/* must be a file -- treat it line by line */
	//mexPrintf("optind is %d \n", optind);
 	//return (readfile (argv[optind]));
 
  if (single == TRUE){  	   	 	/* only one calculation requested? */
    //mexPrintf("single is true \n");
 	output_str_arr = do_calculations(mz, tol);       /* do it, then exit ... */
         }
 else
 	{				/* otherwise run a loop */
 	while (input(comment, &mz))
 		{
		do_calculations(mz, tol);
 		mexPrintf("\n");
 		}
 	}
 //mexPrintf("I am at the end \n");
/*
for (j=0;j<3;j++){
	mexPrintf("In the main(), the formula is: %s \n",output_str_arr[j]);
}
*/
return output_str_arr;
}

void mexFunction(
		 int          nlhs,
		 mxArray      **plhs,
		 int          nrhs,
		 const mxArray **prhs
		 )
{
  char *input_pr[nrhs];
  //int maxFM = 100;
  string* final_output;

  int i, tmp;
  int result_cnt =0;
//in order to create a cell of strings output:
 
  mwSize ndim = 1;
   mxArray *mx_arr;

  for (i=0; i<nrhs; i++){
    input_pr[i] = (char*) mxArrayToString(prhs[i]);
  }
  /*
  while ((tmp = getopt(nrhs, input_pr, "hvpnt:m:c:C:H:N:M:O:D:1:S:F:L:B:P:I:A:")) != EOF){
  mexPrintf("the number is: %d \n", tmp);
  }
  */
  final_output = main_function(nrhs,input_pr);

  string* p_itr = &final_output[0]; 
  while(*p_itr != ""){
    p_itr++;
    result_cnt++;
  }

  const int dim_int [1] = {result_cnt};
  const mwSize *dims = dim_int;
  mx_arr = mxCreateCellArray(ndim,dims);

  //mexPrintf("The size is:%d \n",result_cnt);
  for (i = 0; i < result_cnt; i++){
    mxSetCell(mx_arr,i,mxCreateString(final_output[i].c_str()));
  }
  plhs[0] = mx_arr;
  //mexPrintf("the result is:%s \n",final_output);
 //lear: because main() is now a funciton in matlab, each time the program runs, the value of optind will not be erased. That's why we need to manually set it back to zero
  optind = 0;
  return;
}
