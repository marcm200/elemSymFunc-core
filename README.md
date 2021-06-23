Specialized computer algebra system with interactive text mode designed to work with elementary symmetric functions
in the periodic points of a cycle of length p to construct integral polynomials in two variables - the product
of periodic points and the seed c - to describe the hyperbolic components of exact period p
in the quadratic Mandelbrot case z^2+c (see literature list at the end):

## Result

The final equation fulfills:

A complex number c is in a period-p component if the polynomial has a complex solution VARSYM
with ||VARSYM|| < 2^-p, where VARSYM is the notational variable symbol for the product of periodic points and
is here denoted by the p-th letter (uppercase) in the German alphabet.

## General

The software handles a set of equations (identified by a number id) where the left hand side is only consisting of
elementary symmetric functions A,B,C,... (sigma1, sigma2, sigma3 etc) as well as integral numbers and polynomials in the seed
c, whereas the right hand side consists of polynomial terms
in the roots of a cycle b_k. Elementary ysymmetric functions are denoted by capital letters and can
be manipulated, roots or the seed (lowercase) c not.

After every manipulation of an equation, it is automatically simplified by the following steps:
- division by the integral gcd
- terms only consisting of numbers or polynomials in the seed are moved to the left side
- shared c-factors of the form c^exponent - if occuring on both sides or one side is zero -  are factored out and discarded
- the right-hand side is analyzed and fully symmetric (sub)expressions in the roots are converted to terms of elementary symmetric functions and moved to the left hand side.

The goal is to finally have a set of equations where the rhs equals 0. Then all elementary symmetric function
variables are eliminated in various ways until an equation occurs that only consists
of the product of periodic points' variable (the p-th letter in the alphabet) and
the seed c.

## Quick start

Run the software and enter

<br>`skript(p4a.s`<br>
This constructs the equation for period-4, prints its final form into the log file (PRINT-statement in the skript)
and automatically creates an description file (DESCRIBEON-command) that explains the manipulations done.

Provided are scripts that perform a series of operations to arrive at the equations for
periods 4-8, named p4a.s, p5a.s etc. Time demands are seconds for periods 4-6, 20 min for
period-7 and about a day for period-8 on an AMD Ryzen 9 3900X.

## Commands

### Settings

`PERIOD(p` <br>
Sets the period to p and constructs the elementary symmetric functions (C++ array sigma[]). For longer peiods (>= 9) this may take some time.

`SIMPLIFYON` <br>
switches automatic simplification ON - which is the standard configuration.

`SIMPLIFYOFF` <br>
switches automatic simplification OFF - recommended when using SETSTR (see below)-

`QUIT`<br>

### Input/Output

`LOAD(fn`<br>
destroys all currently constructed equations and replaces them with the ones in file fn
The correct period to recognize the variables has to be set beforehand.

`SAVE(fn`<br>
stores the current list of equations in file 'fn

`LOADEQ(id,fn`<br>
loads the equation stored in file fn into identification number id.

`SAVEEQ(id,fn`<br>
stores equation id into file fn

`PRINT(id`<br>
prints the full equation into the log file for further usage.

### Creating equations

`ADDFSS(tgtid,factor,cexponent,bexponent0[,bexponent1,...]` <br>
creates the fully symmetric sum of the term factor * c^cexponent * b0^bexponent0 * b1^bexponent1*... and simplifies it. The result is stored under identification number 'id.

`CONSTRBINDIV(id,a0,a1,adiff,type,bexponent0,[bexponent1,...]` <br>
constructs a specialized function based on binomials (see the description in e.g. file p6a.txt for equation 10).

`SETSTR(id,text` <br>
sets equation id to: right hand side to 0, left hand side to the string text using a simple parser. Brackets are not allowed, variables are case-sensitive (uppercase C and lowercase c have different meanings) and have a text length of 2000 characters. If longer equations are to be used, they have to be manually split into parts, individually set and added (see ADDEQ command). It is recommended to switch automatic simplificatgion OFF (SIMPLIFYOFF) before the first SETSTR
and switch it back on after the last step (SETSTR or ADDEQ, see e.g. file p8a.s, command ADDEQ(5308,...)

### Manipulating equartions

`ADDEQ(tgtid,id0,factor0,cexponent0,id1,factor1,cexponent1` <br>
adds the two multiples factor0 * c^cexponent0 * eq(id0) + factor1 * c^cexponent1 * eq(id1) and stores the result in tgtid

`MULEQ(tgtid,id0,factor0,cexponent0,id1,factor1,cexponent1` <br>
multiplies the two multiples factor0 * c^cexponent0 * eq(id0) * factor1 * c^cexponent1 * eq(id1) and stores the result in tgtid

`DEFEQ(tgtid,id0,factor0,cexponent0,id1,factor1,cexponent1[,power]` <br>
for equations with empty right side, tgtid is set to factor0 * c^cexponent0 * eq(id0) / ( factor1 * c^cexponent1 * eq(id1)) ^ power if the polynomial division is remainder-free. If power is missing, 1 is taken. This is used to discard (intermediate) results that do not describe components of interest.

### Querying equations

`DEL(id[,id1]`<br>
deletes equation id or, if provided with two parameters, the entire id-range id..id1

`SIMP(id`<br>
analyzses the variables of the left-hand side (A-Z) with its coefficient polynomials. Result is 'linear if the highest-exponent of a variable is 1, simple if it is linear and the coefficient polynomial consists solely of integral numbers and a polynomial in the seed c. A variable is perfect if it is simple and no higher-alphabet-indexed variable occurs. This equation is usually then used to solve for that "highest-indexed" variable.

`COEFF(id,varsym,varexp`<br>
returns the coefficient polynomial of equation id of the variable-exponent combination varsym^varexp

`COEFFS(varsym,varexp`<br>
returns the coefficient polynomials for all available equations

`HIPOW(id,varsym`<br>
returns the highest exponent of varsym in equation id

`HIPOWS(varsym`<br>
returns the highest exponents of a variable symbol in all equations

### Elimination of variables
Works only if the right-hand side is equal to 0.

<br>`SOLVE(varsym,id`<br>
solves the equation id for variable symbol varsym. The right hand side must be zero, and the variable occuring in linear form only

`SOLT(tgtid,id0,varsym`<br>
substitutes the solution for a variable symbol varsym into equation id0, converts it to a polynomial equation by bringing
every term to the same denominator and discards it (as equation's right hand side is 0).

`SOLTALL(varsym,offset`<br>
substitutes the solution for variable symbol varsym into all equations currently present, sets the result
as equastion current_id+offset and deletes the current_id equation

`ELMT(tgtid,id0,varsym,varexp,id1`<br>
calculates the equation coeff1* eq(id0) - coeff0* eq(id1) were coeffN is the coefficient of equation N of
the variable-exponent combination varsym^varexp

`ELMTALL(varsym,varexp,id0,offset`<br>
eliminates the combination varsym^varexp of all equastions (except id0) against id0, sets the result
into current_id+offset and deletes the equation current_id

`RESTABV(tgtid,id0,id1,varsym`<br>
eliminates the variable varsym, computing the resultant of the two equadtions id0 and id1 and stores the result in equation tgtid. The resultant is computed using the Bareiss' 2-step algorithm starting from Sylvester's matrix

`RESTABVALL(tgtid,id0,varsym,offset`<br>
computes the resultant of id0 with all other equations to eliminate variable varsym. The result is stored in current_id+offset and equation current_id is deleted

### Certificate

`VALIDATE(id`<br>
performs two validation checks for equation id (rhs must be 0, only seed c and the highest-index variable denoting the
product of periodic points (say VARSYM) is allowed to occur): The constant term w.r.t VARSYM=0 must equal the polynomial describing the hyperbolic
centers of exact period p. The VARSYM-degree must equal the maximum number of cycles of exact length p for a given seed c.

`DESCRIBEON(fn`<br>
starts the description mode that explains the manipulations and stores it in file fn

`DESCRIBEOFF`<br>
ends the description mode and prints the literature list

`DESCR(te´xt`<br>
prints the text at the current description file position

### External

`XFACTOR2(tgtid,id0` needs the externally available program maxima to factor a polynomial (equation id0, right hand side must be 0). This requres
to adjust the path for maxima in the main.cpp sourcecode at positions MAXIMAPATH. The output of maxima is parsed and the factors
are set into the list of equations starting at id number tgtid.

## Limitations
- in case of an error, the software prints a message and mostly exits "dirty". Saving of intermediate results is recommended.
- currently, the maximum period length is 12 (C++ variable ANZSYMVARS) and can be increased up to 26. Currently (June 2021) I am trying to construct period-9.
- Input length of lines is 2000 characters. Especially for SETSTR, this requires to split the function to set into several parts and adding them. Care must be taken to switch off automatic simplification for those parts.
- factoring is currently performed externally - either manually or by the command XFACTOR2 which requires maxima and was curently tested with version 5.42.2.
- only the variables that describe elementary symmetric functions (uppercase letters A-Z) are accessible. The seed c or the
root variable s b_k are not.
- currently integral numbers up to ~10^4500 are supported (C++ constant MAXBIGINTDIGITS in file bigint.cpp, 9 * 512) but can be increased at compile time. Overflow is monitored. Lowering the value can lead to a speed gain but may render stored files e.g. SAVE or SAVEEQ invalid.

## Further information:
<br>[1] B Blum Smith, S Coskey. The fundamental theorem on symmetric polynomials. History's first whiff of Galois tehory. 2010.
<br>[2] K Conrad. Symmetric functions.
<br>[3] D Giarrusso, Y Fisher. A parameterization of the period 3 hyperbolic component of the Mandelbrot set. Proc Am Math Soc 123(12) 1995.
<br>[4] A Brown. Equations for periodic solutions of a logisitc difference equation. J Austral Math Soc 23 1981.
<br>[5] J Stephenson, T Ridgway. Formulae for cycles in the MandelBrot set. I, II, III. Physics A 1992.
<br>[6] A Healy. Resultants, resolvents and the computation of Galois group.
<br>[7] EH Bareiss. Sylvester's identity and multistep integer-preserving Gaussean elimination. 1968
<br>[8] https://fractalforums.org/fractal-mathematics-and-new-theories/28/explicit-equations-for-the-interior-of-hyperbolic-components/3463/msg29136#new

### Contact

Marc Meidlinger, June 2021<br>
marc.meidlinger@web.de



