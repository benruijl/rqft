
#define PSI "1337"
#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define H "25"
#define GHO "82"
#define GHOBAR "-82"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"
#define QMASSLESS "1,2,3,4,5"
#define QBARMASSLESS "-1,-2,-3,-4,-5"
#define QMASSIVE "6,"
#define QBARMASSIVE "-6,"
#define L "11,12,13"
#define LBAR "-11,-12,-13"
#define BACKGROUND "210"

S vev, pi;

Auto S mass, mUV;
Auto S yukawa;
CTable masses(-10000:10000);
CTable gyq(-10000:10000);
CTable logmasses(-10000:10000);
CTable charges(-10000:10000);

Fill gyq(1) = yukawad; * d
Fill gyq(2) = yukawau; * u
Fill gyq(3) = yukawas; * s
Fill gyq(4) = yukawac; * c
Fill gyq(5) = yukawab; * b
Fill gyq(6) = yukawat; * t
Fill gyq(11) = 0; * e-
Fill gyq(12) = 0; * mu-
Fill gyq(13) = 0; * ta-
Fill gyq(-1) = yukawad; * d
Fill gyq(-2) = yukawau; * u
Fill gyq(-3) = yukawas; * s
Fill gyq(-4) = yukawac; * c
Fill gyq(-5) = yukawab; * b
Fill gyq(-6) = yukawat; * t
Fill gyq(-11) = 0; * e+
Fill gyq(-12) = 0; * mu+
Fill gyq(-13) = 0; * ta+

Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = massb;
Fill masses(6) = masst;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = massb;
Fill masses(-6) = masst;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;

Fill masses(-82) = 0;
Fill masses(82) = 0;
Fill masses(21) = 0;
Fill masses(22) = 0;
Fill masses(210) = 0;
Fill masses(25) = massh;
Fill masses(1337) = 0;

* note: this set needs to be complete for the UV expansion
Set allmasses: massu, massd, massc, masss, masst, massb, masse, massmu, masstau, massh, massw, massz, mUVexp;

Fill charges(1) = -1/3; * d
Fill charges(2) = 2/3; * u
Fill charges(3) = -1/3; * s
Fill charges(4) = 2/3; * c
Fill charges(5) = -1/3; * b
Fill charges(6) = 2/3; * t
Fill charges(11) = -1; * e
Fill charges(-1) = 1/3; * d
Fill charges(-2) = -2/3; * u
Fill charges(-3) = 1/3; * s
Fill charges(-4) = -2/3; * c
Fill charges(-5) = 1/3; * b
Fill charges(-6) = -2/3; * t
Fill charges(-11) = 1; * e
*--#] setup :

*--#[ feynman-rules :

************************************************
* Substitute the Feynman rules for the numerator
************************************************

#procedure FeynmanRulesGlobal()
* extract the global factors from the Feynman rules, including colour

id vx(`BACKGROUND', `GLU', `GLU', `GLU', ?a) =  vx(`GLU', `GLU', `GLU', `GLU',?a); * TODO: check
id vx(x1?{`QBAR'}, `BACKGROUND', x2?{`Q'}, ?a) = vx(x1, `GLU', x2, ?a); * TODO: check

* first, we split up the quartic gluon vertex into distinct colour factors
* we also generate an extra dummy index
Multiply counter(1);
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(idx5?) = counter(idx5 + 1) *(
    +vx(`GLU', `GLU', `GLU', `GLU', 1, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
    +vx(`GLU', `GLU', `GLU', `GLU', 2, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
    +vx(`GLU', `GLU', `GLU', `GLU', 3, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
);
repeat id vx(`H', `GLU', `GLU', `GLU', `GLU', p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?)*counter(idx6?) = counter(idx6 + 1) * (
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
);

repeat id vx(`BACKGROUND', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(idx7?) = counter(idx7 + 1) *(
    +vx(`BACKGROUND', `GLU', `GLU', `GLU', p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx7)
);

id counter(x?) = 1;

* make a copy of the Feynman rules
id prop(?a) = prop(?a)*tmp(prop(?a));
id vx(?a) = vx(?a)*tmp(vx(?a));
repeat id tmp(x1?)*tmp(x2?) = tmp(x1*x2);

* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* projectors
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = 1;
repeat id prop(`BACKGROUND', in, p?, idx1?)*prop(`BACKGROUND', out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = 1;
repeat id prop(`GHO', in, p?, idx1?)*prop(`GHO', out, p?, idx2?) = d_(colA[idx2], colA[idx1])/cOlNA;
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1])/cOlNR/4;
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = 1;
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = d_(colF[idx1], colF[idx2]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = i_;
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = i_ * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = - i_ * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(`PSI', virtual, p?, idx1?, idx2?) = -i_;
id prop(`PSI', in, p?, idx1?) = 1;
id prop(`PSI', out, p?, idx1?) = 1;


if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;


* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * cOlf(colA[idx3], colA[idx2], colA[idx1]);
id vx(`GHOBAR', `BACKGROUND', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * cOlf(colA[idx3], colA[idx2], colA[idx1]);
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_;
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_;
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ghhh * i_;

id vx(`H', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = - i_ * d_(colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );
id vx(`H', `GLU', `GLU', `GLU', p4?, p1?, p2?, p3?, idx4?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );

#do i=3,6
    id vx(<x1?{`PSI',}>,...,<x`i'?{`PSI',}>, p1?, ...,p`i'?, idx1?, ..., idx`i'?) = (-1*i_)^(`i'-2);
#enddo


id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]);

id vx(`BACKGROUND', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]);
id vx(`BACKGROUND', `GLU', `GHOBAR', `GHO', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) = -gs^2 * i_ *
    cOlf(colAdum[idx5], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[idx5]);


id vx(`GLU', `GLU', `GLU', `GLU', 1, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    cOlf(colAdum[idx5], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[idx5]);
id vx(`GLU', `GLU', `GLU', `GLU', 2, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    cOlf(colAdum[idx5], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[idx5]);
id vx(`GLU', `GLU', `GLU', `GLU', 3, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    cOlf(colAdum[idx5], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[idx5]);

id vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    cOlf(colAdum[idx6], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[idx6]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    cOlf(colAdum[idx6], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[idx6]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    cOlf(colAdum[idx6], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[idx6]);

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id tmp(x?) = x;

.sort:feynman-rules-global;

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, ?a, cOli1?) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color(1);
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
repeat id cOlf(cOlj1?,cOlj2?,cOlj3?)*color(x?) = color(x * cOlf(cOlj1,cOlj2,cOlj3));
id cOlNA^n?*color(x?) = color(x * cOlNA^n);
id cOlNR^n?*color(x?) = color(x * cOlNR^n);

B+ color;
.sort:color-prep;
Keep brackets;

* Only evaluate this part per unique color by bracketing
Argument color;
    #call color
    #call simpli

    id  cOlI2R = cOlcR*cOlNR/cOlNA;
    id  cOlNR/cOlNA*cOlcR = cOlI2R;
    id  cOld33(cOlpR1,cOlpR2) = [dabc^2];
    id  cOlNR/cOlNA = nf/cf/2;
    id  cOlcR = cf;
    id  cOlcA = ca;
    id  cOlI2R = nf/2;
	id	cOld44(cOlpA1,cOlpA2) = [d4AA];
	id	cOld44(cOlpR1,cOlpR2) = [d4RR];
	id	cOld44(cOlpR1,cOlpA1) = [d4RA];

* set the SU(3) values
*    id [dabc^2] = 15/18;
*    id [d4AA] = 135;
*    id [d4RR] = 5/12;
*    id [d4RA] = 15/2;
*    id cOlNR = 3;
*    id cOlNA = 8;
*    id cf = 4 / 3;
*    id ca = 3;
*    id nf = 1;
EndArgument;
.sort:color;
S tf;

id color(x?) = x;

id cOlNR^n? = ca^n;
id cOlNA^n? = cf^n;

#endprocedure

#procedure FeynmanRulesMomentum()
* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pdgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);

#if 0
    id prop(1,virtual,p?,13,12) = p.p;
    id prop(1,virtual,p?,17,16) = p.p;
    id prop(21,virtual,p?,7,6) = p.p;

    id vx(?a) = 1;
    id prop(?a) = 1;
#endif

* do the spin sum external particles
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
repeat id prop(`GHO', in, p?, idx1?)*prop(`GHO', out, p?, idx2?) = 1;
repeat id prop(`BACKGROUND', in, p?, idx1?)*prop(`BACKGROUND', out, p?, idx2?) = (p.p*d_(lorentz[idx1],lorentz[idx2])-p(lorentz[idx1])*p(lorentz[idx2]));
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = 1;
id prop(`PHO', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = 1;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(`PSI', virtual, p?, idx1?, idx2?) = 1;
id prop(`PSI', in, p?, idx1?) = 1;
id prop(`PSI', out, p?, idx1?) = 1;


if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = p3(lorentz[idx2]);
id vx(`GHOBAR', `BACKGROUND', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = (p3(lorentz[idx2])-p1(lorentz[idx2]));
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = 1;

id vx(`BACKGROUND', `GLU', `QBAR', `Q', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = g_(lorentz[idx1],lorentz[idx2]);

id vx(`H', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = p3(lorentz[idx2])*p2(lorentz[idx3]) - p2.p3 * d_(lorentz[idx2], lorentz[idx3]);
id vx(`H', `GLU', `GLU', `GLU', p4?, p1?, p2?, p3?, idx4?, idx1?, idx2?, idx3?) =
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2])
;

#do i=3,6
    id vx(<x1?{`PSI',}>,...,<x`i'?{`PSI',}>, p1?, ...,p`i'?, idx1?, ..., idx`i'?) = 1;
#enddo

id D = rat(4-2*ep, 1);
.sort:feynman-rules-vertices-1;

* construct gamma string, drop odd-length gamma traces and symmetrize the trace
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
id gamma(s1?,?a,s1?) = gammatrace(?a)*delta_(mod_(nargs_(?a), 2));
id gammatrace = 4;

.sort:gamma-filter;

* TODO: use momentum conservation to reduce the number of different terms
#do i=1,1
    id once ifnomatch->skip vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) =
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2]);

*    id once ifnomatch->skip vx(`BACKGROUND', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) =
*    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2])
*    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
*    + d_(lorentz[idx1], lorentz[idx3]) * p2(lorentz[idx2]) * [1/1-xi] 
*    
*    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
*    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
*    - d_(lorentz[idx1], lorentz[idx2]) * p3(lorentz[idx3]) * [1/1-xi] 
*    
*    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
*    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
*    ;

    redefine i "0";
    label skip;

    B+ vx;
    .sort:3g;
    Keep brackets;
#enddo

id vx(`GLU', `GLU', `GLU', `GLU', 1, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 2, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 3, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

id vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id D = rat(4-2*ep, 1);
.sort:feynman-rules-vertices-2;

Multiply counter(1);
repeat id gammatrace(?a,p?,?b)*counter(i?) = vec(p,lorentzdummy[i])*gammatrace(?a,lorentzdummy[i],?b)*counter(i+1);
id counter(n?) = 1;

id gammatrace(?a) = gammatracetensor(?a);
id vec(p?,mu?) = p(mu);
id D = rat(4-2*ep, 1);

#ifdef `FOURDIM'
    .sort:fourdim-poly;
    Polyratfun;
    id rat(x1?,x2?) = x1/x2;
#endif

B+ gammatracetensor;
.sort:gamma-to-tensor;
Keep brackets;

#ifdef `FOURDIM'
    #do i=1,10
        id once gammatracetensor(?a) = g_(`i',?a);
        trace4 `i';
        B+ gammatracetensor;
        .sort:trace-`i';
        Keep brackets;
    #enddo
    .sort:fourdim-trace;
    Polyratfun rat;
#else
* at this stage all indices should be inside the gammatracetensor only
* FIXME: not the case anymore!
    #call Gstring(gammatracetensor,1)
    id D = rat(4-2*ep, 1);
#endif

id pzero = 0; * Substitute the 0-momentum by 0
.sort:feynman-rules-final;
#endprocedure

#procedure FeynmanRules()
    #call FeynmanRulesGlobal()
    #call FeynmanRulesMomentum()
#endprocedure

#procedure ExtractMomenta()
* strip momentum tags
    repeat id f?{vx,prop}(?a,p?) = f(?a);

* extract momenta from the Feynman rules
    id prop(x1?,x2?,?a,p1?,x3?,?b) = prop(x1,x2,?a,p1,x3,?b,conf(?a,p1));
    id vx(?a,x1?,p1?,?b,p2?,x3?,?c) = vx(?a,x1,p1,?b,p2,x3,?c,conf(p1,?b,p2));
    argument prop,vx;
        splitarg conf;
        repeat id conf(?a,-p?vector_,?b) = conf(?a,p,?b);
        repeat id conf(?a,p?,?b,p?,?c) = conf(?a,p,?b,?c);
    endargument;
    id f?{prop,vx}(?a,conf(?b)) = f(?a,?b);
#endprocedure
*--#] feynman-rules :