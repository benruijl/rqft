#-
#:ContinuationLines 0
Off statistics;
On nospacesinnumbers;

*--#[ setup :

S D, ep(:`MAXPOLE');
V energyselector,p1,...,p40,ps1,...,ps40,k1,...,k40,c1,...,c40,cs1,...,cs40,fmb1,...,fmb40,fmbs1,...,fmbs40; * force this internal ordering in FORM
Auto V p,k,c;
Auto S lm,ext,E;
* note: a change of spinor dimension here needs to be reflected in Gstring and in the Feynman rules
Auto I mu=D,s=4;
Symbol ge, gs, ghhh, type, in, out, virtual, [1/1-xi];
Auto S x, idx, t, n;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, gammatrace(c), GGstring, NN, vector,g(s),delta(s),tmps(s),T, counter,color, prop, replace;


#include feynman_rules.h


CF f, vx, vxs(s), uvx, vec, vec1;
CF subs, configurations, conf, tder, cmb, cbtofmb, fmbtocb, diag, forestid, der, energy, spatial(s), onshell, uvcutoff;
CF subgraph, uvconf, uvconf1, uvconf2, uvprop, uv, uvtopo, irtopo, intuv, integrateduv, gluonbubble;
CT gammatracetensor(c),opengammastring;

S UVRenormFINITE, massct, ICT, mUVexp, logmu, logmUVexp, logmt, z3, mi1L1, alarmMi1L1;
Fill logmasses(6) = logmt;
Fill logmasses(-6) = logmt;

CF integratedct, rat, num, den, tmp;
Set ts: t0,...,t20;
CT penergy,energyct;
NF energync;
Symbol ca,cf,nf,[dabc^2],[d4RR],[d4RA],[d4AA];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli40;
Set colA: cOlj1,...,cOlj40;
Set colAdum: cOljj1,...,cOljj40;

Polyratfun rat;

*--#] setup :

#include tensorreduce.frm
#include integrateduv.frm

* Load the diagrams
#include- `DIAG'_in.h

.sort:load-input;

* Fill in the Feynman rules for the parts that do not depend on kinematics
#call FeynmanRulesGlobal()

* collect the edges and vertices belonging to the forests
#call ExtractMomenta()
id forestid(?a) = forestid(?a,1);
repeat id forestid(?a,k1?,?b,x?)*f?{prop,vx}(?c,k1?,?d) = forestid(?a,k1,?b,x*f(?c,k1,?d));
id forestid(?a) = forestid(forestid(?a));
argtoextrasymbol tonumber,forestid,1;
#redefine tensorforeststart "`extrasymbols_'"

.sort:tensorforest-splitoff-1;
#redefine tensorforestend "`extrasymbols_'"
#do ext={`tensorforeststart'+1},`tensorforestend'
    L tensorforest{`ext'-`tensorforeststart'} = extrasymbol_(`ext');
#enddo
#define tensorforestcount "{`tensorforestend'-`tensorforeststart'}"


id forestid(x1?,?a,x3?) = forestid(?a)*forest(x1)*x3;

.sort:tensorforest-splitoff-2;
Hide F;

id forestid(?a) = 1;

id cbtofmb(?a) = replace_(?a);
.sort:fmb-replace-2;

* get the momentum dependence in the fmb
#call ExtractMomenta()

* move the vertices and propagators into a subgraph, starting from the most nested ones first
#do i=0,5
    repeat id subgraph(`i',?e,x1?,fmb1?,?a,fmb2?,?b)*f?{prop,vx}(?c,fmb2?,?d) = subgraph(`i',?e,x1*f(?c,fmb2,?d),fmb1,?a,fmb2,?b);
    repeat id subgraph(`i',?e,x1?,fmb1?,?a)*f?{prop,vx}(?b,fmb1?,?c) = subgraph(`i',?e,x1*f(?b,fmb1,?c),fmb1,?a);
#enddo
repeat id subgraph(?a,p?) = subgraph(?a);

.sort:graphs;


* Create the local UV counterterm
#define uvdiagtotalstart "`extrasymbols_'"
#do i = 1,1
* Split off all subgraphs without dependencies to a separate expression
    id subgraph(x1?, x2?) = uvconf1(x1, x2);
    argtoextrasymbol tonumber,uvconf1,2;
    #redefine uvdiagstart "`extrasymbols_'"
    .sort:uvdiag-1;
    #redefine uvdiagend "`extrasymbols_'"
    #do ext={`uvdiagstart'+1},`uvdiagend'
        L uvdiag`ext' = extrasymbol_(`ext');
    #enddo
    .sort:uvdiag-2;

* fill in results from subdiagrams
    #do ext={`uvdiagtotalstart'+1},`uvdiagstart'
        id uvconf1(`ext') = uvdiag`ext';
    #enddo
    .sort:uv-subgrap-fillin;
    Hide tensorforest1,...,tensorforest`tensorforestcount';

* Apply Feynman rules to the UV subgraph
    id opengammastring(?a) = gamma(?a);
    #call FeynmanRulesMomentum()

* linearize the gamma matrices and convert them to tensors
* this makes them suitable for differentiation
    repeat id once gamma(s1?,?a,p?!vector_,?b,s2?) = p(mudummy)*gamma(s1,?a,mudummy,?b,s2);
    id gamma(?a) = opengammastring(?a);

    id uvprop(k?,t1?,p?,m?) = uvprop(k,t1,p,m)*uvconf2(p);

    splitarg uvconf2;
    chainout uvconf2;
    id uvconf2(-p?vector_) = uvconf2(p);
    repeat id uvconf2(p?)*uvconf2(p?) = uvconf2(p);

    id uvconf2(p?) = replace_(p, t * p);

    argument uvprop,1,vxs,uvtopo,irtopo,diag,onshell,intuv;
        Multiply replace_(t, 1);
    endargument;

* Taylor expand the propagators to the right depth
* p carries a t dependence that determines the order
* t1 determines the powers of the UV propagator
* rescale all masses in the numerator coming from the expansion if we are UV expanding, including mUVexp
    if (count(uvprop,1));
        repeat id m?allmasses = tmp(m);
        id tmp(m?) = t*m;
    endif;

    B tmax,t,uvprop;
    .sort:expand-prep;
    Keep brackets;

* expand the propagators without loop momentum dependence
    id uvprop(k?,t1?,0,m?) = uvprop(k,t1,1,m) * (1 - (mUVexp^2*t^2-m^2*t^2) * t1 + (mUVexp^2*t^2-m^2*t^2)^2 * t1^2 + ALARM * t^5);
    id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
    repeat;
        id once ifnomatch->skiptruncation0 uvprop(k?,t1?,p?,m?)*t^x1?*tmax^x2? = uvprop(k,t1,1,m) * t^x1*tmax^x2 * theta_(x2-x1) *
            (1 +
                (-2*p.k-(p.p+mUVexp^2*t^2-m^2*t^2)) * t1 +
                (+4*p.k^2+4*p.k*(p.p+mUVexp^2*t^2-m^2*t^2)+(p.p+mUVexp^2*t^2-m^2*t^2)^2) * t1^2 +
                (-8*p.k^3-12*p.k^2*(p.p+mUVexp^2*t^2-m^2*t^2)) * t1^3 +
                (16*p.k^4) * t1^4 +
                ALARM * t^5);
        id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
        label skiptruncation0;
    endrepeat;

    id t = 1;
    id tmax = 1;

    if (count(ALARM, 1));
        Print "UV Taylor expansion depth exceeded.";
        exit "";
    endif;

    .sort:uv-expansion;

* match the denominator structure to a diagram
    id uvtopo(x?,?a) = uvtopo(x,1,?a);
    id uvprop(?a,m?) = uvprop(?a);
    repeat id t1?ts^n2? = tmp(t1,n2);
    repeat id uvprop(k?,t1?,n1?)*uvprop(k?,t1?,n2?) = uvprop(k,t1,n1+n2); 
    repeat id uvprop(k?,t1?,n1?)*uvtopo(x?,x1?,?a)*tmp(t1?,n2?) = uvprop(k,n1 + n2)*uvtopo(x,x1*t1^n2,?a);
    id uvprop(k?,t1?ts,n?) = uvprop(k, n);
    id uvtopo(?a) = 1;
    if (count(uvtopo, 1, irtopo, 1, tmp, 1));
        Print "Unsubstituted UV topology: %t";
        exit "Critical error";
    endif;

* compute the integrated UV counterterm
    Multiply counter(1);
    repeat id k1?.k2?*counter(n?) = vec(k1,n)*vec(k2,n)*counter(n + 1);
    id opengammastring(?a) = gamma(?a);
    repeat id gamma(s1?,?a,k1?,?b,s2?)*counter(n?) = gamma(s1,?a,n,?b,s2)*vec(k1,n)*counter(n + 1);
    repeat id k1?(mu?)*counter(n?) = vec(k1,n)*vec1(mu,n)*counter(n + 1);

* convert every k^0 into k.p0select, where p0select is effectively (1,0,0,0)
    repeat id penergy(k1?)*counter(n?) = vec(k1,n)*vec(p0select,n)*counter(n + 1);
    id counter(x?) = 1;

    Multiply replace_(vec, vec1); * consider all vectors as external
    repeat id uvprop(k1?,n?)*vec1(k1?,n1?) = uvprop(k1,n)*vec(k1,n1);

    .sort:topo-match;

    #call TensorReduce()

* contract all metrics
    repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);
    id g(n1?,n1?) = rat(4-2*ep,1);
    repeat id vec1(p1?,n1?)*g(n1?,n2?) = vec1(p1,n2);


    .sort:tensor-projection-loop;

    if (count(uvprop,1));
        id k1?.k2? = g(k1,k2);
    endif;

    #call IntegrateUV()

    id D = rat(4-2*ep,1);
    id ep^n1? = rat(ep^n1,1);
    .sort:ibp-reduction;

    id g(k1?,k2?) = k1.k2; * k1,k2 should be external only
    id vec1(p?,n1?)*vec1(mu?,n1?) = p(mu);
    id vec1(k1?,n?)*vec1(k2?,n?) = k1.k2;
    id k1?.p0select = penergy(k1);

    repeat id gamma(s1?,?a,n?,?b,s2?) = gamma(s1,?a,lorentzdummy[n],?b,s2);
    id g(n1?,n2?) = d_(lorentzdummy[n1], lorentzdummy[n2]);
    id vec1(k1?,n?) = k1(lorentzdummy[n]);
    id vec1(mu?,n1?) = d_(mu, lorentzdummy[n1]);
    id gamma(?a) = opengammastring(?a);

* Simplify all open gamma strings so that the Ds are contracted out
    #call Gstring(opengammastring,0)


    if (match(opengammastring(s1?,?a,mu?lorentzdummy,?b,s2?)) || count(vec1,1));
        Print "Warning: dummy index left in gamma string after projecting: %t";
* this happens when there are two gamma strings in a subgraph
* TODO: is there a risk for dummy index collisions?
*        exit "Critical error";
    endif;

* Substitute the masters and expand in ep
    #call SubstituteMasters()

* construct MS-bar contributions (the poles)
    #ifndef `NOMSBARSUBTRACTION'
        argument rat;
            id ep^n? = ep^n*thetap_(-n);
        endargument;
    #endif

    Print +s;
    .sort:uv-subgraph-done;
    UnHide tensorforest1,...,tensorforest`tensorforestcount';
    Hide uvdiag{`uvdiagstart'+1},...,uvdiag`uvdiagend';

* in the next loop the UV subgraph will be added to the tensorforest or to a UV graph that embeds it
    if (count(uvconf1, 1)) redefine i "0";

* fill in the unsubstituted subgraph into the supergraph
    repeat id subgraph(x1?,?a,n?,?b,x2?)*uvconf1(n?,x3?) = subgraph(x1,?a,?b,x2*uvconf1(x3));
    id uvconf1(n?,x?) = uvconf1(x);
    .sort:uv3;
#enddo

Multiply replace_(mUVexp, mUV);
.sort
UnHide F;
.sort
Drop tensorforest1,...,tensorforest`tensorforestcount';

#do i=1,`tensorforestcount'
    id forestid(`i') = tensorforest`i'; * fill in the tensor
#enddo

id opengammastring(?a) = gamma(?a);

#call FeynmanRulesMomentum()

#write<`DIAG'_out.h> "%E",F

Print +s;
.end
