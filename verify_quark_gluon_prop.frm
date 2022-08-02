#-
S nf,cf,n,gs,e,ca,ep,x,[dabc^2],z3,z5,cl2,sqrt3;
Auto V p;
CF rat;

L QuarkNLO = 
    #do i=0,0
        +
        #include quark_nlo/quark_nlo_`i'_out.h
    #enddo
    ;

L GhostNLO = 
    #do i=0,0
        +
        #include ghost_nlo/ghost_nlo_`i'_out.h
    #enddo
    ;

L QuarkNNLO = 
    #do i=0,5
        +
        #include quark_nnlo/quark_nnlo_`i'_out.h
    #enddo
    ;

L GhostNNLO = 
    #do i=0,5
        +
        #include ghost_nnlo/ghost_nnlo_`i'_out.h
    #enddo
    ;

L QuarkNNNLO = 
    #do i=0,77
        +
        #include quark_nnnlo/quark_nnnlo_`i'_out.h
    #enddo
    ;

L GhostNNNLO = 
    #do i=0,77
        +
        #include ghost_nnnlo/ghost_nnnlo_`i'_out.h
    #enddo
    ;

id gs = 1;
Multiply 1/p1.p1;

.sort

L ReferenceQuarkNLO = ep^-1* (-cf);

L ReferenceQuarkNNLO = ep^-2* (
       + cf*ca
       + 1/2*cf^2
    )
    + ep^-1 * (
        - 17/4*cf*ca
        + 3/4*cf^2
        + 1/2*nf*cf
    );

L ReferenceQuarkNNNLO = 
   ep^-3* (
       - 7/4*cf*ca^2
       - cf^2*ca
       - 1/6*cf^3
       + 1/6*nf*cf*ca
   )+

   ep^-2* (
       + 104/9*cf*ca^2
       + 29/12*cf^2*ca
       - 3/4*cf^3
       - 28/9*nf*cf*ca
       - 1/6*nf*cf^2
       + 2/9*nf^2*cf
   )+

   ep^-1 * (
       - 10559/432*cf*ca^2
       + 5/2*cf*ca^2*z3
       + 143/12*cf^2*ca
       - 4*cf^2*ca*z3
       - 1/2*cf^3
       + 1301/216*nf*cf*ca
       - 1/2*nf*cf^2
       - 5/27*nf^2*cf
   );


L ReferenceGhostNLO = ep^-1* 1/2 * (ca);

L ReferenceGhostNNLO = ep^-2* (
       - ca^2
       + 1/4*nf*ca
    )
    + ep^-1 * (
        + 49/48*ca^2
       - 5/24*nf*ca
    );

L ReferenceGhostNNNLO = 
   ep^-3* (
        + 343/144*ca^3
       - 19/18*nf*ca^2
       + 1/9*nf^2*ca
   )+

   ep^-2* (
       - 1873/432*ca^3
       + 349/216*nf*ca^2
       + 1/2*nf*cf*ca
       - 5/54*nf^2*ca
   )+

   ep^-1 * (
       + 229/81*ca^3
       + 1/4*ca^3*z3
       - 5/1296*nf*ca^2
       - 3/2*nf*ca^2*z3
       - 15/8*nf*cf*ca
       + 2*nf*cf*ca*z3
       - 35/324*nf^2*ca
   );

L QuarkCheckNLO = -QuarkNLO - ReferenceQuarkNLO;
L GhostCheckNLO = GhostNLO - ReferenceGhostNLO;
L QuarkCheckNNLO = i_*QuarkNNLO - ReferenceQuarkNNLO;
L GhostCheckNNLO = -i_*GhostNNLO - ReferenceGhostNNLO;
L QuarkCheckNNNLO = QuarkNNNLO - ReferenceQuarkNNNLO;
L GhostCheckNNNLO = -GhostNNNLO - ReferenceGhostNNNLO;

.sort
id rat(x?) = x;

B ep;
Print +s;
.end