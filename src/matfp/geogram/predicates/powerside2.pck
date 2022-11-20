#include "kernel.pckh"

Sign predicate(side2)(
    point(p0), scalar w0,
    point(p1), scalar w1,
    point(p2), scalar w2,
    point(q0), point(q1)  DIM
) 
group[degree=2] w0 w1 w2;
{
    /* note: 1* for ensuring mcc's degree compatibility */
    scalar l1 = 1*sq_dist(p1,p0);
    scalar l2 = 1*sq_dist(p2,p0);

    scalar R1 = 1*w_diff(w0,w1);
    scalar R2 = 1*w_diff(w0,w2);

    scalar L1 = l1 + R1;
    scalar L2 = l2 + R2;

    scalar a10 = 2*dot_at(p1,q0,p0);
    scalar a11 = 2*dot_at(p1,q1,p0);

    scalar a20 = 2*dot_at(p2,q0,p0);
    scalar a21 = 2*dot_at(p2,q1,p0);

    scalar Delta = a11 - a10;

    /*
     *       [ Lambda0 ]   [ -1 ]        [  a11 ]
     * Delta [         ] = [    ] * L1 + [      ]
     *       [ Lambda1 ]   [  1 ]        [ -a10 ]
     */
    scalar DeltaLambda0 = a11 - L1;
    scalar DeltaLambda1 =  L1 - a10;

    scalar r = Delta*L2 - a20*DeltaLambda0 - a21*DeltaLambda1;

    Sign Delta_sign = sign(Delta);
    Sign r_sign     = sign(r);

    generic_predicate_result(Delta_sign*r_sign);

    begin_sos3(p0,p1,p2)
       sos(p0, Sign(Delta_sign*sign(Delta-a21+a20)))
       sos(p1, Sign(Delta_sign*sign(a21-a20)))
       sos(p2, NEGATIVE)
    end_sos
}