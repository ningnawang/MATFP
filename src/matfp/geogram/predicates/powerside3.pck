#include "kernel.pckh"

Sign predicate(side3)(
    point(p0), scalar w0,
    point(p1), scalar w1,
    point(p2), scalar w2,
    point(p3), scalar w3,
    point(q0), point(q1), point(q2)  DIM
) 
group[degree=2] w0 w1 w2 w3;
{

    /* note: 1* for ensuring mcc's degree compatibility */
    scalar l1 = 1*sq_dist(p1,p0);
    scalar l2 = 1*sq_dist(p2,p0);
    scalar l3 = 1*sq_dist(p3,p0);

    scalar R1 = 1*w_diff(w0,w1);
    scalar R2 = 1*w_diff(w0,w2);
    scalar R3 = 1*w_diff(w0,w3);

    scalar L1 = l1 + R1;
    scalar L2 = l2 + R2;
    scalar L3 = l3 + R3;

    scalar a10 = 2*dot_at(p1,q0,p0);
    scalar a11 = 2*dot_at(p1,q1,p0);
    scalar a12 = 2*dot_at(p1,q2,p0);
    scalar a20 = 2*dot_at(p2,q0,p0);
    scalar a21 = 2*dot_at(p2,q1,p0);
    scalar a22 = 2*dot_at(p2,q2,p0);

    scalar a30 = 2*dot_at(p3,q0,p0);
    scalar a31 = 2*dot_at(p3,q1,p0);
    scalar a32 = 2*dot_at(p3,q2,p0);

    /*
     * [ b00 b01 b02 ]           [  1   1   1  ]-1
     * [ b10 b11 b12 ] = Delta * [ a10 a11 a12 ]
     * [ b20 b21 b22 ]           [ a20 a21 a22 ]
     */

    scalar b00 = a11*a22-a12*a21;
    scalar b01 = a21-a22;
    scalar b02 = a12-a11;
    scalar b10 = a12*a20-a10*a22;
    scalar b11 = a22-a20;
    scalar b12 = a10-a12;
    scalar b20 = a10*a21-a11*a20;
    scalar b21 = a20-a21;
    scalar b22 = a11-a10;

    scalar Delta = b00+b10+b20;

    /*
     *       [ Lambda0 ]   [ b01 b02 ]   [ L1 ]   [ b00 ]
     * Delta [ Lambda1 ] = [ b11 b12 ] * [    ] + [ b10 ]
     *       [ Lambda2 ]   [ b21 b22 ]   [ L2 ]   [ b20 ]
     */

    scalar DeltaLambda0 = b01*L1 + b02*L2 + b00;
    scalar DeltaLambda1 = b11*L1 + b12*L2 + b10;
    scalar DeltaLambda2 = b21*L1 + b22*L2 + b20;

    scalar r = Delta*L3 - (
           a30 * DeltaLambda0 + 
           a31 * DeltaLambda1 +
           a32 * DeltaLambda2
    );

    Sign Delta_sign = sign(Delta);
    Sign r_sign     = sign(r);

    generic_predicate_result(Delta_sign*r_sign);

    begin_sos4(p0,p1,p2,p3)
       sos(p0, Sign(Delta_sign*sign(
                 Delta-((b01+b02)*a30+(b11+b12)*a31+(b21+b22)*a32)
       )))
       sos(p1, Sign(Delta_sign*sign((a30*b01)+(a31*b11)+(a32*b21))))
       sos(p2, Sign(Delta_sign*sign((a30*b02)+(a31*b12)+(a32*b22))))
       sos(p3, NEGATIVE)
    end_sos
}