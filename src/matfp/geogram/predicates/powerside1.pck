#include "kernel.pckh"

Sign predicate_power(side1)(
   point(p0), scalar w0,
   point(p1), scalar w1,
   point(q0)  DIM
)
group[degree=2] w0 w1;
{
  /* note: 1* for ensuring mcc's degree compatibility */
  scalar r = 1*sq_dist(p0,p1);
  r += w_diff(w0,w1);
  r -= 2*dot_at(p1,q0,p0);

  generic_predicate_result(sign(r));
  begin_sos2(p0,p1)
     sos(p0,POSITIVE) 
     sos(p1,NEGATIVE) 
  end_sos
}


Sign predicate_power_new(side1)(
   scalar p0x, scalar p0y, scalar p0z, double w0,
   scalar p1x, scalar p1y, scalar p1z, double w1,
   scalar q0x, scalar q0y, scalar q0z  DIM
)
group p0x p1x q0x;
group p0y p1y q0y;
group p0z p1z q0z;
group[degree=2] w0 w1;
{
   /* note: 1* for ensuring mcc's degree compatibility */

   scalar dpx = p0x - p1x;
   scalar dpy = p0y - p1y;
   scalar dpz = p0z - p1z;

   scalar r = square(dpx) + square(dpy) + square(dpz) + (w0-w1);


   scalar dqx = p0x - q0x;
   scalar dqy = p0y - q0y;
   scalar dqz = p0z - q0z;
   r-= 2* ( square(dqx) + square(dqy) + square(dqz) );

   generic_predicate_result(sign(r));
}