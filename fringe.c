/*
  Two functions quadFringe and dipoleFringe are provided for tracking fringe field effects.
  The argument vec is the phase space vector to be tracked.
  The argument inFringe is a flag: -1 for going into the magnet and +1 for going out.
  There are magnet/fringe parameters need to be setup (indicated below).

  Physics approximations:

  assuming the fringes are symmetric
  quadFringe: soft linear fringe maps, leading order nonlinear generator
  dipoleFringe: a leading order nonlinear generator for hard edge

  If highOrder sets to 1, higher order terms are included for sufficient symplecticity
  (not necessarily more physics), but speed will sacrify.
 
* Original code by C.-X. Wang.
* Modified for use in elegant by M. Borland.
*/

#include "mdb.h"
#include "track.h"

void quadFringe(double **coord, long np, double K1, int inFringe, int higherOrder)
{
  /* vec = {x, qx, y, qy, s, delta}
     need paramenters K1, R11, R12, R21, R33, R34, R43, r11, r12, r21, r33, r34, r43
     */
  long ip;
  double ks=inFringe*K1/0.75563;   /* scaling factor based on the field used for computing the fringe matrices, may need better treatment */
  double R11, R12, R21, R33, R34, R43;
  double r11, r12, r21, r33, r34, r43;
  double R22, R44, r22, r44;
  double x, px, y, py, delta;
  double *vec;
  double a, dx, dpx, dy, dpy, ds;

  if (inFringe<0){
    /* Entrance fringe */
    r11=1.00022; r12=-0.0000121256; r21=0.00985779; r33=0.999781; r34=0.0000121256; r43=-0.00987381;
    R11=1.00032; R12=0.0000216824; R21=-0.00948671; R33=0.999684; R34=-0.0000216824; R43=0.00948972;
  } else {
    R11=1.00022, R12=-0.0000121256, R21=0.00985779, R33=0.999781, R34=0.0000121256, R43=-0.00987381;
    r11=1.00032, r12=0.0000216824, r21=-0.00948671, r33=0.999684, r34=-0.0000216824, r43=0.00948972;
  }


  R11 = 1 + ks*(R11-1);       /* assuming fringe matrix is very close to I */
  R12 = ks*R12;
  R21 = ks*R21;
  R33 = 1 + ks*(R33-1);
  R34 = ks*R34;
  R43 = ks*R43;

  r11 = 1 + ks*(r11-1);
  r12 = ks*r12;
  r21 = ks*r21;
  r33 = 1 + ks*(r33-1);
  r34 = ks*r34;
  r43 = ks*r43;

  R22 = (1 + R12*R21)/R11;
  R44 = (1 + R34*R43)/R33;
  r22 = (1 + r12*r21)/r11;
  r44 = (1 + r34*r43)/r33;     /* up till here can be moved out for faster tracking */


  for (ip=0; ip<np; ip++) {
    vec = coord[ip];
    delta = vec[5];
    x  = R11*vec[0] + R12*vec[1];
    px = R21*vec[0] + R22*vec[1];
    y  = R33*vec[2] + R34*vec[3];
    py = R43*vec[2] + R44*vec[3];
    
    a = inFringe*K1/(12*(1 + delta));

    if (higherOrder) {
	dx  = (a*x*(8*(ipow(x,2) + 3*ipow(y,2)) + 4*ipow(a,2)*(5*ipow(x,6) + 21*ipow(x,4)*ipow(y,2) 
		  - 25*ipow(x,2)*ipow(y,4) - ipow(y,6)) + ipow(a,3)*(35*ipow(x,8) + 84*ipow(x,6)*ipow(y,2) 
		  + 498*ipow(x,4)*ipow(y,4) - 108*ipow(x,2)*ipow(y,6) + 3*ipow(y,8)) + 12*a*ipow(ipow(x,2)
		  - ipow(y,2),2)))/8;
	dpx = (a*(12*a*(-4*py*x*y + px*(ipow(x,2) - 5*ipow(y,2)))*(ipow(x,2) - ipow(y,2)) 
		  - 24*(-2*py*x*y + px*(ipow(x,2) + ipow(y,2))) + 4*ipow(a,2)*(-2*py*x*y*(3*ipow(x,4) 
		  + 50*ipow(x,2)*ipow(y,2) - 21*ipow(y,4)) + px*(ipow(x,6) + 75*ipow(x,4)*ipow(y,2) 
		  - 105*ipow(x,2)*ipow(y,4) - 35*ipow(y,6))) + 3*ipow(a,3)*(-8*py*x*y*(ipow(x,6) 
		  - 27*ipow(x,4)*ipow(y,2) + 83*ipow(x,2)*ipow(y,4) + 7*ipow(y,6)) + px*(ipow(x,8) 
		  - 108*ipow(x,6)*ipow(y,2) + 830*ipow(x,4)*ipow(y,4) + 196*ipow(x,2)*ipow(y,6) 
		  + 105*ipow(y,8)))))/8;
	dy  = (a*y*(-8*(3*ipow(x,2) + ipow(y,2)) + 4*ipow(a,2)*(ipow(x,6) + 25*ipow(x,4)*ipow(y,2) 
		  - 21*ipow(x,2)*ipow(y,4) - 5*ipow(y,6)) + ipow(a,3)*(3*ipow(x,8) - 108*ipow(x,6)*ipow(y,2) 
		  + 498*ipow(x,4)*ipow(y,4) + 84*ipow(x,2)*ipow(y,6) + 35*ipow(y,8)) + 12*a*ipow(ipow(x,2) 
		  - ipow(y,2),2)))/8;
	dpy = (a*(-8*px*x*y*(6 - 6*a*(ipow(x,2) - ipow(y,2)) + ipow(a,2)*(21*ipow(x,4) 
		  - 50*ipow(x,2)*ipow(y,2) - 3*ipow(y,4)) + 3*ipow(a,3)*(7*ipow(x,6) + 83*ipow(x,4)*ipow(y,2)
		  - 27*ipow(x,2)*ipow(y,4) + ipow(y,6))) + py*(315*ipow(a,3)*ipow(x,8) 
		  + 28*ipow(a,2)*ipow(x,6)*(5 + 21*a*ipow(y,2)) + 30*a*ipow(x,4)*(2 + 14*a*ipow(y,2) 
		  + 83*ipow(a,2)*ipow(y,4)) + ipow(y,2)*(24 + 12*a*ipow(y,2) - 4*ipow(a,2)*ipow(y,4)
		  + 3*ipow(a,3)*ipow(y,6)) - 12*ipow(x,2)*(-2 + 6*a*ipow(y,2) + 25*ipow(a,2)*ipow(y,4) 
		  + 27*ipow(a,3)*ipow(y,6)))))/8;

      } else {
	dx  = a*(ipow(x,3) + 3*x*ipow(y,2))/3;
	dpx = a*(-px*ipow(x,2) + 2*py*x*y - px*ipow(y,2));
	dy  = a*(-ipow(x,2)*y - ipow(y,3)/3);
	dpy = a*(-py*ipow(x,2) - 2*px*x*y - py*ipow(y,2));
      }
	
    ds  = (a/(1+delta))*(3*py*y*ipow(x,2) - px*ipow(x,3) - 3*px*x*ipow(y,2) + py*ipow(y,3));

    x  += dx;
    px += dpx;
    y  += dy;
    py += dpy;

    vec[0] = r11*x + r12*px;
    vec[1] = r21*x + r22*px;
    vec[2] = r33*y + r34*py;
    vec[3] = r43*y + r44*py;
    vec[4] += ds;
  }
}


void dipoleFringe(double **coord, long np, double h, int inFringe, int higherOrder)
{
  /* vec = {x, qx, y, qy, s, delta} */
  double *vec;
  double x, px, y, py, delta;
  double a, dx, dpx, dy, dpy, ds;
  long ip;
  
  for (ip=0; ip<np; ip++) {
    vec = coord[ip];
    x = vec[0];
    px = vec[1];
    y = vec[2];
    py = vec[3];
    delta = vec[5];
    
    a = inFringe*h/(8*(1 + delta));

    if (higherOrder) {
	dx  = (a*(7*ipow(a,7)*ipow(x,9) + 7*ipow(a,8)*ipow(x,10) + ipow(a,5)*ipow(x,7)*(7 + 27*ipow(a,2)*ipow(y,2)) 
		  + ipow(a,6)*ipow(x,8)*(7 + 30*ipow(a,2)*ipow(y,2)) + ipow(a,4)*ipow(x,6)*(7 + 24*ipow(a,2)*ipow(y,2)
		  + 30*ipow(a,4)*ipow(y,4)) + ipow(a,3)*ipow(x,5)*(7 + 21*ipow(a,2)*ipow(y,2) + 99*ipow(a,4)*ipow(y,4))
		  + 3*a*x*ipow(y,2)*(-7 + 35*ipow(a,2)*ipow(y,2) - 91*ipow(a,4)*ipow(y,4) + 246*ipow(a,6)*ipow(y,6))
		  + ipow(a,2)*ipow(x,4)*(7 + 21*ipow(a,2)*ipow(y,2) - 90*ipow(a,4)*ipow(y,4) + 1140*ipow(a,6)*ipow(y,6))
		  + ipow(x,3)*(7*a + 189*ipow(a,5)*ipow(y,4) - 975*ipow(a,7)*ipow(y,6)) 
		  + 3*ipow(y,2)*(7 - 7*ipow(a,2)*ipow(y,2) + 21*ipow(a,4)*ipow(y,4) - 39*ipow(a,6)*ipow(y,6) 
		  + 82*ipow(a,8)*ipow(y,8)) - 7*ipow(x,2)*(-1 - 6*ipow(a,2)*ipow(y,2) + 21*ipow(a,4)*ipow(y,4) 
		  - 96*ipow(a,6)*ipow(y,6) + 315*ipow(a,8)*ipow(y,8))))/7;
	dpx = (a*(2*py*y*(7 - 7*a*x - 28*ipow(a,2)*ipow(y,2) + 70*x*ipow(a,3)*ipow(y,2) + 56*x*ipow(a,5)*(ipow(x,2)
		  - 6*ipow(y,2))*ipow(y,2) + 84*ipow(a,4)*ipow(y,2)*(-ipow(x,2) + ipow(y,2)) 
		  + 3*x*ipow(a,7)*ipow(y,2)*(ipow(x,4) - 262*ipow(x,2)*ipow(y,2) + 409*ipow(y,4))
		  - 4*ipow(a,6)*(5*ipow(x,4)*ipow(y,2) - 162*ipow(x,2)*ipow(y,4) + 57*ipow(y,6))
		  + 20*ipow(a,8)*(33*ipow(x,4)*ipow(y,4) - 162*ipow(x,2)*ipow(y,6) + 29*ipow(y,8)))
		  + px*(-3*ipow(a,7)*ipow(x,6)*ipow(y,2) - 24*ipow(a,6)*ipow(x,5)*ipow(y,2)*(-1 + 55*ipow(a,2)*ipow(y,2))
		  + 3*ipow(a,5)*ipow(x,4)*ipow(y,2)*(-28 + 655*ipow(a,2)*ipow(y,2))
		  + 24*ipow(a,4)*ipow(x,3)*ipow(y,2)*(7 - 90*ipow(a,2)*ipow(y,2) + 630*ipow(a,4)*ipow(y,4))
		  + 3*a*ipow(y,2)*(-21 + 70*ipow(a,2)*ipow(y,2) - 196*ipow(a,4)*ipow(y,4) + 513*ipow(a,6)*ipow(y,6))
		  - 7*a*ipow(x,2)*(-1 + 30*ipow(a,2)*ipow(y,2) - 240*ipow(a,4)*ipow(y,4) + 1227*ipow(a,6)*ipow(y,6))
		  - 2*x*(7 - 84*ipow(a,2)*ipow(y,2) + 420*ipow(a,4)*ipow(y,4) - 1596*ipow(a,6)*ipow(y,6)
		  + 5220*ipow(a,8)*ipow(y,8)))))/7;
	dy  = -(a*y*(ipow(a,7)*ipow(x,6)*ipow(y,2) + 8*ipow(a,6)*ipow(x,5)*ipow(y,2)*(-1 + 33*ipow(a,2)*ipow(y,2))
		  - 8*ipow(a,4)*ipow(x,3)*ipow(y,2)*(7 - 54*ipow(a,2)*ipow(y,2) + 270*ipow(a,4)*ipow(y,4))
		  + ipow(x,4)*(28*ipow(a,5)*ipow(y,2) - 393*ipow(a,7)*ipow(y,4)) + 3*a*ipow(y,2)*(7 - 14*ipow(a,2)*ipow(y,2)
		  + 28*ipow(a,4)*ipow(y,4) - 57*ipow(a,6)*ipow(y,6)) + a*ipow(x,2)*(-7 + 70*ipow(a,2)*ipow(y,2)
		  - 336*ipow(a,4)*ipow(y,4) + 1227*ipow(a,6)*ipow(y,6)) + 2*x*(7 - 28*ipow(a,2)*ipow(y,2) 
		  + 84*ipow(a,4)*ipow(y,4) - 228*ipow(a,6)*ipow(y,6) + 580*ipow(a,8)*ipow(y,8))))/7;
	dpy = (a*(-6*px*y*(7 - 7*a*x + 14*ipow(a,2)*(ipow(x,2) - ipow(y,2)) + 70*x*ipow(a,3)*ipow(y,2) 
		  + 7*ipow(a,4)*(ipow(x,4) - 14*ipow(x,2)*ipow(y,2) + 9*ipow(y,4)) + 7*ipow(a,5)*(ipow(x,5)
		  + 18*ipow(x,3)*ipow(y,2) - 39*x*ipow(y,4)) + 4*ipow(a,6)*(2*ipow(x,6) - 15*ipow(x,4)*ipow(y,2)
		  + 168*ipow(x,2)*ipow(y,4) - 39*ipow(y,6)) + ipow(a,7)*(9*ipow(x,7) + 66*ipow(x,5)*ipow(y,2)
		  - 975*ipow(x,3)*ipow(y,4) + 984*x*ipow(y,6)) + 10*ipow(a,8)*(ipow(x,8) + 2*ipow(x,6)*ipow(y,2) 
		  + 114*ipow(x,4)*ipow(y,4) - 294*ipow(x,2)*ipow(y,6) + 41*ipow(y,8))) + py*(63*ipow(a,7)*ipow(x,8)
		  + 70*ipow(a,8)*ipow(x,9) + 7*ipow(a,5)*ipow(x,6)*(7 + 27*ipow(a,2)*ipow(y,2)) 
		  + 8*ipow(a,6)*ipow(x,7)*(7 + 30*ipow(a,2)*ipow(y,2)) + 6*ipow(a,4)*ipow(x,5)*(7 + 24*ipow(a,2)*ipow(y,2)
		  + 30*ipow(a,4)*ipow(y,4)) + 5*ipow(a,3)*ipow(x,4)*(7 + 21*ipow(a,2)*ipow(y,2) + 99*ipow(a,4)*ipow(y,4))
		  + 3*a*ipow(x,2)*(7 + 189*ipow(a,4)*ipow(y,4) - 975*ipow(a,6)*ipow(y,6))
		  + 3*a*ipow(y,2)*(-7 + 35*ipow(a,2)*ipow(y,2) - 91*ipow(a,4)*ipow(y,4) + 246*ipow(a,6)*ipow(y,6))
		  + 4*ipow(a,2)*ipow(x,3)*(7 + 21*ipow(a,2)*ipow(y,2) - 90*ipow(a,4)*ipow(y,4) + 1140*ipow(a,6)*ipow(y,6))
		  - 14*x*(-1 - 6*ipow(a,2)*ipow(y,2) + 21*ipow(a,4)*ipow(y,4) - 96*ipow(a,6)*ipow(y,6)
		  + 315*ipow(a,8)*ipow(y,8)))))/7;

      } else {

	dx  = a*(ipow(x,2) + 3*ipow(y,2));
	dpx = 2*a*(-px*x + py*y);
	dy  = -2*a*x*y;
	dpy = 2*a*(py*x - 3*px*y);

      }
	
    ds  = (a/(1+delta))*(2*py*x*y - px*ipow(x,2) - 3*px*ipow(y,2));

    vec[0] += dx;
    vec[1] += dpx;
    vec[2] += dy;
    vec[3] += dpy;
    vec[4] += ds;
  }
}

