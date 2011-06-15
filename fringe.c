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
  double ks=inFringe*K1/(-0.75563);   /* scaling factor based on the field used for computing the fringe matrices, may need better treatment */
  double R11, R12, R21, R33, R34, R43;
  double r11, r12, r21, r33, r34, r43;
  double R22, R44, r22, r44;
  double x, px, y, py, delta;
  double *vec;
  double a, dx, dpx, dy, dpy, ds;

  if (inFringe<0){
    /* Entrance fringe */
    R11=0.999684, R12=-0.0000216824, R21=0.00948972, R33=1.00032, R34=0.0000216824, R43=-0.00948671;
    r11=0.999781, r12=0.0000121256, r21=-0.00987381, r33=1.00022, r34=-0.0000121256, r43=0.00985779; 
   } else {
    R11=0.999781, R12=0.0000121256, R21=-0.00987381, R33=1.00022, R34=-0.0000121256, R43=0.00985779; 
    r11=0.999684, r12=-0.0000216824, r21=0.00948972, r33=1.00032, r34=0.0000216824, r43=-0.00948671;
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


void dipoleFringe(double *vec, double h, long inFringe, long higherOrder)
{
  /* vec = {x, qx, y, qy, s, delta} */
  double x, px, y, py, delta;
  double a, dx, dpx, dy, dpy, ds;

  x = vec[0];
  px = vec[1];
  y = vec[2];
  py = vec[3];
  delta = vec[5];
  dx = dpx = dy = dpy = ds = 0;
  
  a = inFringe*h/(8*(1 + delta));

  if (higherOrder) {
    double xr[11], yr[11], ar[11];
    long i, j;
    xr[0] = yr[0] = ar[0] = 1;
    for (j=0, i=1; i<11; i++, j++) {
      xr[i] = xr[j]*x;
      yr[i] = yr[j]*y;
      ar[i] = ar[j]*a;
    }
    dx  = (a*(7*ar[7]*xr[9] + 7*ar[8]*xr[10] + ar[5]*xr[7]*(7 + 27*ar[2]*yr[2]) 
		  + ar[6]*xr[8]*(7 + 30*ar[2]*yr[2]) + ar[4]*xr[6]*(7 + 24*ar[2]*yr[2]
		  + 30*ar[4]*yr[4]) + ar[3]*xr[5]*(7 + 21*ar[2]*yr[2] + 99*ar[4]*yr[4])
		  + 3*a*x*yr[2]*(-7 + 35*ar[2]*yr[2] - 91*ar[4]*yr[4] + 246*ar[6]*yr[6])
		  + ar[2]*xr[4]*(7 + 21*ar[2]*yr[2] - 90*ar[4]*yr[4] + 1140*ar[6]*yr[6])
		  + xr[3]*(7*a + 189*ar[5]*yr[4] - 975*ar[7]*yr[6]) 
		  + 3*yr[2]*(7 - 7*ar[2]*yr[2] + 21*ar[4]*yr[4] - 39*ar[6]*yr[6] 
		  + 82*ar[8]*yr[8]) - 7*xr[2]*(-1 - 6*ar[2]*yr[2] + 21*ar[4]*yr[4] 
		  - 96*ar[6]*yr[6] + 315*ar[8]*yr[8])))/7;
    dpx = (a*(2*py*y*(7 - 7*a*x - 28*ar[2]*yr[2] + 70*x*ar[3]*yr[2] + 56*x*ar[5]*(xr[2]
		  - 6*yr[2])*yr[2] + 84*ar[4]*yr[2]*(-xr[2] + yr[2]) 
		  + 3*x*ar[7]*yr[2]*(xr[4] - 262*xr[2]*yr[2] + 409*yr[4])
		  - 4*ar[6]*(5*xr[4]*yr[2] - 162*xr[2]*yr[4] + 57*yr[6])
		  + 20*ar[8]*(33*xr[4]*yr[4] - 162*xr[2]*yr[6] + 29*yr[8]))
		  + px*(-3*ar[7]*xr[6]*yr[2] - 24*ar[6]*xr[5]*yr[2]*(-1 + 55*ar[2]*yr[2])
		  + 3*ar[5]*xr[4]*yr[2]*(-28 + 655*ar[2]*yr[2])
		  + 24*ar[4]*xr[3]*yr[2]*(7 - 90*ar[2]*yr[2] + 630*ar[4]*yr[4])
		  + 3*a*yr[2]*(-21 + 70*ar[2]*yr[2] - 196*ar[4]*yr[4] + 513*ar[6]*yr[6])
		  - 7*a*xr[2]*(-1 + 30*ar[2]*yr[2] - 240*ar[4]*yr[4] + 1227*ar[6]*yr[6])
		  - 2*x*(7 - 84*ar[2]*yr[2] + 420*ar[4]*yr[4] - 1596*ar[6]*yr[6]
		  + 5220*ar[8]*yr[8]))))/7;
    dy  = -(a*y*(ar[7]*xr[6]*yr[2] + 8*ar[6]*xr[5]*yr[2]*(-1 + 33*ar[2]*yr[2])
		  - 8*ar[4]*xr[3]*yr[2]*(7 - 54*ar[2]*yr[2] + 270*ar[4]*yr[4])
		  + xr[4]*(28*ar[5]*yr[2] - 393*ar[7]*yr[4]) + 3*a*yr[2]*(7 - 14*ar[2]*yr[2]
		  + 28*ar[4]*yr[4] - 57*ar[6]*yr[6]) + a*xr[2]*(-7 + 70*ar[2]*yr[2]
		  - 336*ar[4]*yr[4] + 1227*ar[6]*yr[6]) + 2*x*(7 - 28*ar[2]*yr[2] 
		  + 84*ar[4]*yr[4] - 228*ar[6]*yr[6] + 580*ar[8]*yr[8])))/7;
    dpy = (a*(-6*px*y*(7 - 7*a*x + 14*ar[2]*(xr[2] - yr[2]) + 70*x*ar[3]*yr[2] 
		  + 7*ar[4]*(xr[4] - 14*xr[2]*yr[2] + 9*yr[4]) + 7*ar[5]*(xr[5]
		  + 18*xr[3]*yr[2] - 39*x*yr[4]) + 4*ar[6]*(2*xr[6] - 15*xr[4]*yr[2]
		  + 168*xr[2]*yr[4] - 39*yr[6]) + ar[7]*(9*xr[7] + 66*xr[5]*yr[2]
		  - 975*xr[3]*yr[4] + 984*x*yr[6]) + 10*ar[8]*(xr[8] + 2*xr[6]*yr[2] 
		  + 114*xr[4]*yr[4] - 294*xr[2]*yr[6] + 41*yr[8])) + py*(63*ar[7]*xr[8]
		  + 70*ar[8]*xr[9] + 7*ar[5]*xr[6]*(7 + 27*ar[2]*yr[2]) 
		  + 8*ar[6]*xr[7]*(7 + 30*ar[2]*yr[2]) + 6*ar[4]*xr[5]*(7 + 24*ar[2]*yr[2]
		  + 30*ar[4]*yr[4]) + 5*ar[3]*xr[4]*(7 + 21*ar[2]*yr[2] + 99*ar[4]*yr[4])
		  + 3*a*xr[2]*(7 + 189*ar[4]*yr[4] - 975*ar[6]*yr[6])
		  + 3*a*yr[2]*(-7 + 35*ar[2]*yr[2] - 91*ar[4]*yr[4] + 246*ar[6]*yr[6])
		  + 4*ar[2]*xr[3]*(7 + 21*ar[2]*yr[2] - 90*ar[4]*yr[4] + 1140*ar[6]*yr[6])
		  - 14*x*(-1 - 6*ar[2]*yr[2] + 21*ar[4]*yr[4] - 96*ar[6]*yr[6]
		  + 315*ar[8]*yr[8]))))/7;
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

