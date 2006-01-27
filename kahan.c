/* file: Kahan.c
 * This is the implementation of Kahan's algorithm. 
 * It has been tested with GNU compiler without optimization flag. 
 * Other compiler or flags could get unexpected results.
 * Yusong Wang, 2005
 */

double Kahan (long length, double a[], double *error)
{
  
  double sum = 0.0, C = 0.0, Y, T;
  long i; 

  /* Kahan's summation formula */
  for (i=0; i<length; i++) {
     Y = a[i] - C;
     T = sum + Y;
     C = (T-sum) - Y;
     sum = T;
  }

  *error = -C;
  return sum;
}


double KahanPlus (double oldSum, double b, double *error)
{
  double sum, Y;

  Y = b + *error;
  sum = oldSum + Y;
  *error = Y - (sum-oldSum);
 
  return sum; 
}
