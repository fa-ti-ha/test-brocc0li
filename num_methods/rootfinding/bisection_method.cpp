#include <bits/stdc++.h>
using namespace std;
int degree;
vector<double> coeff;
double f(double x)
{
  double val = 0;
  double d = degree;
  int n = coeff.size();
  for (int i = 0; i < n; i++)
  {
    //  cout<<degree<<endl;
    double t = coeff[i] * (pow(x, d));
    val += t;
    d--;
  }

  return val;
}
int main()
{

  cout << "please enter the degree of polynomial equation :";
  cin >> degree;
  coeff.resize(degree + 1);
  cout << "please enter the coefficient of polynomial :";
  for (int i = 0; i < degree + 1; i++)
    cin >> coeff[i];
  double xmax = 5;
  double a, b, c, root = 0;
  double it = 0;
  // cout << xmax << endl;
  for (double i = -xmax; i <= xmax; i += 0.5)
  {
    a = i;
    b = i + 0.5;
    double fa = f(a), fb = f(b);
    // cout << a << " " << b << " " << fa * fb << endl;
    if (fa * fb < 0)
    {

      double e = 0.0001;
      do
      {
        it++;
        // cout<<"ok"<<endl;
        c = (a + b) / 2.0;
        
        //c= a - f(a) * ((b - a) / (f(b) - f(a)));
        // cout<<"a="<<a<<" f(a)="<<f(a)<<" b="<<b<<" f(b)="<<f(b)<<" c="<<c<<" f(c)="<<f(c)<<endl;

        if (f(c) * f(a) < 0)
          b = c;
        else
          a = c;
      } while (fabs(f(c)) > e && fabs(b - a) > e);
      cout << "iteration is:" << it << endl
           << endl;
      cout << "the root " << root << " is: " << c << endl;
      cout << "search interval [" << i << "," << i + 0.5 << "]" << endl;
    }
  }
}
