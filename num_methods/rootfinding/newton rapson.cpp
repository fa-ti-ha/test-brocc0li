#include <bits/stdc++.h>
using namespace std;
int degree;
vector<double> coeff;
void print()
{
  int d = degree;
  bool m = true;
  for (int i = 0; i < degree + 1; i++)
  {
    if (m)
    {
      cout << coeff[i] << "X^" << d;
      m = false;
      d--;
      continue;
    }
    if (coeff[i] == 0)
    {
      d--;
      continue;
    }

    if (coeff[i] > 0)
      cout << "+";
    if (i != degree)
      cout << coeff[i] << "X^" << d;
    else
      cout << coeff[i];
    d--;
  }
  cout << "=0" << endl;
}
double f(double x)
{

  double val = 0;
  double deg = coeff.size() - 1;
  for (int i = 0; i < coeff.size(); i++)
  {
    val += coeff[i] * pow(x, deg);
    deg--;
  }
  return val;
}
double df(double x)
{
  double val = 0;
  double deg = coeff.size() - 1;
  for (int i = 0; i < coeff.size() - 1; i++)
  {
    val += coeff[i] * deg * pow(x, deg - 1);
    deg--;
  }
  return val;
}
int main()
{
  cout << "enter degree:";
  cin >> degree;
  cout << "enter the coefficient:";
  coeff.resize(degree + 1);
  for (int i = 0; i < degree + 1; i++)
  {
    cin >> coeff[i];
  }
  cout << "equation is: ";
  print();
  double xmax = sqrt((coeff[1] / coeff[0]) * (coeff[1] / coeff[0]) - 2 * (coeff[2] / coeff[0]));

  double num = 0, e = 0.00001;

  for (double i = -xmax; i <= xmax; i += .65)
  {
    double df1, x0 = i, x1 = i + 0.65, x2, fx2;
    double fx0 = f(x0), fx1 = f(x1);
    if (fx0 * fx1 < 0)
    {
      num++;

      int it = 0;
      do
      {
        df1 = df(x1);
        fx1 = f(x1);
        x2 = x1 - (fx1 / df1);
        fx2 = f(x2);
        it++;
        if (abs(x2 - x1) < e || abs(fx2 - fx1) < e)
        {
          break;
        }
        cout << "search interval is [" << x0 << "," << x1 << "]" << endl;
        cout << "the root " << num << " is: " << x2 << endl;
        cout << "iteration is:" << it << endl
             << endl;
        x1 = x2;

      } while (1);
    }
  }

  return 0;
}
