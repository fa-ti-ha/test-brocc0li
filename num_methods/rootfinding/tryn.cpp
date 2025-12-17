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
double df(double x)
{
    double val = 0, d = degree;
    for (int i = 0; i < degree; i++)
    {
        val += coeff[i] * d * pow(x, d - 1);
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
    double root = 0, step = 0.05;
    double it = 0;
    double e = 0.0001;

    double x0, x1, x2, f0, f1, f2, d0, d1, d2;
    for (double i = -xmax; i <= xmax; i += step)
    {
        x0 = i;
        x1 = i + step;

        f0 = f(x0);
        f1 = f(x1);

        if (f0 * f1 < 0)
        {
            root++;
            it++;
            do
            {
                d1 = df(x1);
                f1 = f(x1);
                x2 = x1 - f1 / d1;
                f2 = f(x2);
                if (fabs(f2 - f1) <= e && fabs(x2 - x1) <= e)
                {
                    break;
                }
                x1 = x2;

            } while (1);
            cout << endl;

            cout << "iteration is:" << it << endl;
            cout << "the root " << root << " is: " << x2 << endl;
            cout << "search interval [" << i << "," << i + step << "]" << endl;
        }
    }
}
