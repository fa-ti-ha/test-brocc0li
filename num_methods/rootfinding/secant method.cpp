#include <bits/stdc++.h>
using namespace std;
vector<double> coeff;
void print(vector<double> coeff)
{
    int power = coeff.size() - 1;
    bool m = false;
    for (int i = 0; i < coeff.size(); i++)
    {
        if (coeff[i] == 0)
        {
            power--;
            continue;
        }
        if (power == 0)
        {
            if (coeff[i] > 0)
                cout << "+";
            // else cout<<"-";
            cout << coeff[i];
            continue;
        }
        if (!m)
        {
            m = true;
            cout << coeff[i] << "X^" << power;
        }
        else
        {
            if (coeff[i] > 0)
                cout << "+";
            // else cout<<"-";
            cout << coeff[i] << "X^" << power;
        }
        power--;
    }
    cout << "=0" << endl;
}
double f(double x)
{
    double val = 0;
    double power = coeff.size() - 1;
    for (int i = 0; i < coeff.size(); i++)
    {
        val += coeff[i] * pow(x, power);
        power--;
    }
    return val;
}
int main()
{
    int degree;
    cout << "enter the degree: ";
    cin >> degree;
    cout << "enter the coefficient :";
    for (int i = 0; i <= degree; i++)
    {
        int x;
        cin >> x;
        coeff.push_back(x);
    }
    cout << "equation is: ";
    print(coeff);
    cout << endl;
    double xmax = 0, e = 0.001;
    for (int i = 0; i < coeff.size(); i++)
    {
        double temp = coeff[i] / coeff[0];
        xmax = max(xmax, temp);
    }
    xmax++;
    double c = -xmax;
    while (c <= xmax)
    {
        double x0 = c;
        double x1 = c + 0.45;
        double fx0 = f(x0), fx1 = f(x1);
        if (fx0 * fx1 < 0)
        {
            cout << "search interval is: [" << x0 << "," << x1 << "]" << endl;
            int it = 0;
            do
            {
                double x2 = x1 - f(x1) * ((x1 - x0) / (fx1 - fx0));
                it++;
                double fx2 = f(x2);
                if (abs(x1 - x2) < e && fx2 < e)
                {
                    cout << "root is =" << x2 << endl
                         << "iteration is = " << it << endl
                         << endl;
                    break;
                }
                x0 = x1;
                x1 = x2;
                fx0 = fx1;
                fx1 = fx2;

            } while (1);
        }
        c = c + 0.45;
    }

    return 0;
}
