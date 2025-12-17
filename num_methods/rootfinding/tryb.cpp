#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef double db;

ll degree;
vector<db> cof;

db f(double x)
{

    double sum = 0;
    ll n = cof.size();
    for (ll i = 0; i < n; i++)
    {
        sum += (cof[i] * pow(x, (n - 1 - i)));
    }
    return sum;
}

int main()
{
    cout << "enter degree" << endl;
    cin >> degree;
    cof.resize(degree + 1);

    double xmax, root = 0, it, step = .5, a, b, c, e = .0001;
    for (auto &u : cof)
        cin >> u;
    xmax = sqrt((cof[1] / cof[0]) * (cof[1] / cof[0]) - 2 * (cof[2] / cof[0]));
    it = 0;
    for (double i = -xmax; i <= xmax - step; i += step)
    {
        a = i;
        b = i + step;

        if (f(a) * f(b) < 0)
        {

            c = a + b;
            c /= 2;

            while (fabs(a - b) > e)
            {
                c = a + b;
                c /= 2;
                if (f(a) * f(c) < 0)
                    b = c;
                else if (f(b) * f(c) < 0)
                    a = c;
                it++;
            }
            cout << "bracket : [" << i << " , " << i + step << "]" << endl;
            cout << "root : " << c << endl;
            cout << "it steps : " << it << endl;
        }
    }
}