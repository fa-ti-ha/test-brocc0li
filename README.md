# Numerical Methods

## Table of Contents


- [Non Linear Methods](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)



- [Linear Methods](#linear-methods)
  
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
   

- [Interpolation Methods](#solution-of-interpolation)
  - [Newton's Forward Interpolation Method](#newtons-forward-interpolation-method)
    - [Theory](#newtons-forward-interpolation-theory)
    - [Code](#newtons-forward-interpolation-code)
    - [Input](#newtons-forward-interpolation-input)
    - [Output](#newtons-forward-interpolation-output)
  - [Newton's Backward Interpolation Method](#newtons-backward-interpolation-method)
    - [Theory](#newtons-backward-interpolation-theory)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)


 
- [Integration Methods](#solution-of-numerical-integrations)
  - [Simpson's One-Third Rule](#simpsons-one-third-rule)
    - [Theory](#simpsons-one-third-rule-theory)
    - [Code](#simpsons-one-third-rule-code)
    - [Input](#simpsons-one-third-rule-input)
    - [Output](#simpsons-one-third-rule-output)
  - [Simpson's Three-Eighths Rule](#simpsons-three-eighths-rule)
    - [Theory](#simpsons-three-eighth-rule-theory)
    - [Code](#simpsons-three-eighth-rule-code)
    - [Input](#simpsons-three-eighth-rule-input)
    - [Output](#simpsons-three-eighth-rule-output)
  

---

## Non Linear Methods

### Bisection Method

### Bisection Method Theory

The bisection method is a numerical approach used to locate a real root of a nonlinear equation f(x) = 0. It relies on the fact that a continuous function must cross the x-axis if it changes sign over an interval. Starting from two points that lie on opposite sides of the root, the method gradually narrows the interval until the root is sufficiently approximated. Because the root always remains within the chosen interval, this method is considered very stable and reliable.


At each step, the current interval is divided into two equal halves. One half is discarded based on the sign of the function, and the remaining half continues to contain the root. This repeated narrowing continues until the approximation meets the required accuracy.

#### Basic Formula


Midpoint Calculation:
```
c = (a + b) / 2

```



## Integration Method

### Simpson's One-Third Rule

### Simpson's One-Third Rule Theory


### Simpson's One-Third Rule Code
```cpp
#include <bits/stdc++.h>
using namespace std;

void print(vector<double> &coeff)
{
    int n = coeff.size();
    int pow = n - 1;
    bool m = true;
    for (int i = 0; i < n; i++)
    {
        if (coeff[i] == 0)
        {
            pow--;
            continue;
        }
        if (i == n - 1)
        {
            if (coeff[i] < 0)
                cout << coeff[i] << "=0";
            else
                cout << "+" << coeff[i] << "=0";
        }
        else
        {
            if (m)
            {
                cout << coeff[i] << "X^" << pow;
                m = false;
            }
            else
            {
                if (coeff[i] > 0)
                    cout << "+" << coeff[i] << "X^" << pow;
                else
                    cout << coeff[i] << "X^" << pow;
            }
            pow--;
        }
    }
    cout << endl;
}

double f(double x, vector<double> &coeff)
{
    double val = 0;
    int n = coeff.size();
    int p = n - 1;
    for (int i = 0; i < n; i++)
    {
        val += coeff[i] * pow(x, p);
        p--;
    }
    return val;
}

double simp1_3rd(double u, double l, int interval, vector<double> &coeff)
{
    if (interval % 2 != 0)
        interval++; // ensure even
    double h = (u - l) / interval;
    double ans = f(u, coeff) + f(l, coeff);
    for (int i = 1; i < interval; i++)
    {
        double x = l + i * h;
        double y = f(x, coeff);
        if (i % 2 == 0)
            ans += 2 * y;
        else
            ans += 4 * y;
    }
    ans = ans * (h / 3.0);
    return ans;
}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int test;
    cin >> test;

    for (int t = 1; t <= test; t++)
    {
        cout << "Testcase: " << t << endl;

        int n;
        cout << "Enter the degree: ";
        cin >> n;
        cout << "Enter equation coefficients:" << endl;
        vector<double> coeff(n + 1);
        for (int i = 0; i <= n; i++)
            cin >> coeff[i];

        double u, l;
        cout << "Enter upper limit: ";
        cin >> u;
        cout << "Enter lower limit: ";
        cin >> l;

        int interval;
        cout << "Enter the interval: ";
        cin >> interval;

        double p;
        cout << "Enter the value of p: ";
        cin >> p;

        cout << "Polynomial: ";
        print(coeff);

        double result = simp1_3rd(u, l, interval, coeff);
        cout << "Integral of f(x) from " << l << " to " << u << " is: " << result << endl
             << endl;
    }
}
```

### Simpson's One-Third Rule Input
```

5

2
1 -3 2
2
0
4
1

3
2 0 -1 1
5
1
10
2

1
4 -2
6
0
3
1

0
5
0
1
2
1

2
1 0 -1
3
-1
2
2

```

### Simpson's One-Third Rule Output
```
Testcase: 1
Enter the degree: Enter equation coefficients:
Enter upper limit: Enter lower limit: Enter the interval: Enter the value of p: Polynomial: 1X^2-3X^1+2=0
Integral of f(x) from 0 to 2 is: 0.666667

Testcase: 2
Enter the degree: Enter equation coefficients:
Enter upper limit: Enter lower limit: Enter the interval: Enter the value of p: Polynomial: 2X^3-1X^1+1=0
Integral of f(x) from 1 to 5 is: 304

Testcase: 3
Enter the degree: Enter equation coefficients:
Enter upper limit: Enter lower limit: Enter the interval: Enter the value of p: Polynomial: 4X^1-2=0
Integral of f(x) from 0 to 6 is: 60

Testcase: 4
Enter the degree: Enter equation coefficients:
Enter upper limit: Enter lower limit: Enter the interval: Enter the value of p: Polynomial: +5=0
Integral of f(x) from 1 to 0 is: -5

Testcase: 5
Enter the degree: Enter equation coefficients:
Enter upper limit: Enter lower limit: Enter the interval: Enter the value of p: Polynomial: 1X^2-1=0
Integral of f(x) from -1 to 3 is: 5.33333
```

---

### Simpson's Three-Eighth Rule 

### Simpson's Three-Eighths Rule Theory
 

### Simpson's Three-Eighth Rule Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x)
{
  return 1.0 / (1 + x * x);
}

void printFunction()
{
  cout << "f(x) = 1 / (1 + x^2)" << endl;
}

double simp_3_8th(double u, double l, int interval)
{

  if (interval % 3 != 0)
  {
    interval += (3 - interval % 3);
  }

  double h = (u - l) / interval;
  double ans = f(u) + f(l);

  for (int i = 1; i < interval; i++)
  {
    double x = l + i * h;
    double y = f(x);

    if (i % 3 == 0)
      ans += 2 * y;
    else
      ans += 3 * y;
  }

  ans = ans * (3 * h / 8.0);
  return ans;
}

int main()
{
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);

  int test;
  cin >> test;

  for (int t = 1; t <= test; t++)
  {
    cout << "Testcase: " << t << endl;

    double u, l;
    cout << "Enter upper limit: ";
    cin >> u;

    cout << "Enter lower limit: ";
    cin >> l;

    int interval;
    cout << "Enter the interval: ";
    cin >> interval;

    printFunction();

    double result = simp_3_8th(u, l, interval);
    cout << "Integral of f(x) from " << l << " to " << u
         << " is: " << result << endl
         << endl;
  }

  return 0;
}

```

### Simpson's Three-Eighth Rule Input
```

5

1
0
3

2
0
6

3
0
9

4
0
12

5
0
15

```

### Simpson's Three-Eighth Rule Output
```
Testcase: 1
Enter upper limit: Enter lower limit: Enter the interval: f(x) = 1 / (1 + x^2)
Integral of f(x) from 0 to 1 is: 0.784615

Testcase: 2
Enter upper limit: Enter lower limit: Enter the interval: f(x) = 1 / (1 + x^2)
Integral of f(x) from 0 to 2 is: 1.10638

Testcase: 3
Enter upper limit: Enter lower limit: Enter the interval: f(x) = 1 / (1 + x^2)
Integral of f(x) from 0 to 3 is: 1.2483

Testcase: 4
Enter upper limit: Enter lower limit: Enter the interval: f(x) = 1 / (1 + x^2)
Integral of f(x) from 0 to 4 is: 1.32508

Testcase: 5
Enter upper limit: Enter lower limit: Enter the interval: f(x) = 1 / (1 + x^2)
Integral of f(x) from 0 to 5 is: 1.37267


```

---
