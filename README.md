# Numerical Methods

## Table of Contents

* [Solution of Linear Equations](#solution-of-linear-equations)

  * [Gauss Elimination Method](#gauss-elimination-method)
  * [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
  * [LU Decomposition Method](#lu-decomposition-method)
  * [Matrix Inversion](#matrix-inversion)

* [Solution of Non-Linear Equations](#solution-of-non-linear-equations)

  * [Bisection Method](#bisection-method)
  * [False Position Method](#false-position-method)
  * [Secant Method](#secant-method)
  * [Newton Raphson Method](#newton-raphson-method)

* [Solution of Interpolation](#solution-of-interpolation)

  * [Newton's Forward Interpolation Method](#newtons-forward-interpolation-method)
  * [Newton's Backward Interpolation Method](#newtons-backward-interpolation-method)
  * [Newton's Divided Difference Method](#divided-difference-method)

* [Solution of Ordinary Differential Equations (ODE)](#solution-of-ordinary-differential-equations-ode)

  * [Runge Kutta Method](#runge-kutta-method)

* [Solution of Numerical Differentiation](#solution-of-numerical-differentiation)

  * [Numerical Differentiation by Forward Interpolation Method](#numerical-differentiation-by-forward-interpolation-method)
  * [Numerical Differentiation by Backward Interpolation Method](#numerical-differentiation-by-backward-interpolation-method)

* [Solution of Numerical Integrations](#solution-of-numerical-integrations)

  * [Simpson's One-Third Rule](#simpsons-one-third-rule)
  * [Simpson's Three-Eighths Rule](#simpsons-three-eighths-rule)

* [Solution of Curve Fitting Model](#solution-of-curve-fitting-model)

  * [Least Square Regression Method for Linear Equations](#least-square-regression-method-for-linear-equations)
  * [Least Square Regression Method for Transcendental Equations](#least-square-regression-method-for-transcendental-equations)
  * [Least Square Regression Method for Polynomial Equations](#least-square-regression-method-for-polynomial-equations)

---

## Solution of Linear Equations

---

### Gauss Elimination Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Gauss Jordan Elimination Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### LU Decomposition Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Matrix Inversion
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

## Solution of Non-Linear Equations

---

### Bisection Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### False Position Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Secant Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Newton Raphson Method
---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

## Solution of Interpolation

---

### Newton's Forward Interpolation Method

---


#### Theory


Newton’s Forward Interpolation Method is a numerical technique used to estimate unknown values of a function from equally spaced tabulated data. 
It constructs an interpolation polynomial using finite forward differences of the function values. 
It is most accurate when the value to be interpolated lies near the beginning of the data.


#### Basic Idea

Given a set of equally spaced data points:

x₀, x₁, x₂, ..., xₙ
and corresponding function values:

f(x₀), f(x₁), f(x₂), ..., f(xₙ)

We calculate **forward differences** (Δy, Δ²y, Δ³y, ...) and use them to form the interpolation polynomial.


#### Forward Difference Table


The forward difference table organizes the differences systematically:

| x   | f(x)   | Δf(x)  | Δ²f(x) | Δ³f(x) | ... |
|-----|--------|--------|--------|--------|-----|
| x₀  | f₀     | Δf₀    | Δ²f₀   | Δ³f₀   | ... |
| x₁  | f₁     | Δf₁    | Δ²f₁   | ...    |     |
| x₂  | f₂     | Δf₂    | ...    |        |     |
| x₃  | f₃     | ...    |        |        |     |
| ... | ...    |        |        |        |     |

Where:

Δf₀ = f₁ - f₀  
Δ²f₀ = Δf₁ - Δf₀  
Δ³f₀ = Δ²f₁ - Δ²f₀  
and so on.


#### Formula

Let h = x₁ - x₀  (equal spacing)
Let u = (x - x₀) / h

Newton Forward Interpolation Polynomial:

P(x) = f₀ + uΔf₀ + u(u-1)/2! Δ²f₀ + u(u-1)(u-2)/3! Δ³f₀ + ... + u(u-1)...(u-n+1)/n! Δⁿf₀


#### Steps to Apply

1. Construct the forward difference table from the given data.
2. Compute Δf₀, Δ²f₀, ..., Δⁿf₀.
3. Calculate u = (x - x₀)/h for the required x.
4. Substitute into the interpolation polynomial to find P(x).


#### Conditions of Applicability

- Data points must be equally spaced.
- The value to be interpolated should lie near the beginning of the table.
- Function should be continuous over the interval.
  

#### Advantages

- Simple and systematic for equally spaced data.
- Forward difference table reduces repeated calculations.
- Good accuracy near the beginning of the data set.


#### Limitations

- Not suitable for unequally spaced data.
- Accuracy decreases for values far from x₀.
- Higher-order differences may introduce rounding errors.



#### Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double forward_interpolation(vector<double> &x, vector<double> &y, int n, double xn)
{
    vector<vector<double>> dif(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
        dif[i][0] = y[i];

    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            dif[i][j] = dif[i + 1][j - 1] - dif[i][j - 1];

    cout << "Forward difference table:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
            cout << setw(12) << dif[i][j];
        cout << endl;
    }
    cout << endl;

    double h = x[1] - x[0];
    double u = (xn - x[0]) / h;

    double res = y[0];
    double term = 1.0;
    double fact = 1.0;

    for (int i = 1; i < n; i++)
    {
        term *= (u - (i - 1));
        fact *= i;
        res += (term * dif[0][i]) / fact;
    }
    return res;
}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    cout << fixed << setprecision(3);

    int t;
    cin >> t;
    for (int test = 1; test <= t; test++)
    {
        cout << "Test case : " << test << endl;

        int n;
        cin >> n;
        vector<double> x(n), y(n);

        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            int x1, x2, yf;
            cin >> x1 >> x2 >> yf;
            x[i] = x2;
            sum += yf;
            y[i] = sum;
        }

        int d1;
        double xn;
        cin >> d1 >> xn;

        double res = forward_interpolation(x, y, n, xn);
        cout << "the value of Y at " << d1 << "-" << xn << " : " << res << endl;

        double ex1, ex2, ey;
        cin >> ex1 >> ex2 >> ey;
        sum += ey;
        x.push_back(ex2);
        y.push_back(sum);

        double res2 = forward_interpolation(x, y, n + 1, xn);
        cout << "error is:" << fabs(res2 - res) << endl
             << endl;
    }
    return 0;
}

```

#### Input

```
3
5
30 40 31
40 50 42
50 60 51
60 70 55
70 80 31
40 45
80 90 25

4
10 20 15
20 30 20
30 40 18
40 50 22
20 25
50 60 30

6
0 10 5
10 20 9
20 30 14
30 40 20
40 50 27
50 60 35
10 15
60 70 45
```

#### Output

```
Test case : 1
Forward difference table:
      31.000      42.000       9.000      -5.000     -23.000
      73.000      51.000       4.000     -28.000
     124.000      55.000     -24.000
     179.000      31.000
     210.000

the value of Y at 40-45.000 : 51.461
Forward difference table:
      31.000      42.000       9.000      -5.000     -23.000      69.000
      73.000      51.000       4.000     -28.000      46.000
     124.000      55.000     -24.000      18.000
     179.000      31.000      -6.000
     210.000      25.000
     235.000

error is:1.887

Test case : 2
Forward difference table:
      15.000      20.000      -2.000       6.000
      35.000      18.000       4.000
      53.000      22.000
      75.000

the value of Y at 20-25.000 : 25.625
Forward difference table:
      15.000      20.000      -2.000       6.000      -2.000
      35.000      18.000       4.000       4.000
      53.000      22.000       8.000
      75.000      30.000
     105.000

error is:0.078

Test case : 3
Forward difference table:
       5.000       9.000       5.000       1.000       0.000       0.000
      14.000      14.000       6.000       1.000       0.000
      28.000      20.000       7.000       1.000
      48.000      27.000       8.000
      75.000      35.000
     110.000

the value of Y at 10-15.000 : 8.938
Forward difference table:
       5.000       9.000       5.000       1.000       0.000       0.000       1.000
      14.000      14.000       6.000       1.000       0.000       1.000
      28.000      20.000       7.000       1.000       1.000
      48.000      27.000       8.000       2.000
      75.000      35.000      10.000
     110.000      45.000
     155.000

error is:0.021


```

---

### Newton's Backward Interpolation Method

---

#### Theory

Newton’s Backward Interpolation Method is a numerical technique used to estimate unknown values of a function from equally spaced tabulated data. 
It constructs an interpolation polynomial using finite backward differences of the function values. 
It is most accurate when the value to be interpolated lies near the end of the data.


#### Basic Idea

Given a set of equally spaced data points:

x₀, x₁, x₂, ..., xₙ
and corresponding function values:

f(x₀), f(x₁), f(x₂), ..., f(xₙ)

We calculate backward differences (∇y, ∇²y, ∇³y, ...) and use them to form the interpolation polynomial.


#### Backward Difference Table

The backward difference table organizes the differences systematically:

| x   | f(x)   | ∇f(x)  | ∇²f(x) | ∇³f(x) | ... |
|-----|--------|--------|--------|--------|-----|
| x₀  | f₀     |        |        |        |     |
| x₁  | f₁     | ∇f₁    |        |        |     |
| x₂  | f₂     | ∇f₂    | ∇²f₂   |        |     |
| x₃  | f₃     | ∇f₃    | ∇²f₃   | ∇³f₃   |     |
| ... | ...    | ...    | ...    | ...    | ... |
| xₙ  | fₙ      | ∇fₙ    | ∇²fₙ   | ∇³fₙ    | ... |

Where:

∇fₙ = fₙ - fₙ₋₁  
∇²fₙ = ∇fₙ - ∇fₙ₋₁  
∇³fₙ = ∇²fₙ - ∇²fₙ₋₁  
and so on.


#### Formula

Let h = x₁ - x₀  (equal spacing)
Let u = (x - xₙ) / h

Newton Backward Interpolation Polynomial:

P(x) = fₙ + u∇fₙ + u(u+1)/2! ∇²fₙ + u(u+1)(u+2)/3! ∇³fₙ + ... + u(u+1)...(u+n-1)/n! ∇ⁿfₙ


#### Steps to Apply

1. Construct the backward difference table from the given data.
2. Compute ∇fₙ, ∇²fₙ, ..., ∇ⁿfₙ.
3. Calculate u = (x - xₙ)/h for the required x.
4. Substitute into the interpolation polynomial to find P(x).


#### Conditions of Applicability

- Data points must be equally spaced.
- The value to be interpolated should lie near the end of the table.
- Function should be continuous over the interval.


#### Advantages

- Simple and systematic for equally spaced data.
- Backward difference table reduces repeated calculations.
- Good accuracy near the end of the data set.


#### Limitations

- Not suitable for unequally spaced data.
- Accuracy decreases for values far from xₙ.
- Higher-order differences may introduce rounding errors.


#### Code

```cpp
#include<bits/stdc++.h>
using namespace std;

double error(vector<double>&x,vector<double>&y,double val)
{
    int n=x.size();
    vector<vector<double>>dif(n,vector<double>(n));
    for(int i=0;i<n;i++) dif[i][0]=y[i];
    for(int j=1;j<n;j++)
        for(int i=0;i<n-j;i++)
            dif[i][j]=(dif[i+1][j-1]-dif[i][j-1])/(x[i+j]-x[i]);

    double del=1.0;
    for(int i=1;i<n;i++)
        del*= (val-x[i-1]);

    return dif[0][n-1]*del;
}

int main()
{
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);

    cout<<fixed<<setprecision(3);

    int t;
    cin>>t;
    for(int test=1;test<=t;test++)
    {
        cout<<"Test case : "<<test<<endl;

        int n;
        cin>>n;
        vector<double>x(n),y(n);
        for(int i=0;i<n;i++)
            cin>>x[i]>>y[i];

        vector<vector<double>>dif(n,vector<double>(n,0));
        for(int i=0;i<n;i++)
            dif[i][0]=y[i];

        for(int j=1;j<n;j++)
            for(int i=n-1;i>=j;i--)
                dif[i][j]=dif[i][j-1]-dif[i-1][j-1];

        cout<<"Backward difference table:"<<endl;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<=i;j++)
                cout<<setw(12)<<dif[i][j];
            cout<<endl;
        }
        cout<<endl;

        double xx;
        cin>>xx;

        double h=x[n-1]-x[n-2];
        double v=(xx-x[n-1])/h;

        double res=y[n-1];
        double term=1.0;
        double fact=1.0;

        for(int i=1;i<n;i++)
        {
            term*=(v+i-1);
            fact*=i;
            res+=(term*dif[n-1][i])/fact;
        }

        cout<<"answer is : "<<res<<endl;

        double nx,ny;
        cin>>nx>>ny;
        x.push_back(nx);
        y.push_back(ny);

        cout<<"error is : "<<fabs(error(x,y,xx))<<endl<<endl;
    }
    return 0;
}

```

#### Input

```
3
5
10 5
20 9
30 14
40 20
50 27
45
60 35

4
1 2
2 4
3 9
4 16
3
5 25

6
0 1
1 1
2 2
3 6
4 24
5 120
4
6 720

```

#### Output

```
Test case : 1
Backward difference table:
       5.000
       9.000       4.000
      14.000       5.000       1.000
      20.000       6.000       1.000       0.000
      27.000       7.000       1.000       0.000       0.000

answer is : 23.375
error is : 0.000

Test case : 2
Backward difference table:
       2.000
       4.000       2.000
       9.000       5.000       3.000
      16.000       7.000       2.000      -1.000

answer is : 9.000
error is : 0.000

Test case : 3
Backward difference table:
       1.000
       1.000       0.000
       2.000       1.000       1.000
       6.000       4.000       3.000       2.000
      24.000      18.000      14.000      11.000       9.000
     120.000      96.000      78.000      64.000      53.000      44.000

answer is : 24.000
error is : 0.000

```

---

### Newton's Divided Difference Method

---

#### Theory

Newton’s Divided Difference Interpolation Method is a numerical technique used to estimate unknown values of a function from a set of **unequally spaced data points**. 
It constructs an interpolation polynomial using **divided differences** of the function values. 
This method generalizes Newton's forward and backward methods and works for both equally and unequally spaced data.


#### Basic Idea

Given a set of data points:

x₀, x₁, x₂, ..., xₙ
and corresponding function values:

f(x₀), f(x₁), f(x₂), ..., f(xₙ)

We calculate **divided differences** (f[xᵢ, xⱼ], f[xᵢ, xⱼ, xₖ], ...) recursively and use them to form the interpolation polynomial. 
Divided differences generalize forward/backward differences to unequally spaced points.


#### Divided Difference Table

The divided difference table organizes the differences systematically:

| x   | f(x)   | 1st Divided Difference | 2nd Divided Difference | 3rd Divided Difference | ... |
|-----|--------|----------------------|----------------------|----------------------|-----|
| x₀  | f₀     | f[x₀,x₁]             | f[x₀,x₁,x₂]          | f[x₀,x₁,x₂,x₃]       | ... |
| x₁  | f₁     | f[x₁,x₂]             | f[x₁,x₂,x₃]          | ...                  |     |
| x₂  | f₂     | f[x₂,x₃]             | ...                   |                      |     |
| x₃  | f₃     | ...                   |                      |                      |     |
| ... | ...    |                      |                      |                      | ... |

Where:

1st Divided Difference: f[xᵢ, xᵢ₊₁] = (f(xᵢ₊₁) - f(xᵢ)) / (xᵢ₊₁ - xᵢ)  
2nd Divided Difference: f[xᵢ, xᵢ₊₁, xᵢ₊₂] = (f[xᵢ₊₁, xᵢ₊₂] - f[xᵢ, xᵢ₊₁]) / (xᵢ₊₂ - xᵢ)  
3rd Divided Difference: f[xᵢ, xᵢ₊₁, xᵢ₊₂, xᵢ₊₃] = (f[xᵢ₊₁, xᵢ₊₂, xᵢ₊₃] - f[xᵢ, xᵢ₊₁, xᵢ₊₂]) / (xᵢ₊₃ - xᵢ)  
and so on.


#### Formula

Newton Divided Difference Polynomial:

P(x) = f(x₀) 
       + (x - x₀)f[x₀,x₁] 
       + (x - x₀)(x - x₁)f[x₀,x₁,x₂] 
       + (x - x₀)(x - x₁)(x - x₂)f[x₀,x₁,x₂,x₃] 
       + ... 
       + (x - x₀)(x - x₁)...(x - xₙ₋₁)f[x₀,x₁,...,xₙ]


#### Steps to Apply

1. Arrange the given data points in a table.
2. Construct the divided difference table recursively.
3. Use the top row of divided differences to construct the interpolation polynomial.
4. Substitute the required value of x into the polynomial to find P(x).


#### Conditions of Applicability

- Data points can be equally or unequally spaced.
- Function should be continuous over the interval.


#### Advantages

- Works for unequally spaced data points.
- Systematic and can be extended to higher orders easily.
- Provides an explicit polynomial for interpolation.


#### Limitations

- Computationally more intensive for large datasets.
- Accuracy may decrease for very high-order polynomials due to rounding errors.



#### Code

```cpp
#include<bits/stdc++.h>
using namespace std;
double error(vector<double>&x,vector<double>&y,double val)
{
    int n=x.size();
    vector<vector<double>>dif(n,vector<double>(n));
    for(int i=0; i<n; i++) dif[i][0]=y[i];
    for(int j=1; j<n; j++)
    {
        for(int i=0; i<n-j; i++)
        {
            dif[i][j]=(dif[i+1][j-1]-dif[i][j-1])/(x[i+j]-x[i]);
        }
    }
    double del=1.0;
    for(int i=1; i<n; i++)
    {
        del=del*(val-x[i-1]);
    }
    double e=dif[0][n-1]*del;
    return e;

}
double ddi(vector<double>x,vector<double>y,double val)
{
    int n=x.size();
    vector<vector<double>>tb(n,vector<double>(n,0));
    for(int i=0; i<n; i++) tb[i][0]=y[i];
    for(int j=1; j<n; j++)
    {
        for(int i=0; i<n-j; i++)
        {
            tb[i][j]=(tb[i+1][j-1]-tb[i][j-1])/(x[i+j]-x[i]);
        }
    }
    cout<<"divided difference table is :"<<endl;
    for(int j=0; j<n; j++)
    {
        for(int i=0; i<n-j; i++)
        {
          if(i<n-j-1)  cout<<tb[i][j]<<setw(12);
          else cout<<tb[i][j];
        }
        cout<<endl;
    }
    cout<<endl;
    double res=y[0],del=1.0;

    for(int i=1; i<n; i++)
    {
        del=del*(val-x[i-1]);
        res+=del*tb[0][i];
    }
    return res;
}
int main()
{
//cout<<"ok"<<endl;
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
cout << fixed << setprecision(3);
    int t;
    cin>>t;
    for(int test=1; test<=t; test++)
    {

        cout<<"Test case : "<<test<<endl;
        int n;
        cin>>n;
        vector<double>x(n),y(n);
        for(int i=0; i<n; i++)
        {
            cin>>x[i]>>y[i];
        }
//cout<<"enter the x:";
        double val;
        cin>>val;
        double res1=ddi(x,y,val);
        cout<<"value of y at x= "<<val<<" is :"<<res1<<endl;
        cout<<"error is :"<<error(x,y,val)<<"%"<<endl<<endl;
    }
    return 0;
}

```

#### Input

```
3
5
1 1
2.5 2
4 4
6 5
7.5 7
3.3

4
0 0
1 2
3 10
4.5 20
2.7

6
0 1
0.5 2
1.7 5
3 6
4.2 10
6 20
3.5

```

#### Output

```
Test case : 1
divided difference table is :
1.000       2.000       4.000       5.000       7.000
0.667       1.333       0.500       1.333
0.222      -0.238       0.238
-0.092       0.095
0.029

value of y at x= 3.300 is :3.161
error is :0.100%

Test case : 2
divided difference table is :
0.000       2.000      10.000      20.000
2.000       4.000       6.667
0.667       0.762
0.021

value of y at x= 2.700 is :8.431
error is :-0.029%

Test case : 3
divided difference table is :
1.000       2.000       5.000       6.000      10.000      20.000
2.000       2.500       0.769       3.333       5.556
0.294      -0.692       1.026       0.741
-0.329       0.464      -0.066
0.189      -0.096
-0.048

value of y at x= 3.500 is :6.973
error is :0.315%


```

---

## Solution of Numerical Differentiation

---

### Numerical Differentiation by Forward Interpolation Method

---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Numerical Differentiation by Backward Interpolation Method

---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

## Solution of Ordinary Differential Equations (ODE)

---

### Runge Kutta Method

---


#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

## Solution of Numerical Integrations

---

### Simpson's One-Third Rule

---


#### Theory


Simpson’s 1/3 Rule is a numerical integration method used to approximate the definite integral of a function when an exact analytical solution is difficult or impossible to obtain. It provides higher accuracy than the Trapezoidal Rule by approximating the integrand using parabolic arcs instead of straight lines.

#### Basic Idea

In Simpson’s 1/3 Rule, 
the interval [a, b] is divided into an even number of equal sub-intervals of width:

h = (b − a) / n, where n is even
a is the lower limit and b is the upper limit.

The function values at these equally spaced points are used to construct quadratic polynomials over pairs of intervals. The area under each parabola is then calculated to approximate the total area under the curve.

#### Mathematical Formula

Let:

x₀ = a  
x₁ = a + h  
x₂ = a + 2h  
...  
xₙ = b  

and

yᵢ = f(xᵢ)

Then Simpson’s 1/3 Rule is given by:

∫ₐᵇ f(x) dx ≈ (h / 3) [ y₀ + yₙ + 4(y₁ + y₃ + ... + yₙ₋₁) + 2(y₂ + y₄ + ... + yₙ₋₂) ]

#### Conditions of Applicability

- The number of sub-intervals must be even  
- The data points must be equally spaced  
- The function should be smooth and continuous over the interval  

#### Advantages

- Higher accuracy than the Trapezoidal Rule  
- Simple and easy to apply  
- Requires fewer intervals for good accuracy  

#### Limitations

- Cannot be applied when the number of intervals is odd  
- Not suitable for unequally spaced data  
- Accuracy decreases for highly oscillatory functions




#### Code

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

#### Input

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

#### Output

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

### Simpson's Three-Eighths Rule

---

#### Theory 



Simpson’s 3/8 Rule is a numerical integration method used to approximate the definite integral of a function when an exact analytical solution is difficult or impossible to obtain. It is an extension of Simpson’s 1/3 Rule and uses cubic polynomials (third-degree) to approximate the integrand, providing higher accuracy for certain functions.


#### Basic Idea

In Simpson’s 3/8 Rule, the interval [a, b] is divided into a multiple of 3 equal sub-intervals of width:

h = (b − a) / n,  where n is a multiple of 3

The function values at these equally spaced points are used to construct cubic polynomials over sets of three intervals. The area under each cubic curve is calculated to approximate the total integral.


#### Mathematical Formula

Let:

x₀ = a
x₁ = a + h
x₂ = a + 2h
x₃ = a + 3h
...
xₙ = b

and

yᵢ = f(xᵢ)

Then Simpson’s 3/8 Rule is given by:

∫ₐᵇ f(x) dx ≈ (3h / 8) [ y₀ + yₙ + 3(y₁ + y₂ + y₄ + y₅ + ... + yₙ₋₁ + yₙ₋₂) + 2(y₃ + y₆ + ... + yₙ₋₃) ]


#### Conditions of Applicability

- The number of sub-intervals must be a multiple of 3
- The data points must be equally spaced
- The function should be smooth and continuous over the interval


#### Advantages

- More accurate than the Trapezoidal Rule and 1/3 Rule for functions requiring cubic approximation
- Can handle curves with higher-order behavior
- Simple to apply for equally spaced data


#### Limitations

- Cannot be applied if the number of intervals is not a multiple of 3
- Not suitable for unequally spaced data
- Accuracy decreases for highly oscillatory functions




#### Code

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

#### Input

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

#### Output

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

## Solution of Curve Fitting Model

---

### Least Square Regression Method for Linear Equations

--- 

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Least Square Regression Method for Transcendental Equations

---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```

---

### Least Square Regression Method for Polynomial Equations

---

#### Theory

*Add theory here*

#### Code

```cpp
// add your code here
```

#### Input

```
add input here
```

#### Output

```
add output here
```
