#include<bits/stdc++.h>
using namespace std;
void print(vector<vector<double>>a)
{
    for(int i=0; i<a.size(); i++)
    {
        for(int j=0; j<a.size(); j++)
        {
            cout<<a[i][j]<<" ";

        }
        cout<<endl;

    }
    cout<<endl;
    return ;
}
int main()
{
//cout<<"ok"<<endl;
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
   int t;
   cin>>t;
   for(int test=1;test<=t;test++){

   cout<<"Test case : "<<test<<endl;
    int n;
   //cout<<"enter the no of eqaution:";
    cin>>n;
    //cout<<"enter the augmented matrix:"<<endl;
    vector<vector<double>>a(n,vector<double>(n)),L(n,vector<double>(n,0)),U(n,vector<double>(n,0));
    vector<double>b(n,0),x(n,0),y(n,0);
    int r=0;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n+1; j++)
        {
            if(j==n)
            {
                cin>>b[r];
                r++;
            }
            else cin>>a[i][j];

        }
    }
    bool f=false;
    for(int i=0; i<n; i++)
    {
        for(int j=i; j<n; j++)
        {
            double sum=0;
            for(int k=0; k<i; k++)
            {
                sum+=L[i][k]*U[k][j];
            }
            U[i][j]=a[i][j]-sum;
        }
        for(int j=i; j<n; j++)
        {
            double sum=0;
            for(int k=0; k<i; k++)
            {
                sum+=L[j][k]*U[k][i];
            }
            if(U[i][i]!=0)
            {
                L[j][i]=(a[j][i]-sum)/U[i][i];

            }
            else
            {
                f=true;
            }

        }

    }
    cout<<"L matrix is : "<<endl;
    print(L);
    cout<<"U matrix is :"<<endl;
    print(U);
    for(int i=0; i<n; i++)
    {
        double sum=0;
        for(int j=0; j<i; j++)
        {
            sum+=L[i][j]*y[j];
        }
        y[i]=b[i]-sum;
    }
    cout<<"Y matrix is :"<<endl;
    for(int i=0; i<n; i++)
    {
        cout<<y[i]<<endl;
    }
    cout<<endl;
    if(f)
    {

        if(y[n-1]==0)
        {
            cout<<"INFINITE SOLUTION "<<endl;
        }
        else cout<<"NO solution"<<endl;

    }
    else
    {
        cout<<"Unique solution "<<endl;
        for(int i=n-1; i>=0; i--)
        {
            double sum=0;
            for(int j=i+1; j<n; j++)
            {
                sum+=U[i][j]*x[j];
            }
            x[i]=(y[i]-sum)/U[i][i];
        }
        cout<<"The solution is (x1,x2,x3....) :";
        for(int i=0; i<n; i++)cout<<x[i]<<" ";
        cout<<endl;
    }
    cout<<endl<<endl;
   }
    return 0;
}
