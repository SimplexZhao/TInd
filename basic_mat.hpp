//A class contains Mat Rotate, Multi & Inv.
#include<iomanip>
#include<iostream>
using namespace std;
class MatOperation
{
	public:    
	MatOperation(){};
    virtual ~MatOperation(){};
	
	// Generate R matrix, left=R*right
	// left:image right:world
	// flag==2, φ-ω-κ | flag==1, x-y-z| flag==0, z-y-x
	static void GetRotationMat(double* mat, double a1, double a2, double a3, int flag=1)
	{
		if(flag==0)
		{
			double r=a1;
			double p=a2;
			double y=a3;
			mat[0] = cos(r)*cos(y);
			mat[1] = -cos(r)*sin(y);
			mat[2] = -sin(r);
			mat[3] = cos(p)*sin(y)+sin(r)*sin(p)*cos(y);
			mat[4] = cos(p)*cos(y)-sin(r)*sin(p)*sin(y);
			mat[5] = cos(r)*sin(p);
			mat[6] = sin(r)*cos(y)*cos(p)-sin(p)*sin(y);
			mat[7] = -cos(p)*sin(r)*sin(y)-sin(p)*cos(y);
			mat[8] = cos(r)*cos(p);
		}
		if(flag==1)
		{
			double r=a1;
			double p=a2;
			double y=a3;
			mat[0] = cos(p)*cos(y);
			mat[1] = -cos(p)*sin(y);
			mat[2] = -sin(p);
			mat[3] = cos(r)*sin(y)+sin(p)*sin(r)*cos(y);
			mat[4] = cos(r)*cos(y)-sin(p)*sin(r)*sin(y);
			mat[5] = cos(p)*sin(r);
			mat[6] = sin(p)*cos(y)*cos(r)-sin(r)*sin(y);
			mat[7] = -cos(r)*sin(p)*sin(y)-sin(r)*cos(y);
			mat[8] = cos(p)*cos(r);
		}
		if(flag==2)
		{
			double phi=a1;  // a1 b1 c1
			double omega=a2;// a2 b2 c2
			double kappa=a3;// a3 b3 c3
			mat[0] = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
			mat[1] = cos(omega)*sin(kappa);
			mat[2] = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
			mat[3] = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
			mat[4] = cos(omega)*cos(kappa);
			mat[5] = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
			mat[6] = -sin(phi)*cos(omega);
			mat[7] = -sin(omega);
			mat[8] = cos(phi)*cos(omega);
		}
	}

	// get the vector v in unit length
	static double getunit(double* v, int len=3)
	{
		double m=0;
		for(int i=0;i<len;++i)m+=v[i]*v[i];
		m=sqrt(m);
		if(m<1e-8)return m;
		for(int i=0;i<len;++i)v[i]/=m;
		return m;
	}

	// Generate a square unit matrix
	static void OnesMatrix(double* mat, int width)
	{
		memset(mat, 0, sizeof(double)*width*width);
		int idx=0;
		for(int i=0;i<width;++i)
		{
			mat[idx]=1;
			idx+=width+1;
		}
	}

	// 3D vector's cross-product
	static void Xprod_3d(double *v1, double *v2, double* res)
	{
		res[0]=v1[1]*v2[2]-v1[2]*v2[1];
		res[1]=-v1[0]*v2[2]+v1[2]*v2[0];
		res[2]=v1[0]*v2[1]-v1[1]*v2[0];
	}

	static double Dotprod_nd(double*v1, double*v2, int len)
	{
		if(!(len>0))return -0.2333;
		double result=0;
		for(int i=0;i<len;++i)
		{
			result+=v1[i]*v2[i];
		}
		return result;
	}

	// result=v1-v2
	static void Subtraction(double* result, double*v1, double*v2, int len)
	{
		if(!(len>0))return;
		for(int i=0;i<len;++i)
		{
			result[i]=v1[i]-v2[i];
		}
	}

	// Show a m*n matrix; for mx1 vectors n=1 by default
	static void ShowMat(double *p, int precision, int m, int n=1)
	{
		if(!m*n>0)return;
		int iter=0;
		cout<<endl<<"[";
		cout.flags(ios::fixed);
		cout.precision(precision);
		for(int i=0;i<m-1;++i)
		{
			for(int j=0;j<n-1;++j)
			{
				cout<<p[iter]<<", ";
				++iter;
			}
			cout<<p[iter]<<";"<<endl;
			++iter;
		}
		for(int j=0;j<n-1;++j)
		{
			cout<<p[iter]<<", ";
			++iter;
		}
		cout<<p[iter]<<"]"<<endl;
	}

	static void WriteMat(double *p, int m, int n=1, int floatprecision=10)
	{
		if(!m*n>0)return;
		int iter=0;
		ofstream ofs("mat.csv");
		ofs.flags(ios::fixed);
		ofs.precision(floatprecision);
		for(int i=0;i<m-1;++i)
		{
			for(int j=0;j<n-1;++j)
			{
				ofs<<p[iter]<<", ";
				++iter;
			}
			ofs<<p[iter]/*<<";"*/<<endl;
			++iter;
		}
		for(int j=0;j<n-1;++j)
		{
			ofs<<p[iter]<<", ";
			++iter;
		}
		ofs<<p[iter]/*<<";"*/<<endl;
		ofs.close();
	}

	static void Zero(double *p, int n)
	{
		if (n == 0)
			return;
		while (n--) *p++ = 0;
	}
	
	//src m*n, unable to rotate self.
	static bool Rotate(double *src, double *dst, int m, int n)
	{
		if (src == 0 || dst == 0)return 0;
		double *ps = src, *pd = dst;
		for (int i = 0; i < m; i++)
		{
			pd = dst + i;
			for (int j = 0; j < n; j++)
			{
				*pd = *ps;
				pd += m;
				++ps;
			}
		}
		pd = ps = 0; return 1;
	}
	///////////////////////////////////////////////
	//设A为m x l阶矩阵，B为l x n阶矩阵，C为m x n阶矩阵，计算 C=A x B的子程序为：
	static bool MatrixMulti(double * A, double * B, double * C, int M, int N, int L)
	{
		if (A == 0 || B == 0 || C == 0) return 0;
		int i, j, k;
		Zero(C, M*N);

		for (i = 0; i<M; i++)
		{
			for (j = 0; j<N; j++)
			{
				//			for (k=0;k<L;k++) C[i*N+j] += A[i*L+k]*B[k*N+j];
				for (k = 0; k<L; k++) *(C + i*N + j) += *(A + i*L + k)* *(B + k*N + j);	//等效！！！
			}
		}
		return 1;
	}

	//10*10 以及更小对称正定矩阵求逆 a为n*n阶对称正定矩阵，n为矩阵阶数
	static int MatrixInversion(double *a, int n)
	{
		int i, j, k, m;
		double w, g, *b;

		b = new double[n];

		for (k = 0; k <= n - 1; k++)
		{
			w = a[0] + 1.0e-25;
			m = n - k - 1;
			for (i = 1; i < n; i++)
			{
				g = a[i * n];
				b[i] = g / w;

				if (i <= m)
					b[i] = -b[i];

				int tmpOffset1 = (i - 1) * n - 1;
				int tmpOffset2 = i * n;

				for (j = 1; j <= i; j++)
					a[tmpOffset1 + j] = a[tmpOffset2 + j] + g * b[j];
			}
			a[n * n - 1] = 1.0 / w;

			for (i = 1; i <= n - 1; i++)
				a[(n - 1) * n + i - 1] = b[i];
		}

		for (i = 0; i <= n-2; i++)
			for (j = i + 1; j <= n-1; j++)
				a[i * n + j] = a[j * n + i];
		delete[] b;	b = 0;
		return(2);
	}
	// 10*10 以及更小矩阵求逆
	static bool MatrixInversion_General(double *a, int n)
        {
            int        *is, *js, i, j, k, l, u, v;
            double     d, p;
        
            is = new int [n];   // = malloc(n *sizeof(int));
            js = new int [n];   // = malloc(n *sizeof(int));
        
            for(k = 0; k <= n - 1; k++)
            {
                d = 0.0;
            
                for(i = k; i <= n - 1; i++)
                    for(j = k; j <= n - 1; j++)
                    {
                        l = i * n + j;
                        p = fabs(a[l]);
                    
                        if(p > d)
                        {
                            d = p;
                            is[k] = i;
                            js[k] = j;
                        }
                    }
                
                    if(d + 1.0 == 1.0)
                    {
                        delete is;
                        delete js;
                        return (0);
                    }
                
                    if(is[k] != k)
                        for(j = 0; j <= n - 1; j++)
                        {
                            u = k * n + j;
                            v = is[k] * n + j;
                            p = a[u];
                            a[u] = a[v];
                            a[v] = p;
                        }
                    
                        if(js[k] != k)
                            for(i = 0; i <= n - 1; i++)
                            {
                                u = i * n + k;
                                v = i * n + js[k];
                                p = a[u];
                                a[u] = a[v];
                                a[v] = p;
                            }
                            l = k * n + k;
                            a[l] = 1.0 / a[l];
                        
                            for(j = 0; j <= n - 1; j++)
                                if(j != k)
                                {
                                    u = k * n + j;
                                    a[u] = a[u] * a[l];
                                }
                            
							for(i = 0; i <= n - 1; i++)
								if(i != k)
									for(j = 0; j <= n - 1; j++)
										if(j != k)
										{
											u = i * n + j;
											a[u] = a[u] - a[i * n + k] * a[k * n + j];
										}
									
									for(i = 0; i <= n - 1; i++)
										if(i != k)
										{
											u = i * n + k;
											a[u] = -a[u] * a[l];
										}
            }
        
            for(k = n - 1; k >= 0; k--)
            {
                if(js[k] != k)
                    for(j = 0; j <= n - 1; j++)
                    {
                        u = k * n + j;
                        v = js[k] * n + j;
                        p = a[u];
                        a[u] = a[v];
                        a[v] = p;
                    }
                
                    if(is[k] != k)
                        for(i = 0; i <= n - 1; i++)
                        {
                            u = i * n + k;
                            v = i * n + is[k];
                            p = a[u];
                            a[u] = a[v];
                            a[v] = p;
                        }
            }
            delete is; delete js;
            return (1);
    }


//---------------------------------------------------

};