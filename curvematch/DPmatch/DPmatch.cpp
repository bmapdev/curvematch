#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <cblas.h>
#include <iostream>
#include <vector>
#include "DPmatch.h"
using namespace std;

//void copy(int T, float *src, float *dest);

DPmatch::DPmatch()
{
    v_pshape1 = NULL;
    v_pshape1 = NULL;
    
}

DPmatch::DPmatch(int dim, int siz)
{
    v_pshape1 = new shape(dim,siz);
    v_pshape2 = new shape(dim,siz);
}


DPmatch::~DPmatch()
{
    free(v_pshape1);
    free(v_pshape2);
}


void DPmatch::Initialize(int dim, int siz, float *v11, float *v12, float *v13, float *v21, float *v22, float *v23)
{
    copy(siz,v11,(v_pshape1)->m_v11);
    copy(siz,v12,(v_pshape1)->m_v12);

    if (dim == 3)
        copy(siz,v13,(v_pshape1)->m_v13);
    
    copy(siz,v21,(v_pshape2)->m_v11);
    copy(siz,v22,(v_pshape2)->m_v12);
    if (dim == 3)
        copy(siz,v23,(v_pshape2)->m_v13);

}


bool DPmatch::Read_v1_v2(char *v_szFileName)
{
	// Open the Dir_Fn file to read the array of Direction functions

	FILE *vl_fpDir_Fn = NULL;
	vl_fpDir_Fn = fopen(v_szFileName,"rb");
	int vl_iNum_Shapes = 0;
	if( vl_fpDir_Fn == NULL)
	{
		printf("\nError: %s File not Found\n",v_szFileName);
		return false;
	}
	int vl_iT = 0;
	size_t vl_readsize = 0;
    int n = 0;
    vl_readsize = fread(&vl_iNum_Shapes,sizeof(vl_iNum_Shapes),1,vl_fpDir_Fn);
    vl_readsize = fread(&n,sizeof(n),1,vl_fpDir_Fn);

	vl_readsize = fread(&vl_iT,sizeof(vl_iT),1,vl_fpDir_Fn);

	v_pshape1 = new shape(n,vl_iT);	
	v_pshape2 = new shape(n,vl_iT);

	float *vect = new float[vl_iT];		
    
    float *v11 = new float[vl_iT];		
	float *v12 = new float[vl_iT];		
	float *v13 = new float[vl_iT];		    
    float *v21 = new float[vl_iT];		
	float *v22 = new float[vl_iT];		
	float *v23 = new float[vl_iT];		    
	
	vl_readsize = fread(v11,sizeof(float),vl_iT,vl_fpDir_Fn);
    
	vl_readsize = fread(v12,sizeof(float),vl_iT,vl_fpDir_Fn);

    if(n == 3)
        vl_readsize = fread(v13,sizeof(float),vl_iT,vl_fpDir_Fn);

	vl_readsize = fread(v21,sizeof(float),vl_iT,vl_fpDir_Fn);
	vl_readsize = fread(v22,sizeof(float),vl_iT,vl_fpDir_Fn);

    if(n == 3)
        vl_readsize = fread(v23,sizeof(float),vl_iT,vl_fpDir_Fn);

    Initialize(n, vl_iT, v11, v12, v13, v21, v22, v23);
    
	fclose(vl_fpDir_Fn);
	
	delete [] vect;

	return true;
}




float DPmatch::MatchQ(float *gamma)
{
//    int Nbrs[][2] = {{0 ,0},{0, 1},{ 1, 0},{ 1, 2},{ 2, 1},{ 0, 2},{ 2, 0},{ 0, 3},{ 2, 3},{ 3 ,2},{ 3 ,2},{ 0, 4},{ 1, 4}, \
    { 2, 4},{ 3, 4},{ 4, 3},{ 4, 2},{ 4, 1},{ 4, 0} };
    int Nbrs[][2] = {{1 ,1},{1, 2},{ 2, 1},{ 2, 3},{ 3, 2},{ 1, 3},{ 3, 1},{1, 4},{ 3, 4},{4 ,3},{ 4,1},{ 1, 5},{ 2, 5}, \
    { 3, 5},{ 4, 5},{5, 4},{ 5, 3},{5, 2},{ 5, 1},{1,6},{5,6},{6,5},{6,1},{1,7},{2,7},{3,7},{4,7},{5,7},{6,7}, \
		{4,10}, {4,30}, {4,40}, {4,50}, {5,50},{6,50},{8,50},{10,50},{15,50},{20,50},{25,50},{30,50},{35,50},{40,50}    };
    const int NBR_SIZ = 43;	
//	int Nbrs[][2] = {{1, 1}, {1, 2}, {2, 1}, {2, 3}, {3, 2},{ 1, 3},{ 3, 1}, {1, 4}, {3, 4}, {4, 3}, {4, 1}};
//	const int NBR_SIZ = 11;
    int NumPlot = 15;
    float shfx = 0.2;
    float shfy = 0.2;
    int N = 0;
    N = v_pshape1->m_iT;
    int i = 0;
    int j = 0;
    int Num = 0;
    float CandE[NBR_SIZ];
    int k = 0,l = 0;
    float minCandE = 10000;
    int minCandE_idx = 0;
    float **Path_x = NULL;
    float **Path_y = NULL;
    float cost = 0;
    Path_x = (float **)malloc(N*sizeof(float *));  // This is Path(i,j,1) in Match.m
    Path_y = (float **)malloc(N*sizeof(float *));  // This is Path(i,j,2) in Match.m
	float **Energy = NULL;

    float *x = NULL; 
    float *y = NULL;
    float *xnew = NULL; 
    float *ynew = NULL;
    x = (float *)malloc(N*sizeof(float));
    y = (float *)malloc(N*sizeof(float));
    float *xx1 = (float *) malloc(N*sizeof(float) );	

    xnew = (float *)malloc(N*sizeof(float));
    ynew = (float *)malloc(N*sizeof(float));

    int cnt = 0;
    Energy = (float **)malloc(N*sizeof(float *));

    for(i = 0;i < N;i ++)
    {
            Energy[i] = (float *)calloc(N,sizeof(float));	
            Path_x[i] = (float *)calloc(N,sizeof(float));	
            Path_y[i] = (float *)calloc(N,sizeof(float));	
            //Forming energies associated with different paths
            Energy[0][i] = 5;	
            Energy[i][0] = 5;
    }

/*		for(i = 1; i < N ; i ++)
    {
            Energy[0][i] = 5;	
            Energy[i][0] = 5;
    }*/
    Energy[0][0] = 0;
    xx1[0] = 0;
    for(i = 1 ; i < N; i ++)
    {
            for(j = 1; j < N ; j ++)
            {
                    minCandE = 10000;
                    for(Num = 0; Num < NBR_SIZ; Num ++)
                    {
                            k = i - Nbrs[Num][0];
                            l = j - Nbrs[Num][1];
                            if(k >= 0 && l >= 0)
                            {
                                    CandE[Num] = Energy[k][l] + Match_CostQ(k,l,i,j);
                            }
                            else
                            {
                                    CandE[Num] = 10000;
                            }
                            if(CandE[Num] < minCandE )
                            {
                                    minCandE = CandE[Num];
                                    minCandE_idx = Num;	
                            }
                            Energy[i][j] = minCandE;

                            Path_x[i][j] = i - Nbrs[minCandE_idx][0];
                            Path_y[i][j] = j - Nbrs[minCandE_idx][1];		
                    }
            }	
            xx1[i] = (float ) i/(N - 1);				
    }

    x[0] = N-1; y[0] = N-1;
    while ( x[cnt] > 0 )
    {
            i = (int) y[cnt];
            j = (int) x[cnt];
            y[cnt + 1] = Path_x[i] [j]  ;
            x[cnt + 1] = Path_y[i] [j] ;
            cnt ++;
    }
//    printf("\n%d\n",cnt);

    for (i = 0; i < cnt; i ++)
    {
            xnew[i] = (x[cnt-i-1] -x[cnt-1] )/(x[0] - x[cnt - 1]);
            ynew[i] = (y[cnt-i-1] -y[cnt-1] )/(y[0] - y[cnt - 1]);
    }

    xnew[cnt-1] = 1;
    ynew[cnt-1] = 1;
    xnew[0] = 0;
    ynew[0] = 0;

    linint(xnew, ynew, cnt, xx1, gamma, N);
//		cost = Cost_Group_Action_By_Gamma(v_pshape1->m_pfPhi, v_pshape1->m_pfTheta,  \
                    v_pshape2->m_pfPhi, v_pshape2->m_pfTheta, gamma, N,a,b);

    cost = 0;				
    for(i = 0;i < N;i ++)
    {
            free(Energy[i]);
            free(Path_x[i]);
            free(Path_y[i]);
    }

    free(Path_x);	
    free(Path_y);
    free(Energy);

    free(x);
    free(y);
    free(xx1);
    free(xnew);
    free(ynew);




    return cost;
}


float DPmatch::Match_CostQ(int k,int l,int i,int j)
{
	int x;
	float y;
	int y1,y2;
	float slope = 0;
	float E,E1,E2;
	float f;
	float vec11, vec12, vec21, vec22,vec23;
	
	E1 = 0; E2 = 0;
	E = 0;
	slope =(float ) ( i - k)/(j - l);
//	printf("slope1 %f\n",slope);
	if (slope == 0)
		printf("\nslope zero\n");
	for(x = l; x <= j ; x ++)
	{
            y = k + (x - l) * slope;
            y1 = (int )floorf(y);
            y2 = (int )ceilf(y);	
            f = y - y1;
            vec21 = (f*v_pshape2->m_v11[y2] + (1 - f)*v_pshape2->m_v11[y1])*sqrt(slope);	
            vec22 = (f*v_pshape2->m_v12[y2] + (1 - f)*v_pshape2->m_v12[y1])*sqrt(slope);
            if(v_pshape1->m_n == 3)
            {
                vec23 = (f*v_pshape2->m_v13[y2] + (1 - f)*v_pshape2->m_v13[y1])*sqrt(slope);
                E2 = E2 + (v_pshape1->m_v11[x] - vec21)*(v_pshape1->m_v11[x] - vec21) + (v_pshape1->m_v12[x] - vec22)*(v_pshape1->m_v12[x] - vec22) + (v_pshape1->m_v13[x] - vec23)*(v_pshape1->m_v13[x] - vec23);
            }
            else
            {
                E2 = E2 + (v_pshape1->m_v11[x] - vec21)*(v_pshape1->m_v11[x] - vec21) + (v_pshape1->m_v12[x] - vec22)*(v_pshape1->m_v12[x] - vec22);
            }

	}
	E = E2/(v_pshape1->m_iT);
	//cout << E << " ";	
	return E;
        
}



bool DPmatch::WriteGammaV(FILE *vl_fpDirFn ,float *gamma,int N)
{
	size_t bytes_written;

	if( vl_fpDirFn == NULL)
	{
		printf("\nError: File cannot be opened for appending");
		return false;
	}
	bytes_written = fwrite((void*)&N,sizeof(int),1,vl_fpDirFn);
	bytes_written = fwrite((void*)gamma,sizeof(float),N,vl_fpDirFn);
	if(bytes_written == N)
		return true;
	else
		return false;
}

void DPmatch::linint(float *xnew, float *ynew, int cnt, float *xx, float *yy, int n)
{
	int i = 0;
	int idx = 0;
	//Assume xnew and xx are sorted. 
	//Find the interval where xx[0] is located
	float m = 0;
	for (i = 0; i < n; i ++)
	{
		while(idx < cnt - 1)
		{
			if(xx[i] >= xnew[idx] && xx[i] <= xnew[idx+1] )
			{
				yy[i] = ynew[idx] + (xx[i] - xnew[idx]) *( ynew[idx+1] - ynew[idx] ) / ( xnew[idx+1] - xnew[idx] );
				break;
			}
			else if ( xx[i] <= xnew[idx+1] )
			{
				// Need to extrapolate
				// Use the slope 	of points corresponding to idx and idx + 1
				m = (ynew[idx+1]  - ynew[idx])/(xnew[idx+1]  - xnew[idx]);
				yy[i] = ynew[idx] + m * ( xx[i] - xnew[idx]);
				break;
			}
			else if(xx[i] >= xnew[cnt-1] )
			{
				//This is beyond the boundary of xnew
				// Need to extrapolate
				// Use the slope 	of points corresponding to cnt -1 and cnt -2
				m = (ynew[cnt-1]  - ynew[cnt-2])/(xnew[cnt-1]  - xnew[cnt-2]);
				yy[i] = ynew[cnt - 1] + m * ( xx[i] - xnew[cnt-1]);			
				break;
			}
			else 
				idx ++;
			}
		}
}

void DPmatch::spline(float *x, float *y, int n, float yp1, float ypn, float *y2)
{
	int i,k;
	float p,qn,sig,un,*u;
	
	u = (float *)malloc(n*sizeof(float));
	if (yp1 > 0.99e20) 
		y2[1]=u[1]=0.0;
	else 
	{ 
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}

	
	for (i=2;i<=n-1;i++) 
	{ 
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;

	}

	if (ypn > 0.99e30) 
		qn=un=0.0; 
	else 
	{ 
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}

	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--) 
	{
		y2[k]=y2[k]*y2[k+1]+u[k]; 
	}
	free(u);
}

void DPmatch::splint(float *xa, float *ya, float *y2a, int n, float x, float *y)
{
	int klo,khi,k;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) 
	{
		k=(khi+klo) >> 1;
		if (xa[k] > x) 
			khi=k;
		else
			klo=k;
		h=xa[khi]-xa[klo];
		if (h == 0.0) 
			printf("\n Bad xa input to routine splint\n ");
		a=( xa[khi]-x)/h;
		b=(x-xa[klo])/h;
		*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

	}
}


void DPmatch::copy(int T, float *src, float *dest)
{
	for(int i = 0; i < T; i++)
	{
		dest[i] = src[i];
	}
	
}
