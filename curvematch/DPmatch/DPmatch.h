//#include "params.h"
#include "shape.h"
#include <stdio.h>
#include <vector>

class DPmatch
{
    
public:
    
    shape *v_pshape1;
    shape *v_pshape2;
    
    DPmatch();
    DPmatch(int dim, int siz);
    ~DPmatch();
    bool Read_v1_v2(char *);
    
    void Initialize(int dim, int siz, float *v11, float *v12, float *v13, float *v21, float *v22, float *v23);
    
    float MatchQ(float * gamma);
    
    float Match_CostQ(int ,int ,int ,int);
    
    bool WriteGammaV(FILE * ,float *,int);
    
    void linint(float *xnew, float *ynew, int cnt, float *xx, float *yy, int n);
    void copy(int T, float *src, float *dest);
    
    void splint(float *xa, float *ya, float *y2a, int n, float x, float *y);
    
    void spline(float *x, float *y, int n, float yp1, float ypn, float *y2);
};
