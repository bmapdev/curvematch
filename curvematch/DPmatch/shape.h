#if !defined(SHAPE_H)
#define SHAPE_H
#include <math.h>

const float PI =  (float)3.14159265358979;

class shape
{
	public:
		//Member variables
		float *m_v11;
		float *m_v12;
        float *m_v13;
                
		int m_iT;
		int m_n;
		//Member functions
		shape();//default constructor
		shape(int, int);
		~shape();
};
#endif
