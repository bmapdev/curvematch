#include <stdlib.h>
#include "shape.h"
#include <math.h>

shape::shape()
{
		
}

shape::shape(int n, int v_iT)
{
    m_v11 = NULL;
    m_v12 = NULL;
    m_v13 = NULL;
    
    m_v11 = (float *)malloc(v_iT*sizeof(float));
    m_v12 = (float *)malloc(v_iT*sizeof(float));
    if(n == 3)
        m_v13 = (float *)malloc(v_iT*sizeof(float));
	m_iT = v_iT;
    m_n = n;
}

shape::~shape()
{
	free(m_v11);
	free(m_v12);
	free(m_v13);    
}