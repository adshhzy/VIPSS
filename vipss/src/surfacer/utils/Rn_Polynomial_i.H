/* -------------------------------------------------------------------------
 * DESCR   :	Accessors for 1 dimension
 * ------------------------------------------------------------------------- */
inline double 
RNPolynomial::operator[](int i) const 
{ 
    return Array<double>::operator[](i); 
}

inline double &
RNPolynomial::operator[](int i) 
{ 
    return Array<double>::operator[](i); 
}

inline double  
RNPolynomial::Coef(int which) const 
{ 
    ASSERT( m_iDim == 1 );
    return (*this)[which]; 
}

inline double &
RNPolynomial::Coef(int which) 
{ 
    ASSERT( m_iDim == 1 );
    return (*this)[which]; 
}

/* -------------------------------------------------------------------------
 * DESCR   :	Accessors for 2 dimensions
 * ------------------------------------------------------------------------- */
inline double &
RNPolynomial::Coef(int which_s, int which_t) 
{ 
    ASSERT( m_iDim == 2 );
    return (*this)[which_t * (Degree()+1) + which_s]; 
}

inline double 
RNPolynomial::Coef(int which_s, int which_t) const
{ 
    ASSERT( m_iDim == 2 );
    return (*this)[which_t * (Degree()+1) + which_s]; 
}
