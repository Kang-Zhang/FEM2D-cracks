#ifndef __VECTOR_HPP__FED
#define __VECTOR_HPP__FED

#include <cmath>

namespace VectLib
{
    
template <int N, typename T>
class Vector {
public:

    Vector operator + (const Vector & v2) const
    {
	Vector res;
	for (int i = 0 ; i < N ; i ++)
	    res._data[i] = _data[i] + v2._data[i];
	return res;
    }

    Vector operator - (const Vector & v2) const
    {
	Vector res;
	for (int i = 0 ; i < N ; i ++)
	    res._data[i] = _data[i] - v2._data[i];
	return res;
    }

    friend Vector operator * (const T & x, const Vector & vec) {
	Vector res;
	for (int i = 0 ; i < N ; i ++)
	    res._data[i] = x * vec._data[i];
	return res;
    }

    friend Vector operator / (const Vector & vec, const T & x) {
	Vector res;
	for (int i = 0 ; i < N ; i ++)
	    res._data[i] = vec._data[i] / x;
	return res;
    }

    Vector & operator += (const Vector & vec) {
	for (int i = 0 ; i < N ; i ++)
	    _data[i] += vec._data[i];
	return * this;
    }

    Vector & operator -= (const Vector & vec) {
	for (int i = 0 ; i < N ; i ++)
	    _data[i] -= vec._data[i];
	return * this;
    }

    T operator * (const Vector & vec) const {
	T res = 0;
	for (int i = 0 ; i < N ; i ++)
	    res += _data[i] * vec._data[i];
	return res;
    }

    static T normSq (const Vector & vec) {
	T res = 0;
	for (int i = 0 ; i < N ; i ++)
	    res += vec._data[i] * vec._data[i];
	return res;
    }

    static Vector normalize (const Vector & vec, const T & eps) {
	T len = sqrt (normSq (vec));
	if (len >= eps)
	    return vec / len;
	else
	    return vec;
    }

 /*
    friend Vector operator*(const T & scalar,const Vector & vec)
	{
	    Vector ans;
	    for(int i = 0 ; i < N ; i++)
		ans[i] = scalar * vec._data[i];
	    return ans;
	}
 */

/*
    Vector operator * (const T & x) const
    {
	Vector res;
	for (int i = 0 ; i < N ; i ++)
	    res._data[i] = x * _data[i];
	return res;
    }
*/
    
    Vector () {;}
    
protected:
    T _data [N];
};

/*    
template <int N, class T>
Vector <N,T> operator * (const T & x, const Vector <N,T> & vec)
{
    Vector <N,T> res;
    for (int i = 0 ; i < N ; i ++)
	res._data[i] = x * vec._data[i];
    return res;
}
*/

template <typename T>
class Vector2 : public Vector <2, T>
{
public:
    Vector2 (const T & val1, const T & val2)
    {
	(this-> _data) [0] = val1;
	(this-> _data) [1] = val2;
    }
    Vector2 (const Vector <2,T> & v) : Vector <2,T> (v) {;}
    Vector2 () : Vector <2,T> () {;}
    T & x () { return this-> _data[0]; }
    T & y () { return this-> _data[1]; }
    const T & x () const { return this-> _data[0]; }
    const T & y () const { return this-> _data[1]; }
};

template <typename T>
class Vector3 : public Vector <3, T>
{
public:
    Vector3 (const T & val1, const T & val2, const T & val3)
    {
	this-> _data [0] = val1;
	this-> _data [1] = val2;
	this-> _data [2] = val3;
    }
    Vector3 (const Vector <3,T> & v) : Vector <3,T> (v) {;}
    Vector3 () : Vector <3,T> () {;}
    Vector3 (const Vector2 <T> & v) : Vector2 <T> (v) { this-> _data [2] = 0; }
    T & x_ () { return this-> _data[0]; }
    T & y_ () { return this-> _data[1]; }
    T & z_ () { return this-> _data[2]; }
    const T & x () const { return this-> _data[0]; }
    const T & y () const { return this-> _data[1]; }
    const T & z () const { return this-> _data[2]; }
};

typedef Vector2 <double> V2D;
typedef Vector3 <double> V3D;

template <int N, typename T>
T dist_to_line_seg (const Vector <N,T> & p, const Vector <N,T> & p1, const Vector <N,T> & p2)
// returns the distance of a point p to a line segment p1,p2
{
    Vector <N,T> v12 = p2 - p1;
    T u = (p-p1) * v12 / Vector<N,T>::normSq (v12);
//    std::cerr << "u = " << u << "\n";
    if (! std::isfinite (u))
	// p1 & p2 must nearly coinside
	return std::sqrt (Vector<N,T>::normSq (p - p1));
    if (u <= 0)
	return std::sqrt (Vector<N,T>::normSq (p - p1));
    if (u >= 1)
	return std::sqrt (Vector<N,T>::normSq (p - p2));
    return std::sqrt (Vector<N,T>::normSq (p - u * v12 - p1));
}

};
    
// ifndef __VECTOR_HPP_FED
#endif
