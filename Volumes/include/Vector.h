#pragma once

#include "math.h"

class Vector {
public:
	// initializes all components to zero
	Vector();

	// initializes with specified components
	Vector(double x, double y, double z = 0.0);

	// copy constructor
	Vector(const Vector& v);

	// access
	double& operator[](int index);
	const double& operator[](int index) const;

	// math
	Vector operator*(double s) const;
	Vector operator/(double s) const;
	Vector operator+(const Vector& v) const;
	Vector operator-(const Vector& v) const;
	Vector operator-() const;
	double operator*(const Vector& v) const;

	Vector& operator*=(double s);
	Vector& operator/=(double s);
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);


	Vector operator^(const Vector& v) const;

	// returns Euclidean length
	double norm() const;

	// returns Euclidean length squared
	double norm2() const;

	// normalizes vector
	void normalize();

	// returns unit vector in the direction of this vector
	Vector unit() const;

	// members
	double x, y, z;
};

Vector operator*(double s, const Vector& v);
double dot(const Vector& u, const Vector& v);
Vector cross(const Vector& u, const Vector& v);

inline Vector::Vector():
x(0.0),
y(0.0),
z(0.0)
{

}

inline Vector::Vector(double x_, double y_, double z_):
x(x_),
y(y_),
z(z_)
{

}

inline Vector::Vector(const Vector& v):
x(v.x),
y(v.y),
z(v.z)
{

}

inline double& Vector::operator[](int index)
{
	return (&x)[index];
}

inline const double& Vector::operator[](int index) const
{
	return (&x)[index];
}

inline Vector Vector::operator*(double s) const
{
	return Vector(x*s, y*s, z*s);
}

inline Vector Vector::operator/(double s) const
{
	return (*this)*(1.0/s);
}

inline Vector Vector::operator+(const Vector& v) const
{
	return Vector(x + v.x, y + v.y, z + v.z);
}

inline double Vector::operator*(const Vector& v) const
{
	return x*v.x + y*v.y + z*v.z;
}

inline Vector Vector::operator^(const Vector& v) const
{
	return Vector(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
}

inline Vector Vector::operator-(const Vector& v) const
{
	return Vector(x - v.x, y - v.y, z - v.z);
}

inline Vector Vector::operator-() const
{
	return Vector(-x, -y, -z);
}

inline Vector& Vector::operator*=(double s)
{
	x *= s;
	y *= s;
	z *= s;

	return *this;
}

inline Vector& Vector::operator/=(double s)
{
	(*this) *= (1.0/s);

	return *this;
}

inline Vector& Vector::operator+=(const Vector& v)
{
	x += v.x;
	y += v.y;
	z += v.z;

	return *this;
}

inline Vector& Vector::operator-=(const Vector& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;

	return *this;
}

inline double Vector::norm() const
{
	return sqrt(norm2());
}

inline double Vector::norm2() const
{
	return dot(*this, *this);
}

inline void Vector::normalize()
{
	(*this) /= norm();
}

inline Vector Vector::unit() const
{
   return (*this) / norm();
}

inline Vector operator*(double s, const Vector& v)
{
	return v*s;
}

inline double dot(const Vector& u, const Vector& v)
{
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline Vector cross(const Vector& u, const Vector& v)
{
	return Vector(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}
