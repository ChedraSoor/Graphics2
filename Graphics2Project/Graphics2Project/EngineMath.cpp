/**
* @file EngineMath.cpp
*
*/

#include "EngineMath.h"
#include "math.h"

//////////////////////////////////////////////////////////////////////////
// Common math functions
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// General Utility functions
//////////////////////////////////////////////////////////////////////////

// Are two floating point numbers equal to each other
// Floating Point Error Safe
//
// IN:		a		The first number
//			b		The second number
//
// RETURN: TRUE iff |a-b| < Tolerance
//
// NOTE:	EPSILON is tolerance
bool IsEqual(float a, float b)
{
	// NOTE: Do not modify.
	return fabs(a - b) < EPSILON;
}

// Is a floating point value equal to zero
// Floating Point Error Safe
//
// IN:		a		The number to check
//
// RETURN:	TRUE iff |a| < Tolerance
//
// NOTE:	Tolerance set by EPSILON
bool IsZero(float a)
{
	// NOTE: Do not modify
	return (fabs(a))<EPSILON;
}

// RETURN: MAX of two numbers
float Max(float a, float b)
{
	// NOTE: Do not modify.
	return (a > b) ? a : b;
}

// RETURN: MIN of two numbers
float Min(float a, float b)
{
	// NOTE: Do not modify.
	return (a < b) ? a : b;
}

// RETURN: Converts input to radian measure
float Degrees_To_Radians(float Deg)
{
	// NOTE: Do not modify.
	return Deg * PI / 180.0f;
}

// RETURN: Converts input to degree measure
float Radians_To_Degrees(float Rad)
{
	// NOTE: Do not modify.
	return Rad * 180.0f / PI;
}
////////////////////////////////////////////////////////////////////////
// Linear Algebra Functions Day 1
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Vector Functions
//////////////////////////////////////////////////////////////////////////

// Check if two TVECTOR's are equal to each other
//
// IN:		v		First Vector
//			w		Second Vector
//
// RETURN:  True if v==w, False otherwise
//
// NOTE:	Use's all four components
//			Should be floating point error safe.
bool Vector_IsEqual(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	return IsEqual(v.w, w.w) && IsEqual(v.x, w.x) && IsEqual(v.y, w.y) && IsEqual(v.z, w.z);
	
}

// ADD two TVECTOR's togother
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v + w
//
// NOTE:	Use's all four components
TVECTOR Vector_Add(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR vector;

	vector.w = v.w + w.w;
	vector.x = v.x + w.x;
	vector.y = v.y + w.y;
	vector.z = v.z + w.z;
	
	v = vector;

	return v;
}

// SUBTRACT one TVECTOR from another
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v - w
//
// NOTE:	Use's all four components
TVECTOR Vector_Sub(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR vector;

	vector.w = v.w - w.w;
	vector.x = v.x - w.x;
	vector.y = v.y - w.y;
	vector.z = v.z - w.z;

	v = vector;

	return v;
}

// MULTIPLY all four components of a TVECTOR by a scalar
//
// IN:		v		The vector to scale
//			s		The value to scale by
//
// RETURN:  s * v
TVECTOR Vector_Scalar_Multiply(TVECTOR v, float s)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR vector;

	vector.w = v.w * s;
	vector.x = v.x * s;
	vector.y = v.y * s;
	vector.z = v.z * s;

	v = vector;

	return v;
}

// NEGATE all the components of a TVECTOR
//
// IN:		v		The vector to negate
//
// RETURN:	-1 * v
//
// NOTE:	Use's all four components
TVECTOR Vector_Negate(TVECTOR v)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR vector;

		vector.w = -1 * v.w;
		vector.x = -1 * v.x;
		vector.y = -1 * v.y;
		vector.z = -1 * v.z;

	//v = vector;
	return vector;
}

// Perform a Dot Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (DOT) w
//
// NOTE:	Use's all four components
float Vector_Dot(TVECTOR v, TVECTOR w)
{
	float prod = (v.x * w.x) + (v.y * w.y) + (v.z * w.z) + (v.w * w.w);
	// TODO LAB 1: Replace with your implementation.
	return  prod;
}

// Perform a Cross Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (CROSS) w
//
// NOTE:	The w-component of each vector is not used.
//			The resultant vector will have a w-component of zero.
TVECTOR Vector_Cross(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR value;
	value.x = ((v.y * w.z) - (v.z * w.y));
	value.y = -((v.x * w.z) - (v.z * w.x));
	value.z = ((v.x * w.y) - (v.y * w.x));
	value.w = 0;

	return value;
}

// Find the squared length of a TVECTOR
//
// IN:		v		The vector to find the squared length of
//
// RETURN:	Squared Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_LengthSq(TVECTOR v)
{
	float square = Vector_Dot(v, v);
	return square;
}

// Find the length of a TVECTOR
//
// IN:		v		The vector to find the length of
//
// RETURN:	Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_Length(TVECTOR v)
{
	float length = sqrtf(Vector_LengthSq(v));
	// TODO LAB 1: Replace with your implementation.
	return length;
}

// Normalize a TVECTOR
//
// IN:		v		The vector to normalize
//
// RETURN:	Normalized version of v
//
// NOTE:	Use's all four components
TVECTOR Vector_Normalize(TVECTOR v)
{
	// TODO LAB 1: Replace with your implementation.
	float normal;
	normal = Vector_Length(v);
	if (IsZero(normal))
	{
		v.w = 0;
		v.x = 0;
		v.y = 0;
		v.z = 0;

	}
	else
	{
		v.x = (v.x / normal);
		v.y = (v.y / normal);
		v.z = (v.z / normal);
		v.w = (v.w / normal);
	}
	return v;
}

// Makes a TVECTOR's w-component normalized
//
// IN:		v		The vector (point object) to homogenise
//
// RETURN:	The homogenised vector (point)
//
// NOTE:	If the w-component of the vector is 0 then the
//			function will return a zero vector with a w-component
//			of 0.
TVECTOR Vector_Homogenise(TVECTOR v)
{
	// TODO LAB 1: Replace with your implementation.
	if (IsZero(v.w))
	{
		v.w = 0;
		v.x = 0;
		v.y = 0;
		v.z = 0;

	}
	else
	{
		v.x = (v.x / v.w);
		v.y = (v.y / v.w);
		v.z = (v.z / v.w);
		v.w = (v.w / v.w);
	}
	return v;
}

// Get a TVECTOR made from the maximun components of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A maximized vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Maximize(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR max;

	max.x = Max(v.x, w.x);
	max.y = Max(v.y, w.y);
	max.z = Max(v.z, w.z);
	max.w = Max(v.w, w.w);

	return max;

}

// Get a TVECTOR made from the minimum components of two TVECTOR's
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A minimum vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Minimize(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR min;

	min.x = Min(v.x, w.x);
	min.y = Min(v.y, w.y);
	min.z = Min(v.z, w.z);
	min.w = Min(v.w, w.w);

	return min;
	
}

// Get a TVECTOR made from the average of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A vector made from the average of two vectors
//
// NOTE:	Use's all four components

TVECTOR Vector_Average(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	TVECTOR average;

	average.x = (v.x + w.x) / 2;
	average.y = (v.y + w.y) / 2;
	average.z = (v.z + w.z) / 2;
	average.w = (v.w + w.w) / 2;

	return average;
}

// Find the angle between two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:  The angle in degrees between the two vectors
//
// NOTE:	If either vector is a zero vector then the return
//			value will be 0.
float Vector_AngleBetween(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.
	float value;
		float angle;

	 value = Vector_Dot(v, w) / (Vector_Length(v)*Vector_Length(w));
	 angle = acosf(value);
	return Radians_To_Degrees(angle);
}

// Get the distance one TVECTOR points in the direction of another
// TVECTOR
//
// IN:		v		The first vector
//			w		The direction of the component
//
// RETURN:	The distance that v points in the direction of w.
//
// NOTE:	If w or v is a zero vector then the return value is zero.
float Vector_Component(TVECTOR v, TVECTOR w)
{
	// TODO LAB 1: Replace with your implementation.'
	float component;
	if (w.x == 0 && w.y == 0 && w.z == 0 && w.w == 0)
	{
		component = 0;
		return component;
	}
	else
	{
		Vector_Normalize(w);
		component = Vector_Dot(v, w) / Vector_Length(w);
		return component;
	}
	
}

// Get the TVECTOR that represents v projected on w.
//
// IN:		v		The first vector
//			w		The direction of the projection
//
// RETURN:	The projection of v onto w
//
// NOTE:	If w or v is a zero vector then the return value is zero.
TVECTOR Vector_Project(TVECTOR v, TVECTOR w)
{
	auto comp = Vector_Component(v, w);

	TVECTOR mat = Vector_Normalize(w);
	TVECTOR project;

	project.x = mat.x * comp;
	project.y = mat.y * comp;
	project.z = mat.z * comp;
	project.w = mat.w * comp;

	return project;
}

////////////////////////////////////////////////////////////////////////
// Functions Lab  #2
///////////////////////////////////////////////////////////////////////


// Get the reflection of v across w
//
// IN:		v		The vector to reflect
//			w		The "axis" to reflect across
//
// RETURN:	v reflected across w
//
// NOTE:	If w is a zero vector then return -v.
TVECTOR Vector_Reflect(TVECTOR v, TVECTOR w)
{
	// TODO LAB 2: Replace with your implementation.
	TVECTOR vec;
	if (w.x == 0 && w.y == 0 && w.z == 0 && w.w == 0)
	{
		return Vector_Negate(v);
	}
	else
	{
		// R = (2 * Vector_Dot(v,w)) * w - v ===> (2 * dot) * W - V; 
		w = Vector_Normalize(w);
		float dot = Vector_Dot(v, w);
		vec.x = (2 * dot) * w.x - v.x;
		vec.y = (2 * dot) * w.y - v.y;
		vec.z = (2 * dot) * w.z - v.z;
		vec.w = (2 * dot) * w.w - v.w;
		v = vec;
		return v;
	}
}

//////////////////////////////////////////////////////////////////////////
// Matrix Functions
//////////////////////////////////////////////////////////////////////////

// Get a [0] matrix
//
// RETURN: A 0 4x4 matrix
TMATRIX Matrix_Zero(void)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m;
	m._e11 = 11;
	m._e14 = 14;
	m._e23 = 23;

	m._e11 = 0;
	m._e12 = 0;
	m._e13 = 0;
	m._e14 = 0;

	m._e21 = 0;
	m._e22 = 0;
	m._e23 = 0;
	m._e24 = 0;

	m._e31 = 0;
	m._e32 = 0;
	m._e33 = 0;
	m._e34 = 0;

	m._e41 = 0;
	m._e42 = 0;
	m._e43 = 0;
	m._e44 = 0;

	
	return m;
}

// Get a [I] matrix
//
// RETURN: A 4x4 Identity matrix
TMATRIX Matrix_Identity(void)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e22 = 1;
	m._e33 = 1;
	m._e44 = 1;

	return m;
}

// Get a translation matrix
//
// IN:		x		Amount of translation in the x direction
//			y		Amount of translation in the y direction
//			z		Amount of translation in the z direction
//
// RETURN:	The translation matrix
TMATRIX Matrix_Create_Translation(float x, float y, float z)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e22 = 1;
	m._e33 = 1;
	m._e44 = 1;

	m._e41 = x;
	m._e42 = y;
	m._e43 = z;

	return m;
}

// Create a scale matrix
//
// IN:		x		Amount to scale in the x direction
//			y		Amount to scale in the y direction
//			z		Amount to scale in the z direction
//
// RETURN:	The scale matrix
TMATRIX Matrix_Create_Scale(float x, float y, float z)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e22 = 1;
	m._e33 = 1;
	m._e44 = 1;

	m._e11 = m._e11 * x;
	m._e22 = m._e22 * y;
	m._e33 = m._e33 * z;
	m._e44 = 1;

	return m;
}

// Get a rotation matrix for rotation about the x-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A X-Rotation Matrix
TMATRIX Matrix_Create_Rotation_X(float Deg)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e44 = 1;

	float rad = Degrees_To_Radians(Deg);
	float cos = cosf(rad);
	float sin = sinf(rad);

	m._e22 = cos;
	m._e23 = -sin;
	m._e32 = sin;
	m._e33 = cos;

	return m;
}

// Get a rotation matrix for rotation about the y-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Y-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Y(float Deg)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e22 = 1;
	m._e44 = 1;

	float rad = Degrees_To_Radians(Deg);
	float cos = cosf(rad);
	float sin = sinf(rad);

	m._e11 = cos;
	m._e13 = sin;
	m._e31 = -sin;
	m._e33 = cos;

	return m;
}

// Get a rotation matrix for rotation about the z-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Z-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Z(float Deg)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, };
	m._e33 = 1;
	m._e44 = 1;

	float rad = Degrees_To_Radians(Deg);
	float cos = cosf(rad);
	float sin = sinf(rad);

	m._e11 = cos;
	m._e12 = -sin;
	m._e21 = sin;
	m._e22 = cos;
	return m;
}

// ADD two matrices together
//
// IN:		m		The first matrix
//			n		The second matrix
//
// RETURN: m + n
TMATRIX Matrix_Matrix_Add(TMATRIX m, TMATRIX n)
{
	// TODO LAB 2: Replace with your implementation.	
	m._e11 = m._e11 + n._e11;
	m._e12 = m._e12 + n._e12;
	m._e13 = m._e13 + n._e13;
	m._e14 = m._e14 + n._e14;

	m._e21 = m._e21 + n._e21;
	m._e22 = m._e22 + n._e22;
	m._e23 = m._e23 + n._e23;
	m._e24 = m._e24 + n._e24;

	m._e31 = m._e31 + n._e31;
	m._e32 = m._e32 + n._e32;
	m._e33 = m._e33 + n._e33;
	m._e34 = m._e34 + n._e34;

	m._e41 = m._e41 + n._e41;
	m._e42 = m._e42 + n._e42;
	m._e43 = m._e43 + n._e43;
	m._e44 = m._e44 + n._e44;
	return m;
}

// SUBTRACT two matrices
//
// IN:		m		The first matrix (left hand side)
//			n		The second matrix (right hand side)
//
// RETURN: m - n
TMATRIX Matrix_Matrix_Sub(TMATRIX m, TMATRIX n)
{
	// TODO LAB 2: Replace with your implementation.
	m._e11 = m._e11 - n._e11;
	m._e12 = m._e12 - n._e12;
	m._e13 = m._e13 - n._e13;
	m._e14 = m._e14 - n._e14;
					
	m._e21 = m._e21 - n._e21;
	m._e22 = m._e22 - n._e22;
	m._e23 = m._e23 - n._e23;
	m._e24 = m._e24 - n._e24;
					
	m._e31 = m._e31 - n._e31;
	m._e32 = m._e32 - n._e32;
	m._e33 = m._e33 - n._e33;
	m._e34 = m._e34 - n._e34;
					
	m._e41 = m._e41 - n._e41;
	m._e42 = m._e42 - n._e42;
	m._e43 = m._e43 - n._e43;
	m._e44 = m._e44 - n._e44;


	return m;
}

// Multiply a matrix by a scalar
//
// IN:		m		The matrix to be scaled (right hand side)
//			s		The value to scale by   (left hand side)
//
// RETURN:	The matrix formed by s*[m]
TMATRIX Matrix_Scalar_Multiply(TMATRIX m, float s)
{
	// TODO LAB 2: Replace with your implementation.
		m._e11 = m._e11 * s;
		m._e12 = m._e12	* s;	 
		m._e13 = m._e13	* s; 
		m._e14 = m._e14	* s; 
				  
		m._e21	= m._e21 * s;
		m._e22 = m._e22	* s;
		m._e23 = m._e23	* s;
		m._e24 = m._e24	* s;
				 
		m._e31	 =	m._e31 * s;
		m._e32 = m._e32 * s;
		m._e33 = m._e33 * s;
		m._e34 = m._e34 * s;
				  

		m._e41 = m._e41	* s;
		m._e42 = m._e42	* s;
		m._e43 = m._e43	* s;
		m._e44 = m._e44	* s;

	return m;	  
}

// Negate a matrix
//
// IN:		m		The matrix to negate
//
// RETURN:  The negation of m
TMATRIX Matrix_Negate(TMATRIX m)
{
	// TODO LAB 2: Replace with your implementation.
			m._e11 = - 1 * m._e11;
			m._e12 = -1 *m._e12;
			m._e13 = -1 *m._e13;
			m._e14 = -1 *m._e14;
					
			m._e21 = -1 *m._e21;
			m._e22 = -1 *m._e22;
			m._e23 = -1 *m._e23;
			m._e24 = -1 *m._e24;
					
			m._e31 = -1 *m._e31;
			m._e32 = -1 *m._e32;
			m._e33 = -1 *m._e33;
			m._e34 = -1 *m._e34;
					
			m._e41 = -1 *m._e41;
			m._e42 = -1 *m._e42;
			m._e43 = -1 *m._e43;
			m._e44 = -1 *m._e44;

	return m;
}

// Transpose a matrix
//
// IN:		m		The matrix to transpose
//
// RETURN:	The transpose of m
TMATRIX Matrix_Transpose(TMATRIX m)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m2 = m;
	m._e11 = m2._e11;
	m._e12 = m2._e21;
	m._e13 = m2._e31;
	m._e14 = m2._e41;
			 
	m._e21 = m2._e12;
	m._e22 = m2._e22;
	m._e23 = m2._e32;
	m._e24 = m2._e42;
			 
	m._e31 = m2._e13;
	m._e32 = m2._e23;
	m._e33 = m2._e33;
	m._e34 = m2._e43;
			 
	m._e41 = m2._e14;
	m._e42 = m2._e24;
	m._e43 = m2._e34;
	m._e44 = m2._e44;

	return m;
}

// Multipy a matrix and a vector
//
// IN:		m		The matrix (left hand side)
//			v		The vector (right hand side)
//
// RETURN:	[m]*v
TVECTOR Matrix_Vector_Multiply(TMATRIX m, TVECTOR v)
{
	// TODO LAB 2: Replace with your implementation.
	TVECTOR vec = v;

	vec.x = m._e11 * v.x + m._e12 * v.y + m._e13 * v.z + m._e14 * v.w;
	vec.y = m._e21 * v.x + m._e22 * v.y + m._e23 * v.z + m._e24 * v.w;
	vec.z = m._e31 * v.x + m._e32 * v.y + m._e33 * v.z + m._e34 * v.w;
	vec.w = m._e41 * v.x + m._e42 * v.y + m._e43 * v.z + m._e44 * v.w;

	v = vec;

	return v;
}

// Multipy a vector and a matrix
//
// IN:		v		The vector ( left hand side)
//			m		The matrix (right hand side)
//
// RETURN:	v*[m]
TVECTOR Vector_Matrix_Multiply(TVECTOR v, TMATRIX m)
{
	// TODO LAB 2: Replace with your implementation.
	TVECTOR matrix;
	m = Matrix_Transpose(m);
	matrix.x = v.x * m._e11 + v.y * m._e12 + v.z * m._e13 + v.w * m._e14;

	matrix.y = v.x * m._e21 + v.y * m._e22 + v.z * m._e23 + v.w * m._e24;

	matrix.z = v.x * m._e31 + v.y * m._e32 + v.z * m._e33 + v.w * m._e34;

	matrix.w = v.x * m._e41 + v.y * m._e42 + v.z * m._e43 + v.w * m._e44;

	v = matrix;

	return v;
}
// Multiply a matrix by a matrix
//
// IN:		m		First Matrix (left hand side)
//			n		Second Matrix (right hand side)
//
// RETURN:	[m]*[n]
TMATRIX Matrix_Matrix_Multiply(TMATRIX m, TMATRIX n)
{
	TMATRIX matrix = m;
	matrix._e11 = m._e11*n._e11 + m._e12*n._e21 + m._e13*n._e31 + m._e14*n._e41;
	matrix._e12 = m._e11*n._e12 + m._e12*n._e22 + m._e13*n._e32 + m._e14*n._e42;
	matrix._e13 = m._e11*n._e13 + m._e12*n._e23 + m._e13*n._e33 + m._e14*n._e43;
	matrix._e14 = m._e11*n._e14 + m._e12*n._e24 + m._e13*n._e34 + m._e14*n._e44;

	matrix._e21 = m._e21*n._e11 + m._e22*n._e21 + m._e23*n._e31 + m._e24*n._e41;
	matrix._e22 = m._e21*n._e12 + m._e22*n._e22 + m._e23*n._e32 + m._e24*n._e42;
	matrix._e23 = m._e21*n._e13 + m._e22*n._e23 + m._e23*n._e33 + m._e24*n._e43;
	matrix._e24 = m._e21*n._e14 + m._e22*n._e24 + m._e23*n._e34 + m._e24*n._e44;

	matrix._e31 = m._e31*n._e11 + m._e32*n._e21 + m._e33*n._e31 + m._e34*n._e41;
	matrix._e32 = m._e31*n._e12 + m._e32*n._e22 + m._e33*n._e32 + m._e34*n._e42;
	matrix._e33 = m._e31*n._e13 + m._e32*n._e23 + m._e33*n._e33 + m._e34*n._e43;
	matrix._e34 = m._e31*n._e14 + m._e32*n._e24 + m._e33*n._e34 + m._e34*n._e44;

	matrix._e41 = m._e41*n._e11 + m._e42*n._e21 + m._e43*n._e31 + m._e44*n._e41;
	matrix._e42 = m._e41*n._e12 + m._e42*n._e22 + m._e43*n._e32 + m._e44*n._e42;
	matrix._e43 = m._e41*n._e13 + m._e42*n._e23 + m._e43*n._e33 + m._e44*n._e43;
	matrix._e44 = m._e41*n._e14 + m._e42*n._e24 + m._e43*n._e34 + m._e44*n._e44;

	m = matrix;
	return m;

	// TODO LAB 2: Replace with your implementation.
	
}

////////////////////////////////////////////////////////////////////////
// Matrix Functions Lab # 3
///////////////////////////////////////////////////////////////////////

// HELPER FUNCTION  *** NOT GRADED, ONLY SUGGESTED ***
// USE THIS FUNCTION TO FIND THE DETERMINANT OF A 3*3
// MATRIX. IT CAN BE USED IN THE MATRIX DETERMINANT
// AND MATRIX INVERSE FUNCTIONS BELOW
// 
// RETURN:	The determinant of a 3x3 matrix
float Matrix_Determinant(float e_11,float e_12,float e_13,
						 float e_21,float e_22,float e_23,
						 float e_31,float e_32,float e_33)
{
	float determinant = 0;
	float dA = 0, dB = 0, dC = 0;
	dA = e_11 * (e_22 * e_33 - e_23 * e_32);
	dB = e_12 * (e_21 * e_33 - e_23 * e_31);
	dC = e_13 * (e_21 * e_32 - e_22 * e_31);
	determinant = dA - dB + dC;

	return determinant;
}

// Get the determinant of a matrix
//
// IN:		m		The ONE!
//
// RETURN:	It's deterinant
float Matrix_Determinant(TMATRIX m)
{
	// TODO LAB 3: Replace with your implementation.'
	float determinant = 0;
	float dA = 0, dB = 0, dC = 0, dD = 0;
	dA = m._e11 * Matrix_Determinant(m._e22, m._e23, m._e24,
									 m._e32, m._e33, m._e34,
									 m._e42, m._e43, m._e44);

	dB = m._e12 * Matrix_Determinant(m._e21, m._e23, m._e24,
									 m._e31, m._e33, m._e34,
									 m._e41, m._e43, m._e44);

	dC = m._e13 * Matrix_Determinant(m._e21, m._e22, m._e24,
									 m._e31, m._e32, m._e34,
									 m._e41, m._e42, m._e44);

	dD = m._e14 * Matrix_Determinant(m._e21, m._e22, m._e23,
									 m._e31, m._e32, m._e33,
									 m._e41, m._e42, m._e43);

	determinant = dA - dB + dC - dD;

	return determinant;
}

// Get the inverse of a matrix
//
// IN:		m		The matrix to inverse
//
// RETURN:	The Inverse of [m]
//
// NOTE: Returns the matrix itself if m is not invertable.
TMATRIX Matrix_Inverse(TMATRIX m)
{
	// TODO LAB 3: Replace with your implementation.
	TMATRIX matrix;
	float determinant = Matrix_Determinant(m);
	if (IsZero(determinant))
	{
		return m;
	}
	else
	{
		matrix._e11 = 1 * Matrix_Determinant(m._e22, m._e23, m._e24,
											 m._e32, m._e33, m._e34,
											 m._e42, m._e43, m._e44);
		
		matrix._e12 = -1 * Matrix_Determinant(m._e21, m._e23, m._e24,
											  m._e31, m._e33, m._e34,
											  m._e41, m._e43, m._e44);

		matrix._e13 = 1 * Matrix_Determinant(m._e21, m._e22, m._e24,
											 m._e31, m._e32, m._e34,
											 m._e41, m._e42, m._e44);

		matrix._e14 = -1 * Matrix_Determinant(m._e21, m._e22, m._e23,
											  m._e31, m._e32, m._e33,
											  m._e41, m._e42, m._e43);

		matrix._e21 = -1 * Matrix_Determinant(m._e12, m._e13, m._e14,
											  m._e32, m._e33, m._e34,
											  m._e42, m._e43, m._e44);

		matrix._e22 = 1 * Matrix_Determinant(m._e11, m._e13, m._e14,
											 m._e31, m._e33, m._e34,
											 m._e41, m._e43, m._e44);

		matrix._e23 = -1 * Matrix_Determinant(m._e11, m._e12, m._e14,
											  m._e31, m._e32, m._e34,
											  m._e41, m._e42, m._e44);

		matrix._e24 =  1 * Matrix_Determinant(m._e11, m._e12, m._e13,
											  m._e31, m._e32, m._e33,
											  m._e41, m._e42, m._e43);

		matrix._e31 = 1 * Matrix_Determinant(m._e12, m._e13, m._e14,
											 m._e22, m._e23, m._e24,
											 m._e42, m._e43, m._e44);

		matrix._e32 = -1 * Matrix_Determinant(m._e11, m._e13, m._e14,
											  m._e21, m._e23, m._e24,
											  m._e41, m._e43, m._e44);

		matrix._e33 = 1 * Matrix_Determinant(m._e11, m._e12, m._e14,
											 m._e21, m._e22, m._e24,
											 m._e41, m._e42, m._e44);

		matrix._e34 = -1 * Matrix_Determinant(m._e11, m._e12, m._e13,
											  m._e21, m._e22, m._e23,
										      m._e41, m._e42, m._e43);

		matrix._e41 = -1 * Matrix_Determinant(m._e12, m._e13, m._e14,
											  m._e22, m._e23, m._e24,
											  m._e32, m._e33, m._e34);

		matrix._e42 = 1 * Matrix_Determinant(m._e11, m._e13, m._e14,
										     m._e21, m._e23, m._e24,
									         m._e31, m._e33, m._e34);

		matrix._e43 = -1 * Matrix_Determinant(m._e11, m._e12, m._e14,
											  m._e21, m._e22, m._e24,
											  m._e31, m._e32, m._e34);

		matrix._e44 = 1 * Matrix_Determinant(m._e11, m._e12, m._e13,
											 m._e21, m._e22, m._e23,
											 m._e31, m._e32, m._e33);

		m._e11 = matrix._e11 / determinant;
		m._e12 = matrix._e21 / determinant;
		m._e13 =  matrix._e31 /determinant;
		m._e14 =  matrix._e41 /determinant;
				  	 
		m._e21 =  matrix._e12 /determinant;
		m._e22 =  matrix._e22 /determinant;
		m._e23 =  matrix._e32 /determinant;
		m._e24 =  matrix._e42 /determinant;
				
		m._e31 =  matrix._e13 /determinant;
		m._e32 =  matrix._e23 /determinant;
		m._e33 =  matrix._e33 /determinant;
		m._e34 =  matrix._e43 /determinant;
				
		m._e41 =  matrix._e14 /determinant;
		m._e42 =  matrix._e24 /determinant;
		m._e43 =  matrix._e34 /determinant;
		m._e44 =  matrix._e44 /determinant;

		return m;
	}
}

TMATRIX Matrix_Projection_M(TMATRIX m, float widthS, float heightS, float nZ, float fZ)
{
	float aspect = 1;
	if (widthS > heightS)
	{
		aspect = widthS / heightS;
	}
	else if (heightS > widthS)
	{
		aspect = heightS / widthS;
	}


	m._e22 = 1 / tan(Degrees_To_Radians(90 * 0.5f));

	m._e11 = m._e22 * aspect;

	m._e33 = fZ / (fZ - nZ);

	m._e34 = 1;

	m._e43 = -(fZ * nZ) / (fZ - nZ);

	return m;
}
