// CQuaternion.h: interface for the CQuaternion class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CQuaternion_H__4408006D_97B6_408E_B02D_C04FC8B65274__INCLUDED_)
#define AFX_CQuaternion_H__4408006D_97B6_408E_B02D_C04FC8B65274__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Matrix3X3.h"

class CQuaternion  
{
public:
	inline CQuaternion (
		double fW = 1.0,
		double fX = 0.0, double fY = 0.0, double fZ = 0.0)
	{
		m_w = fW;
		m_x = fX;
		m_y = fY;
		m_z = fZ;
	}
	inline CQuaternion (const CQuaternion& rkQ)
	{
		m_w = rkQ.m_w;
		m_x = rkQ.m_x;
		m_y = rkQ.m_y;
		m_z = rkQ.m_z;
	}

	/// construct a quaternion from euler angle
	void FromEulerAngle(double pitch, double yaw, double roll);
	bool ToEulerAngleXZY(double &yAngle, double &xAngle, double &zAngle);

	bool ToEulerAnglesZXY(double &zAngle, double &xAngle, double &yAngle) const;

	/// Construct a quaternion from a rotation matrix
	inline CQuaternion(const CMatrix3x3& rot)
	{
		this->FromRotationMatrix(rot);
	}
	/// Construct a quaternion from an angle/axis
	inline CQuaternion(const double& rfAngle, const CVector3D& rkAxis)
	{
		this->FromAngleAxis(rfAngle, rkAxis);
	}

	/// Construct a quaternion from 3 orthonormal local axes
	inline CQuaternion(const CVector3D& xaxis, const CVector3D& yaxis, const CVector3D& zaxis)
	{
		this->FromAxes(xaxis, yaxis, zaxis);
	}
	/// Construct a quaternion from 3 orthonormal local axes
	inline CQuaternion(const CVector3D* akAxis)
	{
		this->FromAxes(akAxis);
	}
	/// Construct a quaternion from 4 manual w/x/y/z values
	inline CQuaternion(double* valptr)
	{
		memcpy(&m_w, valptr, sizeof(double)*4);
	}
	
	/// Array accessor operator
	inline double operator [] ( const size_t i ) const
	{
		assert( i < 4 );
		
		return *(&m_w+i);
	}
	
	/// Array accessor operator
	inline double& operator [] ( const size_t i )
	{
		assert( i < 4 );
		
		return *(&m_w+i);
	}
	
	/// Pointer accessor for direct copying
	inline double* ptr()
	{
		return &m_w;
	}
	
	/// Pointer accessor for direct copying
	inline const double* ptr() const
	{
		return &m_w;
	}
	
	void FromRotationMatrix (const CMatrix3x3& kRot);
	void ToRotationMatrix (CMatrix3x3& kRot) const;
	void FromAngleAxis (const double& rfAngle, const CVector3D& rkAxis);
	void ToAngleAxis (double& rfAngle, CVector3D& rkAxis) const;
	
	void FromAxes (const CVector3D* akAxis);
	void FromAxes (const CVector3D& xAxis, const CVector3D& yAxis, const CVector3D& zAxis);
	void ToAxes (CVector3D* akAxis) const;
	void ToAxes (CVector3D& xAxis, CVector3D& yAxis, CVector3D& zAxis) const;
	/// Get the local x-axis
	CVector3D xAxis(void) const;
	/// Get the local y-axis
	CVector3D yAxis(void) const;
	/// Get the local z-axis
	CVector3D zAxis(void) const;

	// ---------------- get the quaternion which convert from fromAxis to toAxis 
	void FromTwoAxes(const CVector3D &fromAxis, const CVector3D &toAxis);
	
	inline CQuaternion& operator= (const CQuaternion& rkQ)
	{
		m_w = rkQ.m_w;
		m_x = rkQ.m_x;
		m_y = rkQ.m_y;
		m_z = rkQ.m_z;
		return *this;
	}
	
    void operator*= (double fScalar);

	CQuaternion operator+ (const CQuaternion& rkQ) const;
	CQuaternion operator- (const CQuaternion& rkQ) const;
	CQuaternion operator* (const CQuaternion& rkQ) const;
	CQuaternion operator* (double fScalar) const;
	friend CQuaternion operator* (double fScalar,
		const CQuaternion& rkQ);
	CQuaternion operator- () const;
	inline bool operator== (const CQuaternion& rhs) const
	{
		return (rhs.m_x == m_x) && (rhs.m_y == m_y) &&
			(rhs.m_z == m_z) && (rhs.m_w == m_w);
	}
	inline bool operator!= (const CQuaternion& rhs) const
	{
		return !operator==(rhs);
	}

	
	// functions of a quaternion
	double Dot (const CQuaternion& rkQ) const;  // dot product
	double Norm () const;  // squared-length
	/// Normalises this quaternion, and returns the previous length
	double Normalise(void); 
	CQuaternion Inverse () const;  // apply to non-zero quaternion
	CQuaternion UnitInverse () const;  // apply to unit-length quaternion
	CQuaternion Exp () const;
	CQuaternion Log () const;
	CQuaternion Conjugate()const;

	CQuaternion ClosestTo(const CQuaternion &quat) const
	{
		double cos = Dot(quat);
		if(cos < 0.0)
			return -(*this);

		return *this;			
	}
	
	// rotation of a vector by a quaternion
	CVector3D operator* (const CVector3D& rkVector) const;
	
	/** Calculate the local roll element of this quaternion.
	@param reprojectAxis By default the method returns the 'intuitive' result
	that is, if you projected the local Y of the quaterion onto the X and
	Y axes, the angle between them is returned. If set to false though, the
	result is the actual yaw that will be used to implement the quaternion,
	which is the shortest possible path to get to the same orientation and 
	may involve less axial rotation. 
	*/
	double getRoll(bool reprojectAxis = true) const;
	/** Calculate the local pitch element of this quaternion
	@param reprojectAxis By default the method returns the 'intuitive' result
	that is, if you projected the local Z of the quaterion onto the X and
	Y axes, the angle between them is returned. If set to true though, the
	result is the actual yaw that will be used to implement the quaternion,
	which is the shortest possible path to get to the same orientation and 
	may involve less axial rotation. 
	*/
	double getPitch(bool reprojectAxis = true) const;
	/** Calculate the local yaw element of this quaternion
	@param reprojectAxis By default the method returns the 'intuitive' result
	that is, if you projected the local Z of the quaterion onto the X and
	Z axes, the angle between them is returned. If set to true though, the
	result is the actual yaw that will be used to implement the quaternion,
	which is the shortest possible path to get to the same orientation and 
	may involve less axial rotation. 
	*/
	double getYaw(bool reprojectAxis = true) const;		
	/// Equality with tolerance (tolerance is max angle difference)
	bool equals(const CQuaternion& rhs, const double& tolerance) const;
	
	// spherical linear interpolation
	static CQuaternion Slerp (double fT, const CQuaternion& rkP,
		const CQuaternion& rkQ, bool shortestPath = false);

	//notice that the input quaternion may be changed
    static CQuaternion CubicSlerp(double fT, CQuaternion &q0, CQuaternion &q1,
									  CQuaternion &q2, CQuaternion &q3);
	
	static CQuaternion SlerpExtraSpins (double fT,
		const CQuaternion& rkP, const CQuaternion& rkQ,
		int iExtraSpins);
	
	// setup for spherical quadratic interpolation
	static void Intermediate (const CQuaternion& rkQ0,
		const CQuaternion& rkQ1, const CQuaternion& rkQ2,
		CQuaternion& rka, CQuaternion& rkB);
	
	// spherical quadratic interpolation
	static CQuaternion Squad (double fT, const CQuaternion& rkP,
		const CQuaternion& rkA, const CQuaternion& rkB,
		const CQuaternion& rkQ, bool shortestPath = false);
	
	// normalised linear interpolation - faster but less accurate (non-constant rotation velocity)
	static CQuaternion nlerp(double fT, const CQuaternion& rkP, 
		const CQuaternion& rkQ, bool shortestPath = false);
	
	// cutoff for sine near zero
	static const double ms_fEpsilon;
	
	// special values
	static const CQuaternion ZERO;
	static const CQuaternion IDENTITY;
	
	double m_w, m_x, m_y, m_z;
	
	/** Function for writing to a stream. Outputs "CQuaternion(w, x, y, z)" with w,x,y,z
	being the member values of the quaternion.
	*/
	inline friend std::ostream& operator <<
		( std::ostream& o, const CQuaternion& q )
	{
		o << "CQuaternion(" << q.m_w << ", " << q.m_x << ", " << q.m_y << ", " << q.m_z << ")";
		return o;
	}

	friend double operator % (const CQuaternion &quat1, const CQuaternion &quat2)
	{
		return quat1.Dot(quat2);
	}

	friend CQuaternion DoubleQuater(const CQuaternion &p, const CQuaternion &q)
	{
		return 2.0 * (p.Dot(q)) * q - p;
	}

	friend CQuaternion Bisect(const CQuaternion &quat1, const CQuaternion &quat2)
	{
		CQuaternion quat = quat1 + quat2;
		quat.Normalise();

		return quat;
	}
	
};

#endif // !defined(AFX_CQuaternion_H__4408006D_97B6_408E_B02D_C04FC8B65274__INCLUDED_)
