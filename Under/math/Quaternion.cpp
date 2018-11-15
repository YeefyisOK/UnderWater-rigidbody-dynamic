// ReQuaternion.cpp: implementation of the ReQuaternion class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "Quaternion.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const double CQuaternion::ms_fEpsilon = 1e-04;
const CQuaternion CQuaternion::ZERO(0.0,0.0,0.0,0.0);
const CQuaternion CQuaternion::IDENTITY(1.0,0.0,0.0,0.0);
//      |                                                                         |
//      | w^2 + x^2 - y^2 - z^2     2 * x * y - 2 * w * z   2 * x * z + 2 * w * y |
//      |                                                                         |
//  R = | 2 * x * y + 2 * w * z     w^2 - x^2 + y^2 - z^2   2 * y * z - 2 * w * x |
//      |                                                                         |
//      | 2 * x * z - 2 * w * y     2 * y * z + 2 * w * x   w^2 - x^2 - y^2 + z^2 |

//-----------------------------------------------------------------------
void CQuaternion::FromRotationMatrix (const CMatrix3x3& kRot)
{
	// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
	// article "ReQuaternion Calculus and Fast Animation".
	
	double fTrace = kRot[0][0] + kRot[1][1] + kRot[2][2];
	double fRoot;
	
	if ( fTrace > 0.0 )
	{
		// |w| > 1/2, may as well choose w > 1/2
		fRoot = sqrt(fTrace + 1.0);  // 2w
		m_w = 0.5*fRoot;
		fRoot = 0.5/fRoot;  // 1/(4w)
		m_x = (kRot[2][1] - kRot[1][2]) * fRoot;
		m_y = (kRot[0][2] - kRot[2][0])*fRoot;
		m_z = (kRot[1][0] - kRot[0][1])*fRoot;
	}
	else
	{
		// |w| <= 1/2
		static size_t s_iNext[3] = { 1, 2, 0 };
		size_t i = 0;
		if ( kRot[1][1] > kRot[0][0] )
			i = 1;
		if ( kRot[2][2] > kRot[i][i] )
			i = 2;
		size_t j = s_iNext[i];
		size_t k = s_iNext[j];
		
		fRoot = sqrt(kRot[i][i]-kRot[j][j]-kRot[k][k] + 1.0);
		double* apkQuat[3] = { &m_x, &m_y, &m_z };
		*apkQuat[i] = 0.5*fRoot;
		fRoot = 0.5/fRoot;
		m_w = (kRot[k][j] - kRot[j][k]) * fRoot;
		*apkQuat[j] = (kRot[j][i] + kRot[i][j]) * fRoot;
		*apkQuat[k] = (kRot[k][i] + kRot[i][k]) * fRoot;
	}
}

// the quaternion should be normalized
//-----------------------------------------------------------------------
void CQuaternion::ToRotationMatrix (CMatrix3x3& kRot) const
{
	double fTx  = m_x + m_x;
	double fTy  = m_y + m_y;
	double fTz  = m_z + m_z;
	double fTwx = fTx * m_w;
	double fTwy = fTy * m_w;
	double fTwz = fTz * m_w;
	double fTxx = fTx * m_x;
	double fTxy = fTy * m_x;
	double fTxz = fTz * m_x;
	double fTyy = fTy * m_y;
	double fTyz = fTz * m_y;
	double fTzz = fTz * m_z;

	kRot[0][0] = 1.0f - (fTyy + fTzz);
	kRot[0][1] = fTxy - fTwz;
	kRot[0][2] = fTxz + fTwy;
	kRot[1][0] = fTxy + fTwz;
	kRot[1][1] = 1.0f - (fTxx + fTzz);
	kRot[1][2] = fTyz - fTwx;
	kRot[2][0] = fTxz - fTwy;
	kRot[2][1] = fTyz + fTwx;
	kRot[2][2] = 1.0f - (fTxx + fTyy);
}

//-----------------------------------------------------------------------
//angle is in radian 
void CQuaternion::FromAngleAxis (const double& rfAngle,
								const CVector3D& rkAxis)
{
	// assert:  axis[] is unit length
	//
	// The quaternion representing the rotation is
	//   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
	
	double fHalfAngle ( 0.5*rfAngle );
	double fSin = sin(fHalfAngle);
	m_w = cos(fHalfAngle);

	m_x = fSin * rkAxis.X();
	m_y = fSin * rkAxis.Y();
	m_z = fSin * rkAxis.Z();
}

//-----------------------------------------------------------------------
void CQuaternion::ToAngleAxis (double& rfAngle, CVector3D& rkAxis) const
{
	// The quaternion representing the rotation is
	//   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
	
	double fSqrLength = m_x*m_x+m_y*m_y+m_z*m_z;
	if ( fSqrLength > 0.0 )
	{
		rfAngle = 2.0 * acos(m_w);
		double fInvLength = 1.0 / sqrt(fSqrLength);
		rkAxis.Set(m_x*fInvLength, m_y*fInvLength, m_z*fInvLength);
	}
	else
	{
		// angle is 0 (mod 2*pi), so any axis will do
		rfAngle = 0.0;
		rkAxis.Set(1.0, 0.0, 0.0);
	}
}

//-----------------------------------------------------------------------
//from three orthogonal axes, each axes form a column
void CQuaternion::FromAxes (const CVector3D* akAxis)
{
	CMatrix3x3 kRot;
	
	for (size_t iCol = 0; iCol < 3; iCol++)
	{
		kRot[0][iCol] = akAxis[iCol].X();
		kRot[1][iCol] = akAxis[iCol].Y();
		kRot[2][iCol] = akAxis[iCol].Z();
	}
	
	FromRotationMatrix(kRot);
}
//-----------------------------------------------------------------------
void CQuaternion::FromAxes (const CVector3D& xaxis, const CVector3D& yaxis, const CVector3D& zaxis)
{
	CMatrix3x3 kRot;
	
	kRot[0][0] = xaxis.X();
	kRot[1][0] = xaxis.Y();
	kRot[2][0] = xaxis.Z();
	
	kRot[0][1] = yaxis.X();
	kRot[1][1] = yaxis.Y();
	kRot[2][1] = yaxis.Z();
	
	kRot[0][2] = zaxis.X();
	kRot[1][2] = zaxis.Y();
	kRot[2][2] = zaxis.Z();
	
	FromRotationMatrix(kRot);
	
}
//-----------------------------------------------------------------------
void CQuaternion::ToAxes (CVector3D* akAxis) const
{
	CMatrix3x3 kRot;
	
	ToRotationMatrix(kRot);
	
	for (size_t iCol = 0; iCol < 3; iCol++)
	{
		akAxis[iCol].Set(kRot[0][iCol], kRot[1][iCol], kRot[2][iCol]);
	}
}
//-----------------------------------------------------------------------
CVector3D CQuaternion::xAxis(void) const
{
	//double fTx  = 2.0*x;
	double fTy  = 2.0*m_y;
	double fTz  = 2.0*m_z;
	double fTwy = fTy*m_w;
	double fTwz = fTz*m_w;
	double fTxy = fTy*m_x;
	double fTxz = fTz*m_x;
	double fTyy = fTy*m_y;
	double fTzz = fTz*m_z;
	
	return CVector3D(1.0-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
}
//-----------------------------------------------------------------------
CVector3D CQuaternion::yAxis(void) const
{
	double fTx  = 2.0*m_x;
	double fTy  = 2.0*m_y;
	double fTz  = 2.0*m_z;
	double fTwx = fTx*m_w;
	double fTwz = fTz*m_w;
	double fTxx = fTx*m_x;
	double fTxy = fTy*m_x;
	double fTyz = fTz*m_y;
	double fTzz = fTz*m_z;
	
	return CVector3D(fTxy-fTwz, 1.0-(fTxx+fTzz), fTyz+fTwx);
}
//-----------------------------------------------------------------------
CVector3D CQuaternion::zAxis(void) const
{
	double fTx  = 2.0*m_x;
	double fTy  = 2.0*m_y;
	double fTz  = 2.0*m_z;
	double fTwx = fTx*m_w;
	double fTwy = fTy*m_w;
	double fTxx = fTx*m_x;
	double fTxz = fTz*m_x;
	double fTyy = fTy*m_y;
	double fTyz = fTz*m_y;
	
	return CVector3D(fTxz+fTwy, fTyz-fTwx, 1.0-(fTxx+fTyy));
}

//-----------------------------------------------------------------------
void CQuaternion::ToAxes (CVector3D& xaxis, CVector3D& yaxis, CVector3D& zaxis) const
{
	CMatrix3x3 kRot;
	
	ToRotationMatrix(kRot);
	
	xaxis.Set(kRot[0][0], kRot[1][0], kRot[2][0]);	
	yaxis.Set(kRot[0][1], kRot[1][1], kRot[2][1]);	
	zaxis.Set(kRot[0][2], kRot[1][2], kRot[2][2]);
}

// ---------------- get the quaternion which convert from fromAxis to toAxis 
void CQuaternion::FromTwoAxes(const CVector3D &fromAxis, const CVector3D &toAxis)
{
	CVector3D from = fromAxis;
	CVector3D to = toAxis;
	from.Normalize();
	to.Normalize();

	double c = from % to;
	if (fabs(c - 1) < ms_fEpsilon)  //no rotaion
	{
		m_w = 0.0;
		m_x = 0.0;
		m_y = 0.0;
		m_z = 1.0;
	}
	else
	{
		CVector3D rotateAxis = from * to;
		rotateAxis.Normalize();

		double dHalfAngle = acos(c) * 0.5;

		c = cos(dHalfAngle);
		double s = sin(dHalfAngle);

		m_w = c;
		m_x = rotateAxis[0] * s;
		m_y = rotateAxis[1] * s;
		m_z = rotateAxis[2] * s;
	}
	
}

//-----------------------------------------------------------------------
CQuaternion CQuaternion::operator+ (const CQuaternion& rkQ) const
{
	return CQuaternion(m_w+rkQ.m_w,m_x+rkQ.m_x,m_y+rkQ.m_y,m_z+rkQ.m_z);
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::operator- (const CQuaternion& rkQ) const
{
	return CQuaternion(m_w-rkQ.m_w,m_x-rkQ.m_x,m_y-rkQ.m_y,m_z-rkQ.m_z);
}

//Multiplying q1 with q2 applies the rotation q2 to q1
//-----------------------------------------------------------------------
CQuaternion CQuaternion::operator* (const CQuaternion& rkQ) const
{
	// NOTE:  Multiplication is not generally commutative, so in most
	// cases p*q != q*p.
	
	return CQuaternion
        (
		m_w * rkQ.m_w - m_x * rkQ.m_x - m_y * rkQ.m_y - m_z * rkQ.m_z,
		m_w * rkQ.m_x + m_x * rkQ.m_w + m_y * rkQ.m_z - m_z * rkQ.m_y,
		m_w * rkQ.m_y + m_y * rkQ.m_w + m_z * rkQ.m_x - m_x * rkQ.m_z,
		m_w * rkQ.m_z + m_z * rkQ.m_w + m_x * rkQ.m_y - m_y * rkQ.m_x
        );
}
//-----------------------------------------------------------------------
inline CQuaternion CQuaternion::operator* (double fScalar) const
{
	return CQuaternion(fScalar*m_w,fScalar*m_x,fScalar*m_y,fScalar*m_z);
}
//-----------------------------------------------------------------------
inline CQuaternion operator* (double fScalar, const CQuaternion& rkQ)
{
	return CQuaternion(fScalar*rkQ.m_w,fScalar*rkQ.m_x,fScalar*rkQ.m_y,
		fScalar*rkQ.m_z);
}
//-----------------------------------------------------------------------
inline void CQuaternion::operator*= (double fScalar)
{
	m_w *= fScalar;
	m_x *= fScalar;
	m_y *= fScalar;
	m_z *= fScalar;
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::operator- () const
{
	return CQuaternion(-m_w,-m_x,-m_y,-m_z);
}
//-----------------------------------------------------------------------
inline double CQuaternion::Dot (const CQuaternion& rkQ) const
{
	return m_w*rkQ.m_w+m_x*rkQ.m_x+m_y*rkQ.m_y+m_z*rkQ.m_z;
}
//-----------------------------------------------------------------------
inline double CQuaternion::Norm () const
{
	return m_w*m_w+m_x*m_x+m_y*m_y+m_z*m_z;
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::Inverse () const
{
	double fNorm = m_w*m_w+m_x*m_x+m_y*m_y+m_z*m_z;
	if ( fNorm > 0.0 )
	{
		double fInvNorm = 1.0/fNorm;
		return CQuaternion(m_w*fInvNorm,-m_x*fInvNorm,-m_y*fInvNorm,-m_z*fInvNorm);
	}
	else
	{
		// return an invalid result to flag the error
		return ZERO;
	}
}

//The conjugate of a quaternion is the same as the inverse, as long as the quaternion is unit length
inline CQuaternion CQuaternion::Conjugate()const
{
	return CQuaternion(m_w, -m_x, -m_y, -m_z);
}
//-----------------------------------------------------------------------
inline CQuaternion CQuaternion::UnitInverse () const
{
	// assert:  'this' is unit length
	return CQuaternion(m_w,-m_x,-m_y,-m_z);
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::Exp () const
{
	// If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
	// exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
	// use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.
	
	double fAngle = sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
	double fSin = sin(fAngle);
	
	CQuaternion kResult;
	kResult.m_w = cos(fAngle);
	
	if ( fabs(fSin) >= ms_fEpsilon )
	{
		double fCoeff = fSin/fAngle;
		kResult.m_x = fCoeff*m_x;
		kResult.m_y = fCoeff*m_y;
		kResult.m_z = fCoeff*m_z;
	}
	else
	{
		kResult.m_x = m_x;
		kResult.m_y = m_y;
		kResult.m_z = m_z;
	}
	
	return kResult;
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::Log () const
{
	// If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
	// log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
	// sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.
	
	CQuaternion kResult;
	kResult.m_w = 0.0;
	
	if ( fabs(m_w) < 1.0 )
	{
		double fAngle = acos(m_w) ;
		double fSin = sin(fAngle);
		if ( fabs(fSin) >= ms_fEpsilon )
		{
			double fCoeff = fAngle/fSin;
			kResult.m_x = fCoeff*m_x;
			kResult.m_y = fCoeff*m_y;
			kResult.m_z = fCoeff*m_z;
			return kResult;
		}
	}
	
	kResult.m_x = m_x;
	kResult.m_y = m_y;
	kResult.m_z = m_z;
	
	return kResult;
}
//-----------------------------------------------------------------------
CVector3D CQuaternion::operator* (const CVector3D& v) const
{
	// nVidia SDK implementation
	CVector3D uv, uuv;
	CVector3D qvec(m_x, m_y, m_z);
	uv = qvec * v;
	uuv = qvec * uv;
	uv *= (2.0f * m_w);
	uuv *= 2.0f;
	
	return v + uv + uuv;
	
}
//-----------------------------------------------------------------------
bool CQuaternion::equals(const CQuaternion& rhs, const double& tolerance) const
{
	double fCos = Dot(rhs);
	double angle = acos(fCos);
	
	return (fabs(angle) <= tolerance || fabs(angle - PI) < tolerance);
	
	
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::Slerp (double fT, const CQuaternion& rkP,
							  const CQuaternion& rkQ, bool shortestPath)
{
	//                       sin( (1 - u) *theta )        sin( u * theta)
	//   slerp(q1, q2, u) = -----------------------q1 +  -----------------q2,
	//                            sin( theta )              sin( theta )
	//
	//   where cos (theta) = q1 dot q2
	//
	double fCos = rkP.Dot(rkQ);
	CQuaternion rkT;
	
	// Do we need to invert rotation?
	if (fCos < 0.0f && shortestPath)
	{
		fCos = -fCos;
		rkT = -rkQ;
	}
	else
	{
		rkT = rkQ;
	}
	
	if (fabs(fCos) < 1 - ms_fEpsilon)
	{
		// Standard case (slerp)
		double fSin = sqrt(1 - fCos * fCos);
		double fAngle = atan2(fSin, fCos);
		double fInvSin = 1.0f / fSin;
		double fCoeff0 = sin((1.0f - fT) * fAngle) * fInvSin;
		double fCoeff1 = sin(fT * fAngle) * fInvSin;
		return fCoeff0 * rkP + fCoeff1 * rkT;
	}
	else
	{
		// There are two situations:
		// 1. "rkP" and "rkQ" are very close (fCos ~= +1), so we can do a linear
		//    interpolation safely.
		// 2. "rkP" and "rkQ" are almost inverse of each other (fCos ~= -1), there
		//    are an infinite number of possibilities interpolation. but we haven't
		//    have method to fix this case, so just use linear interpolation here.
		CQuaternion t = (1.0f - fT) * rkP + fT * rkT;
		// taking the complement requires renormalisation
		t.Normalise();
		return t;
	}
}

//notice that the input quaternion may be changed
CQuaternion CQuaternion::CubicSlerp(double fT, CQuaternion &q0, CQuaternion &q1,
									  CQuaternion &q2, CQuaternion &q3)
{
	q1 = q1.ClosestTo(q0);
	q2 = q2.ClosestTo(q1);
	q3 = q3.ClosestTo(q2);


    double temp = 1.0 / 3;
	CQuaternion an = CQuaternion::Slerp(temp, q1, Bisect(DoubleQuater(q0,q1),q2));
	CQuaternion bn = CQuaternion::Slerp(temp, q2, Bisect(DoubleQuater(q1,q2),q3));
	CQuaternion bn1 = DoubleQuater(bn,q2);

	an = an.ClosestTo(q1);
	bn1 = bn1.ClosestTo(an);
	q2 = q2.ClosestTo(bn1);

	//q1, an, bn1, q2
	CQuaternion p0 = CQuaternion::Slerp(fT, q1, an, true);
	CQuaternion p1 = CQuaternion::Slerp(fT, an, bn1, true);
	CQuaternion p2 = CQuaternion::Slerp(fT, bn1, q2, true);

	CQuaternion p01 = CQuaternion::Slerp(fT, p0, p1, true);
	CQuaternion p12 = CQuaternion::Slerp(fT, p1, p2, true);

	return CQuaternion::Slerp(fT, p01, p12, true);
}

//-----------------------------------------------------------------------
CQuaternion CQuaternion::SlerpExtraSpins (double fT,
										const CQuaternion& rkP, const CQuaternion& rkQ, int iExtraSpins)
{
	double fCos = rkP.Dot(rkQ);
	double fAngle = acos(fCos);
	
	if ( fabs(fAngle) < ms_fEpsilon )
		return rkP;
	
	double fSin = sin(fAngle);
	double fPhase = PI*iExtraSpins*fT;
	double fInvSin = 1.0/fSin;
	double fCoeff0 = sin((1.0-fT)*fAngle - fPhase)*fInvSin;
	double fCoeff1 = sin(fT*fAngle + fPhase)*fInvSin;
	return fCoeff0*rkP + fCoeff1*rkQ;
}
//-----------------------------------------------------------------------
void CQuaternion::Intermediate (const CQuaternion& rkQ0,
							   const CQuaternion& rkQ1, const CQuaternion& rkQ2,
							   CQuaternion& rkA, CQuaternion& rkB)
{
	// assert:  q0, q1, q2 are unit quaternions
	
	CQuaternion kQ0inv = rkQ0.UnitInverse();
	CQuaternion kQ1inv = rkQ1.UnitInverse();
	CQuaternion rkP0 = kQ0inv*rkQ1;
	CQuaternion rkP1 = kQ1inv*rkQ2;
	CQuaternion kArg = 0.25*(rkP0.Log()-rkP1.Log());
	CQuaternion kMinusArg = -kArg;
	
	rkA = rkQ1*kArg.Exp();
	rkB = rkQ1*kMinusArg.Exp();
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::Squad (double fT,
							  const CQuaternion& rkP, const CQuaternion& rkA,
							  const CQuaternion& rkB, const CQuaternion& rkQ, bool shortestPath)
{
	double fSlerpT = 2.0*fT*(1.0-fT);
	CQuaternion kSlerpP = Slerp(fT, rkP, rkQ, shortestPath);
	CQuaternion kSlerpQ = Slerp(fT, rkA, rkB);
	return Slerp(fSlerpT, kSlerpP ,kSlerpQ);
}
//-----------------------------------------------------------------------
double CQuaternion::Normalise(void)
{
	double len = Norm();
	double factor = 1.0f / sqrt(len);
	(*this) *= factor;
	return len;
}

//the angle is in radian
void CQuaternion::FromEulerAngle(double pitch, double yaw, double roll)
{
	//Basically we create 3 Quaternions, one for pitch, one for yaw, one for roll
	//and multiply those together
	//the calculation below does the same, just shorter

	double p = pitch * 0.5;
	double sinp = sin(p);
	double cosp = cos(p);

	double y = yaw * 0.5;
	double siny = sin(y);
	double cosy = cos(y);

	double r = roll * 0.5;
	double sinr = sin(r);
	double cosr = cos(r);

	m_x = sinr * cosp * cosy - cosr * sinp * siny;
	y = cosr * sinp * cosy + sinr * cosp * siny;
	m_z = cosr * cosp * siny - sinr * sinp * cosy;
	m_w = cosr * cosp * cosy + sinr * sinp * siny;

	Normalise();
}

bool CQuaternion::ToEulerAngleXZY(double &xAngle, double &zAngle, double &yAngle)
{
	// rot =  cy*cz          -sz              cz*sy
	//        sx*sy+cx*cy*sz  cx*cz          -cy*sx+cx*sy*sz
	//       -cx*sy+cy*sx*sz  cz*sx           cx*cy+sx*sy*sz

	//     = 1-2*y^2-2*z^2    2xy-2wz         2xz+2wy
	//       2xy+2zw          1-2*x^2-2*z^2   2yz-2wx
	//       2xz-2yw          2yz+2wx         1-2*x^2-2*y^2

	double fx = 2.0 * m_x;
	double fy = 2.0 * m_y;
	double fz = 2.0 * m_z;

	double data10 = fx * m_y + fz * m_w;
	double data12 = fy * m_z - fx * m_w;

	double sz = -(fx * m_y - fz * m_w);
	if (sz > 1.0 - ms_fEpsilon)
	{
		// WARNING.  Not a unique solution.
		zAngle = HALF_PI;
		yAngle = 0.0;  // any angle works

		double fRpY = atan2(-data10, data12);
		xAngle = fRpY - yAngle;
		return false;
	}
	else if (sz < -1.0 + ms_fEpsilon)
	{
		// WARNING.  Not a unique solution.
		zAngle = -HALF_PI;
		yAngle = 0.0;  // any angle works

		double fRmY = atan2(-data10, data12);		
		xAngle = yAngle - fRmY;
		return false;
	}

	zAngle = asin(sz);
	double data21 = fx * m_z - fy * m_w;
	double data11 = 1 - fx * m_x - fz * m_z;
	xAngle = atan2(data21, data11);

	double data02 = fx * m_z + fy * m_w;
	double data00 = 1 - fy * m_y - fz * m_z;
	yAngle = atan2(data02, data00);
	return true;
}

bool CQuaternion::ToEulerAnglesZXY(double &zAngle, double &xAngle, double &yAngle) const
{
	if (fabs(m_w) < ms_fEpsilon)
	{
		zAngle = 0.0;
		xAngle = 0.0;
		yAngle = 0.0;
		return true;
	}
	// rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
	//        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
	//       -cx*sy           sx              cx*cy

	//     = 1-2*y^2-2*z^2    2xy-2wz         2xz+2wy
	//       2xy+2zw          1-2*x^2-2*z^2   2yz-2wx
	//       2xz-2yw          2yz+2wx         1-2*x^2-2*y^2

	double fx = 2.0 * m_x;
	double fy = 2.0 * m_y;
	double fz = 2.0 * m_z;
	
	double data00 = 1.0 - fy * m_y - fz * m_z;
	double data02 = fx * m_z + fy * m_w;

	double data21 = fy * m_z + fx * m_w;
	if (data21 > 1 - ms_fEpsilon)   // x angle is 90 degree
	{
		// WARNING.  Not a unique solution.
		xAngle = HALF_PI;
		yAngle = 0.0;  // any angle works
		
		double fRpY = atan2(data02, data00);		
		zAngle = fRpY - yAngle;
		return false;
	}
	else if (data21 < -1.0 + ms_fEpsilon)  // x angle is -90 degree
	{
		// WARNING.  Not a unique solution.
		xAngle = -HALF_PI;
		yAngle = 0.0;              // any angle works
		
		double fRmY = atan2(data02, data00);		
		zAngle = yAngle - fRmY;
		return false;
	}
	
	xAngle = asin(data21);

	double data01 = fx * m_y - fz * m_w;
	double data11 = 1 - fx * m_x - fz * m_z;
	zAngle = atan2(-data01, data11);

	double data20 = fx * m_z - fy * m_w;
	double data22 = 1 - fx * m_x - fy * m_y;
	yAngle = atan2(-data20, data22);
	return true;
	
}

//-----------------------------------------------------------------------
double CQuaternion::getRoll(bool reprojectAxis) const
{
	if (reprojectAxis)
	{
		// roll = atan2(localx.y, localx.x)
		// pick parts of xAxis() implementation that we need
		double fTx  = 2.0*m_x;
		double fTy  = 2.0*m_y;
		double fTz  = 2.0*m_z;
		double fTwz = fTz*m_w;
		double fTxy = fTy*m_x;
		double fTyy = fTy*m_y;
		double fTzz = fTz*m_z;
		
		// CVector3D(1.0-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
		
		return atan2(fTxy+fTwz, 1.0-(fTyy+fTzz));
		
	}
	else
	{
		return atan2(2*(m_x*m_y + m_w*m_z), m_w*m_w + m_x*m_x - m_y*m_y - m_z*m_z);
	}
}
//-----------------------------------------------------------------------
double CQuaternion::getPitch(bool reprojectAxis) const
{
	if (reprojectAxis)
	{
		// pitch = atan2(localy.z, localy.y)
		// pick parts of yAxis() implementation that we need
		double fTx  = 2.0*m_x;
		double fTy  = 2.0*m_y;
		double fTz  = 2.0*m_z;
		double fTwx = fTx*m_w;
		double fTxx = fTx*m_x;
		double fTyz = fTz*m_y;
		double fTzz = fTz*m_z;
		
		// CVector3D(fTxy-fTwz, 1.0-(fTxx+fTzz), fTyz+fTwx);
		return atan2(fTyz+fTwx, 1.0-(fTxx+fTzz));
	}
	else
	{
		// internal version
		return atan2(2*(m_y*m_z + m_w*m_x), m_w*m_w - m_x*m_x - m_y*m_y + m_z*m_z);
	}
}
//-----------------------------------------------------------------------
double CQuaternion::getYaw(bool reprojectAxis) const
{
	if (reprojectAxis)
	{
		// yaw = atan2(localz.x, localz.z)
		// pick parts of zAxis() implementation that we need
		double fTx  = 2.0*m_x;
		double fTy  = 2.0*m_y;
		double fTz  = 2.0*m_z;
		double fTwy = fTy*m_w;
		double fTxx = fTx*m_x;
		double fTxz = fTz*m_x;
		double fTyy = fTy*m_y;
		
		// CVector3D(fTxz+fTwy, fTyz-fTwx, 1.0-(fTxx+fTyy));
		
		return atan2(fTxz+fTwy, 1.0-(fTxx+fTyy));
		
	}
	else
	{
		// internal version
		return asin(-2*(m_x*m_z - m_w*m_y));
	}
}
//-----------------------------------------------------------------------
CQuaternion CQuaternion::nlerp(double fT, const CQuaternion& rkP,
							 const CQuaternion& rkQ, bool shortestPath)
{
	CQuaternion result;
	double fCos = rkP.Dot(rkQ);
	if (fCos < 0.0f && shortestPath)
	{
		result = rkP + fT * ((-rkQ) - rkP);
	}
	else
	{
		result = rkP + fT * (rkQ - rkP);
	}
	result.Normalise();
	return result;
}
