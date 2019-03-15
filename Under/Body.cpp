#include"Body.h"
#include<iostream>
using namespace std;

int Body::idnum = 0;
Body::Body(PICnew *m_picnew, Matrix3d R, Vector3d y) {

	this->m_picnew = m_picnew;
	onepointST *aonepoint = new onepointST();
	int facenum = m_picnew->faceandnormal.size();
	this->faceNum = facenum;
	for (int i = 0;i < facenum;i++) {
		aonepoint->id=idnum ;//静态变量，只与类有关
		idnum++;
		aonepoint->normal =R* m_picnew->faceandnormal[i].faceNormal+y;
		//这个面上三点坐标
		Vector3d p0 = m_picnew->vertexandnormal[m_picnew->faceandnormal[i].vertexIndex[0]].coordinate;
		Vector3d p1 = m_picnew->vertexandnormal[m_picnew->faceandnormal[i].vertexIndex[1]].coordinate;
		Vector3d p2 = m_picnew->vertexandnormal[m_picnew->faceandnormal[i].vertexIndex[2]].coordinate;

		double a = sqrt((p0 - p1)(0)*(p0 - p1)(0) + (p0 - p1)(1)*(p0 - p1)(1) + (p0 - p1)(2)*(p0 - p1)(2));
		double b = sqrt((p0 - p2)(0)*(p0 - p2)(0) + (p0 - p2)(1)*(p0 - p2)(1) + (p0 - p2)(2)*(p0 - p2)(2));
		double c = sqrt((p2 - p1)(0)*(p2 - p1)(0) + (p2 - p1)(1)*(p2 - p1)(1) + (p2 - p1)(2)*(p2 - p1)(2));
		double p = (a + b + c) / 2;
		aonepoint->area = sqrt(p*(p - a)*(p - b)*(p - c) );
		aonepoint->midpoint = R *( (p0 + p1 + p2) / 3) + y;
		this->v_onepoint.push_back(*aonepoint);

	}
	Vector3d Zero;//先设置个初始速度，加了重力之后，记得改回来！！！！！
	Zero.setZero();
	Vector3d one(0, -1, 0);
	Vector3d zero(0, 0, 0);
	VectorXd temp(6);
	temp.block(0, 0, 3, 1) = zero;
	temp.block(3, 0, 3, 1) = one;
	this->epsilon = temp;
	Matrix4d tempg;
	tempg.block(0, 0, 3, 3) = R;
	tempg.block(0, 3, 3, 1) = y;
	g = tempg;

	this->K = computeKB();
}

void Body::Subexpressions(double &w0, double &w1, double &w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2) {
	double temp0 = w0 + w1;
	f1 = temp0 + w2;
	double temp1 = w0 * w0;
	double temp2 = temp1 + w1 * temp0;
	f2 = temp2 + w2 * f1;
	f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
	g0 = f2 + w0 * (f1 + w0);
	g1 = f2 + w1 * (f1 + w1);
	g2 = f2 + w2 * (f1 + w2);
}
MatrixXd Body::computeKB() {
	MatrixXd res(6, 6);
	res.setZero(6, 6);
	Matrix3d J = computeJ();//计算结束得到mass，注意并行问题
	res.block(0, 0, 3, 3) = J;
	Matrix3d identity;
	identity.setIdentity(3, 3);
	res.block(3, 3, 3, 3) = bodyMass * identity;//标量乘以矩阵
	return res;

}
Matrix3d Body::computeJ() {
	Matrix3d inertia;
	const double mult[10] = { 1.0 / 6, 1.0 / 24, 1.0 / 24, 1.0 / 24,
		1.0 / 60, 1.0 / 60, 1.0 / 60, 1.0 / 120, 1.0 / 120, 1.0 / 120 };
	double intg[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx 

	double f1x = 0, f2x = 0, f3x = 0;
	double f1y = 0, f2y = 0, f3y = 0;
	double f1z = 0, f2z = 0, f3z = 0;
	double g0x = 0, g1x = 0, g2x = 0;
	double g0y = 0, g1y = 0, g2y = 0;
	double g0z = 0, g1z = 0, g2z = 0;
	
	for (int t = 0; t < faceNum; t++) { // get vertices of triangle t 
		/*
		i0 = index[3*t]; 	i1 = index[3*t+1];		i2 = index[3*t+2];
		x0 = p[i0].x; 		y0 = p[i0].y; 		z0 = p[i0].z;
		x1 = p[i1].x;       y1 = p[i1].y;       z1 = p[i1].z;
		x2 = p[i2].x;       y2 = p[i2].y;       z2 = p[i2].z;
		*/
		int index0 = m_picnew->faceandnormal[t].vertexIndex[0];
		int index1 = m_picnew->faceandnormal[t].vertexIndex[1];
		int index2 = m_picnew->faceandnormal[t].vertexIndex[2];
		double x0 = m_picnew->vertexandnormal[index0].coordinate(0);
		double y0 = m_picnew->vertexandnormal[index0].coordinate(1);
		double z0 = m_picnew->vertexandnormal[index0].coordinate(2);

		double x1 = m_picnew->vertexandnormal[index1].coordinate(0);
		double y1 = m_picnew->vertexandnormal[index1].coordinate(1);
		double z1 = m_picnew->vertexandnormal[index1].coordinate(2);

		double x2 = m_picnew->vertexandnormal[index2].coordinate(0);
		double y2 = m_picnew->vertexandnormal[index2].coordinate(1);
		double z2 = m_picnew->vertexandnormal[index2].coordinate(2);
		// get edges and cross product of edges 
		double a1 = x1 - x0; double b1 = y1 - y0; double c1 = z1 - z0;
		double a2 = x2 - x0; double b2 = y2 - y0; double c2 = z2 - z0;
		double d0 = b1 * c2 - b2 * c1; double d1 = a2 * c1 - a1 * c2; double d2 = a1 * b2 - a2 * b1;
		// compute integral terms
		Subexpressions(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
		Subexpressions(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
		Subexpressions(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
		// update integrals 
		intg[0] += d0 * f1x;
		intg[1] += d0 * f2x; intg[2] += d1 * f2y; intg[3] += d2 * f2z;
		intg[4] += d0 * f3x; intg[5] += d1 * f3y; intg[6] += d2 * f3z;
		intg[7] += d0 * (y0*g0x + y1 * g1x + y2 * g2x);
		intg[8] += d1 * (z0*g0y + z1 * g1y + z2 * g2y);
		intg[9] += d2 * (x0*g0z + x1 * g1z + x2 * g2z);
	}
	for (int i = 0;i < 10;i++) {
		intg[i] *= mult[i];
	}
	volume = intg[0];
	cout << "volume" << volume << endl;
	for (int i = 0;i < 10;i++) {
		intg[i] *= bodyDensity;
	}
	bodyMass = intg[0];
	// 质心
	Vector3d cm(intg[1] / bodyMass, intg[2] / bodyMass, intg[3] / bodyMass);
	// 相对于质心的惯性张量   0x 1y 2z
	inertia(0, 0) = intg[5] + intg[6] - bodyMass * (cm(1)*cm(1) + cm(2)*cm(2));//yz 12
	inertia(1, 1) = intg[4] + intg[6] - bodyMass * (cm(2)*cm(2) + cm(0)*cm(0));//zx 20
	inertia(2, 2) = intg[4] + intg[5] - bodyMass * (cm(0)*cm(0) + cm(1)*cm(1));//xy 01
	inertia(0, 1) = inertia(1, 0) = -(intg[7] - bodyMass * cm(0)*cm(1));
	inertia(1, 2) = inertia(2, 1) = -(intg[8] - bodyMass * cm(1)*cm(2));
	inertia(0, 2) = inertia(2, 0) = -(intg[9] - bodyMass * cm(2)*cm(0));
	Matrix3d R = g.block(0, 0, 3, 3);
	Vector3d y = g.block(0, 3, 3, 1);
	this->masscenter = R * cm + y;

	return inertia;
}
Matrix3d Body::so3_cay(Vector3d tempw) {
	Matrix3d daOmega = so3_ad(tempw);
	double wlength2 = tempw(0)*tempw(0) + tempw(1)*tempw(1) + tempw(2)*tempw(2);
	Matrix3d res = Matrix3d::Identity();
	double xishu = 4 / (4 + wlength2);
	res = res + xishu * (daOmega + daOmega * daOmega / 2);
	//cout << "so3cay" << res;
	return res;
}
Matrix4d Body::se3_cay(VectorXd tempepsilon) {
	Vector3d tempw = tempepsilon.block(0, 0, 3, 1);//w
	Vector3d tempv = tempepsilon.block(3, 0, 3, 1);//v
	Matrix3d so3cay = so3_cay(delta_t*tempw);
	double wlength2 = tempw(0)*tempw(0) + tempw(1)*tempw(1) + tempw(2)*tempw(2);
	Matrix3d daOmega = so3_ad(tempw);
	double xishu = 2 / (4 + wlength2);
	Matrix3d tempmat = Matrix3d::Identity();
	tempmat = 2 * tempmat + daOmega;
	Matrix3d B = xishu * tempmat;
	Matrix4d res;
	res.block(0, 0, 3, 3) = so3cay;
	res.block(0, 3, 3, 1) = delta_t * B*tempv;
	//cout << "delta_t * B*v" << delta_t * B*tempv << endl;
	res(3, 0) = 0;
	res(3, 1) = 0;
	res(3, 2) = 0;
	res(3, 3) = 1;
	//cout << "se3_cay" << res << endl;
	return res;
}

MatrixXd Body::se3_ad(VectorXd tempepsilon) {
	Vector3d w = tempepsilon.block(0, 0, 3, 1);//w
	Vector3d v = tempepsilon.block(3, 0, 3, 1);//v
	MatrixXd res(6, 6);
	res.block(0, 0, 3, 3) = so3_ad(w);
	res.block(0, 3, 3, 3) = Matrix3d::Zero();
	res.block(3, 0, 3, 3) = so3_ad(v);
	res.block(3, 3, 3, 3) = so3_ad(w);
	return res;
}
VectorXd Body::se3_DEP(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk) {
	//Vector3d tempy = gk.block(0, 3, 3, 1);
	//Matrix3d Y = so3_ad(tempy);
	//MatrixXd ctln1 = se3_Ctln(delta_t * epsilon_now);
	//MatrixXd ctln2 = se3_Ctln(-delta_t * epsilon_last);
	//return ctln1 * K*epsilon_now -
	//	ctln2 * K*epsilon_last - delta_t * tf;//
	Vector3d tempy = gk.block(0, 3, 3, 1);
	Matrix3d Y = so3_ad(tempy);
	VectorXd tf = tsfs2tf(Y);
	MatrixXd ctln1 = se3_Ctln(delta_t * epsilon_now);
	MatrixXd ctln2 = se3_Ctln(-delta_t * epsilon_last);
	return ctln1 * K*epsilon_now -
		ctln2 * K*epsilon_last - delta_t * tf;//
}
VectorXd Body::tsfs2tf(Matrix3d Y) {
	Matrix3d R = g.block(0, 0, 3, 3);
	Matrix3d Rt = R.transpose();
	Matrix3d zero = Matrix3d::Zero();
	//zero.setZero(3, 3);
	Matrix3d negRtY = zero - Rt * Y;
	MatrixXd trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
	return trans * tsfs;
}
VectorXd Body::Unconstr_Dyn(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk) {
	VectorXd fepsilonk_est = se3_DEP(epsilon_now, epsilon_last, gk);//f epsilonk的估计值
	MatrixXd Jacobian(6, 6);
	Vector3d w = epsilon_now.block(0, 0, 3, 1);
	Vector3d v = epsilon_now.block(3, 0, 3, 1);
	Jacobian(0, 0) = K(0, 0) + delta_t * 0.5*w(2)*K(1, 0) - delta_t * 0.5*w(1)*K(2, 0);
	Jacobian(0, 1) = K(0, 1) + delta_t * 0.5*w(2)*K(1, 1) - delta_t * 0.5*w(1)*K(2, 1);
	Jacobian(0, 2) = K(0, 2) + delta_t * 0.5*w(2)*K(1, 2) - delta_t * 0.5*w(1)*K(2, 2);
	Jacobian(0, 3) = K(0, 3) + delta_t * 0.5*w(2)*K(1, 3) - delta_t * 0.5*w(1)*K(2, 3);
	Jacobian(0, 4) = K(0, 4) + delta_t * 0.5*w(2)*K(1, 4) - delta_t * 0.5*w(1)*K(2, 4);
	Jacobian(0, 5) = K(0, 5) + delta_t * 0.5*w(2)*K(1, 5) - delta_t * 0.5*w(1)*K(2, 5);

	Jacobian(1, 0) = -delta_t * 0.5*w(2)*K(0, 0) + K(1, 0) + delta_t * 0.5*w(0)*K(2, 0);
	Jacobian(1, 1) = -delta_t * 0.5*w(2)*K(0, 1) + K(1, 1) + delta_t * 0.5*w(0)*K(2, 1);
	Jacobian(1, 2) = -delta_t * 0.5*w(2)*K(0, 2) + K(1, 2) + delta_t * 0.5*w(0)*K(2, 2);
	Jacobian(1, 3) = -delta_t * 0.5*w(2)*K(0, 3) + K(1, 3) + delta_t * 0.5*w(0)*K(2, 3);
	Jacobian(1, 4) = -delta_t * 0.5*w(2)*K(0, 4) + K(1, 4) + delta_t * 0.5*w(0)*K(2, 4);
	Jacobian(1, 5) = -delta_t * 0.5*w(2)*K(0, 5) + K(1, 5) + delta_t * 0.5*w(0)*K(2, 5);

	Jacobian(2, 0) = delta_t * 0.5*w(1)*K(0, 0) - delta_t * 0.5*w(0)*K(1, 0) + K(2, 0);
	Jacobian(2, 1) = delta_t * 0.5*w(1)*K(0, 1) - delta_t * 0.5*w(0)*K(1, 1) + K(2, 1);
	Jacobian(2, 2) = delta_t * 0.5*w(1)*K(0, 2) - delta_t * 0.5*w(0)*K(1, 2) + K(2, 2);
	Jacobian(2, 3) = delta_t * 0.5*w(1)*K(0, 3) - delta_t * 0.5*w(0)*K(1, 3) + K(2, 3);
	Jacobian(2, 4) = delta_t * 0.5*w(1)*K(0, 4) - delta_t * 0.5*w(0)*K(1, 4) + K(2, 4);
	Jacobian(2, 5) = delta_t * 0.5*w(1)*K(0, 5) - delta_t * 0.5*w(0)*K(1, 5) + K(2, 5);

	Jacobian(3, 0) = delta_t * 0.5*v(2)*K(1, 0) - delta_t * 0.5*v(1)*K(2, 0) + K(3, 0) +
		delta_t * 0.5*w(2)*K(4, 0) - delta_t * 0.5*w(1)*K(5, 0);
	Jacobian(3, 1) = delta_t * 0.5*v(2)*K(1, 1) - delta_t * 0.5*v(1)*K(2, 1) + K(3, 1) +
		delta_t * 0.5*w(2)*K(4, 1) - delta_t * 0.5*w(1)*K(5, 1);
	Jacobian(3, 2) = delta_t * 0.5*v(2)*K(1, 2) - delta_t * 0.5*v(1)*K(2, 2) + K(3, 2) +
		delta_t * 0.5*w(2)*K(4, 2) - delta_t * 0.5*w(1)*K(5, 2);
	Jacobian(3, 3) = delta_t * 0.5*v(2)*K(1, 3) - delta_t * 0.5*v(1)*K(2, 3) + K(3, 3) +
		delta_t * 0.5*w(2)*K(4, 3) - delta_t * 0.5*w(1)*K(5, 3);
	Jacobian(3, 4) = delta_t * 0.5*v(2)*K(1, 4) - delta_t * 0.5*v(1)*K(2, 4) + K(3, 4) +
		delta_t * 0.5*w(2)*K(4, 4) - delta_t * 0.5*w(1)*K(5, 4);
	Jacobian(3, 5) = delta_t * 0.5*v(2)*K(1, 5) - delta_t * 0.5*v(1)*K(2, 5) + K(3, 5) +
		delta_t * 0.5*w(2)*K(4, 5) - delta_t * 0.5*w(1)*K(5, 5);

	Jacobian(4, 0) = -delta_t * 0.5*v(2)*K(0, 0) + delta_t * 0.5*v(0)*K(2, 0) -
		delta_t * 0.5*w(2)*K(3, 0) + K(4, 0) + delta_t * 0.5*w(0)*K(5, 0);
	Jacobian(4, 1) = -delta_t * 0.5*v(2)*K(0, 1) + delta_t * 0.5*v(0)*K(2, 1) -
		delta_t * 0.5*w(2)*K(3, 1) + K(4, 1) + delta_t * 0.5*w(0)*K(5, 1);
	Jacobian(4, 2) = -delta_t * 0.5*v(2)*K(0, 2) + delta_t * 0.5*v(0)*K(2, 2) -
		delta_t * 0.5*w(2)*K(3, 2) + K(4, 2) + delta_t * 0.5*w(0)*K(5, 2);
	Jacobian(4, 3) = -delta_t * 0.5*v(2)*K(0, 3) + delta_t * 0.5*v(0)*K(2, 3) -
		delta_t * 0.5*w(2)*K(3, 3) + K(4, 3) + delta_t * 0.5*w(0)*K(5, 3);
	Jacobian(4, 4) = -delta_t * 0.5*v(2)*K(0, 4) + delta_t * 0.5*v(0)*K(2, 4) -
		delta_t * 0.5*w(2)*K(3, 4) + K(4, 4) + delta_t * 0.5*w(0)*K(5, 4);
	Jacobian(4, 5) = -delta_t * 0.5*v(2)*K(0, 5) + delta_t * 0.5*v(0)*K(2, 5) -
		delta_t * 0.5*w(2)*K(3, 5) + K(4, 5) + delta_t * 0.5*w(0)*K(5, 5);

	Jacobian(5, 0) = delta_t * 0.5*v(1)*K(0, 0) - delta_t * 0.5*v(0)*K(1, 0) +
		delta_t * 0.5*w(1)*K(3, 0) - delta_t * 0.5*w(0)*K(4, 0) + K(5, 0);
	Jacobian(5, 1) = delta_t * 0.5*v(1)*K(0, 1) - delta_t * 0.5*v(0)*K(1, 1) +
		delta_t * 0.5*w(1)*K(3, 1) - delta_t * 0.5*w(0)*K(4, 1) + K(5, 1);
	Jacobian(5, 2) = delta_t * 0.5*v(1)*K(0, 2) - delta_t * 0.5*v(0)*K(1, 2) +
		delta_t * 0.5*w(1)*K(3, 2) - delta_t * 0.5*w(0)*K(4, 2) + K(5, 2);
	Jacobian(5, 3) = delta_t * 0.5*v(1)*K(0, 3) - delta_t * 0.5*v(0)*K(1, 3) +
		delta_t * 0.5*w(1)*K(3, 3) - delta_t * 0.5*w(0)*K(4, 3) + K(5, 3);
	Jacobian(5, 4) = delta_t * 0.5*v(1)*K(0, 4) - delta_t * 0.5*v(0)*K(1, 4) +
		delta_t * 0.5*w(1)*K(3, 4) - delta_t * 0.5*w(0)*K(4, 4) + K(5, 4);
	Jacobian(5, 5) = delta_t * 0.5*v(1)*K(0, 5) - delta_t * 0.5*v(0)*K(1, 5) +
		delta_t * 0.5*w(1)*K(3, 5) - delta_t * 0.5*w(0)*K(4, 5) + K(5, 5);
	VectorXd delta_epsilion = Jacobian.inverse() * fepsilonk_est;
	//cout << "delta_epsilion是是是=" << delta_epsilion << endl;
	return epsilon_now - delta_epsilion;
}
Matrix3d Body::so3_ad(Vector3d omega) {
	Matrix3d res;
	res(0, 0) = 0;
	res(0, 1) = -omega(2);
	res(0, 2) = omega(1);
	res(1, 0) = omega(2);
	res(1, 1) = 0;
	res(1, 2) = -omega(0);
	res(2, 0) = -omega(1);
	res(2, 1) = omega(0);
	res(2, 2) = 0;
	return res;
}
MatrixXd Body::se3_Ctln(VectorXd tempepsilon) {
	MatrixXd identity(6, 6);
	identity.setIdentity();
	MatrixXd temp = 0.5f*se3_ad(tempepsilon);
	return identity - temp;
}
void Body::nextTime() {
	//lp_ = computelp_();
	//cout << "oldlp" << lp << endl;
	//cout << "nextlp_:" << lp_<< endl;
	//lp=computeNextlp();//lp
	//cout << "lp:" << lp << endl;
	//VectorXd tempwv= computeNextwv();
	//cout << "tempwv" << tempwv << endl;
	//cout << "w:"<<w(0)<<" " << w(1) << " " << w(2) << endl;
	//cout << "v:"<<v(0) << " " << v(1) << " " << v(2) << endl;
	Vector3d tempy = g.block(0, 3, 3, 1);
	Matrix3d Y = so3_ad(tempy);
	VectorXd tf = tsfs2tf(Y);
	cout << "R" << g.block(0,0,3,3)<<endl<<"y"<<tempy << endl;
	cout << "tf物体坐标系的力矩和力" << tf << endl;
	VectorXd epsilon_last = epsilon;
	cout << "epsilon:" << epsilon << endl;
	VectorXd delta_epsilon = delta_t * K.inverse()*(se3_ad(delta_t*epsilon_last)*K*epsilon_last + tf);
	VectorXd epsilon_now = epsilon_last + delta_epsilon;
	//cout << "epsilon_now 预估" << epsilon_now << endl;
	VectorXd res = se3_DEP(epsilon_now, epsilon_last, g);
	//cout << "偏差res=" << res << endl;
	//牛顿迭代法求解方程组
	int i = 0;
	double cancha = 1e-14;
	do {//设置残差值
		epsilon_now = Unconstr_Dyn(epsilon_now, epsilon_last, g);
		//cout << "epsilon_now 迭代后的值" << epsilon_now << endl;
		res = se3_DEP(epsilon_now, epsilon_last, g);
		//cout << "偏差res=" << res << endl;
		i++;
	} while ((res(0) > cancha || res(1) > cancha || res(2) > cancha ||
		res(3) > cancha || res(4) > cancha || res(5) > cancha ||
		res(0) < -cancha || res(1) < -cancha || res(2) < -cancha ||
		res(3) < -cancha || res(4) < -cancha || res(5) < -cancha
		) && i < 50);
	//cout << "迭代了多少次？" << i << endl;
	//cout << "偏差res=" << res << endl;
	Matrix4d se3cay = se3_cay(delta_t * epsilon_now);
	//cout << "se3cay==" << se3cay << endl;
	g = g * se3cay;
	//cout << "g" << g << endl;
	//R = g.block(0, 0, 3, 3);
	//y = g.block(0, 3, 3, 1);
	cout << "迭代后的epsilon_now:" << epsilon_now << endl;
	this->epsilon = epsilon_now;
	//Vector3d y_ = computey_();
	//cout << "y_" << y_ << endl;
	//y=computeNexty(y_);//y
	//cout << "y!!!!!!!!!!!!!!:" << y << endl;
	//R=computeNextR();//R
	//cout << "R:" << R << endl;
	//cout << "w:" << epsilon.block(0, 0, 3, 1) << endl;
	//cout << "v:" << epsilon.block(3, 0, 3, 1) << endl;
	//w = tempwv.block(0, 0, 3, 1);//w
	//v = tempwv.block(3, 0, 3, 1);//v
}
float* Body::GetRotAndTransData() {//没有trans
	Matrix3d R = g.block(0, 0, 3, 3);
	static float data[16];//!!!!
	data[0] = R(0, 0);
	data[1] = R(1, 0);
	data[2] = R(2, 0);
	data[3] = 0;

	data[4] = R(0, 1);
	data[5] = R(1, 1);
	data[6] = R(2, 1);
	data[7] = 0;

	data[8] = R(0, 2);
	data[9] = R(1, 2);
	data[10] = R(2, 2);
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;
	return data;
}
void Body::computetsfs(VectorXd traction) {
	Vector3d fs(0, 0, 0);
	Vector3d ts(0, 0, 0);
	for (int i = 0;i < faceNum;i++) {
		Vector3d faceifs=traction.block(3 * i, 0, 3, 1);//第i个面上的力
		fs += faceifs;
		Vector3d r = v_onepoint[i].midpoint - masscenter;
		Vector3d faceits = faceifs.cross(r);//第i面上计算得到的力矩
		//cout << "力矩在第" << i << "个面是：" << faceits << endl;
		ts += faceits;
	}
	VectorXd temptsfs(6);
	temptsfs.block(0, 0, 3, 1) = ts;
	temptsfs.block(3, 0, 3, 1) = fs;
	tsfs = temptsfs;//赋值成员变量
	cout << "tsfs" << tsfs << endl;
}
