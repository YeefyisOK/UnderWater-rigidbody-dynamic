#include "DynamicFormula2.h"
double DynamicFormula2::mu=0.5;

//计算每个面上的应力
VectorXd DynamicFormula2::computetraction(vector<Body*> m_body){
	int n = m_body[0]->idnum;//获取总面片数
	MatrixXd coefficient(3*n, 3*n);
	coefficient.setZero();
	VectorXd traction(3 * n);
	MatrixXd H(3*n, 3*n);
	H.setZero();
	VectorXd u(3*n);		
	for (int p = 0;p < m_body.size();p++) {
		for (int q = 0;q < m_body[p]->v_onepoint.size();q++) {
			Vector3d y = m_body[p]->v_onepoint[q].midpoint;// 得到y
			int b= m_body[p]->v_onepoint[q].id;//下标符号
			Vector3d ny = m_body[p]->v_onepoint[q].normal;//得到y的法向
			Vector3d w = m_body[p]->epsilon.block(0, 0, 3, 1);
			Vector3d v = m_body[p]->epsilon.block(3, 0, 3, 1);
			u.block(3 * b,0,3,1) = w.cross(m_body[p]->v_onepoint[q].midpoint) + v;
			for (int i = 0;i < m_body.size();i++) {
				for (int j = 0;j < m_body[i]->v_onepoint.size();j++) {
					Vector3d x = m_body[i]->v_onepoint[j].midpoint;//得到sourcepoint x
					int a = m_body[i]->v_onepoint[j].id;//下标符号
					Vector3d nx = m_body[i]->v_onepoint[j].normal;//得到x的法向
					if (x != y) {
						Matrix3d KS=	computeKij(x, nx, y, ny)*m_body[p]->v_onepoint[q].area;
						coefficient.block(3*a, 3*b, 3, 3) = KS;
						Matrix3d HS = computeHij(x, nx, y, ny)*m_body[p]->v_onepoint[q].area;
						//cout << "HS"<<endl << HS << endl;
						H.block(3*a, 3*b, 3, 3) = HS;												
					}					
				}
			}
		}
	}
	for (int i = 0;i < n;i++) {
		Matrix3d HSum;
		HSum.setZero();
		Matrix3d coefficientSum;
		coefficientSum.setZero();
		for (int j = 0;j < n&&i!=j;j++) {
			HSum += H.block(3 * i, 3 * j, 3, 3);
			coefficientSum += coefficient.block(3 * i, 3 * j, 3, 3);
		}
		Matrix3d identity;
		identity.setIdentity();
		H.block(3 * i, 3 * i, 3, 3) = identity-HSum;//对角线上的奇异元素是这样计算吗
		coefficient.block(3 * i, 3 * i, 3, 3) = identity-coefficientSum;
	}
	//cout <<"coefficient"<<endl<< coefficient << endl;
	//解线性方程组
	VectorXd b = H * u;
	//cout << "H:" << endl << H << endl;
	//cout << "u:" << endl << u << endl;
	//cout << "b:" << endl << b << endl;
	traction = coefficient.fullPivHouseholderQr().solve(b);//sigma=strength NAN!
	cout << "traction:" << endl << traction << endl;

	return traction;
}

Matrix3d DynamicFormula2::computeKij(Vector3d x, Vector3d nx, Vector3d y, Vector3d ny) {
	Matrix3d res;
	res.setZero();
	double r =sqrt( (x(0) - y(0))* (x(0) - y(0)) +
		(x(1) - y(1))*(x(1) - y(1)) + (x(2) - y(2))*(x(2) - y(2)) );//r的平方
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			double rknk = 0;
			for (int k = 0;k < 3;k++) {
				//cout << "yk-xk" << y(k) - x(k) << endl;
				//cout << "nx(k)" << nx(k) << endl;
				rknk += (y(k) - x(k)) / r * nx(k);
			}
			res(i, j) = 3 / (4 * 3.1415926*r*r)*(y(i) - x(i))/r*(y(j) - x(j))/r*rknk;
			//cout << "i:" << i << "j:" << j << res(i, j);
		}
	}
	return res;
}
double DynamicFormula2::dirac(int a, int b) {
	if (a == b) {
		return 1;
	}
	else {
		return 0;
	}
}

Matrix3d DynamicFormula2::computeHij(Vector3d x, Vector3d nx, Vector3d y, Vector3d ny) {
	Matrix3d H;
	H.setZero();
	double r = sqrt((x(0) - y(0))* (x(0) - y(0)) +
		(x(1) - y(1))*(x(1) - y(1)) + (x(2) - y(2))*(x(2) - y(2)));//r的平方
	Vector3d n = (y - x) / r;//y到x的单位向量
	double rlnl=0;
	for (int l = 0;l < 3;l++) {
		rlnl += (y(l) - x(l)) / r * ny(l);
	}
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3 ;j++) {
			for (int k = 0;k < 3;k++) {
				H(i, j) += mu / (4 * 3.1415926*r*r*r)*
				(
					(3 * (dirac(i, j)* (y(k) - x(k)) / r + dirac(j, k) * (y(i) - x(i)) / r)
						- 30 * (y(i) - x(i)) / r * (y(j) - x(j)) / r * (y(k) - x(k)) / r
					)  *rlnl
						+ 3 * (n(i)*(y(j) - x(j)) / r * (y(k) - x(k)) / r + n(k) *(y(i) - x(i)) / r * (y(j) - x(j)) / r
								 )
						+ 2 * dirac(i, k)*n(j)
				)*nx(k);
				//cout << "H i =" << i << "j=" << j << H(i, j) << endl;
			}
		}
	}
	return H;
}


