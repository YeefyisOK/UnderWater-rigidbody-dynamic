#include "PICnew.h"
//PICnew::PICnew() {
//
//}
PICnew::PICnew(PIC *m_PIC) {
	for (int i = 0;i < m_PIC->F.size();i++) {
		FaceandNormalST *m_fan = new FaceandNormalST();
		m_fan->vertexIndex[0] = m_PIC->F[i].V[0];
		m_fan->vertexIndex[1] = m_PIC->F[i].V[1];
		m_fan->vertexIndex[2] = m_PIC->F[i].V[2];
		//计算面法向
		Vector3d temp0(m_PIC->V[m_fan->vertexIndex[0]].X, m_PIC->V[m_fan->vertexIndex[0]].Y, m_PIC->V[m_fan->vertexIndex[0]].Z);
		Vector3d temp1(m_PIC->V[m_fan->vertexIndex[1]].X, m_PIC->V[m_fan->vertexIndex[1]].Y, m_PIC->V[m_fan->vertexIndex[1]].Z);
		Vector3d temp2(m_PIC->V[m_fan->vertexIndex[2]].X, m_PIC->V[m_fan->vertexIndex[2]].Y, m_PIC->V[m_fan->vertexIndex[2]].Z);
		Vector3d normal = (temp1 - temp0).cross(temp2 - temp0);
		double length = normal.norm();
			//sqrt(
			//normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2)
		//);
		m_fan->faceNormal = normal/length;
		cout << i << "面法向为" << normal / length << endl;
		this->faceandnormal.push_back(*m_fan);
	}
	struct numandnormalST
	{
		int num;
		Vector3d normal;
	};//遍历face，记录每个顶点的法向
	int n = m_PIC->V.size();
	//numandnormal vn[n];
	vector <numandnormalST> vn;//顶点法向   0   1个   n1；1  4个   n1+n
	numandnormalST *numandnormal = new numandnormalST();
	for (int i = 0;i < n;i++) {//顶点数量2+n3+n4;
		vn.push_back(*numandnormal);
	}
	for (int i = 0;i < faceandnormal.size();i++) {//遍历所有的面
		vn[faceandnormal[i].vertexIndex[0]].num++;
		vn[faceandnormal[i].vertexIndex[0]].normal += faceandnormal[i].faceNormal;
		vn[faceandnormal[i].vertexIndex[1]].num++;
		vn[faceandnormal[i].vertexIndex[1]].normal += faceandnormal[i].faceNormal;
		vn[faceandnormal[i].vertexIndex[2]].num++;
		vn[faceandnormal[i].vertexIndex[2]].normal += faceandnormal[i].faceNormal;
			
		/*cout << "vn[faceandnormal[i].vertexIndex[0]]" 
			<< vn[faceandnormal[i].vertexIndex[0]].normal << endl;*/
	}
	for (int i = 0;i < vn.size();i++) {
		Vector3d temp(m_PIC->V[i].X, m_PIC->V[i].Y, m_PIC->V[i].Z);
		/*if (vn[i].normal(0) == 0 && vn[i].normal(1) == 0 && vn[i].normal(2) == 0) {
			Vector3d pianyi(0, 0.01, 0);
			Vector3d guiyihua = pianyi / pianyi.norm();
			vn[i].normal = guiyihua;
		}*/
		VertexandNormalST *vandn = new VertexandNormalST();
		vandn->coordinate = temp;
		vandn->vertexNormal = vn[i].normal / vn[i].num;
		cout << "第i个顶点的法向是" << vn[i].normal << endl;
		this->vertexandnormal.push_back(*vandn);
	}

}