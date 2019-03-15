#include<windows.h>
#include <fstream>
#include<sstream> 
#include<GL/glut.h>
#include<stdio.h>
#include<string>
#include"Kirchhoff.h"
#include "PIC.h"
#include"PICnew.h"
#include"Body.h"
#include"DynamicFormula2.h"
#include<iostream>
using namespace std;
using namespace Eigen;
/*
定时器：glutTimerFunc(unsigned int millis, void (*func)(int value), int value);
重绘标志：glutPosRedisplay();
定时器第一个参数是每隔millis毫秒便调用func函数，并且创一个value参数进去
因为一个定时器只被调用一次，所以需要多次调用定时器
*/
//obj读取   ../source/yuanpan.obj

int id = 0;
long imagewidth = 600;
long imageheight = 800;
int modelNum = 2;//模型数量 sanlengzhui  2个的话改成2
string name[2] = { "H:\\MeshData\\sanlengzhui.obj", "H:\\MeshData\\cube.obj"};//ell0  myproplab
//模型数组
vector<PIC*> v_pic;
vector<PICnew*>  v_picnew;
vector<Body*> v_body;
Vector3d y_0(0, 0, 0);
Vector3d y_1(3, 0, 0);
Vector3d y[2] = { y_0,y_1 };
//窗口的大小
GLdouble windowWidth;
GLdouble windowHeight;
Vector3d omega(0, 0, 0);
Vector3d velocity(0, 0, 0);
Matrix3d R = Matrix3d::Identity();//设置为单位阵 在init()改不是单位阵
double m_fluidDensity = 0.98;
double m_bodyDensity = 1.75;
double delta_t = 0.04;

//DynamicFormula m_DF(omega, velocity, R, y, delta_t);
bool mouseLeftDown;
bool mouseRightDown;
double mouseX, mouseY;
double cameraDistance = 0;
double cameraAngleX;
double cameraAngleY;
void ReadPIC()//把所有obj文件读取到v_pic中,v_pic转化v_picnew
{
	for (int i = 0;i < modelNum;i++) {
		ifstream ifs(name[i]);//cube bunny Eight
		string s;
		Mian *f;
		Vertex *v;
		FaXiangLiang *vn;
		WenLi  *vt;
		PIC *m_pic = new PIC();
		while (getline(ifs, s))
		{
			if (s.length() < 2)continue;
			if (s[0] == 'v') {
				if (s[1] == 't') {//vt 0.581151 0.979929 纹理
					istringstream in(s);
					vt = new WenLi();
					string head;
					in >> head >> vt->TU >> vt->TV;
					m_pic->VT.push_back(*vt);
				}
				else	if (s[1] == 'n') {//vn 0.637005 -0.0421857 0.769705 法向量
					istringstream in(s);
					vn = new FaXiangLiang();
					string head;
					in >> head >> vn->NX >> vn->NY >> vn->NZ;
					m_pic->VN.push_back(*vn);
				}
				else {//v -53.0413 158.84 -135.806 点
					istringstream in(s);
					v = new Vertex();
					string head;
					in >> head >> v->X >> v->Y >> v->Z;
					m_pic->V.push_back(*v);
				}
			}
			else if (s[0] == 'f') {//f 2443//2656 2442//2656 2444//2656 面
				for (int k = s.size() - 1; k >= 0; k--) {
					if (s[k] == '/')s[k] = ' ';
				}
				istringstream in(s);
				f = new Mian();
				string head;
				in >> head;
				int i = 0;
				while (i < 3)
				{//
					if (m_pic->V.size() != 0)
					{
						in >> f->V[i];
						f->V[i] -= 1;
					}
					if (m_pic->VT.size() != 0)
					{
						in >> f->T[i];
						f->T[i] -= 1;
					}
					if (m_pic->VN.size() != 0)
					{
						in >> f->N[i];
						f->N[i] -= 1;
					}
					i++;
				}
				m_pic->F.push_back(*f);
			}
		}
		v_pic.push_back(m_pic);

	}
	for (int i = 0;i < v_pic.size();i++) {
		PICnew *a_picnew = new PICnew( v_pic[i]);
		v_picnew.push_back(a_picnew);
	}

}
void GLDraw(int k)//k表示是第几个模型
{
	glColor3f(0.0, 1.0, 0.0);     //绿
	PIC *m_pic = v_pic[k];
	for (int i = 0; i < m_pic->F.size(); i++)
	{
		glBegin(GL_TRIANGLE_FAN);                            // 绘制三角形GL_TRIANGLES;GL_LINE_LOOP;GL_LINES;GL_POINTS
		if (m_pic->VT.size() != 0)glTexCoord2f(m_pic->VT[m_pic->F[i].T[0]].TU, m_pic->VT[m_pic->F[i].T[0]].TV);  //纹理    
		if (m_pic->VN.size() != 0)glNormal3f(m_pic->VN[m_pic->F[i].N[0]].NX, m_pic->VN[m_pic->F[i].N[0]].NY, m_pic->VN[m_pic->F[i].N[0]].NZ);//法向量
		glVertex3f(m_pic->V[m_pic->F[i].V[0]].X, m_pic->V[m_pic->F[i].V[0]].Y, m_pic->V[m_pic->F[i].V[0]].Z);        // 上顶点

		if (m_pic->VT.size() != 0)glTexCoord2f(m_pic->VT[m_pic->F[i].T[1]].TU, m_pic->VT[m_pic->F[i].T[1]].TV);  //纹理
		if (m_pic->VN.size() != 0)glNormal3f(m_pic->VN[m_pic->F[i].N[1]].NX, m_pic->VN[m_pic->F[i].N[1]].NY, m_pic->VN[m_pic->F[i].N[1]].NZ);//法向量
		glVertex3f(m_pic->V[m_pic->F[i].V[1]].X, m_pic->V[m_pic->F[i].V[1]].Y, m_pic->V[m_pic->F[i].V[1]].Z);        // 左下

		if (m_pic->VT.size() != 0)glTexCoord2f(m_pic->VT[m_pic->F[i].T[2]].TU, m_pic->VT[m_pic->F[i].T[2]].TV);  //纹理
		if (m_pic->VN.size() != 0)glNormal3f(m_pic->VN[m_pic->F[i].N[2]].NX, m_pic->VN[m_pic->F[i].N[2]].NY, m_pic->VN[m_pic->F[i].N[2]].NZ);//法向量
		glVertex3f(m_pic->V[m_pic->F[i].V[2]].X, m_pic->V[m_pic->F[i].V[2]].Y, m_pic->V[m_pic->F[i].V[2]].Z);        // 右下
		glEnd();
		// 三角形绘制结束    
		/*if(m_pic.VN.size()!=0){
		glBegin(GL_LINES);                            // 绘制三角形
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X/YU,m_pic.V[m_pic.F[i].V[0]].Y/YU, m_pic.V[m_pic.F[i].V[0]].Z/YU);        // 上顶点
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X/YU+m_pic.VN[m_pic.F[i].N[0]].NX
		,m_pic.V[m_pic.F[i].V[0]].Y/YU+m_pic.VN[m_pic.F[i].N[0]].NY
		, m_pic.V[m_pic.F[i].V[0]].Z/YU+m_pic.VN[m_pic.F[i].N[0]].NZ);                    // 左下
		glEnd();                                // 三角形绘制结束
		}*/
	}
}
void write() {
	glPushMatrix();
	glPopMatrix();
}
void init() {
	/*
	Vector3d temp = R.row(0);
	R.row(0) = R.row(1);
	R.row(1) = temp;
	m_DF.setR(R);*/
	ReadPIC();
	cout << "y0" << y[0] << endl;
	cout << "y1" << y[1] << endl;
	for (int i = 0;i < modelNum;i++) {
		PICnew *picnew = v_picnew[i];
		Body *m_body = new Body(v_picnew[i], R, y[i]);
		v_body.push_back(m_body);
	}
	DynamicFormula2 m_DF2 ;
	VectorXd tractions= m_DF2.computetraction(v_body);//内部是3*1子矩阵
	int j = 0;
	for (int i = 0;i < modelNum;i++) {//对模型数量
		int facenum = v_body[i]->faceNum;
		//Vector3d 
		VectorXd abodytraction(facenum * 3);
		abodytraction=tractions.block(3*i,0, facenum*3,1);//一个物体面片受到的外力
		v_body[i]->computetsfs(abodytraction);
	}
	//基尔霍夫张量
	//CKirchhoff m_K(m_pic, m_bodyDensity, m_fluidDensity);
	//MatrixXd K = m_K.computeK();//初始时得到K矩阵
	//VectorXd m_tsfs = m_K.computetsfs();
	//cout << "tsfs" << m_tsfs << endl;
	//m_DF.tsfs = m_tsfs;

	//m_DF.q = R;
	//m_DF.setK(K);
	//m_DF.lp = m_DF.computelp();//计算初始的lp
	glClearColor(0.0, 0.0, 0.0, 0.0);  //背景色
	//深度测试的相关设置 
	glClearDepth(1.0);                    //设置深度缓存的初始值 
	glDepthFunc(GL_LEQUAL);           //深度测试的方法 
	glEnable(GL_DEPTH_TEST);          //启用深度测试
	GLfloat direction[] = { 0, -1, 0, 0 }; // 平行光源, GL_POSITION属性的最后一个参数为0
	GLfloat ambient[] = { 0.3f, 0.3f, 0.3f, 1.0f };  // 环境强度
	GLfloat diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };  // 散射强度
	GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f }; // 镜面强度
	glLightfv(GL_LIGHT0, GL_POSITION, direction);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

	GLfloat direction1[] = { 0, 0, 1, 0 };
	glLightfv(GL_LIGHT1, GL_POSITION, direction1);

	glEnable(GL_LIGHTING);   //开关:使用光
	glEnable(GL_LIGHT0);     //打开0#灯
	glEnable(GL_LIGHT1);     //打开0#灯
	glShadeModel(GL_SMOOTH);//中多边形内部各点的光照采用双线性插值的方法得到
	//设置多边形材质
	GLfloat mat_ambient[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat mat_diffuse[] = { 0.1, 0.5, 0.8, 1.0 };
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat shininess = 50.0f;

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);
	//使能GL_COLOR_MATERIAL，而后利用glColorMaterial()函数指明glColor*()函数将影响到的材质属性
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);

}

void drawScene()           //绘制
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 10, 0, 0, 0, 0, 1, 0);//4, 0, -2,在实际场景中，应放在所有模型变换之后，但在OpenGL中，由于顺序恰好相反，应放在所有模型变换之前
	glTranslatef(0, -cameraDistance, 0);
	//多物体绘制
	for (int i = 0;i < modelNum;i++) {
		glPushMatrix();
		float *pData = v_body[i]->GetRotAndTransData();
		//cout << pData[0] << " " << pData[4] << " " << pData[8] << " " << pData[12] << endl;
		//cout << pData[1] << " " << pData[5] << " " << pData[9] << " " << pData[13] << endl;
		//cout << pData[2] << " " << pData[6] << " " << pData[10] << " " << pData[14] << endl;
		//cout << pData[3] << " " << pData[7] << " " << pData[11] <<" " << pData[15] << endl;
		Vector3d tempy = v_body[i]->g.block(0, 3, 3, 1);
		cout << "第" << i << "个cube的位置是" << tempy << endl;
		glTranslatef(tempy(0), tempy(1), tempy(2));//与实际相反
		glMultMatrixf(pData);
		//glRotated(m_DF.theta, m_DF.delta_q.x(), m_DF.delta_q.y(), m_DF.delta_q.z());
		GLDraw(i);
		glPopMatrix();
	}
}
//窗口重绘时的响应函数 
void display() {
	//用前面设置的背景色和深度值清除颜色缓存和深度缓存 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//初始化模型视图矩阵 
	drawScene();            //几何图元的绘制 
	glFlush();              //绘制结束，Flush 当前渲染流水线 
	glutSwapBuffers();       //交换前后缓存（只用于双缓冲的模式） 
}
//窗口大小发生变化时的响应函数 
void reshape(int width, int height) {
	glViewport(0, 0, width, height);      //设置视窗大小 
	//设置视景体大小 
	glMatrixMode(GL_PROJECTION);
	double ratio = (double)width / height;
	glLoadIdentity();
	//gluPerspective(60, ratio, 1, 1000);
	int viewsize = 25;

	if (width <= height)
		glOrtho(-viewsize, viewsize, -viewsize * (GLdouble)height / (GLdouble)width, viewsize * (GLdouble)height / (GLdouble)width, 0, 20.0);
	else
		glOrtho(-viewsize * (GLdouble)width / (GLdouble)height, viewsize*(GLdouble)width / (GLdouble)height, -viewsize, viewsize, 0, 20.0);

	//glOrtho(-25, 25, -25, 25, -10, 10);
}

void TimerFunction(int value)
{
	DynamicFormula2 m_DF2;
	VectorXd tractions = m_DF2.computetraction(v_body);//内部是3*1子矩阵
	int j = 0;
	for (int i = 0;i < modelNum;i++) {//对模型数量
		int facenum = v_body[i]->faceNum;
		//Vector3d 
		//Vector3d sumf(0, 0, 0);//traction求和
		VectorXd abodytraction = tractions.block(3 * i, 0, facenum*3, 1);//一个物体面片受到的外力
		v_body[i]->computetsfs(abodytraction);
		v_body[i]->nextTime();
	}
	/*
	glutPostRedisplay 标记当前窗口需要重新绘制。通过glutMainLoop下一次循环时，
	窗口显示将被回调以重新显示窗口的正常面板。多次调用glutPostRedisplay，
	在下一个显示回调只产生单一的重新显示回调
	*/
	glutPostRedisplay(); //标志重新绘制
	glutTimerFunc(delta_t * 1000, TimerFunction, 1);
}
void mouseCB(int button, int state, int x, int y)
{
	mouseX = x;
	mouseY = y;
	if (button == GLUT_RIGHT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			mouseRightDown = true;
		}
		else if (state == GLUT_UP)
			mouseRightDown = false;
	}
}

void mouseMotionCB(int x, int y)
{
	if (mouseRightDown)
	{
		//cameraDistance = 0;
		cameraDistance += (y - mouseY) * 0.02f;
		mouseY = y;
	}
	glutPostRedisplay();
}
int main(int argc, char* argv[]) {
	glutInit(&argc, argv);    //GLUT 库的初始化 
	//显示模式初始化：颜色格式——GLUT_RGBA 
	//           单缓冲或双缓冲——GLUT_SINGLE 或者 GLUT_double //          是否使用深度缓存——GLUT_DEPTH  
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);          //窗口起始位置 
	glutInitWindowSize(800, 500);             //窗口大小 
	glutCreateWindow("UnderWaterRidgebody");       //创建窗口并指定窗口的名称 
	glutMouseFunc(mouseCB);
	glutMotionFunc(mouseMotionCB);

	init();                                 //OpenGL 的初始化设置 

	//设置消息响应的回调函数 
	glutDisplayFunc(display);               //设置窗口重绘时的响应函数 
	glutReshapeFunc(reshape);              //设置窗口大小发生变化时的响应函数 

	//定时器  每delta_t*1000毫秒触发一次
	glutTimerFunc(delta_t * 1000, TimerFunction, 1);

	glutMainLoop();                   //消息循环：获取消息，分发消息，响应消息 
	return 0;
} 
