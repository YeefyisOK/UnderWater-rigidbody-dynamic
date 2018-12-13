#include<windows.h>
#include <fstream>
#include<sstream> 
#include<GL/glut.h>
#include<stdio.h>
#include"Kirchhoff.h"
#include "PIC.h"
#include"DynamicFormula.h"
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
string name = "H:\\MeshData\\cube.obj";//tuoyuan.obj yuanpan20 bunnyclose  myprop2
PIC m_pic;
void drawScene();
//窗口的大小
GLfloat windowWidth;
GLfloat windowHeight;
Vector3f omega(	1, 0, 0);
Vector3f velocity(0, 0, 0);			
Matrix3f R = Matrix3f::Identity();//设置为单位阵 在init()改不是单位阵
Vector3f y(0,0,0);
Vector3f ts(0,0,0);
Vector3f fs(0,0,-50);
MatrixXf K;
float delta_t=1;

DynamicFormula m_DF(omega,velocity,R,y,ts,fs,K,delta_t);
bool mouseLeftDown;
bool mouseRightDown;
float mouseX, mouseY;
float cameraDistance;
float cameraAngleX;
float cameraAngleY;
void ReadPIC()
{
	ifstream ifs(name);//cube bunny Eight
	string s;
	Mian *f;
	Vertex *v;
	FaXiangLiang *vn;
	WenLi  *vt;
	while (getline(ifs, s))
	{
		if (s.length()<2)continue;
		if (s[0] == 'v'){
			if (s[1] == 't'){//vt 0.581151 0.979929 纹理
				istringstream in(s);
				vt = new WenLi();
				string head;
				in >> head >> vt->TU >> vt->TV;
				m_pic.VT.push_back(*vt);
			}
			else	if (s[1] == 'n'){//vn 0.637005 -0.0421857 0.769705 法向量
				istringstream in(s);
				vn = new FaXiangLiang();
				string head;
				in >> head >> vn->NX >> vn->NY >> vn->NZ;
				m_pic.VN.push_back(*vn);
			}
			else{//v -53.0413 158.84 -135.806 点
				istringstream in(s);
				v = new Vertex();
				string head;
				in >> head >> v->X >> v->Y >> v->Z;
				m_pic.V.push_back(*v);
			}
		}
		else if (s[0] == 'f'){//f 2443//2656 2442//2656 2444//2656 面
			for (int k = s.size() - 1; k >= 0; k--){
				if (s[k] == '/')s[k] = ' ';
			}
			istringstream in(s);
			f = new Mian();
			string head;
			in >> head;
			int i = 0;
			while (i<3)
			{//
				if (m_pic.V.size() != 0)
				{
					in >> f->V[i];
					f->V[i] -= 1;
				}
				if (m_pic.VT.size() != 0)
				{
					in >> f->T[i];
					f->T[i] -= 1;
				}
				if (m_pic.VN.size() != 0)
				{
					in >> f->N[i];
					f->N[i] -= 1;
				}
				i++;
			}
			m_pic.F.push_back(*f);
		}
	}
}
void GLDraw()
{													
	glColor3f(0.0, 1.0, 0.0);     //绿
	for (int i = 0; i<m_pic.F.size(); i++)
	{
		glBegin(GL_TRIANGLE_FAN);                            // 绘制三角形GL_TRIANGLES;GL_LINE_LOOP;GL_LINES;GL_POINTS
		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[0]].TU, m_pic.VT[m_pic.F[i].T[0]].TV);  //纹理    
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[0]].NX, m_pic.VN[m_pic.F[i].N[0]].NY, m_pic.VN[m_pic.F[i].N[0]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X , m_pic.V[m_pic.F[i].V[0]].Y , m_pic.V[m_pic.F[i].V[0]].Z );        // 上顶点

		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[1]].TU, m_pic.VT[m_pic.F[i].T[1]].TV);  //纹理
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[1]].NX, m_pic.VN[m_pic.F[i].N[1]].NY, m_pic.VN[m_pic.F[i].N[1]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[1]].X , m_pic.V[m_pic.F[i].V[1]].Y , m_pic.V[m_pic.F[i].V[1]].Z );        // 左下

		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[2]].TU, m_pic.VT[m_pic.F[i].T[2]].TV);  //纹理
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[2]].NX, m_pic.VN[m_pic.F[i].N[2]].NY, m_pic.VN[m_pic.F[i].N[2]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[2]].X , m_pic.V[m_pic.F[i].V[2]].Y , m_pic.V[m_pic.F[i].V[2]].Z );        // 右下
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
void init() {
	/*
	Vector3f temp = R.row(0);
	R.row(0) = R.row(1);
	R.row(1) = temp;
	m_DF.setR(R);*/
	ReadPIC();
	//基尔霍夫张量
	CKirchhoff m_K(m_pic);
	K = m_K.computeK();//初始时得到K矩阵
	m_DF.q = R;
	m_DF.setK(K);
	m_DF.lp=	m_DF.computelp();//计算初始的lp

	glClearColor(0.0, 0.0, 0.0, 0.0);  //背景色
	//深度测试的相关设置 
	glClearDepth(1.0);                    //设置深度缓存的初始值 
	glDepthFunc(GL_LEQUAL);           //深度测试的方法 
	glEnable(GL_DEPTH_TEST);          //启用深度测试
	GLfloat direction[] = { -3.0f, -3.4f, -8.8f, 0.0f }; // 平行光源, GL_POSITION属性的最后一个参数为0
	GLfloat ambient[] = { 0.3f, 0.3f, 0.3f, 1.0f };  // 环境强度
	GLfloat diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };  // 散射强度
	GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f }; // 镜面强度
	glLightfv(GL_LIGHT0, GL_POSITION, direction);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

	glEnable(GL_LIGHTING);   //开关:使用光
	glEnable(GL_LIGHT0);     //打开0#灯
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
//窗口重绘时的响应函数 
void display() {
	//用前面设置的背景色和深度值清除颜色缓存和深度缓存 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//初始化模型视图矩阵 
	drawScene();            //几何图元的绘制 
	glFlush();              //绘制结束，Flush 当前渲染流水线 
	glutSwapBuffers();       //交换前后缓存（只用于双缓冲的模式） 
}

void drawScene()           //绘制
{
	if (mouseRightDown ) {
		glTranslatef(0, 0, -cameraDistance*0.1);
	}
	/*
	glTranslated(m_DF.temp_deltay(0), m_DF.temp_deltay(1), m_DF.temp_deltay(2));*/
	glMultMatrixf(m_DF.GetRotationData());
	//glRotated(m_DF.delta_q.w(), m_DF.delta_q.x(), m_DF.delta_q.y(), m_DF.delta_q.z());
	GLDraw();
}
//窗口大小发生变化时的响应函数 
void reshape(int width, int height) {
	glViewport(0, 0, width, height);      //设置视窗大小 
	//设置视景体大小 
	glMatrixMode(GL_PROJECTION);
	float ratio = (float)width / height;
	glLoadIdentity();
	//gluPerspective(60, ratio, 1, 1000);
	int viewsize=25;
	
	if (width <= height)
		glOrtho(-viewsize, viewsize, -viewsize * (GLfloat)height / (GLfloat)width, viewsize * (GLfloat)height / (GLfloat)width, -10.0, 10.0);
	else
		glOrtho(-viewsize *(GLfloat)width / (GLfloat)height, viewsize*(GLfloat)width / (GLfloat)height	, -viewsize, viewsize, -10.0, 10.0);

	//glOrtho(-25, 25, -25, 25, -10, 10);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(2, 0, 0, 0, 0, 0, 0, 0, 1);//4, 0, -2,

}

void TimerFunction(int value)
{
	m_DF.nextTime();
	//设置tsfs 即下一时刻的受力情况
	/*
	glutPostRedisplay 标记当前窗口需要重新绘制。通过glutMainLoop下一次循环时，
	窗口显示将被回调以重新显示窗口的正常面板。多次调用glutPostRedisplay，
	在下一个显示回调只产生单一的重新显示回调
	*/
	glutPostRedisplay(); //标志重新绘制
	glutTimerFunc(delta_t*1000, TimerFunction, 1);
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
		cameraDistance = 0;
		cameraDistance += (y - mouseY) * 0.2f;
		mouseY = y;
	}
	glutPostRedisplay();
}
int main(int argc, char* argv[]) {
	glutInit(&argc, argv);    //GLUT 库的初始化 
	//显示模式初始化：颜色格式——GLUT_RGBA 
	//           单缓冲或双缓冲——GLUT_SINGLE 或者 GLUT_float //          是否使用深度缓存——GLUT_DEPTH  
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE| GLUT_DEPTH);
	glutInitWindowPosition(100, 100);          //窗口起始位置 
	glutInitWindowSize(800, 500);             //窗口大小 
	glutCreateWindow("UnderWaterRidgebody");       //创建窗口并指定窗口的名称 
	glutMouseFunc(mouseCB);
	glutMotionFunc(mouseMotionCB);

	init();                                 //OpenGL 的初始化设置 

	//设置消息响应的回调函数 
	glutDisplayFunc(display);               //设置窗口重绘时的响应函数 
	glutReshapeFunc(reshape);              //设置窗口大小发生变化时的响应函数 

	//定时器  每500毫秒触发一次
	glutTimerFunc(delta_t*1000, TimerFunction, 1);

	glutMainLoop();                   //消息循环：获取消息，分发消息，响应消息 
	return 0;
} //
