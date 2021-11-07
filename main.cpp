#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/freeglut.h>
#endif

using namespace std;

const float pi = 3.14159;

int winid;
int winw = 400;
int winh = 400;
float cpos[] = {0.0, 1.5, 0.0};
float cangle[] = {0, 0};
float l1pos[] = {0, 5, 0, 1};
float l2pos[] = {-5, 10, -5, 1};
float amb[] = {0.2, 0.2, 0.2, 1.0};
float dif[] = {0.8, 0.8, 0.8, 1.0};
float spec[] = {0.85, 0.85, 0.85, 1.0};
float ambh[] = {0.3, 0.0, 0.0, 1.0};
float difh[] = {0.8, 0.0, 0.0, 1.0};
float spech[] = {1.0, 0.0, 0.0, 1.0};
float l1amb[] = {0.2, 0.2, 0.2, 1.0};
float l1dif[] = {0.8, 0.8, 0.8, 1.0};
float l1spec[] = {0.85, 0.85, 0.85, 1.0};
float c1amb[] = {0.0, 0.0, 0.2, 1.0};
float c1dif[] = {0.0, 0.0, 0.8, 1.0};
float c1spec[] = {0.0, 0.0, 0.85, 1.0};
float c2amb[] = {0.0, 0.2, 0.2, 1.0};
float c2dif[] = {0.0, 0.8, 0.8, 1.0};
float c2spec[] = {0.0, 0.85, 0.85, 1.0};
float c3amb[] = {0.0, 0.2, 0.0, 1.0};
float c3dif[] = {0.0, 0.8, 0.0, 1.0};
float c3spec[] = {0.0, 0.85, 0.0, 1.0};
float c4amb[] = {0.2, 0.2, 0.2, 1.0};
float c4dif[] = {0.8, 0.8, 0.8, 1.0};
float c4spec[] = {0.85, 0.85, 0.85, 1.0};
int currentobject = -1;
float t[3145728];

bool keyw, keys, keya, keyd, keyq, keye, keyt, keyg, keyf, keyh, keyr, keyy, keyi, keyk, keyj, keyl, keyu, keyo, keyz, keyx, keyc, keyv, keyb, keyn, keyW, keyS, keyA, keyD, keyQ, keyE, keyT, keyG, keyF, keyH, keyR, keyY, keyI, keyK, keyJ, keyL, keyU, keyO, keyZ, keyX, keyC, keyV, keyB, keyN = false;
GLuint tex[2];
string texnames[2] = {"grass.ppm", "blocks1.ppm"};

class Point{
	float x, y, z;
	float norm[3];
	public:
	void set(float a, float b, float c){
		x = a;
		y = b;
		z = c;
	}
	float getx(){
		return x;
	}
	float gety(){
		return y;
	}
	float getz(){
		return z;
	}
	void setnormx(float a){
		norm[0] = a;
	}
	float getnormx(){
		return norm[0];
	}
	void setnormy(float a){
		norm[1] = a;
	}
	float getnormy(){
		return norm[1];
	}
	void setnormz(float a){
		norm[2] = a;
	}
	float getnormz(){
		return norm[2];
	}
};

class Triangle{
	Point* e1;
	Point* e2;
	Point* e3;
	float uv[3][2];
	int texid;
	float norm[3];
	public:
	void sete1(Point* a){
		e1 = a;
	}
	void sete2(Point* a){
		e2 = a;
	}
	void sete3(Point* a){
		e3 = a;
	}
	Point* gete1(){
		return e1;
	}
	Point* gete2(){
		return e2;
	}
	Point* gete3(){
		return e3;
	}
	void setuv(float a, int i, int j){
		uv[i][j] = a;
	}
	float getuv(int i, int j){
		return uv[i][j];
	}
	void settexid(string a){
		int i = 0;
		while(a.compare(texnames[i]) != 0){
			i++;
		}
		texid = i;
	}
	int gettexid(){
		return texid;
	}
	void setnormx(float a){
		norm[0] = a;
	}
	float getnormx(){
		return norm[0];
	}
	void setnormy(float a){
		norm[1] = a;
	}
	float getnormy(){
		return norm[1];
	}
	void setnormz(float a){
		norm[2] = a;
	}
	float getnormz(){
		return norm[2];
	}
};

class Plane{
	float a, b, c, d;
	public:
	Plane(float h, float i, float j, float k);
	float geta(){
		return a;
	}
	float getb(){
		return b;
	}
	float getc(){
		return c;
	}
	float getd(){
		return d;
	}
};

Plane::Plane(float h, float i, float j, float k){
	a = h;
	b = i;
	c = j;
	d = k;
}

class Model{
	vector<Point> points;
	vector<Triangle> tri;
	void getcoords(string k, float* u1, float* u2, float* u3){
		int h1;
		int h2;
		string a;
		string b;
		string c;
		h1 = k.find_first_of(",");
		h2 = k.find_last_of(",");
		a = k.substr(0, h1);
		b = k.substr(h1 + 2, h2 - (h1 + 2));
		c = k.substr(h2 + 2, k.length() - (h2 + 2));
		*u1 = stof(a);
		*u2 = stof(b);
		*u3 = stof(c);
	}
	void getcoords2(string k, float* u1, float* u2){
		int h;
		string a;
		string b;
		h = k.find_first_of(",");
		a = k.substr(0, h);
		b = k.substr(h + 2, k.length() - (h + 2));
		*u1 = stof(a);
		*u2 = stof(b);
	}
	public:
	void setmodel(string file){
		int h1;
		int h2;
		int vecpos;
		int r;
		float u[3];
		float u2[2];
		bool isinvec = false;
		string a;
		string c[3];
		string uv[3];
		string t;
		ifstream f(file);
		getline(f, a);
		r = stoi(a);
		points.reserve(r * 3);
		for(int i = 0; i < r; i++){
			getline(f, a);
			h1 = a.find_first_of(";");
			h2 = a.find_last_of(";");
			c[0] = a.substr(0, h1);
			c[1] = a.substr(h1 + 2, h2 - (h1 + 2));
			c[2] = a.substr(h2 + 2, a.length() - (h2 + 2));
			getline(f, a);
			h1 = a.find_first_of(";");
			h2 = a.find_last_of(";");
			uv[0] = a.substr(0, h1);
			uv[1] = a.substr(h1 + 2, h2 - (h1 + 2));
			uv[2] = a.substr(h2 + 2, a.length() - (h2 + 2));
			getline(f, t);
			tri.push_back(Triangle());
			for(int j = 0; j < 3; j++){
				getcoords(c[j], &u[0], &u[1], &u[2]);
				getcoords2(uv[j], &u2[0], &u2[1]);
				for(int k = 0; k < points.size(); k++){
					if(u[0] == points[k].getx() && u[1] == points[k].gety() && u[2] == points[k].getz()){
						isinvec = true;
						vecpos = k;
					}
				}
				if(isinvec == false){
					points.push_back(Point());
					points.back().set(u[0], u[1], u[2]);
					if(j == 0){
						tri.back().sete1(&points.back());
					}else if(j == 1){
						tri.back().sete2(&points.back());
					} else {
						tri.back().sete3(&points.back());
					}
				} else {
					if(j == 0){
						tri.back().sete1(&points[vecpos]);
					}else if(j == 1){
						tri.back().sete2(&points[vecpos]);
					} else {
						tri.back().sete3(&points[vecpos]);
					}
				}
				tri.back().setuv(u2[0], j, 0);
				tri.back().setuv(u2[1], j, 1);
				isinvec = false;
			}
			tri.back().settexid(t);
		}
		f.close();
	}
	Triangle* gettri(int i){
		return &(tri[i]);
	}
	void pushtri(Triangle t){
		tri.push_back(t);
	}
	int tril(){
		return tri.size();
	}
	Point* getpoint(int i){
		return &(points[i]);
	}
	void pushpoint(Point p){
		points.push_back(p);
	}
	int pointsl(){
		return points.size();
	}
	void setnorms(){
		float n1[3];
		float n2[3];
		float h1, h2, h3, h4;
		for(int i = 0; i < tri.size(); i++){
			n1[0] = (*tri[i].gete2()).getx() - (*tri[i].gete1()).getx();
			n1[1] = (*tri[i].gete2()).gety() - (*tri[i].gete1()).gety();
			n1[2] = (*tri[i].gete2()).getz() - (*tri[i].gete1()).getz();
			n2[0] = (*tri[i].gete3()).getx() - (*tri[i].gete1()).getx();
			n2[1] = (*tri[i].gete3()).gety() - (*tri[i].gete1()).gety();
			n2[2] = (*tri[i].gete3()).getz() - (*tri[i].gete1()).getz();
			h1 = (n1[1] * n2[2]) - (n1[2] * n2[1]);
			h2 = (n1[2] * n2[0]) - (n1[0] * n2[2]);
			h3 = (n1[0] * n2[1]) - (n1[1] * n2[0]);
			h4 = sqrt(pow(h1, 2.0) + pow(h2, 2.0) + pow(h3, 2.0));
			tri[i].setnormx(h1 / h4);
			tri[i].setnormy(h2 / h4);
			tri[i].setnormz(h3 / h4);
		}
		for(int i = 0; i < points.size(); i++){
			h1 = 0;
			h2 = 0;
			h3 = 0;
			for(int j = 0; j < tri.size(); j++){
				if(tri[j].gete1() == &points[i] || tri[j].gete2() == &points[i] || tri[j].gete3() == &points[i]){
					h1 += tri[j].getnormx();
					h2 += tri[j].getnormy();
					h3 += tri[j].getnormz();
				}
			}
			h4 = sqrt(pow(h1, 2.0) + pow(h2, 2.0) + pow(h3, 2.0));
			points[i].setnormx(h1 / h4);
			points[i].setnormy(h2 / h4);
			points[i].setnormz(h3 / h4);
		}
	}
};

class Object{
	float pos[3] = {0.0, 0.0, 0.0};
	float angle[3] = {0.0, 0.0, 0.0};
	float scale[3] = {1.0, 1.0, 1.0};
	int type;
	Plane* hitbox[6];
	Model* m;
	public:
	void setposx(float a){
		pos[0] = a;
	}
	float getposx(){
		return pos[0];
	}
	void setposy(float a){
		pos[1] = a;
	}
	float getposy(){
		return pos[1];
	}
	void setposz(float a){
		pos[2] = a;
	}
	float getposz(){
		return pos[2];
	}
	void setanglex(float a){
		angle[0] = a;
	}
	float getanglex(){
		return angle[0];
	}
	void setangley(float a){
		angle[1] = a;
	}
	float getangley(){
		return angle[1];
	}
	void setanglez(float a){
		angle[2] = a;
	}
	float getanglez(){
		return angle[2];
	}
	void setscalex(float a){
		scale[0] = a;
	}
	float getscalex(){
		return scale[0];
	}
	void setscaley(float a){
		scale[1] = a;
	}
	float getscaley(){
		return scale[1];
	}
	void setscalez(float a){
		scale[2] = a;
	}
	float getscalez(){
		return scale[2];
	}
	void sethitbox(float a, float b, float c, float d, int e){
		hitbox[e] = new Plane(a, b, c, d);
	}
	Plane* gethitbox(int a){
		return hitbox[a];
	}
	void setmodel(Model* a){
		m = a;
	}
	Model* getmodel(){
		return m;
	}
	void settype(int a){
		type = a;
	}
	int gettype(){
		return type;
	}
};

vector<Model> models;
vector<Object> objects;
vector<int> intersections;
vector<float> intpos;

void translate(float x, float y, float z, float* u, float* v){
	v[0] = u[0] + x;
	v[1] = u[1] + y;
	v[2] = u[2] + z;
}

void rotatex(float a, float* u, float* v){
	v[0] = u[0];
	v[1] = (cos(a) * u[1]) + (-(sin(a)) * u[2]);
	v[2] = (sin(a) * u[1]) + (cos(a) * u[2]);
}

void rotatey(float a, float* u, float* v){
	v[0] = (cos(a) * u[0]) + (sin(a) * u[2]);
	v[1] = u[1];
	v[2] = (-(sin(a)) * u[0]) + (cos(a) * u[2]);
}

void rotatez(float a, float* u, float* v){
	v[0] = (cos(a) * u[0]) + (-(sin(a)) * u[1]);
	v[1] = (sin(a) * u[0]) + (cos(a) * u[1]);
	v[2] = u[2];
}

void scale(float x, float y, float z, float* u, float* v){
	v[0] = u[0] * x;
	v[1] = u[1] * y;
	v[2] = u[2] * z;
}

float planeintersect(Plane* a, float* b, float* c){
	float h = ((*a).geta() * c[0]) + ((*a).getb() * c[1]) + ((*a).getc() * c[2]);
	float u;
	float t = -1.0;
	if(h != 0.0){
		u = -(((*a).geta() * b[0]) + ((*a).getb() * b[1]) + ((*a).getc() * b[2]) + (*a).getd());
		t = u / h;
	}
	return t;
}

void intersect(int a, float* b, float* c, float* h, bool* r){
	bool b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12;
	b1 = -(*objects[a].gethitbox(3)).getd() < b[1] + (h[0] * c[1]) && b[1] + (h[0] * c[1]) < (*objects[a].gethitbox(2)).getd();
	b2 = -(*objects[a].gethitbox(3)).getd() < b[1] + (h[1] * c[1]) && b[1] + (h[1] * c[1]) < (*objects[a].gethitbox(2)).getd();
	b3 = -(*objects[a].gethitbox(5)).getd() < b[2] + (h[0] * c[2]) && b[2] + (h[0] * c[2]) < (*objects[a].gethitbox(4)).getd();
	b4 = -(*objects[a].gethitbox(5)).getd() < b[2] + (h[1] * c[2]) && b[2] + (h[1] * c[2]) < (*objects[a].gethitbox(4)).getd();
	b5 = -(*objects[a].gethitbox(1)).getd() < b[0] + (h[2] * c[0]) && b[0] + (h[2] * c[0]) < (*objects[a].gethitbox(0)).getd();
	b6 = -(*objects[a].gethitbox(0)).getd() < b[0] + (h[3] * c[0]) && b[0] + (h[3] * c[0]) < (*objects[a].gethitbox(0)).getd();
	b7 = -(*objects[a].gethitbox(5)).getd() < b[2] + (h[2] * c[2]) && b[2] + (h[2] * c[2]) < (*objects[a].gethitbox(4)).getd();
	b8 = -(*objects[a].gethitbox(5)).getd() < b[2] + (h[3] * c[2]) && b[2] + (h[3] * c[2]) < (*objects[a].gethitbox(4)).getd();
	b9 = -(*objects[a].gethitbox(1)).getd() < b[0] + (h[4] * c[0]) && b[0] + (h[4] * c[0]) < (*objects[a].gethitbox(0)).getd();
	b10 = -(*objects[a].gethitbox(1)).getd() < b[0] + (h[5] * c[0]) && b[0] + (h[5] * c[0]) < (*objects[a].gethitbox(0)).getd();
	b11 = -(*objects[a].gethitbox(3)).getd() < b[1] + (h[4] * c[1]) && b[1] + (h[4] * c[1]) < (*objects[a].gethitbox(2)).getd();
	b12 = -(*objects[a].gethitbox(3)).getd() < b[1] + (h[5] * c[1]) && b[1] + (h[5] * c[1]) < (*objects[a].gethitbox(2)).getd();
	r[0] = b1 && b3;
	r[1] = b2 && b4;
	r[2] = b5 && b7;
	r[3] = b6 && b8;
	r[4] = b9 && b11;
	r[5] = b10 && b12;
}

void boxintersect(int a, float* b, float* c){
	float h[6];
	bool r[6];
	h[0] = planeintersect(objects[a].gethitbox(0), b, c);
	h[1] = planeintersect(objects[a].gethitbox(1), b, c);
	h[2] = planeintersect(objects[a].gethitbox(2), b, c);
	h[3] = planeintersect(objects[a].gethitbox(3), b, c);
	h[4] = planeintersect(objects[a].gethitbox(4), b, c);
	h[5] = planeintersect(objects[a].gethitbox(5), b, c);
	intersect(a, b, c, h, r);
	if(r[0] || r[1] || r[2] || r[3] || r[4] || r[5] == true){
		intersections.push_back(a);
		float smallest = 500.00;
		for(int i = 0; i < 6; i++){
			if(r[i] == true){
				if(h[i] < smallest){
					smallest = h[i];
				}
			}
		}
		intpos.push_back(smallest);
	}
}

int ray(float x, float y){
	GLdouble px, py, pz;
	GLdouble m[16], p[16];
	GLint v[4];
	glPushMatrix();
	glRotatef(-cangle[1], 1, 0, 0);
	glRotatef(-cangle[0], 0, 1, 0);
	glTranslatef(-cpos[0], -cpos[1], -cpos[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, m);
	glGetDoublev(GL_PROJECTION_MATRIX, p);
	glGetIntegerv(GL_VIEWPORT, v);
	gluUnProject(x, y, 0.0, m, p, v, &px, &py, &pz);
	float pr1[3];
	float pr2[3];
	pr1[0] = px;
	pr1[1] = py;
	pr1[2] = pz;
	gluUnProject(x, y, 1.0, m, p, v, &px, &py, &pz);
	float h = sqrt(pow(px, 2.0) + pow(py, 2.0) + pow(pz, 2.0));
	pr2[0] = px;
	pr2[1] = py;
	pr2[2] = pz;
	glPopMatrix();
	float temp1[3];
	float temp2[3];
	float temp3[3];
	float temp4[3];
	bool bl;
	intersections.clear();
	intpos.clear();
	for(int i = 0; i < objects.size(); i++){
		temp1[0] = pr1[0];
		temp1[1] = pr1[1];
		temp1[2] = pr1[2];
		temp3[0] = pr2[0] / h;
		temp3[1] = pr2[1] / h;
		temp3[2] = pr2[2] / h;
		translate(-objects[i].getposx(), -objects[i].getposy(), -objects[i].getposz(), temp1, temp2);
		rotatey(-objects[i].getangley(), temp2, temp1);
		rotatey(-objects[i].getangley(), temp3, temp4);
		rotatex(-objects[i].getanglex(), temp1, temp2);
		rotatex(-objects[i].getanglex(), temp4, temp3);
		rotatez(-objects[i].getanglez(), temp2, temp1);
		rotatez(-objects[i].getanglez(), temp3, temp4);
		scale(1 / objects[i].getscalex(), 1 / objects[i].getscaley(), 1 / objects[i].getscalez(), temp1, temp2);
		boxintersect(i, temp2, temp4);
	}
	float smallest = 500.0;
	int r;
	if(intersections.size() > 0){
		for(int i = 0; i < intersections.size(); i++){
			if(intpos[i] < smallest){
				smallest = intpos[i];
				r = i;
			}
		}
		return intersections[r];
	} else {
		return -1;
	}
}

void loadtexture(string k){
	int h;
	string a;
	string b;
	ifstream f(k);
	getline(f, a);
	getline(f, a);
	getline(f, a);
	h = a.find_first_of(",");
	b = a.substr(0, h);
	h = stoi(b);
	getline(f, b);
	for(int i = 0; i < h * h * 3; i++){
		getline(f, b);
		t[i] = (float)stoi(b) / 255.0;
	}
	f.close();
}

void reset(){
	objects.erase(objects.begin() + 1, objects.end());
}

void save(){
	ofstream f;
	f.open("save");
	f << to_string(objects.size()) + "\n";
	for(int i = 1; i < objects.size(); i++){
		f << to_string(objects[i].getposx()) + "\n";
		f << to_string(objects[i].getposy()) + "\n";
		f << to_string(objects[i].getposz()) + "\n";
		f << to_string(objects[i].getanglex()) + "\n";
		f << to_string(objects[i].getangley()) + "\n";
		f << to_string(objects[i].getanglez()) + "\n";
		f << to_string(objects[i].getscalex()) + "\n";
		f << to_string(objects[i].getscaley()) + "\n";
		f << to_string(objects[i].getscalez()) + "\n";
		f << to_string(objects[i].gettype()) + "\n";
	}
	f.close();
}

void load(){
	reset();
	int a;
	string b;
	ifstream f;
	f.open("save");
	getline(f, b);
	a = stoi(b);
	for(int i = 1; i < a; i++){
		objects.push_back(Object());
		objects.back().setmodel(&models[0]);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 0);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 1);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 2);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 4);
		objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 5);
		getline(f, b);
		objects.back().setposx(stof(b));
		getline(f, b);
		objects.back().setposy(stof(b));
		getline(f, b);
		objects.back().setposz(stof(b));
		getline(f, b);
		objects.back().setanglex(stof(b));
		getline(f, b);
		objects.back().setangley(stof(b));
		getline(f, b);
		objects.back().setanglez(stof(b));
		getline(f, b);
		objects.back().setscalex(stof(b));
		getline(f, b);
		objects.back().setscaley(stof(b));
		getline(f, b);
		objects.back().setscalez(stof(b));
		getline(f, b);
		objects.back().settype(stoi(b));
	}
	f.close();
}
		

void init(void){
	glClearColor(0, 0, 0, 0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_TEXTURE_2D);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, amb);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, dif);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, l1amb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, l1dif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, l1spec);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, 1.0, 1.0, 500.0);
	glMatrixMode(GL_MODELVIEW);
	glGenTextures(2, tex);
	for(int i = 0; i < 2; i++){
		loadtexture(texnames[i]);
		glBindTexture(GL_TEXTURE_2D, tex[i]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_FLOAT, t);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}
	models.push_back(Model());
	models.back().setmodel("model");
	models.back().setnorms();
	models.push_back(Model());
	models.back().setmodel("room");
	models.back().setnorms();
	objects.push_back(Object());
	objects.back().setmodel(&models[1]);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 0);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 1);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 2);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 3);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 4);
	objects.back().sethitbox(0.0, 0.0, 0.0, 0.0, 5);
}

void reshape(int w, int h){
	winw = w;
	winh = h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)w / (float)h, 1.0, 500.0);
	glMatrixMode(GL_MODELVIEW);
	glViewport(0, 0, w, h);
}

void keyboard(unsigned char key, int xin, int yin){
	switch(key){
		case 'w':
			keyw = true;
			break;
		case 's':
			keys = true;
			break;
		case 'a':
			keya = true;
			break;
		case 'd':
			keyd = true;
			break;
		case 'q':
			keyq = true;
			break;
		case 'e':
			keye = true;
			break;
		case 't':
			keyt = true;
			break;
		case 'g':
			keyg = true;
			break;
		case 'f':
			keyf = true;
			break;
		case 'h':
			keyh = true;
			break;
		case 'r':
			keyr = true;
			break;
		case 'y':
			keyy = true;
			break;
		case 'i':
			keyi = true;
			break;
		case 'k':
			keyk = true;
			break;
		case 'j':
			keyj = true;
			break;
		case 'l':
			keyl = true;
			break;
		case 'u':
			keyu = true;
			break;
		case 'o':
			keyo = true;
			break;
		case 'z':
			keyz = true;
			break;
		case 'x':
			keyx = true;
			break;
		case 'c':
			keyc = true;
			break;
		case 'v':
			keyv = true;
			break;
		case 'b':
			keyb = true;
			break;
		case 'n':
			keyn = true;
			break;
		case 'W':
			keyW = true;
			break;
		case 'S':
			keyS = true;
			break;
		case 'A':
			keyA = true;
			break;
		case 'D':
			keyD = true;
			break;
		case 'Q':
			keyQ = true;
			break;
		case 'E':
			keyE = true;
			break;
		case 'T':
			keyT = true;
			break;
		case 'G':
			keyG = true;
			break;
		case 'F':
			keyF = true;
			break;
		case 'H':
			keyH = true;
			break;
		case 'R':
			keyR = true;
			break;
		case 'Y':
			keyY = true;
			break;
		case 'I':
			keyI = true;
			break;
		case 'K':
			keyK = true;
			break;
		case 'J':
			keyJ = true;
			break;
		case 'L':
			keyL = true;
			break;
		case 'U':
			keyU = true;
			break;
		case 'O':
			keyO = true;
			break;
		case 'Z':
			keyZ = true;
			break;
		case 'X':
			keyX = true;
			break;
		case 'C':
			keyC = true;
			break;
		case 'V':
			keyV = true;
			break;
		case 'B':
			keyB = true;
			break;
		case 'N':
			keyN = true;
			break;
		case '1':
			save();
			break;
		case '2':
			load();
			break;
		case 'm':
			reset();
			break;
	}
}

void keyboardup(unsigned char key, int xin, int yin){
	switch(key){
		case 'w':
			keyw = false;
			break;
		case 's':
			keys = false;
			break;
		case 'a':
			keya = false;
			break;
		case 'd':
			keyd = false;
			break;
		case 'q':
			keyq = false;
			break;
		case 'e':
			keye = false;
			break;
		case 't':
			keyt = false;
			break;
		case 'g':
			keyg = false;
			break;
		case 'f':
			keyf = false;
			break;
		case 'h':
			keyh = false;
			break;
		case 'r':
			keyr = false;
			break;
		case 'y':
			keyy = false;
			break;
		case 'i':
			keyi = false;
			break;
		case 'k':
			keyk = false;
			break;
		case 'j':
			keyj = false;
			break;
		case 'l':
			keyl = false;
			break;
		case 'u':
			keyu = false;
			break;
		case 'o':
			keyo = false;
			break;
		case 'z':
			keyz = false;
			break;
		case 'x':
			keyx = false;
			break;
		case 'c':
			keyc = false;
			break;
		case 'v':
			keyv = false;
			break;
		case 'b':
			keyb = false;
			break;
		case 'n':
			keyn = false;
			break;
		case 'W':
			keyW = false;
			break;
		case 'S':
			keyS = false;
			break;
		case 'A':
			keyA = false;
			break;
		case 'D':
			keyD = false;
			break;
		case 'Q':
			keyQ = false;
			break;
		case 'E':
			keyE = false;
			break;
		case 'T':
			keyT = false;
			break;
		case 'G':
			keyG = false;
			break;
		case 'F':
			keyF = false;
			break;
		case 'H':
			keyH = false;
			break;
		case 'R':
			keyR = false;
			break;
		case 'Y':
			keyY = false;
			break;
		case 'I':
			keyI = false;
			break;
		case 'K':
			keyK = false;
			break;
		case 'J':
			keyJ = false;
			break;
		case 'L':
			keyL = false;
			break;
		case 'U':
			keyU = false;
			break;
		case 'O':
			keyO = false;
			break;
		case 'Z':
			keyZ = false;
			break;
		case 'X':
			keyX = false;
			break;
		case 'C':
			keyC = false;
			break;
		case 'V':
			keyV = false;
			break;
		case 'B':
			keyB = false;
			break;
		case 'N':
			keyN = false;
			break;
	}
}

void mouse(int b, int s, int x, int y){
	if(b == GLUT_LEFT_BUTTON && s == GLUT_DOWN){
		currentobject = ray((float)x, (float)(winh - y));
	}
	if(b == GLUT_RIGHT_BUTTON && s == GLUT_DOWN){
		int a = ray((float)x, (float)(winh - y));
		if(a != -1){
			objects.erase(objects.begin() + a);
		}
	}
}

void update(int){
	if(keyw == true){
		cpos[0] = cpos[0] - (sin(cangle[0] * (pi / 180.0)) / 2);
		cpos[2] = cpos[2] - (cos(cangle[0] * (pi / 180.0)) / 2);
	}
	if(keys == true){
		cpos[0] = cpos[0] + (sin(cangle[0] * (pi / 180.0)) / 2);
		cpos[2] = cpos[2] + (cos(cangle[0] * (pi / 180.0)) / 2);
	}
	if(keya == true){
		cpos[0] = cpos[0] - (cos(cangle[0] * (pi / 180.0)) / 2);
		cpos[2] = cpos[2] + (sin(cangle[0] * (pi / 180.0)) / 2);
	}
	if(keyd == true){
		cpos[0] = cpos[0] + (cos(cangle[0] * (pi / 180.0)) / 2);
		cpos[2] = cpos[2] - (sin(cangle[0] * (pi / 180.0)) / 2);
	}
	if(keyq == true){
		cpos[1] = cpos[1] - 0.5;
	}
	if(keye == true){
		cpos[1] = cpos[1] + 0.5;
	}
	if(keyi == true){
		if(cangle[1] < 90){
			cangle[1] = cangle[1] + 2;
		}
	}
	if(keyk == true){
		if(cangle[1] > -90){
			cangle[1] = cangle[1] - 2;
		}
	}
	if(keyj == true){
		cangle[0] = cangle[0] + 2;
	}
	if(keyl == true){
		cangle[0] = cangle[0] - 2;
	}
	if(keyt == true){
		objects.push_back(Object());
		objects.back().settype(0);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		objects.back().setmodel(&models[0]);
		keyt = false;
	}
	if(keyg == true){
		objects.push_back(Object());
		objects.back().settype(1);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		keyg = false;
	}
	if(keyf == true){
		objects.push_back(Object());
		objects.back().settype(2);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		keyf = false;
	}
	if(keyh == true){
		objects.push_back(Object());
		objects.back().settype(3);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		keyh = false;
	}
	if(keyr == true){
		objects.push_back(Object());
		objects.back().settype(4);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		keyr = false;
	}
	if(keyy == true){
		objects.push_back(Object());
		objects.back().settype(5);
		objects.back().sethitbox(1.0, 0.0, 0.0, 0.5, 0);
		objects.back().sethitbox(-1.0, 0.0, 0.0, 0.5, 1);
		objects.back().sethitbox(0.0, 1.0, 0.0, 1.0, 2);
		objects.back().sethitbox(0.0, -1.0, 0.0, 0.0, 3);
		objects.back().sethitbox(0.0, 0.0, 1.0, 0.5, 4);
		objects.back().sethitbox(0.0, 0.0, -1.0, 0.5, 5);
		keyy = false;
	}
	if(currentobject != -1){
		if(keyW == true){
			objects[currentobject].setposz(objects[currentobject].getposz() - 0.1);
		}
		if(keyS == true){
			objects[currentobject].setposz(objects[currentobject].getposz() + 0.1);
		}
		if(keyA == true){
			objects[currentobject].setposx(objects[currentobject].getposx() - 0.1);
		}
		if(keyD == true){
			objects[currentobject].setposx(objects[currentobject].getposx() + 0.1);
		}
		if(keyE == true){
			objects[currentobject].setposy(objects[currentobject].getposy() + 0.1);
		}
		if(keyQ == true){
			objects[currentobject].setposy(objects[currentobject].getposy() - 0.1);
		}
		if(keyT == true){
			objects[currentobject].setscalez(objects[currentobject].getscalez() + 0.1);
		}
		if(keyG == true){
			objects[currentobject].setscalez(objects[currentobject].getscalez() - 0.1);
		}
		if(keyF == true){
			objects[currentobject].setscalex(objects[currentobject].getscalex() - 0.1);
		}
		if(keyH == true){
			objects[currentobject].setscalex(objects[currentobject].getscalex() + 0.1);
		}
		if(keyY == true){
			objects[currentobject].setscaley(objects[currentobject].getscaley() + 0.1);
		}
		if(keyR == true){
			objects[currentobject].setscaley(objects[currentobject].getscaley() - 0.1);
		}
		if(keyI == true){
			objects[currentobject].setanglex(objects[currentobject].getanglex() - 1.0);
		}
		if(keyK == true){
			objects[currentobject].setanglex(objects[currentobject].getanglex() + 1.0);
		}
		if(keyJ == true){
			objects[currentobject].setanglez(objects[currentobject].getanglez() + 1.0);
		}
		if(keyL == true){
			objects[currentobject].setanglez(objects[currentobject].getanglez() - 1.0);
		}
		if(keyO == true){
			objects[currentobject].setangley(objects[currentobject].getangley() + 1.0);
		}
		if(keyU == true){
			objects[currentobject].setangley(objects[currentobject].getangley() - 1.0);
		}
	}
	if(keyz == true){
		l1pos[0] = l1pos[0] + 0.1;
	}
	if(keyx == true){
		l1pos[0] = l1pos[0] - 0.1;
	}
	if(keyc == true){
		l1pos[1] = l1pos[1] + 0.1;
	}
	if(keyv == true){
		l1pos[1] = l1pos[1] - 0.1;
	}
	if(keyb == true){
		l1pos[2] = l1pos[2] + 0.1;
	}
	if(keyn == true){
		l1pos[2] = l1pos[2] - 0.1;
	}
	if(keyZ == true){
		l2pos[0] = l2pos[0] + 0.1;
	}
	if(keyX == true){
		l2pos[0] = l2pos[0] - 0.1;
	}
	if(keyC == true){
		l2pos[1] = l2pos[1] + 0.1;
	}
	if(keyV == true){
		l2pos[1] = l2pos[1] - 0.1;
	}
	if(keyB == true){
		l2pos[2] = l2pos[2] + 0.1;
	}
	if(keyN == true){
		l2pos[2] = l2pos[2] - 0.1;
	}
	if(cangle[0] < 0.0){
		cangle[0] += 360.0;
	}else if(cangle[0] >= 360.0){
		cangle[0] -= 360.0;
	}
	glutPostRedisplay();
	glutTimerFunc(1000.0 / 120.0, update, 0);
}

void drawobject2(int i){
	for(int j = 0; j < (*objects[i].getmodel()).tril(); j++){
		glBindTexture(GL_TEXTURE_2D, tex[(*(*objects[i].getmodel()).gettri(j)).gettexid()]);
		glBegin(GL_TRIANGLES);
		glTexCoord2f((*(*objects[i].getmodel()).gettri(j)).getuv(0, 0), (*(*objects[i].getmodel()).gettri(j)).getuv(0, 1));
		glNormal3f((*(*(*objects[i].getmodel()).gettri(j)).gete1()).getnormx(), (*(*(*objects[i].getmodel()).gettri(j)).gete1()).getnormy(), (*(*(*objects[i].getmodel()).gettri(j)).gete1()).getnormz());
		glVertex3f((*(*(*objects[i].getmodel()).gettri(j)).gete1()).getx(), (*(*(*objects[i].getmodel()).gettri(j)).gete1()).gety(), (*(*(*objects[i].getmodel()).gettri(j)).gete1()).getz());
		glTexCoord2f((*(*objects[i].getmodel()).gettri(j)).getuv(1, 0), (*(*objects[i].getmodel()).gettri(j)).getuv(1, 1));
		glNormal3f((*(*(*objects[i].getmodel()).gettri(j)).gete2()).getnormx(), (*(*(*objects[i].getmodel()).gettri(j)).gete2()).getnormy(), (*(*(*objects[i].getmodel()).gettri(j)).gete2()).getnormz());
		glVertex3f((*(*(*objects[i].getmodel()).gettri(j)).gete2()).getx(), (*(*(*objects[i].getmodel()).gettri(j)).gete2()).gety(), (*(*(*objects[i].getmodel()).gettri(j)).gete2()).getz());
		glTexCoord2f((*(*objects[i].getmodel()).gettri(j)).getuv(2, 0), (*(*objects[i].getmodel()).gettri(j)).getuv(2, 1));
		glNormal3f((*(*(*objects[i].getmodel()).gettri(j)).gete3()).getnormx(), (*(*(*objects[i].getmodel()).gettri(j)).gete3()).getnormy(), (*(*(*objects[i].getmodel()).gettri(j)).gete3()).getnormz());
		glVertex3f((*(*(*objects[i].getmodel()).gettri(j)).gete3()).getx(), (*(*(*objects[i].getmodel()).gettri(j)).gete3()).gety(), (*(*(*objects[i].getmodel()).gettri(j)).gete3()).getz());
		glEnd();
	}
} 

void drawobject(int i){
	if(objects[i].gettype() == 0){
		drawobject2(i);
	}else if(objects[i].gettype() == 1){
		glDisable(GL_TEXTURE_2D);
		glMaterialfv(GL_FRONT, GL_AMBIENT, c1amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, c1dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, c1spec);
		glPushMatrix();
		glTranslatef(0.0, 0.5, 0.0);
		glutSolidSphere(0.5, 10, 10);
		glPopMatrix();
		glEnable(GL_TEXTURE_2D);
	}else if(objects[i].gettype() == 2){
		glDisable(GL_TEXTURE_2D);
		glMaterialfv(GL_FRONT, GL_AMBIENT, c2amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, c2dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, c2spec);
		glPushMatrix();
		glTranslatef(0.0, 0.5, 0.0);
		glutSolidCube(1.0);
		glPopMatrix();
		glEnable(GL_TEXTURE_2D);
	}else if(objects[i].gettype() == 3){
		glDisable(GL_TEXTURE_2D);
		glMaterialfv(GL_FRONT, GL_AMBIENT, c3amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, c3dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, c3spec);
		glPushMatrix();
		glRotatef(-90.0, 1.0, 0.0, 0.0);
		glutSolidCylinder(0.5, 1.0, 10, 10);
		glPopMatrix();
		glEnable(GL_TEXTURE_2D);
	}else if(objects[i].gettype() == 4){
		glDisable(GL_TEXTURE_2D);
		glMaterialfv(GL_FRONT, GL_AMBIENT, c1amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, c1dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, c1spec);
		glPushMatrix();
		glRotatef(-90.0, 1.0, 0.0, 0.0);
		glutSolidCone(0.5, 1.0, 10, 10);
		glPopMatrix();
		glEnable(GL_TEXTURE_2D);
	} else {
		glDisable(GL_TEXTURE_2D);
		glMaterialfv(GL_FRONT, GL_AMBIENT, c4amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, c4dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, c4spec);
		glPushMatrix();
		glTranslatef(0.0, 0.5, 0.0);
		glutSolidTeapot(0.5);
		glPopMatrix();
		glEnable(GL_TEXTURE_2D);
	}
	glMaterialfv(GL_FRONT, GL_AMBIENT, l1amb);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, l1dif);
	glMaterialfv(GL_FRONT, GL_SPECULAR, l1spec);
	if(i == currentobject){
		float jr, kr, lr;
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		for(int j = 0; j < 2; j++){
			if(j % 2 == 0){
				jr = 1.0;
			} else {
				jr = -1.0;
			}
			for(int k = 2; k < 4; k++){
				if(k % 2 == 0){
					kr = 1.0;
				} else {
					kr = -1.0;
				}
				for(int l = 4; l < 6; l++){
					if(l % 2 == 0){
						lr = 1.0;
					} else {
						lr = -1.0;
					}
					glVertex3f((*objects[i].gethitbox(j)).getd() * jr, (*objects[i].gethitbox(k)).getd() * kr, (*objects[i].gethitbox(l)).getd() * lr);
				}
			}
		}
		for(int j = 4; j < 6; j++){
			if(j % 2 == 0){
				jr = 1.0;
			} else {
				jr = -1.0;
			}
			for(int k = 0; k < 2; k++){
				if(k % 2 == 0){
					kr = 1.0;
				} else {
					kr = -1.0;
				}
				for(int l = 2; l < 4; l++){
					if(l % 2 == 0){
						lr = 1.0;
					} else {
						lr = -1.0;
					}
					glVertex3f((*objects[i].gethitbox(k)).getd() * kr, (*objects[i].gethitbox(l)).getd() * lr, (*objects[i].gethitbox(j)).getd() * jr);
				}
			}
		}
		for(int j = 2; j < 4; j++){
			if(j % 2 == 0){
				jr = 1.0;
			} else {
				jr = -1.0;
			}
			for(int k = 4; k < 6; k++){
				if(k % 2 == 0){
					kr = 1.0;
				} else {
					kr = -1.0;
				}
				for(int l = 0; l < 2; l++){
					if(l % 2 == 0){
						lr = 1.0;
					} else {
						lr = -1.0;
					}
					glVertex3f((*objects[i].gethitbox(l)).getd() * lr, (*objects[i].gethitbox(j)).getd() * jr, (*objects[i].gethitbox(k)).getd() * kr);
				}
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
	}
}

void display(void){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glRotatef(-cangle[1], 1, 0, 0);
	glRotatef(-cangle[0], 0, 1, 0);
	glTranslatef(-cpos[0], -cpos[1], -cpos[2]);
	glLightfv(GL_LIGHT0, GL_POSITION, l1pos);
	glPushMatrix();
	glTranslatef(l1pos[0], l1pos[1], l1pos[2]);
	glDisable(GL_TEXTURE_2D);
	glutSolidSphere(0.1, 10, 10);
	glEnable(GL_TEXTURE_2D);
	glPopMatrix();
	glLightfv(GL_LIGHT1, GL_POSITION, l2pos);
	glPushMatrix();
	glTranslatef(l2pos[0], l2pos[1], l2pos[2]);
	glDisable(GL_TEXTURE_2D);
	glutSolidSphere(0.1, 10, 10);
	glEnable(GL_TEXTURE_2D);
	glPopMatrix();
	for(int i = 0; i < objects.size(); i++){
		glPushMatrix();
		glTranslatef(objects[i].getposx(), objects[i].getposy(), objects[i].getposz());
		glRotatef(objects[i].getangley(), 0.0, 1.0, 0.0);
		glRotatef(objects[i].getanglex(), 1.0, 0.0, 0.0);
		glRotatef(objects[i].getanglez(), 0.0, 0.0, 1.0);
		glScalef(objects[i].getscalex(), objects[i].getscaley(), objects[i].getscalez());
		drawobject(i);
		glPopMatrix();
	}
	glPopMatrix();
	glutSwapBuffers();
}

int main(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(50, 50);
	winid = glutCreateWindow("Something");
	glutReshapeFunc(reshape);
	glutIgnoreKeyRepeat(1);
	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboardup);
	glutMouseFunc(mouse);
	glutTimerFunc(1000.0 / 120.0, update, 0);
	glutDisplayFunc(display);
	init();
	printf("Controls:\nqweasd move\nrtyfgh create objects or scale\nuiojkl rotate\nzxcvbn move lights\n12 save and load\n\nPressing caps lock switches between controlling the camera and the selected object, or light 1 and light 2. Pressing m clears the scene.\n");
	glutMainLoop();
	return(0);
}
