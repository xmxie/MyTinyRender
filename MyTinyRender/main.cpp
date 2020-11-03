
#include "tgaimage.h"
//#include "base.h"
#include <time.h>
#include "model.h"
#include <fstream>
#include <string>
#include <vector>

#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
using namespace std;
using namespace myGeometry;
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);


const int width = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void triangle(Vec2i to, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color);
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, float* zbuffer, TGAImage &image, TGAColor color);
void bruteline(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);

int main(int argc, char** argv) {
	TGAImage image(width, height, TGAImage::RGB);
	//image.read_tga_file("output.tga");
	//image.show_tga_file();
	////lesson 1 draw a line
	/*
	image.set(52, 41, red);
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	line(99, 0, 0, 99, image, blue);
	*/


	const char* modelPos = "E:/code/tinyrenderer/obj/african_head/african_head.obj";
	Model *model=new Model(modelPos);

	//draw a model with lines
	/*
	for (int i = 0; i < model->nfaces(); i++) {
		std::vector<int> face = model->face(i);
		for (int j = 0; j < 3; j++) {
			Vec3f v0 = model->vert(face[j]);
			Vec3f v1 = model->vert(face[(j + 1) % 3]);
			int x0 = (v0.x + 1.)*width / 2.;
			int y0 = (v0.y + 1.)*height / 2.;
			int x1 = (v1.x + 1.)*width / 2.;
			int y1 = (v1.y + 1.)*height / 2.;
			line(x0, y0, x1, y1, image, white);
		}
	}
	delete model;
	*/

	//rasterize
	//Vec2i p1(1, 1), p2(58, 46), p3(25, 79);
	//triangle(p1, p2, p3,image,red);
	Vec3f light_dir(0, 0, -1);
	
	//with zBuffer
	//float* zbuffer = new float[width*height];
	//memset(zbuffer, -std::numeric_limits<float>::max(), sizeof(float)*width*height);

	float *zbuffer = new float[width*height];
	for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

	for (int i = 0; i < model->nfaces(); i++) {
		std::vector<int> face = model->face(i);
		Vec2i screen_coords[3];
		Vec3f world_coords[3];
		for (int j = 0; j < 3; j++) {
			Vec3f v = model->vert(face[j]);
			screen_coords[j] = Vec2i((v.x + 1.)*width / 2., (v.y + 1.)*height / 2.);
			world_coords[j] = v;
		}
		Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
		n.normalize();
		float intensity = n * light_dir;
		//without zBuffer
		//if (intensity > 0) {
		//	triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
		//}

		////with zBuffer
		if (intensity > 0) {
			triangle(screen_coords[0], screen_coords[1], screen_coords[2],zbuffer,image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
		}

	}
	



	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.show_tga_file();
	image.write_tga_file("output.tga");
	return 0;
}

inline int sign(int x) {
	if (x == 0) return 0;
	return x > 0 ? 1 : -1;
}

void bruteline(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
	for (float t = 0.; t < 1.; t += .01) {
		int x = x0 + (x1 - x0)*t;
		int y = y0 + (y1 - y0)*t;
		image.set(x, y, color);
	}
}

//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//		steep = true;
//	}
//	if (x0 > x1) {
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//	}
//
//	for (int x = x0; x <= x1; x++) {
//		float t = (x - x0) / (float)(x1 - x0);
//		int y = y0 * (1. - t) + y1 * t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
	//Bresenham's line algorithm
	int dx = abs(x1 - x0), dy = abs(y1 - y0);
	int s1 = sign(x1 - x0), s2 = sign(y1 - y0);//处理负斜率
	int x=x0, y=y0;
	int p = dx - 2 * dy;
	bool interchange = false;
	//处理斜率大于1的情况
	if (dx < dy) {
		int tmp = dx; dx = dy; dy = tmp;
		interchange = true;//镜像
	}
	for (int i = 0; i < dx; i++) {
		image.set(x, y, color);
		if (p <= 0) {
			if (interchange) {
				x += s1;
			}
			else {
				y += s2;
			}
			p += 2 * dx;
		}
		if (interchange) {
			y += s2;
		}
		else {
			x += s1;
		}
		p -= 2 * dy;
	}

}

//按照y坐标将顶点从小到大排序并返回最小的y坐标
int sortVertices3(Vec2i& t0, Vec2i& t1, Vec2i& t2) {
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	return t0.y;
}
void lineRaster(Vec2i& t0, Vec2i& t1,vector<vector<int>>& v,const int& offsetY) {
	if (t0.y == t1.y) {
		v[t0.y - offsetY].push_back(min(t0.x, t1.x));
		v[t0.y - offsetY].push_back(max(t0.x, t1.x));
	}
	int dx = abs(t0.x - t1.x), dy = abs(t0.y - t1.y);
	int s1 = sign(t1.x - t0.x), s2 = sign(t1.y - t0.y);//处理负斜率
	int x = t0.x, y = t0.y;
	int p = dx - 2 * dy;
	bool interchange = false;
	//处理斜率大于1的情况
	if (dx < dy) {
		int tmp = dx; dx = dy; dy = tmp;
		interchange = true;//镜像
	}
	for (int i = 0; i <= dx; i++) {
		//在数组内记录两个端点
		if (v[y-offsetY].size() < 2) {
			v[y-offsetY].push_back(x);
			if (v[y-offsetY].size() == 2) {
				sort(v[y-offsetY].begin(), v[y-offsetY].end());
			}
		}
		else {
			if (v[y-offsetY][0] > x) v[y-offsetY][0] = x;
			if (v[y-offsetY][1] < x) v[y-offsetY][1] = x;
		}
		if (p <= 0) {
			if (interchange) {
				x += s1;
			}
			else {
				y += s2;
			}
			p += 2 * dx;
		}
		if (interchange) {
			y += s2;
		}
		else {
			x += s1;
		}
		p -= 2 * dy;
	}
}
void mytriangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
	//1.Sort vertices of the triangle by their y - coordinates;
	int offsetY=sortVertices3(t0, t1, t2);
	//2.Rasterize simultaneously the left and the right sides of the triangle;
	vector<vector<int>> v(t2.y-t0.y+1);//数组用来存储同一行下的两个端点
	lineRaster(t0, t1, v, offsetY);
	lineRaster(t1, t2, v, offsetY);
	lineRaster(t0, t2, v, offsetY);
	for (int i = 0; i < v.size(); i++) {
		int x0 = v[i].front();
		int x1 = v[i].back();
		if (i == 76) {
			int i = 1;
		}
		for(int x=x0;x<=x1;x++){
			image.set(x, i + offsetY, color);
		}
	}
	//3.Draw a horizontal line segment between the left and the right boundary points.
}

//P=(1-u-v)A+uB+vC=A+uAB+vAC,返回(1-u-v,u,v)
template<class T> Vec3f barycentric(const Vec2i& v0,const Vec2i& v1,const Vec2i& v2 , const T& P) {
	Vec3f v = cross(Vec3f((float)(v2[0] - v0[0]), (float)(v1[0] - v0[0]), (float)(v0[0] - P[0])),
					Vec3f((float)(v2[1] - v0[1]), (float)(v1[1] - v0[1]), (float)(v0[1] - P[1])));
	/* 
	   只有三点连线为直线的情况下v.z为0，此时输入的三点不构成三角形   
	*/
	if (abs(v.z) < 1e-2) return Vec3f(-1.0f, 1.0f, 1.0f);
	//除v.z相当于标准化，同时也是对正负号进行处理
	return Vec3f(1.0f-(v.x+v.y)/v.z,v.y/v.z,v.x/v.z);
	
}
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
	if (t0.y == t1.y && t0.y == t2.y) return; // I dont care about degenerate triangles 
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here 
		Vec2i A = t0 + (t2 - t0)*alpha;
		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;
		if (A.x > B.x) std::swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			image.set(j, t0.y + i, color); // attention, due to int casts t0.y+i != A.y 
		}
	}
}
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, float* zbuffer,TGAImage &image, TGAColor color) {
	std::vector<Vec2i> pts = { std::move(t0),std::move(t1),std::move(t2) };
	Vec2f bboxMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxMin[j] = std::max(0.f, std::min(bboxMin[j], (float)pts[i][j]));
			bboxMax[j] = std::min(clamp[j], std::max(bboxMax[j], (float)pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x = bboxMin.x; P.x < bboxMax.x; P.x++) {
		for (P.y = bboxMin.y; P.y < bboxMax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			P.z = 0;
			for (int i = 0; i < 3; i++) P.z += pts[i][2] * bc_screen[i];
			//zbuffer缓存决定是否绘制
			if (zbuffer[int(P.x + P.y*width)] < P.z) {
				zbuffer[int(P.x + P.y*width)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
		
	}
	/*Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::max(0.f, std::min(bboxmin[j], (float)pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], (float)pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			P.z = 0;
			for (int i = 0; i < 3; i++) P.z += pts[i][2] * bc_screen[i];
			if (zbuffer[int(P.x + P.y*width)] < P.z) {
				zbuffer[int(P.x + P.y*width)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
	}*/
}