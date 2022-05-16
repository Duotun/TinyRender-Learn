#include "tgaimage.h"
#include <cmath>
#include <iostream>
#include "model.h"
#include "geometry.h"
#include <vector>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
Model* model = NULL;

const int width = 800;
const int height = 600;


//utilize the barycentric coordinates, cross product way, very trick
Vec3f Barycentric_Cross(Vec2i* pts, Vec2i p)
{
	Vec3f v_cross_x = Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - p.x);
	Vec3f u = cross(v_cross_x, Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - p.y));
	if (std::abs(u.z) < 1.0) return Vec3f(-1.0, 1.0, 1.0);  //z component must be -1 or 1
	//return negative values indicating out of boundary
	return Vec3f(1.0 - (u.x + u.y) / u.z, u.x / u.z, u.y / u.z);   //considering -1 case, so division
}

Vec3f Barycentric_Cross(Vec3f* pts, Vec3f p)
{
	Vec3f v_cross_x = Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - p.x);
	Vec3f u = cross(v_cross_x, Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - p.y)); 
	if (std::abs(u.z) < 1.0) return Vec3f(-1.0, 1.0, 1.0);  //z component must be -1 or 1
	//return negative values indicating out of boundary
	return Vec3f(1.0 - (u.x + u.y) / u.z, u.x / u.z, u.y / u.z);   //considering -1 case, so division
}

//(u, v, w)
Vec3f Barycentric_Efficient(Vec2i* pts, Vec2i p)
{
	//just a method from collision detection book, using dot product only
	Vec3f v0 = Vec3f(pts[1].x - pts[0].x, pts[1].y - pts[0].y, 0);
	Vec3f v1 = Vec3f(pts[2].x - pts[0].x, pts[2].y - pts[0].y, 0);
	Vec3f v2 = Vec3f(p.x - pts[0].x, p.y - pts[0].y, 0);

	float d00 = dot(v0, v0); float d01 = dot(v0,v1);
	float d11 = dot(v1, v1); float d20 = dot(v0, v0);
	float d21 = dot(v1, v1); float denom = d00 * d11 - d01 * d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0f - v - w;
	return Vec3f(u, v, w);
}
//organized method with boundary checking seperate
void triangle_v2(Vec2i pts[3], TGAImage& img, TGAColor color)  //3 is used to fix the array length
{
	// constrain with a bbox
	Vec2i bboxmin(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
	Vec2i bboxmax(-std::numeric_limits<int>::max(), -std::numeric_limits<int>::max());
	for (int i = 0; i < 3; i++)
	{
		bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
		bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));
		bboxmax.x = std::min(img.get_width() - 1, std::max(bboxmax.x, pts[i].x));
		bboxmax.y = std::min(img.get_height() - 1, std::max(bboxmax.y, pts[i].y));
	}
	Vec2i p;
	for (p.x = bboxmin.x; p.x <= bboxmax.x; p.x++)
	{
		for (p.y = bboxmin.y; p.y <= bboxmax.y; p.y++)
		{
			Vec3f bary_res = Barycentric_Cross(pts, p);
			if (bary_res.x < 0 || bary_res.y < 0 || bary_res.z < 0) continue;
			//std::cout << "Here" << std::endl;
			
			img.set(p.x, p.y, color);
		}
	}
}
//three lines and fill
void triangle(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage& img, TGAColor color)
{
	if (v0.y == v1.y && v0.y == v2.y) return; //
	//fill with the line sweeping method
	// draw line segments between the left and right sides
	if (v0.y > v1.y) std::swap(v0, v1);
	if (v0.y > v2.y) std::swap(v0, v2);
	if (v1.y > v2.y) std::swap(v1, v2);   //ascending order, v0, v1, v2

	float threshold = 0.001;  //avoid the zero division
	float height_total = (float)(v2.y - v0.y);
	float height_segment = (float)(v1.y - v0.y + 1.);

	//fill the bottom part and the upper part
	//combine the two filling parts into one loop
	for (int i = v0.y; i <= v2.y; i++)
	{
		// pay attention to the case why the v1.y == v0.y and this is possible for the triangle case!
		bool second_half = i > (v1.y - v0.y) || v1.y == v0.y;
		height_segment = !second_half ? (float)(v1.y - v0.y + 1) : (float)(v2.y - v1.y + 1);
		float alpha = (float)(i - v0.y) / height_total;
		float beta = (float)(i - (!second_half ? v0.y : v1.y)) / height_segment;
		Vec2i a = v0 + (v2 - v0) * alpha;
		Vec2i b = !second_half ? (v0 + (v1 - v0) * beta) : (v1 + (v2 - v1) * beta);

		//make sure a is the left side
		if (a.x > b.x) std::swap(a, b);
		for (int j = a.x; j <= b.x; j++)
		{
			img.set(j, i, color);
		}
	}

}

// line segment according to the position set method
void line(int x0, int y0, int x1, int y1, TGAImage& img, TGAColor color)
{
	bool steep = false;  //if too steep, change the order temporarily for calculate y precisively
	if (std::abs(x0 - x1) < std::abs(y0 -y1))
	{
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1)
	{
		std::swap(x0, x1);
		std::swap(y0, y1);  //reverse the order
	}
	int dx = x1 - x0;
	int dy = y1 - y0;
	int ierror = std::abs(dy) * 2;   //trick for the round up for y
	int error = 0;
	int y = y0;

	for (int x = x0; x <= x1; x++)
	{
		if (steep)
		{
			img.set(y, x, color);   //transpose back
		}
		else
			img.set(x, y, color);

		error += ierror;
		if (error > dx)
		{
			y += (y1 > y0 ? 1 : -1);
			error -= dx * 2;
		}
	}

}

void WireFrameObj(int argc, char ** argv, TGAImage & img)
{
	if (2 == argc)
	{
		model = new Model(argv[1]);
	}
	else
	{
		model = new Model("obj/african_head.obj");
	}

	for (int i = 0; i < model->nfaces(); i++)
	{
		std::vector<int> face = model->face(i);
		for (int j = 0; j < 3; j++) {   //build from the triangular face

			Vec3f v0 = model->vert(face[j]);
			Vec3f v1 = model->vert(face[(j + 1) % 3]);

			int x0 = (v0.x + 1.) * width / 2.;    // here +1 to make the vertices range between 0-1
			int x1 = (v1.x + 1.) * width / 2.;
			int y0 = (v0.y + 1.) * height / 2.;
			int y1 = (v1.y + 1.) * height / 2.;

			line(x0, y0, x1, y1, img, white);
		}
	}
}

Vec3f world2screen(Vec3f v) {
	return Vec3f(int((v.x + 1.0) * width / 2.0 + 0.5), int((v.y + 1.0) * height / 2.0 + 0.5), v.z);
}

void triangle_v3(Vec3f* pts, float* zbuffer, TGAImage& img, TGAColor color)  //3 is used to fix the array length
{
	// constrain with a bbox
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	for (int i = 0; i < 3; i++)
	{
		bboxmin.x = std::max(0.f, std::min(bboxmin.x, pts[i].x));
		bboxmin.y = std::max(0.f, std::min(bboxmin.y, pts[i].y));
		bboxmax.x = std::min(img.get_width() - 1.f, std::max(bboxmax.x, pts[i].x));
		bboxmax.y = std::min(img.get_height() - 1.f, std::max(bboxmax.y, pts[i].y));
	}
	Vec3f p;
	
	// adding z buffer checking here, using barycentric as well
	for (p.x = bboxmin.x; p.x <= bboxmax.x; p.x++)
	{
		for (p.y = bboxmin.y; p.y <= bboxmax.y; p.y++)
		{
			Vec3f bary_res = Barycentric_Cross(pts, p);
			if (bary_res.x < 0 || bary_res.y < 0 || bary_res.z < 0) continue;
			//std::cout << "Here" << std::endl;
			//obtain z values in the barycentric 
			p.z = 0;
			for (int i = 0; i < 3; i++)  p.z += pts[i][2] * bary_res[i];
			
			if (zbuffer[int(p.x + p.y * width)] < p.z)
			{
				//std::cout << "Here?" << color.bgra<< std::endl;
				zbuffer[int(p.x + p.y * width)] = p.z;  //update and draw
				img.set(p.x, p.y, color);
			}
			
		}
	}
}
//just random color
void FlatShadingObj(int argc, char** argv, TGAImage& img)
{
	if (2 == argc)
	{
		model = new Model(argv[1]);
	}
	else
	{
		model = new Model("obj/african_head.obj");
	}
	
	Vec2i screen_coords[3];  //very simple way to project
	for (int i = 0; i < model->nfaces(); i++)
	{
		std::vector<int> face = model->face(i);
		
		for (int j = 0; j < 3; j++) {   //build from the triangular face

			Vec3f worldPos = model->vert(face[j]);
			screen_coords[j] = Vec2i((worldPos.x + 1.0) * width / 2.0, (worldPos.y + 1.0) * height / 2.0);
			//assume no depth here
		}
		triangle_v2(screen_coords, img, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
	}
}

//very simple lambert lighting model
void DiffuseLightObj(int argc, char** argv, TGAImage& img, Vec3f light_dir)
{
	if (2 == argc)
	{
		model = new Model(argv[1]);
	}
	else
	{
		model = new Model("obj/african_head.obj");
	}

	Vec2i screen_coords[3];  //very simple way to project
	Vec3f worldPos[3];
	float light_intensity = 0.0;
	Vec3f normal;
	
	for (int i = 0; i < model->nfaces(); i++)
	{
		std::vector<int> face = model->face(i);

		for (int j = 0; j < 3; j++) {   //build from the triangular face

			worldPos[j] = model->vert(face[j]);
			screen_coords[j] = Vec2i((worldPos[j].x + 1.0) * width / 2.0, (worldPos[j].y + 1.0) * height / 2.0);

		}
		normal = cross((worldPos[2] - worldPos[0]), worldPos[1] - worldPos[0]);
		normal.normalize();
		light_intensity = dot(normal, light_dir);
		if(light_intensity >0)
			triangle_v2(screen_coords, img, TGAColor(light_intensity *255, light_intensity * 255, light_intensity * 255,255));
		//treat negative values as black color
	}
}

//diffuse light obj with zbuffer hidden faces removal
void DiffuseLightObj(int argc, char** argv, float * zbuffer, TGAImage& img, Vec3f light_dir)
{
	if (2 == argc)
	{
		model = new Model(argv[1]);
	}
	else
	{
		model = new Model("obj/african_head.obj");
	}

	Vec3f pts[3];  //very simple way to project
	Vec3f worldPos[3];
	float light_intensity = 0.0;
	Vec3f normal;
	
	for (int i = 0; i < model->nfaces(); i++)
	{
		std::vector<int> face = model->face(i);

		for (int j = 0; j < 3; j++) {   //build from the triangular face

			worldPos[j] = model->vert(face[j]);
			pts[j] = world2screen(worldPos[j]);
		}

		normal = cross((worldPos[2] - worldPos[0]), worldPos[1] - worldPos[0]);
		normal.normalize();
		light_intensity = dot(normal, light_dir);
		
		if (light_intensity > 0)
			triangle_v3(pts,  zbuffer, img, TGAColor(light_intensity * 255, light_intensity * 255, light_intensity * 255, 255));
		//treat negative values as black color
	}
}


int main(int argc, char** argv) {
	TGAImage image(width, height, TGAImage::RGB);
	//image.set(52, 41, red);  //red dot
	
	Vec2i v0[3] = { Vec2i(10, 70), Vec2i(50, 160), Vec2i(70, 80) };
	Vec2i v1[3] = { Vec2i(180, 70), Vec2i(150, 1), Vec2i(70, 180) };
	Vec2i v2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };

	//initialize zbuffer
	float *zbuffer = new float[width * height];
	for (int i = 0; i < width * height; i++) zbuffer[i] = -std::numeric_limits<float>::max();

	//std::cout << "Here?" << std::endl;
	//DiffuseLightObj(argc, argv, zbuffer, image, Vec3f(0, 0, -1));
	
	//FlatShadingObj(argc, argv, image);
	//DiffuseLightObj(argc, argv, image, Vec3f(0, 0, -1));
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

