#include "CSGLBar.h"

#include <cmath>

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"


#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
# define glutSolidCylinder(radio, altura, slices, stacks) gluCylinder(gluNewQuadric(), radio, radio, altura, slices, stacks)
#else
# include <GL/freeglut.h>
# include <GL/glut.h>
#endif


CSGLBar::CSGLBar(Vector3f * endpoint1, Vector3f * endpoint2, const ARGBColor &/*color*/, const double &/*radius*/)
	: CSGLObject(endpoint1)
{
	this->mpPosition = endpoint1;

	this->p1 = endpoint1;
	this->p2 = endpoint2;

	this->mSlices = 20;
	this->mStacks = 20;
}

void
CSGLBar::draw()
{
	glPushMatrix();

	glColor4d(1.0, 0.1, 0.1, 1.0);
	
	double tmp_rot[3];
	
	double vx = this->p2->x - this->p1->x; // x2 - x1
	double vy = this->p2->y - this->p1->y; // y2 - y1
	double vz = this->p2->z - this->p1->z; // z2 - z1
	
	if (vz == 0)
		vz = .0001;

	double tmp_r = sqrt(vx*vx + vy*vy + vz*vz);
	double ax = 57.2957795*acos(vz / tmp_r);
	if (vz < 0.0)
		ax = -ax;
	double rx = -vy*vz;
	double ry = vx*vz;

	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;

	glTranslatef(this->p1->x,
				 this->p1->y,
				 this->p1->z);

	glRotatef(tmp_rot[0], tmp_rot[1], tmp_rot[2], 0.0);

	glutSolidCylinder(0.05, tmp_r, 5, 5);//radius 0.005

	glPopMatrix();
}