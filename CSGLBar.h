#pragma once

// added by Jieling, 30.05.2017
// draw bars connecting neighboring cells

#ifndef CS_GLTOOLS_CSGLBAR_H
#define CS_GLTOOLS_CSGLBAR_H

#include "../CSGLObject.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

class CSGLBar : public CSGLObject
{
public:
	CSGLBar(Vector3f * endpoint1, Vector3f * endpoint2, const ARGBColor &color, const double &radius);

	void draw();
	//void setQuality(int slice, int stacks);

protected:
	Vector3f *p1;
	Vector3f *p2;
	double * mpRadius;
	int mSlices;
	int mStacks;
};

#endif
