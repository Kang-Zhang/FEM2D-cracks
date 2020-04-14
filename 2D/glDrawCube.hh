#ifndef __GL_DRAW_CUBE__
#define __GL_DRAW_CUBE__

void glDrawCube( double cx, double cy, double cz, double d)
// ----------------------------------------------------------------------
// draw a cube centered at cx, cy, cx with edge size d
// ----------------------------------------------------------------------
{
	struct Dtype {
		double x, y, z;
		void draw( void) { glVertex3d( x, y, z); }
	} D[8] = 
		{
			{ cx - d/2, cy - d/2, cz - d/2 } ,
			{ cx + d/2, cy - d/2, cz - d/2 } ,
			{ cx + d/2, cy + d/2, cz - d/2 } ,
			{ cx - d/2, cy + d/2, cz - d/2 } ,
			{ cx - d/2, cy - d/2, cz + d/2 } ,
			{ cx + d/2, cy - d/2, cz + d/2 } ,
			{ cx + d/2, cy + d/2, cz + d/2 } ,
			{ cx - d/2, cy + d/2, cz + d/2 } 
		};
	
	glBegin( GL_QUADS);
	// top face:
	glNormal3d( 0, 0, 1);
	D[4].draw(); D[5].draw(); D[6].draw(); D[7].draw();
	// bottom face:
	glNormal3d( 0, 0, -1);
	D[3].draw(); D[2].draw(); D[1].draw(); D[0].draw();
	// front face:
	glNormal3d( 0, -1, 0);
	D[0].draw(); D[1].draw(); D[5].draw(); D[4].draw();
	// back face:
	glNormal3d( 0, 1, 0);
	D[2].draw(); D[3].draw(); D[7].draw(); D[6].draw();
	// right face:
	glNormal3d( 1, 0, 0);
	D[1].draw(); D[2].draw(); D[6].draw(); D[5].draw();
	// left face:
	glNormal3d( -1, 0, 0);
	D[3].draw(); D[0].draw(); D[4].draw(); D[7].draw();
	glEnd();
}

#endif
