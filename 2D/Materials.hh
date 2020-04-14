#ifndef __MATERIALS_HH__
#define __MATERIALS_HH__

#include <GL/gl.h>

class Material {
private:
	float amb[4];
	float dif[4];
	float spec[4];
	float emis[4];
	float shin[1];
public:
	Material(
		float a_r, float a_g, float a_b, float a_a,
		float d_r, float d_g, float d_b, float d_a,
		float s_r, float s_g, float s_b, float s_a,
		float e_r, float e_g, float e_b, float e_a,
		float sh)
	{
		amb[0] = a_r; amb[1] = a_g; amb[2] = a_b; amb[3] = a_a;
		dif[0] = d_r; dif[1] = d_g; dif[2] = d_b; dif[3] = d_a;
		spec[0] = s_r; spec[1] = s_g; spec[2] = s_b; spec[3] = s_a;
		emis[0] = e_r; emis[1] = e_g; emis[2] = e_b; emis[3] = e_a;
		shin[0] = sh;
	}
		
	void set( void) {
		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);
		
	};
};

Material Material_UnselectedNode(
	0.6, 0.0, 0.0, 1.0,
	1.0, 0.0, 0.0, 1.0,
	1.0, 0.0, 0.0, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_SelectedNode(
	0.4, 0.4, 0.0, 1.0,
	0.75, 0.75, 0.0, 1.0,
	0.75, 0.75, 0.0, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_UnselectedErrorNode(
	0.6, 0.0, 0.6, 1.0,
	1.0, 0.0, 1.0, 1.0,
	1.0, 0.0, 1.0, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_SelectedErrorNode(
	1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_BluePlastic(
	0.0, 0.0, 0.4, 1.0,
	0.0, 0.0, 0.5, 1.0,
	0.1, 0.1, 0.1, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_YellowPlastic(
	0.4, 0.4, 0.0, 1.0,
	0.5, 0.5, 0.0, 1.0,
	0.1, 0.1, 0.1, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_GreenPlastic(
	0.0, 0.4, 0.0, 1.0,
	0.0, 0.4, 0.0, 1.0,
	0.1, 0.1, 0.1, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);

Material Material_RedPlastic(
	0.4, 0.0, 0.0, 1.0,
	0.4, 0.0, 0.0, 1.0,
	0.1, 0.1, 0.1, 1.0,
	0.0, 0.0, 0.0, 1.0,
	0.0);


#endif
