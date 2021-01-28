/***********************************************
*    version 11/11/11
*
**********************************************/

#ifndef DEFS_H
#define DEFS_H

typedef double Vec[3];
typedef double Light;	// first we will consider monolight

typedef struct t_tri
{
	int n[3];	// number of the three points
	Vec Nor;	// normal vector for the triangle
	double I_intersect;
	double I_absorb_tri;
	int count_hit_tri;
}Triangle;

typedef struct t_pt
{
	Vec coordinate;
	Vec normal;
}Point;

typedef struct t_sur
{
	int pnum;
	int tnum;
	Point *vertex;
	Triangle *tri;
	double z_buf;
}Surface;

typedef struct t_obj
{
	Surface *S;
	double refr_index_o;
	double refr_index_i;
	double absorb;	// index of absorption
	double scatter;	// index of scattering
	double I_absorb;
	struct t_obj* subobj;
	struct t_obj* nextobj;
	struct t_obj* belongobj;
	int celltype;//0->leaf 1/2->spongy_o/i 3->s_chl 4->s_vacuole
			//5/6->palisade_o/i 7->p_chl 8->p_vacuole
			//9/10->epidermis_o/i 11/12->epidermis_d_o/i
	int cellname;// 0 means leaf surface
	int chlname;// 0 means cell surface
	int vacuolename;
	int count_hit;
}Object;

typedef struct t_isect
{
	Vec inter_pt;
	Vec inter_nor;
	Object *inter_obj;
	int inter_tri;
	double t;
}Isect;

typedef struct t_ray
{
	Vec P;	// origin point
	Vec D;	// direction
	Light I;	// light intensity or energy in each beam
	Object *P_obj;
	int P_tri;
	int interflag;
	int level;
	double R;
	double T;
	Object *trace_obj;
	double refr_index;
}Ray;

#define Infinity 1.7E+308

#define PI 3.141592653589793

#define discardI 0.0001
Object *curobj,*curbelong,*curobj_refl,*curobj_trans;
int TIR;
Isect *intersect;

double I_discard;
double I_discard_Rf;
double I_discard_Tr;

//#define chl_concentration_MS (2.352941e4) //g m^-3
//#define chl_concentration_BS (2.352941e4) //g m^-3
//#define SAC_chlab 3.27 //m^2 g^-1
//#define SAC_water 0.41
//#define SAC_cellwall (1.47e6)*(4.38e-4)
//#define SAC_cyto 0.41
double chl_con_MS;
double chl_con_BS;
double SAC_chlab;
//double SAC_chlab_BS;
double SAC_water;
double SAC_cellwall;
double SAC_cyto;

#define air_refr_index 1
#define wall_refr_index 1.415
#define cyto_refr_index 1.353
#define chlo_refr_index 1.511
#define vacu_refr_index 1.333

#define DBL_EPSILON (1.835973856023010e-11)
#define DBL_EPSILON_trisect (1e-100)

int num_chl_hit;
double debugI;

#define flag_debug 1
#define flag_debug_precal 0
#define flag_debug_printf 0
#define flag_warning 0
#define flag_randseed 1

double RT_debug[500][25][512][3];//cutnum_x*cutnum_y*max_buff*{x,y,z}
int count_RT_debug[500][25];
int ray_i,ray_j;

//#define xmax 9.324500000000001e-05
//#define xmin 0.0
//#define zmax 5.437164103349100e-06
//#define zmin 0.0
//#define ymax 7.561999999998612e-05
//#define ymin -0.1e-6
//
//#define ms_num 24
//int count_ms;
//#define ms_max_chl_num 1
//int ms_chl_num[ms_num];
//#define num_nonMS 4
//
//Object *p_cell_ms[ms_num];
//Object *p_chl_ms[ms_num][ms_max_chl_num];
//Object *p_vac_ms[ms_num];
//Object *p_leaf;
//Object *p_cell_ns[num_nonMS];
//
//#endif
