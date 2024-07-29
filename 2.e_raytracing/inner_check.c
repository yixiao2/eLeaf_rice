/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

/***********************************************
*    version 2021 May
*    Check meshed chloroplasts are all inside their corresponding MS
*    do this after geo_export, before copy *.ply to /MS/ and /nonMS/ and start ray tracing
*    modified from main_all.c
**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "Defs.h"	// other head files
#include "Global.h"
#include "Trace.h"
#include "precalculate.h"

double chl_con_MS, chl_con_BS;
double SAC_water, SAC_chlab, SAC_cyto, SAC_air;
FILE *fout_file4plot, *fout_absorbevents;
Object *p_leaf;
Object *p_cell_ms[msall_num];
Object *p_chl_ms[msall_num][ms_max_chl_num];
Object *p_vac_ms[msall_num];
Object *p_cell_ns[nonms_num];
double I_discard, I_discard_Rf, I_discard_Tr;
int num_chl_hit, ray_i, ray_j;
double debugI;
int count_ms, count_nonms;
int ms_chl_num[msall_num];
double ms_chl_con[msall_num];

int main(int argc, char **argv)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * *
	 *		Part1. read in the object and precalculate the normal vector
	 * * * * * * * * * * * * * * * * * * * * * * * * * */
	time_t start=time(NULL);
	
	FILE *fp;
	char cell[25];
	char chl[22];
	char vacuole[23];
	char str_num[3];

	int i,j,loop1;
	int dim;
	
	Object *curobj;
	Point *temp_p=NULL;
	Triangle *temp_t=NULL;	
	int vertex_num,tri_num;
	int tmp_bsc_num,tmp_msc_num;
	//char tmp_buf[128];
	
	//////////////////////read in count_chl///////////////
	fp=fopen(*++argv,"r");
	fscanf(fp,"%d %d\n",&tmp_bsc_num,&tmp_msc_num);
	for(loop1=1;loop1<=(msall_num);loop1++)
	{
		fscanf(fp,"%d %lf\n",&(ms_chl_num[loop1-1]),&(ms_chl_con[loop1-1]));
	}
	fclose(fp);
	
//////////////////////////read in MS///////////////////////	
	printf("read in MS\n");
	count_ms=0;	
	for(loop1=1;loop1<=(msall_num);loop1++)
	{	   
		count_ms++;
		/********************
		*    cell
		********************/
		strcpy(cell,"MS/ms");
		sprintf(str_num,"%d",loop1);
		strcat(cell,str_num);
		strcat(cell,".ply");
		
		if((fp=fopen(cell,"r"))==NULL)
		{
			printf("error on open cell: %s\n",cell);
		}
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ms[count_ms-1]=curobj;
				
		/********************
		*    chl
		********************/
		for(i=1;i<=ms_chl_num[count_ms-1];i++)
		{
			strcpy(chl,"MS/ms");
			sprintf(str_num,"%d",loop1);
			strcat(chl,str_num);
			sprintf(str_num,"%d",i);
			strcat(chl,"c");
			strcat(chl,str_num);
			strcat(chl,".ply");
			//printf("chl:%s\n",chl);

			if((fp=fopen(chl,"r"))==NULL)
			{
				printf("error on open %s\n",chl);
			}

			fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
			//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
			//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
			//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
			//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
			//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
			//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
			//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
			//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

			curobj=(Object *)malloc(sizeof(Object));
			curobj->S=(Surface *)malloc(sizeof(Surface));
			curobj->S->pnum=vertex_num;
			curobj->S->tnum=tri_num;
			curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
			curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
			temp_p=curobj->S->vertex;
			temp_t=curobj->S->tri;
			for(j=0;j<vertex_num;j++)
			{
				fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);
				temp_p+=1;
			}
			for(j=0;j<tri_num;j++)
			{
				fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
				temp_t+=1;
			}
			fclose(fp);
			precalculate(curobj);
			p_chl_ms[count_ms-1][i-1]=curobj;
		}
		
		/***********************
		*    vacuole
		***********************/
		strcpy(vacuole,"MS/ms");
		sprintf(str_num,"%d",loop1);
		strcat(vacuole,str_num);
		strcat(vacuole,"v.ply");
		
		if((fp=fopen(vacuole,"r"))==NULL)
		{
			printf("error on open vacuole: %s\n",vacuole);
		}
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_vac_ms[count_ms-1]=curobj;
	}
	printf("\n");

//////////////////////////read in leaf///////////////////////
	printf("read in the leaf\n");
	if((fp=fopen("leaf","r"))==NULL)
	{
		printf("error on open file leaf\n");
	}
	fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
	//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
	//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
	//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
	//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
	//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
	//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
	//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
	//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

	curobj=(Object *)malloc(sizeof(Object));
	curobj->S=(Surface *)malloc(sizeof(Surface));
	curobj->S->pnum=vertex_num;
	curobj->S->tnum=tri_num;
	curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
	curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
	temp_p=curobj->S->vertex;
	temp_t=curobj->S->tri;
	for(j=0;j<vertex_num;j++)
	{
		fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);
		temp_p+=1;
	}
	
	for(j=0;j<tri_num;j++)
	{
		fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
		temp_t+=1;
	}
	fclose(fp);
	precalculate(curobj);
	p_leaf=curobj;
	
	////////////////////////////read non MS///////////////////////
	printf("read in nonMS cell\n");
	count_nonms=0;	
	//EPL_l
	for(loop1=1;loop1<=nonms_num;loop1++)//upper and lower eppidermis, bulliform cells
	{
		count_nonms++;
		//////////read in cell/////////
		strcpy(cell,"nonMS/ns");
		sprintf(str_num,"%d",loop1);
		strcat(cell,str_num);		
		strcat(cell,".ply");
		//printf("cell:%s\n",cell);

		if((fp=fopen(cell,"r"))==NULL)
		{
			printf("error on open %s\n",cell);
		}		
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ns[count_nonms-1]=curobj;
		printf("finish read in %s\n",cell);
	}
	printf("\n");
		
	printf("finish read in \n");
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part2. initialize the object
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	///////////////////////leaf///////////////// 
	curobj=p_leaf;
	//curobj->refr_index_o=0.0001;//maybe this can already help mirror rays back?
	//curobj->refr_index_i=air_refr_index;
	//curobj->absorb=0.0;
	//curobj->I_absorb=0.0;
	curobj->subobj=p_cell_ms[0];
	
	curobj->nextobj=NULL;
	curobj->belongobj=NULL;
	curobj->celltype=0;
	curobj->cellname=0;
	curobj->chlname=0;
	curobj->vacuolename=0;
	curobj->count_hit=0;
	///////////////////////ms///////////////// 
	for(loop1=1;loop1<=msall_num;loop1++)
	{
		//////// p_cell_ms/////////
		curobj=p_cell_ms[loop1-1];
		//curobj->refr_index_o=air_refr_index;
		//curobj->refr_index_i=cyto_refr_index;//cell wall is considered specifically in trace.c
		
		//curobj->absorb=SAC_cyto;
		//curobj->I_absorb=0.0;
		curobj->subobj=p_chl_ms[loop1-1][0];
		if(loop1!=msall_num)
		{
			curobj->nextobj=p_cell_ms[loop1];
		}
		else
		{
			curobj->nextobj=p_cell_ns[0];
		}
		curobj->belongobj=p_leaf;
		curobj->celltype=1;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=0;
		curobj->count_hit=0;
		
		////////////////p_chl_ms////////////////
		for(i=1;i<=ms_chl_num[loop1-1];i++)
		{
			curobj=p_chl_ms[loop1-1][i-1];
			//curobj->refr_index_o=cyto_refr_index; // refraction index
			//curobj->refr_index_i=chlo_refr_index;
		
			//curobj->absorb=SAC_chlab*ms_chl_con[loop1-1]; // same concentration for all chloroplasts
			//curobj->I_absorb=0.0;
			curobj->subobj=NULL;
			if(i==ms_chl_num[loop1-1])
				curobj->nextobj=p_vac_ms[loop1-1];
			else
				curobj->nextobj=p_chl_ms[loop1-1][i];
			curobj->belongobj=p_cell_ms[loop1-1];
			curobj->celltype=2;
			curobj->cellname=loop1;
			curobj->chlname=i;
			curobj->vacuolename=0;
			curobj->count_hit=0;
		}
		////////////////p_vac_ms////////////////
		curobj=p_vac_ms[loop1-1];
		//curobj->refr_index_o=cyto_refr_index; // refraction index
		//curobj->refr_index_i=vacu_refr_index;
		
		//curobj->absorb=SAC_water;
		//curobj->I_absorb=0.0;
		curobj->subobj=NULL;
		curobj->nextobj=NULL;
		curobj->belongobj=p_cell_ms[loop1-1];
		curobj->celltype=4;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=1;
		curobj->count_hit=0;
	}
	printf("finish initial the ms\n");

	//////////////////nonMS///////////////////
	for(loop1=1;loop1<=nonms_num;loop1++)
	{
		curobj=p_cell_ns[loop1-1];
		//curobj->refr_index_o=air_refr_index;
		//curobj->refr_index_i=cyto_refr_index;
		
		//curobj->absorb=SAC_cyto;
		//curobj->I_absorb=0.0;
		curobj->subobj=NULL;

		if(loop1!=nonms_num)
		{
			curobj->nextobj=p_cell_ns[loop1];
		}
		else
		{
			curobj->nextobj=NULL;
		}
		
		curobj->belongobj=p_leaf;
		curobj->celltype=3;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=0;
		curobj->count_hit=0;
	}	
	printf("finish initialization\n");
	

	/*****************************************************************
	*            Check chloroplasts are inside MS?
	*****************************************************************/
	Object *curms, *curchlo;	
	int k_pts,m_tri;
	Ray *curray;
	double tmp_dis[1];//distance from TriIsect, tmp_dis is the pointer of tmp_dis[]
	double *all_dis;
	int size_all_dis;
	int tmp_flag;//flag from TriIsect; flag=1 --> positive dis or negative dis
	int tmp_count_isectpts,tmp2_count_isectpts;
	for(i=1;i<=msall_num;i++)
	{
		curms=p_cell_ms[i-1];
		for(j=1;j<=ms_chl_num[i-1];j++)
		{
			curchlo=p_chl_ms[i-1][j-1];
			// check each pts of chloroplasts inside MS or not?
			for(k_pts=0;k_pts<curchlo->S->pnum;k_pts++)
			{
				curray=(Ray *)malloc(sizeof(Ray));
				curray->P[0]=(curchlo->S->vertex+k_pts)->coordinate[0];
				curray->P[1]=(curchlo->S->vertex+k_pts)->coordinate[1];
				curray->P[2]=(curchlo->S->vertex+k_pts)->coordinate[2];
				curray->D[0]=0.0;
				curray->D[1]=0.0;
				curray->D[2]=-1.0;
				if((curray->P[0]>(xmin+DBL_EPSILON))&&
				(curray->P[0]<(xmax-DBL_EPSILON))&&
				(curray->P[1]>(ymin+DBL_EPSILON))&&
				(curray->P[1]<(ymax-DBL_EPSILON)))
				{
				//// if -- it's not on the boundary of "leaf"
				//// else -- set this pts as orig, dir=[0,0,-1]
				//// for loop each triangle of current MS, check how many intersection
					tmp_count_isectpts=0;
					all_dis=calloc(10, sizeof(double));
					size_all_dis=0;
					for(m_tri=0;m_tri<curms->S->tnum;m_tri++)
					{
						////// use TriIsect(Ray *ray,Object *curobj,int i_tri,double *dis) in Trace_dicot.c
						tmp_flag=TriIsect(curray,curms,m_tri,tmp_dis);
						if(tmp_flag==1&&tmp_dis[0]>(-DBL_EPSILON))
						{
							tmp_count_isectpts=tmp_count_isectpts+1;
							*(all_dis+size_all_dis)=*tmp_dis;
							size_all_dis=size_all_dis+1;
						}
					}
					
					//// Revised June 2021
					//// before decide PTS is inside or not based on tmp_count_isectpts, unique_vec(all tmp_dis) to exclude
					//// cases that intersects at the edge shared by two triangles
					tmp2_count_isectpts=unique_dis(all_dis, size_all_dis);

					//// if tmp_count_isectpts=1, OK; 0 or 2, outside; >=3 unusual
					if(tmp2_count_isectpts==1)
					{
						//printf("MS_%d_CHL_%d.ply PTS %d is inside PAL_MS_%d.ply\n",i,j,k_pts,i);
					}
					if(tmp2_count_isectpts==0||tmp2_count_isectpts==2)
					{
						printf("[Warning] MS_%d_CHL_%d.ply PTS %d is not inside PAL_MS_%d.ply\n",i,j,k_pts,i);
					}
					if(tmp2_count_isectpts>=3)
					{
						printf("[Warning][count_isectpts=%d] MS_%d_CHL_%d.ply PTS %d is not fully inside PAL_MS_%d.ply\n",tmp_count_isectpts,i,j,k_pts,i);
					}
				}
				free(curray);
			}
		}
	}
	time_t end=time(NULL);
	printf("Inner_check finished. If no warning appears, test passed.\n");
	printf("time:%f sec\n",difftime(end,start));
	return 1;
}
