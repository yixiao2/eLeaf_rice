/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

/***********************************************
*    version 2017/06/14
*    ray tracing v2 for rice
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

int main(int argc, char *argv[])
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * *
	 *		Part1. read in the object and precalculate the normal vector
	 * * * * * * * * * * * * * * * * * * * * * * * * * */
	time_t start=time(NULL);
//	srand((unsigned)(time(0)));
	int ray_i_start,ray_i_end;
	int ray_j_start,ray_j_end;

	ray_i_start=atoi(*++argv);
	ray_i_end=atoi(*++argv);
	ray_j_start=atoi(*++argv);
	ray_j_end=atoi(*++argv);

//	SAC_cellwall=1.47e2*atof(*++argv);
	SAC_water=1e2*atof(*++argv);
	SAC_chlab=1e-4*atof(*++argv);
//	SAC_chlab_BS=1e-4*atof(*++argv);
	SAC_cyto=SAC_water;
	
	chl_con_MS=atof(*++argv);//(2.352941e4) //g m^-3
//	chl_con_BS=atof(*++argv);//(2.352941e4) //g m^-3
	
	printf("SAC_water:%e SAC_chlab:%e\n",SAC_water,SAC_chlab);
	printf("chl_con_MS:%e\n",chl_con_MS);
	
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
	
	//////////////////////read in count_chl///////////////
	fp=fopen("count_chl","r");
	for(loop1=1;loop1<=(ms_num);loop1++)
	{
		fscanf(fp,"%d\n",&(ms_chl_num[loop1-1]));
	}
	fclose(fp);
	
//////////////////////////read in MS///////////////////////	
	printf("read in MS\n");
	count_ms=0;	
	for(loop1=1;loop1<=(ms_num);loop1++)
	{	   
		count_ms++;
		/********************
		**cell
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
			temp_t->I_absorb_tri=0.0;
			temp_t->count_hit_tri=0;
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ms[count_ms-1]=curobj;
				
		/********************
		**chl
		********************/
		for(i=1;i<=ms_chl_num[loop1-1];i++)
		{
			strcpy(chl,"MS/ms");
			sprintf(str_num,"%d",loop1);
			strcat(chl,str_num);
			sprintf(str_num,"%d",i);
			strcat(chl,"c");
			strcat(chl,str_num);
			strcat(chl,".ply");
			printf("chl:%s\n",chl);

			if((fp=fopen(chl,"r"))==NULL)
			{
				printf("error on open %s\n",chl);
			}

			fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
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
				temp_t->I_absorb_tri=0.0;
				temp_t->count_hit_tri=0;
				temp_t+=1;
			}
			fclose(fp);
			precalculate(curobj);
			p_chl_ms[count_ms-1][i-1]=curobj;
		}
		
		/***********************
		**vacuole
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
			temp_t->I_absorb_tri=0.0;
			temp_t->count_hit_tri=0;
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
	for(loop1=1;loop1<=num_nonMS;loop1++)//9 upper epidermis
	{
		//////////read in cell/////////
		strcpy(cell,"nonMS/ns");
		sprintf(str_num,"%d",loop1);
		strcat(cell,str_num);		
		strcat(cell,".ply");
		printf("cell:%s\n",cell);

		if((fp=fopen(cell,"r"))==NULL)
		{
			printf("error on open %s\n",cell);
		}		
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
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
			temp_t->I_absorb_tri=0.0;
			temp_t->count_hit_tri=0;
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ns[loop1-1]=curobj;
		printf("finish read in %s\n",cell);
	}
	printf("\n");
		
	printf("finish read in \n");

	////////for debug the new precalculate()//////////////
	curobj=p_chl_ms[1][0];//p_leaf p_cell_ms[0] p_chl_ms[0][0] p_vac_ms[0] p_cell_ns[0]
	tri_num=curobj->S->tnum;
	if(flag_debug_precal==1)
	{
		FILE *f_tmp;
		Triangle *tri_tmp;
		f_tmp=fopen("debug_norvec.txt","wb");
		for(i=0;i<tri_num;i++)
		{
			tri_tmp=curobj->S->tri+i;
			fprintf(f_tmp,"%e %e %e\n",tri_tmp->Nor[0],tri_tmp->Nor[1],tri_tmp->Nor[2]);
		}
		fclose(f_tmp);
	}
	///////output normal vector of object////////////////
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part2. initialize the object
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	///////////////////////leaf///////////////// 
	curobj=p_leaf;
	curobj->refr_index_o=0.0001;//maybe this can already help mirror rays back?
	curobj->refr_index_i=air_refr_index;
	curobj->absorb=0.0;
	curobj->I_absorb=0.0;
	curobj->subobj=p_cell_ms[0];
	
	curobj->nextobj=NULL;
	curobj->belongobj=NULL;
	curobj->celltype=0;
	curobj->cellname=0;
	curobj->chlname=0;
	curobj->vacuolename=0;
	curobj->count_hit=0;
	///////////////////////ms///////////////// 
	for(loop1=1;loop1<=ms_num;loop1++)
	{
		//////// p_cell_ms/////////
		curobj=p_cell_ms[loop1-1];
		curobj->refr_index_o=air_refr_index;
		curobj->refr_index_i=cyto_refr_index;//cell wall is considered specifically in trace.c
		
		curobj->absorb=SAC_cyto;
		curobj->I_absorb=0.0;
		curobj->subobj=p_chl_ms[loop1-1][0];
		if(loop1!=ms_num)
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
			curobj->refr_index_o=cyto_refr_index; // refraction index
			curobj->refr_index_i=chlo_refr_index;
		
			curobj->absorb=SAC_chlab*chl_con_MS;
			curobj->I_absorb=0.0;
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
		curobj->refr_index_o=cyto_refr_index; // refraction index
		curobj->refr_index_i=vacu_refr_index;
		
		curobj->absorb=SAC_water;
		curobj->I_absorb=0.0;
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
	for(loop1=1;loop1<=num_nonMS;loop1++)
	{
		curobj=p_cell_ns[loop1-1];
		curobj->refr_index_o=air_refr_index;
		curobj->refr_index_i=cyto_refr_index;
		
		curobj->absorb=SAC_cyto;
		curobj->I_absorb=0.0;
		curobj->subobj=NULL;

		if(loop1!=num_nonMS)
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
	
	I_discard=0.0;
	I_discard_Rf=0.0;
	I_discard_Tr=0.0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part3. generate light and trace
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	int cutnum_x=500,cutnum_y=25;//50000=223.6^2
	double dx=(xmax-xmin)/(cutnum_x+1),dy=(zmax-zmin)/(cutnum_y+1);
/*		
	Ray *curray;
	curray=(Ray *)malloc(sizeof(Ray));
	curray->P[0]=-2.400000e-05;
	curray->P[1]=-2.313500e-06;
	curray->P[2]=zmax;
	curray->D[0]=-4.114551e-01;
	curray->D[1]=8.075262e-01;
	curray->D[2]=-4.226183e-01;
	curray->I=1.0;
//	curray->interflag=0;
	curray->level=1;
	curray->P_obj=p_leaf;
	curray->interflag=0;
//	curray->P_tri=175;
	printf("\n");
	printf("start ray tracing\n");
	printf("\n");
	
	trace(curray,p_cell_s_o[0]);
*/

	if(flag_randseed==1)
	{
		srand(time(NULL));
	}

	Ray *curray;
	double alpha,beta;
	double sum_I_absorb_all=0.0;
	for(ray_i=ray_i_start;ray_i<ray_i_end;ray_i++) ///// 0 < cutnum
	{
		for(ray_j=ray_j_start;ray_j<ray_j_end;ray_j++)
		{
			curray=(Ray *)malloc(sizeof(Ray));
			// initialize the current light
			debugI=0.0;
			num_chl_hit=0;
			
			sum_I_absorb_all=0.0;
			curray->P[0]=xmin+(ray_i+1)*dx;
			curray->P[1]=ymax-0.05e-6;
			curray->P[2]=zmin+(ray_j+1)*dy;
			//curray->P[2]=zmin+0.1e-6;
		
			curray->D[0]=0.0;
			curray->D[1]=-1.0;
			curray->D[2]=0.0;
			//curray->D[2]=1;
			curray->P_tri=-1;
			
			if(flag_debug==1)
			{
				RT_debug[ray_i][ray_j][0][0]=curray->P[0];
				RT_debug[ray_i][ray_j][0][1]=curray->P[1];
				RT_debug[ray_i][ray_j][0][2]=curray->P[2];
				count_RT_debug[ray_i][ray_j]=1;
			}
			
		
			//alpha=rand()%178+1;
			alpha=(double)rand()/RAND_MAX;//0-1 degree
			printf("alpha:%e\n",alpha);
			curray->D[1]=-cos(alpha/180*PI);
			beta=rand()/360;
			//beta=300;
			curray->D[0]=sin(alpha/180*PI)*sin(beta/180*PI);
			curray->D[2]=sin(alpha/180*PI)*cos(beta/180*PI);
		/**/
			curray->I=1.0; 
			curray->interflag=0;
			curray->level=1;
			curray->P_obj=p_leaf;
			curray->trace_obj=p_leaf->subobj;
			curray->refr_index=air_refr_index;
			
			//printf("#####################\n");
			printf("trace ray i=%d j=%d\n",ray_i,ray_j);
			//printf("P:%e %e %e\n",curray->P[0],curray->P[1],curray->P[2]);
			//printf("D:%e %e %e\n",curray->D[0],curray->D[1],curray->D[2]);
			//printf("#####################\n");
			trace(curray,p_leaf->subobj);
		//	trace(curray,p_cell_e_o[0]);
			printf("\n");
			
			sum_p();
		}
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part4. output file
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	int k;

	FILE *fout_tri;
	FILE *fout_cell;
	//////////////////////////////////////////////////////////////////
	if((fout_tri=fopen(*++argv,"wb"))==NULL)
	{
		printf("error on open file result tri\n");
	}
	if((fout_cell=fopen(*++argv,"wb"))==NULL)
	{
		printf("error on open file result cell\n");
	}

	for(i=1;i<=ms_num;i++)
	{
		curobj=p_cell_ms[i-1];
		sum_I_absorb_all+=curobj->I_absorb;
		for(k=0;k<curobj->S->tnum;k++)
		{
		//	sum_obj_absorb+=(curobj->S->tri+k)->I_absorb_tri;
			fprintf(fout_tri,"%d %d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri);
		}
		fprintf(fout_cell,"%d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb);
		
		for(j=1;j<=ms_chl_num[i-1];j++)
		{
			curobj=p_chl_ms[i-1][j-1];
			sum_I_absorb_all+=curobj->I_absorb;
			for(k=0;k<curobj->S->tnum;k++)
			{
			//	sum_obj_absorb+=(curobj->S->tri+k)->I_absorb_tri;
				fprintf(fout_tri,"%d %d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri);
			}
			fprintf(fout_cell,"%d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb);
		}
		
		curobj=p_vac_ms[i-1];
		sum_I_absorb_all+=curobj->I_absorb;
		for(k=0;k<curobj->S->tnum;k++)
		{
		//	sum_obj_absorb+=(curobj->S->tri+k)->I_absorb_tri;
			fprintf(fout_tri,"%d %d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri);
		}
		fprintf(fout_cell,"%d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb);
		
	}
	
	//////////// EP ///////////////////
	for(i=1;i<=num_nonMS;i++)
	{
		curobj=p_cell_ns[i-1];
		sum_I_absorb_all+=curobj->I_absorb;
		for(k=0;k<curobj->S->tnum;k++)
		{
		//	sum_obj_absorb+=(curobj->S->tri+k)->I_absorb_tri;
			fprintf(fout_tri,"%d %d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri);
		}
		fprintf(fout_cell,"%d %d %d %d %d %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb);	
	}
	
	fclose(fout_tri);
	fclose(fout_cell);
	
	FILE *fout;
	if((fout=fopen(*++argv,"wb"))==NULL)
	{
		printf("error open file result_sum\n");
	}
	fprintf(fout,"Absorption %e\n",sum_I_absorb_all);
	fprintf(fout,"Reflectance %e\n",I_discard_Rf);
	fprintf(fout,"Transmitance %e\n",I_discard_Tr);
	fprintf(fout,"Discard %e\n",I_discard);
	fprintf(fout,"Sum_I %e\n",sum_I_absorb_all+I_discard_Rf+I_discard_Tr+I_discard);
	
	time_t end=time(NULL);
	printf("time:%f sec\n",difftime(end,start));
	fprintf(fout,"Time %f\n",difftime(end,start));
	fclose(fout);
	
	///////////debug_RT///////////////////
	if(flag_debug==1)
	{
		FILE *f_debug;
		f_debug=fopen("debug_RT.txt","wb");
		for(ray_i=ray_i_start;ray_i<ray_i_end;ray_i++)
		{
			for(ray_j=ray_j_start;ray_j<ray_j_end;ray_j++)
			{
				if(count_RT_debug[ray_i][ray_j]<513)
					fprintf(f_debug,"%d\n",count_RT_debug[ray_i][ray_j]);
				else
					fprintf(f_debug,"%d\n",512);
				for(k=0;k<512&&k<count_RT_debug[ray_i][ray_j];k++)
				{
					fprintf(f_debug,"%e %e %e\n",RT_debug[ray_i][ray_j][k][0],RT_debug[ray_i][ray_j][k][1],RT_debug[ray_i][ray_j][k][2]);
				}
				fprintf(f_debug,"\n");
			}
		}
		fclose(f_debug);
	}
	return 1;
}
