/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

/***********************************************
*    version 11/11/11
*    version for rice 2013/7/10
*    new version for rice 2017/6/19
**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Defs.h"
#include "Global.h"

int TriIsect(Ray *ray,Object *curobj,int i_tri,double *dis)
{
	//double dis;
	//Triangle *curTri=curobj->S->tri+i_tri;
	Vec v0,v1,v2,dir,orig;	// three points of the triangle
	int i;
	for(i=0;i<3;i++)
	{
		v0[i]=(curobj->S->vertex+(curobj->S->tri+i_tri)->n[0])->coordinate[i];
		v1[i]=(curobj->S->vertex+(curobj->S->tri+i_tri)->n[1])->coordinate[i];
		v2[i]=(curobj->S->vertex+(curobj->S->tri+i_tri)->n[2])->coordinate[i];
		dir[i]=ray->D[i];
		orig[i]=ray->P[i];
	}

	Vec E1,E2,P;
	VecSub(v0,v1,E1);
	VecSub(v0,v2,E2);
	VecCross(dir,E2,P);
	
	double det;
	det=VecDot(E1,P);
	
	Vec T;
	if(det>0)
	{
		VecSub(v0,orig,T);
	}
	else
	{
		VecSub(orig,v0,T);
		det=-det;
	}
	
	// if determinant is near zero, ray lies in plane of triangle
	if(det<DBL_EPSILON_trisect)
	{
		dis[0]=-Infinity;//dis<0 means no intersection
		//printf("det:%e\n",det);
		return 0;
	}
	
	// calculate u and make sure u<=1
	double u;
	u=VecDot(T,P);
	if(u<0.0||u>det)
	{
		dis[0]=-Infinity;
		//printf("u:%f\n",u);
		return 0;
	}
	
	// Q
	Vec Q;
	VecCross(T,E1,Q);
	
	// calculate v and make sure u+v<=1
	double v;
	v=VecDot(dir,Q);
	if(v<0.0||(u+v)>det)
	{
		//printf("u:%f v:%f det:%f\n",u,v,det);
		dis[0]=-Infinity;
		return 0;
	}
	
	// calculate t, scale parameters, ray intersects triangle
	double t;
	t=VecDot(E2,Q);
	
	double fInvDet=1.0/det;
	dis[0]=t*fInvDet;
	
	//if(dis[0]<0.0)
	//printf("dis:%f\n",dis[0]);
	//u=u*fInvDet;
	//v=v*fInvDet;

	return 1;
}

void ObjIsect(Ray *ray,Object *obj,Isect *intersect,int outin)//outin=1 means ray is outside, outin=0 means ray is inside
{
	Object *tmp_obj=obj;
	int i,inter_i,count,flag_intersect,inter_i_neg,tmp_flag;
	int inter_i_front,inter_i_back;
	i=0;
	count=0;
	flag_intersect=0;
	double t,temp_dis[1],t_neg,tmp_dbl;
	double t_front,t_back;
	t=Infinity;
	t_neg=-Infinity;
	t_front=Infinity;
	t_back=Infinity;
	while(i<tmp_obj->S->tnum)
	{
		if(tmp_obj==ray->P_obj&&i==ray->P_tri)
		{
			i++;
		}
		tmp_flag=TriIsect(ray,tmp_obj,i,temp_dis);
		
		if(temp_dis[0]>0.0)
		{
			flag_intersect=1;
			if(temp_dis[0]<t)
			{
				t=temp_dis[0];
				inter_i=i;
			}
			tmp_dbl=VecDot(ray->D,(tmp_obj->S->tri+inter_i)->Nor);
			if(tmp_dbl>0.0)
			{
				if(temp_dis[0]<t_front)
				{
					t_front=temp_dis[0];
					inter_i_front=i;
				}
			}
			else
			{
				if(temp_dis[0]<t_back)
				{
					t_back=temp_dis[0];
					inter_i_back=i;
				}
			}
		}
		else if(temp_dis[0]<0.0)
		{
			if(temp_dis[0]>t_neg)
			{
				t_neg=temp_dis[0];
				inter_i_neg=i;
			}
		}
		i++;
		if(tmp_obj==ray->P_obj&&i==ray->P_tri)//for i=ray->P_tri
		{
			i++;
		}
	}
	
	if(outin==0)
	{
		if(flag_intersect==1)
		{
			tmp_dbl=VecDot(ray->D,(tmp_obj->S->tri+inter_i)->Nor);
			if(tmp_dbl>0.0)
			{
				intersect->inter_tri=inter_i;
				intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i)->Nor[0];
				intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i)->Nor[1];
				intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i)->Nor[2];
				intersect->inter_pt[0]=ray->P[0]+t*ray->D[0];
				intersect->inter_pt[1]=ray->P[1]+t*ray->D[1];
				intersect->inter_pt[2]=ray->P[2]+t*ray->D[2];
				intersect->inter_obj=tmp_obj;
				intersect->t=t;
			}
			else
			{
				if(flag_warning==1)
				{
					printf("Warning:Trace03: ray went outside, but intersect with wrong triangle considering normal vector\n");
				}
				
				if(t_front!=Infinity)
				{
					intersect->inter_tri=inter_i_front;
					intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i_front)->Nor[0];
					intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i_front)->Nor[1];
					intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i_front)->Nor[2];
					intersect->inter_pt[0]=ray->P[0]+t_front*ray->D[0];
					intersect->inter_pt[1]=ray->P[1]+t_front*ray->D[1];
					intersect->inter_pt[2]=ray->P[2]+t_front*ray->D[2];
					intersect->inter_obj=tmp_obj;
					intersect->t=t_front;
				}
				else
				{
					if(flag_warning==1)
					{
						printf("Warning:Trace04: ray went outside, cannot find intersect consider normal vector\n");
					}
					if(t_back!=Infinity)
					{
						intersect->inter_tri=inter_i_back;
						intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i_back)->Nor[0];
						intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i_back)->Nor[1];
						intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i_back)->Nor[2];
						intersect->inter_pt[0]=ray->P[0]+t_back*ray->D[0];
						intersect->inter_pt[1]=ray->P[1]+t_back*ray->D[1];
						intersect->inter_pt[2]=ray->P[2]+t_back*ray->D[2];
						intersect->inter_obj=tmp_obj;
						intersect->t=t_back;
					
					//mirror ray accordingly
					//ray->P --> ray->P-((tmp_obj->S->tri+inter_i)->Nor)*(2*t_back*(-tmp_dbl))
					//ray->P[0]=ray->P[0]-intersect->inter_nor[0]*(2*t_back*(-tmp_dbl));
					//ray->P[1]=ray->P[1]-intersect->inter_nor[1]*(2*t_back*(-tmp_dbl));
					//ray->P[2]=ray->P[2]-intersect->inter_nor[2]*(2*t_back*(-tmp_dbl));
					//ray->D
						ray->D[0]=ray->D[0]+intersect->inter_nor[0]*(2*(-tmp_dbl));
						ray->D[1]=ray->D[1]+intersect->inter_nor[1]*(2*(-tmp_dbl));
						ray->D[2]=ray->D[2]+intersect->inter_nor[2]*(2*(-tmp_dbl));
					}
					else
					{
						printf("You must be kidding...\n");
						exit(1);
					}
				}
			}
		}
		else
		{
			if(t_neg==-Infinity)
				printf("error here\n");
			intersect->inter_tri=inter_i_neg;
			intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i_neg)->Nor[0];
			intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i_neg)->Nor[1];
			intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i_neg)->Nor[2];
			intersect->inter_pt[0]=ray->P[0]+t_neg*ray->D[0];
			intersect->inter_pt[1]=ray->P[1]+t_neg*ray->D[1];
			intersect->inter_pt[2]=ray->P[2]+t_neg*ray->D[2];
			intersect->inter_obj=tmp_obj;
			intersect->t=0.0;
			if(flag_warning==1)
			{
				printf("Warning:Trace01: ray is inside, but has no positive intersection\n");
				printf("error here 1;tfront %e;tback %e;tneg %e\n",t_front,t_back,t_neg);
				printf("ray->P:%e %e %e\n",ray->P[0],ray->P[1],ray->P[2]);
				printf("ray->D:%e %e %e\n",ray->D[0],ray->D[1],ray->D[2]);
			}
		}
	}
	else//outin==1
	{
		if(flag_intersect==0)
		{
			intersect->t=-1;
		}
		else
		{
			tmp_dbl=-VecDot(ray->D,(tmp_obj->S->tri+inter_i)->Nor);
			if(tmp_dbl>0.0)
			{
				intersect->inter_tri=inter_i;
				intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i)->Nor[0];
				intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i)->Nor[1];
				intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i)->Nor[2];
				intersect->inter_pt[0]=ray->P[0]+t*ray->D[0];
				intersect->inter_pt[1]=ray->P[1]+t*ray->D[1];
				intersect->inter_pt[2]=ray->P[2]+t*ray->D[2];
				intersect->inter_obj=tmp_obj;
				intersect->t=t;
			}
			else
			{
				if(t_neg!=-Infinity)
				{
					if(flag_warning==1)
					{
						printf("error here 2;tfront %e;tback %e;tneg %e\n",t_front,t_back,t_neg);
						printf("Warning:Trace02: ray is outside, but intersect with wrong triangle considering normal vector\n");
					}
					intersect->inter_tri=inter_i_neg;
					intersect->inter_nor[0]=(tmp_obj->S->tri+inter_i_neg)->Nor[0];
					intersect->inter_nor[1]=(tmp_obj->S->tri+inter_i_neg)->Nor[1];
					intersect->inter_nor[2]=(tmp_obj->S->tri+inter_i_neg)->Nor[2];
					intersect->inter_pt[0]=ray->P[0]+t_neg*ray->D[0];
					intersect->inter_pt[1]=ray->P[1]+t_neg*ray->D[1];
					intersect->inter_pt[2]=ray->P[2]+t_neg*ray->D[2];
					intersect->inter_obj=tmp_obj;
					intersect->t=0.0;
				}
				else
				{
					flag_intersect=0;
					intersect->t=-1;
				}
			}
		}
	}
	if(intersect->inter_obj==NULL&&outin==0)
	{
		printf("debug\n");
	}
}



void compute_ray_next(Ray *ray,Isect *intersect0,Isect *intersect1,int *level_contact, double *level_contact_refr,int count_contact,Ray *ray_next,int flag_ray_next)
{
	//calculate r1/t1 r2/t2... rn/tn
	int loop_i,TIR=0,flag_stop=0,tmp_level;
	double tmp_n1,tmp_n2,rs,rp,reflectance[10],transmittance[10];
	Vec D1,norm,tmp_norm1,tmp_norm2;
	double sin_t[10],cos_t[10],sin_t1,cos_t1,sin_t2,cos_t2,t,tmp_rnd;
	
	D1[0]=ray->D[0];
	D1[1]=ray->D[1];
	D1[2]=ray->D[2];
	norm[0]=intersect0->inter_nor[0];
	norm[1]=intersect0->inter_nor[1];
	norm[2]=intersect0->inter_nor[2];
	cos_t1=VecDot(D1,norm);
	if(cos_t1<0.0)
	{
		cos_t1=-cos_t1;
		norm[0]=-norm[0];
		norm[1]=-norm[1];
		norm[2]=-norm[2];
	}
	if(cos_t1>=1.0)
		cos_t1=1.0;
	sin_t1=sqrt(1-cos_t1*cos_t1);
	cos_t[0]=cos_t1;
	sin_t[0]=sin_t1;
	for(loop_i=0;loop_i<=count_contact-1&&TIR==0;loop_i++)
	{
		tmp_n1=level_contact_refr[loop_i];
		tmp_n2=level_contact_refr[loop_i+1];
		cos_t1=cos_t[loop_i];
		sin_t1=sin_t[loop_i];
		if((1-cos_t1*cos_t1)*(tmp_n1*tmp_n1/tmp_n2/tmp_n2)>1)
		{
			TIR=1;
		}
		else
		{
			sin_t2=sin_t1*tmp_n1/tmp_n2;
			cos_t2=sqrt(1-sin_t2*sin_t2);
			rs=((cos_t1-(tmp_n2/tmp_n1)*cos_t2)/(cos_t1+(tmp_n2/tmp_n1)*cos_t2))*
				((cos_t1-(tmp_n2/tmp_n1)*cos_t2)/(cos_t1+(tmp_n2/tmp_n1)*cos_t2));
			rp=((cos_t2-(tmp_n2/tmp_n1)*cos_t1)/(cos_t2+(tmp_n2/tmp_n1)*cos_t1))*
				((cos_t2-(tmp_n2/tmp_n1)*cos_t1)/(cos_t2+(tmp_n2/tmp_n1)*cos_t1));
			reflectance[loop_i]=(rs+rp)/2;
			transmittance[loop_i]=1-reflectance[loop_i];
			cos_t[loop_i+1]=cos_t2;
			sin_t[loop_i+1]=sin_t2;
		}
	}
	
	if(TIR!=1)//not total reflection, rand until reflected at boundary 1 or transmitted at boundary n
	{
		tmp_level=0;
		do
		{
			tmp_rnd=(double)(rand()*1.0/RAND_MAX);
			if(tmp_rnd<reflectance[tmp_level])
				tmp_level--;
			else
				tmp_level++;
			
			if(tmp_level<0||tmp_level>count_contact-1)
			{
				flag_stop=1;
				if(tmp_level<0)
					TIR=1;
			}
		}while(flag_stop==0);
	}
	
	if(TIR==1)//here TIR became flag_reflectance
	{
		ray_next->P[0]=intersect0->inter_pt[0];
		ray_next->P[1]=intersect0->inter_pt[1];
		ray_next->P[2]=intersect0->inter_pt[2];
		t=cos_t[0];
		ray_next->D[0]=D1[0]-2*norm[0]*t;
		ray_next->D[1]=D1[1]-2*norm[1]*t;
		ray_next->D[2]=D1[2]-2*norm[2]*t;
		ray_next->P_obj=intersect0->inter_obj;
		ray_next->P_tri=intersect0->inter_tri;
		ray_next->trace_obj=ray->trace_obj;
		ray_next->refr_index=ray->refr_index;
	}
	else//transmitted ray
	{
		ray_next->P[0]=intersect1->inter_pt[0];
		ray_next->P[1]=intersect1->inter_pt[1];
		ray_next->P[2]=intersect1->inter_pt[2];
		
		VecCross(ray->D,norm,tmp_norm1);
		Normalize(tmp_norm1);
		VecCross(norm,tmp_norm1,tmp_norm2);
		cos_t2=cos_t[count_contact];
		sin_t2=sin_t[count_contact];
		ray_next->D[0]=norm[0]*cos_t2+tmp_norm2[0]*sin_t2;
		ray_next->D[1]=norm[1]*cos_t2+tmp_norm2[1]*sin_t2;
		ray_next->D[2]=norm[2]*cos_t2+tmp_norm2[2]*sin_t2;
		ray_next->P_obj=intersect1->inter_obj;
		ray_next->P_tri=intersect1->inter_tri;
		if(flag_ray_next==1)
		{
			ray_next->trace_obj=intersect1->inter_obj->belongobj->subobj;
			ray_next->refr_index=intersect1->inter_obj->refr_index_o;
		}
		else
		{
			ray_next->trace_obj=intersect1->inter_obj->subobj;
			ray_next->refr_index=intersect1->inter_obj->refr_index_i;
		}
	}
	
}

void trace(Ray *ray, Object *obj)
{
	intersect=(Isect *)malloc(sizeof(Isect));
	if(intersect==NULL)
		exit(1);
	intersect->t=-1;
	Isect *tmp_intersect;
	tmp_intersect=(Isect *)malloc(sizeof(Isect));
	if(tmp_intersect==NULL)
		exit(1);
	tmp_intersect->t=-1;
	Isect *intersect_cld;
	intersect_cld=(Isect *)malloc(sizeof(Isect));//used when searching children node
	if(intersect_cld==NULL)
		exit(1);
	intersect_cld->t=-1;
	Isect *intersect_mth;
	intersect_mth=(Isect *)malloc(sizeof(Isect));//used when searching mother node
	if(intersect_mth==NULL)
		exit(1);
	intersect_mth->t=-1;
	Isect *intersect_bk;
	intersect_bk=(Isect *)malloc(sizeof(Isect));//backup
	if(intersect_bk==NULL)
		exit(1);
	Ray *ray_next;//Ray *ray_refl,*ray_trans;
	ray_next=(Ray *)malloc(sizeof(Ray));

	int tmp_level,flag_stop,flag_ray_next=0,loop_i,count_contact,tmp_count_debug;
	int tmp_level_contact[10];
	for(loop_i=0;loop_i<10;loop_i++)
	{
		tmp_level_contact[loop_i]=-1;//level=0=leaf;=1=IAS;=2=cyto;=3=chlo/vacu
	}
	double tmp_level_contact_refr[10];//record the refractive index of each level, since chlo has different refr to vacu
	
	if(obj==NULL)
	{
		intersect->t=-1;
//		printf("NULL object. No intersection\n");
	}
	else
	{
//		printf("TRACE()trinor:%e %e %e\n",(obj->S->tri/*+i*size_tri*/)->Nor[0],
//						(obj->S->tri/*+i*size_tri*/)->Nor[1],
//						(obj->S->tri/*+i*size_tri*/)->Nor[2]);
		curobj=obj;
		if(obj->belongobj==NULL)//leaf? possible?
		{
			curobj=obj;
		}
		else
		{
			curobj=(obj->belongobj)->subobj;
		}
		
		while(curobj!=NULL)//function: search all children nodes of the mother node of curobj
		{
			ObjIsect(ray,curobj,tmp_intersect,1);
			if(tmp_intersect->t<0)
			{
				//curobj=curobj->nextobj;
			}
			else
			{
				if(intersect->t<0)
				{
					IsectCopy(tmp_intersect,intersect);
				}
				else
				{
					if(tmp_intersect->t<intersect->t)
					{
						IsectCopy(tmp_intersect,intersect);
					}
				}
			}
			curobj=curobj->nextobj;
		}
	}
///	printf("finish intersection\n");
///	printf("raylevel:%d\n",ray->level);
///	printf("intersect=NULL? %d\n",intersect==NULL);
	if(intersect->t>=0.0)//search children node of intersect point //search chlldren's children until 1.NULL children node; or 2.not contact anymore
	{
		//calculate absorption
		ray_next->I=ray->I*exp(-(ray->P_obj->absorb)*(intersect->t));//this is OK?
		ray->P_obj->I_absorb+=(ray->I-ray_next->I);
		if(ray->P_obj->celltype!=0)
			(ray->P_obj->S->tri+ray->P_tri)->I_absorb_tri+=(ray->I-ray_next->I);
		
		IsectCopy(intersect,intersect_bk);//backup for refl&refr calculation on contact boundary
		
		flag_stop=0;
		tmp_level=ray->level;
		intersect_cld->inter_obj=intersect->inter_obj;
		tmp_level_contact[0]=ray->level;
		tmp_level_contact_refr[0]=ray->refr_index;
		tmp_level++;
		//tmp_level_contact[1]=tmp_level;
		//tmp_level_contact[1]=intersect->inter_obj->refr_index_i;
		//count_contact=2;
		count_contact=1;
		if(intersect->inter_obj->celltype==1||intersect->inter_obj->celltype==3)//cell wall
		{
			tmp_level_contact[count_contact]=tmp_level;
			tmp_level_contact_refr[count_contact]=wall_refr_index;
			count_contact++;
		}
		do
		{
			curobj=intersect_cld->inter_obj->subobj;
			intersect_cld->t=-1;
			while(curobj!=NULL)
			{
				ObjIsect(ray,curobj,tmp_intersect,1);
				if(tmp_intersect->t>0)
				{
					if(intersect_cld->t<0||tmp_intersect->t<intersect_cld->t)
					{
						IsectCopy(tmp_intersect,intersect_cld);
					}
				}
				curobj=curobj->nextobj;
			}//intersect_cld record closest child node
			flag_stop=1;
			if((intersect_cld->t-intersect->t)>-DBL_EPSILON&&(intersect_cld->t-intersect->t)<DBL_EPSILON)
			{
				flag_stop=0;//should stop search child if only get a contact child here
				tmp_level++;
				if(intersect_cld->inter_obj->celltype==1||intersect_cld->inter_obj->celltype==3)//cell wall
				{
					tmp_level_contact[count_contact]=tmp_level;
					tmp_level_contact_refr[count_contact]=wall_refr_index;
					count_contact++;
				}
				IsectCopy(intersect_cld,intersect);
				//tmp_level_contact[count_contact]=tmp_level;
				//tmp_level_contact[count_contact]=intersect_cld->inter_obj->refr_index_i;
				//count_contact++;
			}	
		}while(flag_stop==0);//intersect_cld can be empty if ray->level=2
		tmp_level_contact[count_contact]=tmp_level;
		tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
		//calculate reflected and transmitted ray
		compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact,tmp_level_contact_refr,count_contact,ray_next,flag_ray_next);
		if(flag_debug_printf==1)
		{
			printf("ray_next->P:\n%e\n%e\n%e\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
			printf("ray_next->D:\n%e\n%e\n%e\n\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
			if(ray_next->trace_obj!=NULL)
				printf("ray_next->trace_obj->celltype:%d\n\n",ray_next->trace_obj->celltype);
			else
				printf("ray_next->trace_obj->celltype:NULL\n\n");
		}
		if(flag_debug==1)
		{
			tmp_count_debug=count_RT_debug[ray_i][ray_j];
			if(tmp_count_debug<512)
			{
				RT_debug[ray_i][ray_j][tmp_count_debug][0]=ray_next->P[0];
				RT_debug[ray_i][ray_j][tmp_count_debug][1]=ray_next->P[1];
				RT_debug[ray_i][ray_j][tmp_count_debug][2]=ray_next->P[2];
			}
			count_RT_debug[ray_i][ray_j]++;
		}
		free(ray);
		free(intersect);
		free(intersect_cld);
		free(intersect_mth);
		free(tmp_intersect);
		free(intersect_bk);
		ray=NULL;
		intersect=NULL;
		intersect_cld=NULL;
		intersect_mth=NULL;
		tmp_intersect=NULL;
		intersect_bk=NULL;
		if(ray_next->I<discardI)
		{
			I_discard+=ray_next->I;
			free(ray_next);
			ray_next=NULL;
		}
		else
		{
			trace(ray_next,ray_next->trace_obj);
		}
	}
	else //no intersection: search all nodes at the same layer with mother node
	{
		if(obj==NULL)
			curobj=ray->P_obj;
		else
			curobj=obj->belongobj;
		tmp_level_contact[0]=ray->level;
		tmp_level_contact_refr[0]=ray->refr_index;
		tmp_level=ray->level-1;
		//tmp_level_contact[1]=tmp_level;
		//count_contact=2;
		count_contact=1;
		ObjIsect(ray,curobj,intersect,0);//intersect with mother node
		IsectCopy(intersect,intersect_bk);//backup for boundary condition
		flag_ray_next=1;//indicate special case for a ray go out of object
		
		if(intersect->inter_obj->celltype==1||intersect->inter_obj->celltype==3)//cell wall
		{
			tmp_level_contact[count_contact]=tmp_level;
			tmp_level_contact_refr[count_contact]=wall_refr_index;
			count_contact++;
		}
		
		if(flag_debug==1)
		{
			tmp_count_debug=count_RT_debug[ray_i][ray_j];
			if(tmp_count_debug<512)
			{
				RT_debug[ray_i][ray_j][tmp_count_debug][0]=intersect->inter_pt[0];
				RT_debug[ray_i][ray_j][tmp_count_debug][1]=intersect->inter_pt[1];
				RT_debug[ray_i][ray_j][tmp_count_debug][2]=intersect->inter_pt[2];
			}
			count_RT_debug[ray_i][ray_j]++;
		}
		
		//search mother's mother etc
		intersect_mth->inter_obj=curobj;//yes, obj can not be leaf here! but curobj can be
		flag_stop=0;
		
		//calculate absorption
		ray_next->I=ray->I*exp(-(ray->P_obj->absorb)*(intersect->t));//this is OK?
		ray->P_obj->I_absorb+=(ray->I-ray_next->I);
		if(ray->P_obj->celltype!=0)
			(ray->P_obj->S->tri+ray->P_tri)->I_absorb_tri+=(ray->I-ray_next->I);
		
		do//search mother's mother etc
		{
			curobj=intersect_mth->inter_obj->belongobj;
			if(curobj==NULL)//intersect->inter_obj=leaf
			{
				flag_stop=1;
			}
			else
			{
				flag_stop=1;
				ObjIsect(ray,curobj,intersect_mth,0);
				if((intersect_mth->t-intersect->t)>-DBL_EPSILON&&(intersect_mth->t-intersect->t)<DBL_EPSILON)
				{
					flag_stop=0;//find a mother's mother contact with mother...
					IsectCopy(intersect_mth,intersect);
					tmp_level--;
						if(intersect_mth->inter_obj->celltype==1||intersect_mth->inter_obj->celltype==3)//cell wall
						{
							tmp_level_contact[count_contact]=tmp_level;
							tmp_level_contact_refr[count_contact]=wall_refr_index;
							count_contact++;
						}
					//tmp_level_contact[count_contact]=tmp_level;
					//count_contact++;
				}
			}
		}while(flag_stop==0);
		//then search at same layer
		if(intersect->inter_obj->belongobj==NULL)//touch the leaf,tmp_level==0; boundary conditions should also belong here.
		{
			//compute leaf refl&refr and mirror boundary
			switch(intersect->inter_tri)
			{
				case 2:
				case 3:
				case 0:
				case 1:
				case 6:
				case 7:
				case 10:
				case 11:
					tmp_level_contact_refr[1]=0.0001;//make sure the total reflection
					compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact,tmp_level_contact_refr,count_contact,ray_next,flag_ray_next);//use intersect_bk to calculate ray_next
					if(flag_debug_printf==1)
					{
						printf("ray_next->P:\n%e\n%e\n%e\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
						printf("ray_next->D:\n%e\n%e\n%e\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
						if(ray_next->trace_obj!=NULL)
							printf("ray_next->trace_obj->celltype:%d\n\n",ray_next->trace_obj->celltype);
						else
							printf("ray_next->trace_obj->celltype:NULL\n\n");
					}
					free(ray);
					free(intersect);
					free(intersect_cld);
					free(intersect_mth);
					free(tmp_intersect);
					free(intersect_bk);
					ray=NULL;
					intersect=NULL;
					intersect_cld=NULL;
					intersect_mth=NULL;
					tmp_intersect=NULL;
					intersect_bk=NULL;
					if(ray_next->I<discardI)
					{
						I_discard+=ray_next->I;
						free(ray_next);
						ray_next=NULL;
					}
					else
					{
						//ray_next->trace_obj=ray->trace_obj;
						trace(ray_next,ray_next->trace_obj);
					}
					break;
				case 8:
				case 9:
				I_discard_Tr+=ray_next->I;
				free(ray);
				free(intersect);
				free(intersect_cld);
				free(intersect_mth);
				free(tmp_intersect);
				free(intersect_bk);
				free(ray_next);
				ray=NULL;
				intersect=NULL;
				intersect_cld=NULL;
				intersect_mth=NULL;
				tmp_intersect=NULL;
				intersect_bk=NULL;
				ray_next=NULL;
				break;
				case 4:
				case 5:
				I_discard_Rf+=ray_next->I;
				free(ray);
				free(intersect);
				free(intersect_cld);
				free(intersect_mth);
				free(tmp_intersect);
				free(intersect_bk);
				free(ray_next);
				ray=NULL;
				intersect=NULL;
				intersect_cld=NULL;
				intersect_mth=NULL;
				tmp_intersect=NULL;
				intersect_bk=NULL;
				ray_next=NULL;
				break;
			}
		}
		else//then search at same layer
		{
			//here if no more contact, then should finished with refr_index_o
			tmp_level_contact[count_contact]=tmp_level;
			tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
			//otherwise, contact with neighbour or child of neighbour, then should finished with refr_index_i
			curobj=intersect->inter_obj->belongobj->subobj;
			flag_stop=0;
			while(curobj!=NULL&&flag_stop==0)
			{
				if(curobj==intersect->inter_obj)
				{
					//skip
				}
				else
				{
					ObjIsect(ray,curobj,tmp_intersect,1);
					if((tmp_intersect->t-intersect->t)>-DBL_EPSILON&&(tmp_intersect->t-intersect->t)<DBL_EPSILON)
					{
						IsectCopy(tmp_intersect,intersect);
						tmp_level++;
							if(tmp_intersect->inter_obj->celltype==1||tmp_intersect->inter_obj->celltype==3)//cell wall
							{
								tmp_level_contact[count_contact]=tmp_level;//this will overwrite the value above
								tmp_level_contact_refr[count_contact]=wall_refr_index;
								count_contact++;
							}
						//tmp_level_contact[count_contact]=tmp_level;
						//count_contact++;
						flag_stop=1;//impossible exist multi obj contact at the same position
					}
				}
				curobj=curobj->nextobj;
			}
			
			if(flag_stop==1)//find a contact neighbour, then search it's child and child's child etc
			{
				flag_ray_next=0;
				intersect_cld->inter_obj=intersect->inter_obj;
				do
				{
					curobj=intersect_cld->inter_obj->subobj;
					if(curobj==NULL)
					{
						flag_stop=1;
					}
					else
					{
						flag_stop=1;
						while(curobj!=NULL&&flag_stop==1)
						{
							ObjIsect(ray,curobj,intersect_cld,1);
							if((intersect_cld->t-intersect->t)>-DBL_EPSILON&&(intersect_cld->t-intersect->t)<DBL_EPSILON)
							{
								IsectCopy(intersect_cld,intersect);
								tmp_level++;
								if(intersect_cld->inter_obj->celltype==1||intersect_cld->inter_obj->celltype==3)//cell wall
								{
									tmp_level_contact[count_contact]=tmp_level;
									tmp_level_contact_refr[count_contact]=tmp_level;
									count_contact++;
								}
								//tmp_level_contact[count]=tmp_level;
								//count++;
								flag_stop=0;
							}
						curobj=curobj->nextobj;
						}
					}
				}while(flag_stop==0);
				tmp_level_contact[count_contact]=tmp_level;
				tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
			}
			
			//compute ray_next
			compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact,tmp_level_contact_refr,count_contact,ray_next,flag_ray_next);
			if(flag_debug_printf==1)
			{
				printf("ray_next->P:\n%e\n%e\n%e\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
				printf("ray_next->D:\n%e\n%e\n%e\n\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
				if(ray_next->trace_obj!=NULL)
					printf("ray_next->trace_obj->celltype:%d\n\n",ray_next->trace_obj->celltype);
				else
					printf("ray_next->trace_obj->celltype:NULL\n\n");
			}
			if(ray_next->I<discardI)
			{
				I_discard+=ray_next->I;
				free(ray_next);
				free(ray);
				free(intersect);
				free(intersect_cld);
				free(intersect_mth);
				free(tmp_intersect);
				free(intersect_bk);
				ray=NULL;
				intersect=NULL;
				intersect_cld=NULL;
				intersect_mth=NULL;
				tmp_intersect=NULL;
				intersect_bk=NULL;
				ray_next=NULL;
			}
			else
			{
				//if(flag_ray_next!=1)//should only change to transmitted ray
				//{
				//	ray_next->trace_obj=intersect->inter_obj->belongobj->subobj;
				//	ray_next->refr_index=intersect->inter_obj->refr_index_o;
				//}
				trace(ray_next,ray_next->trace_obj);
			}
		}	
	}
}
