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
*    big update during the development of eLeaf_dicot 2021 June
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
	
	double det,abs_t1,abs_t2,tmp_scale;
	det=VecDot(E1,P);
	
	Vec T;
	if(det>0)
	{
		VecSub(v0,orig,T);
	}
	else
	{
		//sign of det is affected by the sequence of S->tri and dir
		// is not affected by orig, i.e. does not reflect sign of dis
		VecSub(orig,v0,T);
		det=-det;
	}
	
	// if determinant is near zero, ray and triangle are parallel
	abs_t1=VecDot(E1,E1);
	abs_t2=VecDot(P,P);
	tmp_scale=sqrt(abs_t1)*sqrt(abs_t2);
	if(det<(DBL_EPSILON*tmp_scale))
	{
		dis[0]=-Infinity;//dis<0 means no intersection
		//printf("det:%e\n",det);
		return 0;
	}
	
	// calculate u and make sure u<=det
	double u;
	u=VecDot(T,P);
	abs_t1=fabs(u);
	abs_t2=fabs(det);
	tmp_scale=fmax(abs_t1,abs_t2);
	//if(u<0.0||u>det) //[revise June 2021: numerical robustness] (u-det)>-max(fabs(u),fabs(det))*DBL_EPSILON
	if(u<0.0||(u-det)>(tmp_scale*DBL_EPSILON))
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
	abs_t1=fabs(u+v);
	abs_t2=fabs(det);
	tmp_scale=fmax(abs_t1,abs_t2);
	//if(v<0.0||(u+v)>det) //[revise June 2021: numerical robustness] (u+v-det)>max(fabs(u+v),fabs(det))*DBL_EPSILON
	if(v<0.0||(u+v-det)>(tmp_scale*DBL_EPSILON))
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

void ObjIsect(Ray *ray,Object *obj,Isect *intersect,int outin_ray2obj)
//outin_ray2obj=0 means ray is outside of objs, outin_ray2obj=1 means ray is inside (i.e. detect intersection with mother node)
//outin_ray2obj doesn't affect the calculation of *intersect, but affects the double check of temp_dis by VecDoc(ray->D,tri->Nor)
{
	Object *tmp_obj=obj;
	int i,flag_intersect,tmp_flag;//,count,tmp_flag;
	int inter_i_front=-1,inter_i_back=-1,inter_i_neg=-1;//inter_i
	i=0;
	//count=0;
	flag_intersect=0;
	double temp_dis[1],t_neg,tmp_dbl;
	double t_front,t_back;//t
	//t=Infinity;
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
		
		if(tmp_flag==1&&temp_dis[0]>-DIS_EPSILON)//previously temp_dis[0]>0.0
		//NOTICE: although ray->P_tri is skipped, however, it's still possible to detect a tri next to ray->P_tri if ray->P is on the edge for example
		{
			flag_intersect=1;
			//if(temp_dis[0]<t)
			//{
			//	t=temp_dis[0];
			//	inter_i=i;
			//}
			tmp_dbl=VecDot(ray->D,(tmp_obj->S->tri+i)->Nor);
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
		else// previously if(temp_dis[0]<0.0)
		{
			if(temp_dis[0]>t_neg)
			{
				//TriIsect() can return both positive and negative dis
				t_neg=temp_dis[0];
				inter_i_neg=i;
			}
		}
		i++;
		if(tmp_obj==ray->P_obj&&i==ray->P_tri)//overlap with line134-line137, consider i==ray->P_tri==tmp_obj->S->pnum
		{
			i++;
		}
	}
	
	if(outin_ray2obj==1)//ray is inside
	{
		if(flag_intersect==1)
		{
			if(inter_i_front!=-1)
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
				if(flag_errormsg==1)
				{
					printf("[ErrorMsg] ObjIsect()-02: ray is inside, but didn't find an intersection\n");
					printf("    tfront %e;tback %e;tneg %e\n",t_front,t_back,t_neg);
					printf("    trifront %d;triback %d;trineg %d\n",inter_i_front,inter_i_back,inter_i_neg); 
					printf("    ray->P:%e %e %e\n",ray->P[0],ray->P[1],ray->P[2]);
					printf("    ray->D:%e %e %e\n",ray->D[0],ray->D[1],ray->D[2]);
					printf("    Intersect obj info. Suggest plot ray and obj in MATLAB.\n");
					printf("    cell type:%d\n",obj->celltype);
					printf("    cell name:%d\n",obj->cellname);
					printf("    chl name:%d\n",obj->chlname);
					printf("    vac name :%d\n",obj->vacuolename);
				}
				exit(1);
			}
		}
		else//no intersection, and ray is inside...
		{
			if(flag_errormsg==1)
			{
				printf("[ErrorMsg] ObjIsect()-03: ray is inside, but didn't find an intersection\n");
				printf("    tfront %e;tback %e;tneg %e\n",t_front,t_back,t_neg);
				printf("    trifront %d;triback %d;trineg %d\n",inter_i_front,inter_i_back,inter_i_neg); 
				printf("    ray->P:%e %e %e\n",ray->P[0],ray->P[1],ray->P[2]);
				printf("    ray->D:%e %e %e\n",ray->D[0],ray->D[1],ray->D[2]);
				printf("    Intersect Obj info. Potential singular point in Mesh. \n    Suggest repeat the bug under RT debug mode, \n    then plot ray, objs and intersection in MATLAB.\n");
				printf("    cell type:%d\n",obj->celltype);
				printf("    cell name:%d\n",obj->cellname);
				printf("    chl name:%d\n",obj->chlname);
				printf("    vac name :%d\n",obj->vacuolename);
			}
			exit(1);
		}
	}
	else//outin_ray2obj==0; i.e. ray is outside this obj
	{
		if(flag_intersect==0)
		{
			intersect->t=-Infinity;
		}
		else
		{
			//tmp_dbl=-VecDot(ray->D,(tmp_obj->S->tri+inter_i)->Nor);
			//if(tmp_dbl>0.0)
			//choose inter_i_back and t_back here
			if(inter_i_back!=-1)
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
			}
			else
			{
				if(flag_warningmsg==1)
				{
					printf("[WarningMsg] ObjIsect()-04: ray is outside, but only found intersection with triangle facing away\n");
					printf("    tfront %e;tback %e;tneg %e\n",t_front,t_back,t_neg);
					printf("    trifront %d;triback %d;trineg %d\n",inter_i_front,inter_i_back,inter_i_neg); 
					printf("    ray->P:%e %e %e\n",ray->P[0],ray->P[1],ray->P[2]);
					printf("    ray->D:%e %e %e\n",ray->D[0],ray->D[1],ray->D[2]);
					printf("    Intersect Obj info. Suggest plot ray and obj in MATLAB.\n");
					printf("    cell type:%d\n",obj->celltype);
					printf("    cell name:%d\n",obj->cellname);
					printf("    chl name:%d\n",obj->chlname);
					printf("    vac name :%d\n",obj->vacuolename);
					//exit(1);//This is possible, when ray->P is on the contact surface of A and B
					//ray is inside A, algorithm found an intersection of ray with A
					//then it is searching neighbour with the same dis, but accidentally found B.
				}
				//exit(1);
				flag_intersect=0;
				intersect->t=-Infinity;
			}
		}
	}
}

void compute_ray_next(Ray *ray,Isect *intersect0,Isect *intersect1, double *level_contact_refr,int count_contact,Ray *ray_next,int tmp_outin_ray2obj,int tmp_outin_ray2obj_tr)
{
	//calculate r1/t1 r2/t2... rn/tn
	int loop_i,TIR=0,flag_stop=0,tmp_level,flag_dir_downward;
	double tmp_n1,tmp_n2,rs,rp,reflectance[10];//,transmittance[10];
	Vec D1,norm,tmp_norm1,tmp_norm2;
	double sin_t[10],cos_t[10],sin_t1,cos_t1,sin_t2,cos_t2,tmp_dbl,tmp_rnd,tmp_r12;
	
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
	for(loop_i=0;loop_i<(count_contact-1)&&TIR==0;loop_i++)//number of media=count_contact, number of interface=count_contact-1
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
			reflectance[loop_i]=(rs+rp)/2;//size of reflectacne=count_contact-1
			//transmittance[loop_i]=1-reflectance[loop_i];
			cos_t[loop_i+1]=cos_t2;
			sin_t[loop_i+1]=sin_t2;
		}
	}
	
	if(TIR!=1)//not total reflection, rand until reflected at boundary 1 or transmitted at boundary n
	{
		tmp_level=0;
		flag_dir_downward=1;//downward is TRUE
		do
		{
			tmp_rnd=(double)(rand()*1.0/RAND_MAX);
			if(flag_dir_downward==1)
				tmp_r12=reflectance[tmp_level];
			else
				tmp_r12=reflectance[tmp_level-1];
			if(tmp_rnd<tmp_r12)
			{
				tmp_level--;
				flag_dir_downward=1-flag_dir_downward;
			}
			else
			{
				tmp_level++;
			}
			
			if(tmp_level<=0||tmp_level>=count_contact-1)
			{
				flag_stop=1;
				if(tmp_level<=0)
					TIR=1;//here TIR became flag_reflectance
			}
		}while(flag_stop==0);
	}
	
	if(TIR==1)//here TIR became flag_reflectance
	{
		if(flag_RTdebug_printf==1)
		{
			printf("    Ray - reflection\n");
		}
		ray_next->P[0]=intersect0->inter_pt[0];
		ray_next->P[1]=intersect0->inter_pt[1];
		ray_next->P[2]=intersect0->inter_pt[2];
		tmp_dbl=cos_t[0];
		ray_next->D[0]=D1[0]-2*norm[0]*tmp_dbl;
		ray_next->D[1]=D1[1]-2*norm[1]*tmp_dbl;
		ray_next->D[2]=D1[2]-2*norm[2]*tmp_dbl;
		ray_next->P_obj=intersect0->inter_obj;
		ray_next->P_tri=intersect0->inter_tri;
		if(tmp_outin_ray2obj==0)
		{
			ray_next->outin_P=0;
			ray_next->trace_obj=ray->trace_obj;
		}
		else
		{
			ray_next->outin_P=1;
			ray_next->trace_obj=ray->trace_obj;
		}
	}
	else//refracted ray
	{
		if(flag_RTdebug_printf==1)
		{
			printf("    Ray - refraction\n");
		}
		ray_next->P[0]=intersect1->inter_pt[0];
		ray_next->P[1]=intersect1->inter_pt[1];
		ray_next->P[2]=intersect1->inter_pt[2];
		
		VecCross(ray->D,norm,tmp_norm1);
		Normalize(tmp_norm1);
		VecCross(norm,tmp_norm1,tmp_norm2);
		cos_t2=cos_t[count_contact-1];
		sin_t2=sin_t[count_contact-1];
		ray_next->D[0]=norm[0]*cos_t2+tmp_norm2[0]*sin_t2;
		ray_next->D[1]=norm[1]*cos_t2+tmp_norm2[1]*sin_t2;
		ray_next->D[2]=norm[2]*cos_t2+tmp_norm2[2]*sin_t2;
		ray_next->P_obj=intersect1->inter_obj;
		ray_next->P_tri=intersect1->inter_tri;
		if(tmp_outin_ray2obj_tr==0)
		{
			ray_next->outin_P=1;
			ray_next->trace_obj=intersect1->inter_obj->belongobj->subobj;
		}
		else
		{
			ray_next->outin_P=0;
			ray_next->trace_obj=intersect1->inter_obj->subobj;
		}
	}
	
}

void trace(Ray *ray, Object *obj)
{
	Object *curobj;
	Isect *intersect;	
	intersect=(Isect *)malloc(sizeof(Isect));
	if(intersect==NULL)
		exit(1);
	intersect->t=-Infinity;
	Isect *tmp_intersect;
	tmp_intersect=(Isect *)malloc(sizeof(Isect));
	if(tmp_intersect==NULL)
		exit(1);
	tmp_intersect->t=-Infinity;
	Isect *intersect_cld;
	intersect_cld=(Isect *)malloc(sizeof(Isect));//used when searching children node
	if(intersect_cld==NULL)
		exit(1);
	intersect_cld->t=-Infinity;
	Isect *intersect_mth;
	intersect_mth=(Isect *)malloc(sizeof(Isect));//used when searching mother node
	if(intersect_mth==NULL)
		exit(1);
	intersect_mth->t=-Infinity;
	Isect *intersect_bk;
	intersect_bk=(Isect *)malloc(sizeof(Isect));//backup
	if(intersect_bk==NULL)
		exit(1);
	Ray *ray_next;//Ray *ray_refl,*ray_trans;
	ray_next=(Ray *)malloc(sizeof(Ray));

	int flag_stop,flag_stop2,flag_outin_ray2obj=0,flag_outin_ray2obj_tr,flag_contact,loop_i,count_contact;//,tmp_level,tmp_count_debug;
	//int tmp_level_contact[10];//buff size 10
	double tmp_level_contact_refr[10];//record the refractive index of each level, since chlo has different refr to vacu
	for(loop_i=0;loop_i<10;loop_i++)
	{
		//tmp_level_contact[loop_i]=-1;//level=0=hit leaf;=1=IAS;=2=cyto;=3=chlo/vacu
		tmp_level_contact_refr[loop_i]=-1.0;
	}
	
	if(obj==NULL)
	{
		intersect->t=-Infinity;//TO IMPROVE - full initial?
		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-Intersect: Trace_dicto.c -- Line 531\n");
			printf("    NULL object. No intersection\n");
		}
	}
	else
	{
		/****************************************************/
		/*1. search intersection with obj->belongobj->subobj*/
		/****************************************************/
		// assume ray is outside of all objs
		// TO IMPROVE: force obj to be the first node of a sequential of nodes at the same layer?
		curobj=obj;
		if(obj->belongobj==NULL)//in case obj=leaf
		{
			curobj=obj;
		}
		else
		{
			curobj=(obj->belongobj)->subobj;
		}
		
		while(curobj!=NULL)
		{
			ObjIsect(ray,curobj,tmp_intersect,flag_outin_ray2obj);//NOTICE: outin_ray2obj=0 means ray is outside for all child nodes here
			if(tmp_intersect->t > -DIS_EPSILON)
			{
				// for rice MSC with lobes, it's possible tmp_intersect->inter_obj==ray->P_obj.
				// so don't skip curobj==ray->P_obj here.
				if(tmp_intersect->t > DIS_EPSILON)
				{
					if(intersect->t < -DIS_EPSILON)//basically, this means intersect->t is still the initial value -Infinity
					{
						IsectCopy(tmp_intersect,intersect);
					}
					else
					{
						if(tmp_intersect->t < intersect->t)
						{
							IsectCopy(tmp_intersect,intersect);
						}
					}
				}
				else
				{
					if(tmp_intersect->inter_obj==ray->P_obj && tmp_intersect->inter_tri==ray->P_tri)
					{
						//possible if ray->P_obj is one of the child node
					}
					else
					{
						if(flag_warningmsg==1)
						{
							printf("WarningMsg-Trace()-01: |tmp_intersect->t| < DIS_EPSILON and tmp_intersect->inter_obj!=ray->P_obj at the tri level\n");
							printf("    tmp_intersect->t=%e\n",tmp_intersect->t);
							printf("    tmp_intersect->inter_tri=%d\n",tmp_intersect->inter_tri);
							printf("    tmp_intersect->inter_pt:[%e %e %e]\n",tmp_intersect->inter_pt[0],tmp_intersect->inter_pt[1],tmp_intersect->inter_pt[2]);
							printf("    tmp_intersect->inter_obj->celltype=%d\n",tmp_intersect->inter_obj->celltype);
							printf("    tmp_intersect->inter_obj->cellname=%d\n",tmp_intersect->inter_obj->cellname);
							printf("    tmp_intersect->inter_obj->chlname=%d\n",tmp_intersect->inter_obj->chlname);
							printf("    tmp_intersect->inter_obj->vacuolename=%d\n",tmp_intersect->inter_obj->vacuolename);
						}
						//VerySmallIntersectDistance
						if(tmp_intersect->t > -DBL_EPSILON)
						{
							if(intersect->t < -DIS_EPSILON)//basically, this means intersect->t is still the initial value -Infinity
							{
								IsectCopy(tmp_intersect,intersect);
							}
							else
							{
								if(tmp_intersect->t < intersect->t)
								{
									IsectCopy(tmp_intersect,intersect);
								}
							}
						}
						//else, no intersection
						//exit(1);
					}
				}
			}
			curobj=curobj->nextobj;
		}
	}
	if(flag_RTdebug_printf==1)
	{
		printf("\nDEBUGINFO-Intersect: Trace_dicto.c -- Line 594\n");
		if(intersect==NULL||(intersect->t+Infinity)<DBL_EPSILON)
		{
			printf("    intersect=NULL\n");
		}
		else
		{
			printf("    intersect->t=%e\n",intersect->t);
			printf("    intersect->inter_tri=%d\n",intersect->inter_tri);
			printf("    intersect->inter_pt:[%e %e %e]\n",intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
			printf("    intersect->inter_obj->celltype=%d\n",intersect->inter_obj->celltype);
			printf("    intersect->inter_obj->cellname=%d\n",intersect->inter_obj->cellname);
			printf("    intersect->inter_obj->chlname=%d\n",intersect->inter_obj->chlname);
			printf("    intersect->inter_obj->vacuolename=%d\n",intersect->inter_obj->vacuolename);
		}
		//printf("raylevel:%d\n",ray->level);
	}
	if(intersect->t > -DIS_EPSILON)
	/***********************************************************************/
	/*1-1.[Intersection found] search child and child's child              */
	/*    until 1.NULL i.e. no children node; or 2.no contact with any objs*/
	/***********************************************************************/
	// >-DIS_EPSILON refers line 138 and line 539
        // no intersection if intersect->t==-Infinity
	{
		//calculate absorption - update 2021 Jun
		if(ray->outin_P==1)
		{
			////P of ray is on the inner surface, i.e. ray is from the mother node
			ray_next->I=ray->I*exp(-(ray->P_obj->sac_i)*(intersect->t));
			ray->P_obj->I_absorb_i+=(ray->I-ray_next->I);
			(ray->P_obj->S->tri+ray->P_tri)->I_absorb_tri_i+=(ray->I-ray_next->I);
			//NEW for RT training - added 2021 April
			if(flag_opt_absorbevents==1)
			{
				fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,intersect->t,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,0);
			}
		}
		else
		{
			////P of ray is on the outter surface, i.e. ray is reflected from a neighbour node
			ray_next->I=ray->I*exp(-(ray->P_obj->sac_o)*(intersect->t));
			intersect->inter_obj->I_absorb_o+=(ray->I-ray_next->I);
			(intersect->inter_obj->S->tri+intersect->inter_tri)->I_absorb_tri_o+=(ray->I-ray_next->I);
			//NEW for RT training - added 2021 April
			if(flag_opt_absorbevents==1)
			{
				fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,intersect->t,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,0);
			}
		}	

		// update refr index stack and ray_level
		count_contact=0;
		tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
		count_contact++;
		if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
		{
			//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
			tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
			count_contact++;
		}
		tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
		count_contact++;

		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-SearchContact: Trace_dicto.c -- Line 663\n");
			printf("    Search contact obj start - child node\n");
		}
		IsectCopy(intersect,intersect_bk);//backup for refl&refr calculation on contact boundary
		intersect_cld->inter_obj=intersect->inter_obj;
		flag_stop=0;
		flag_outin_ray2obj=0;
		do
		{
			curobj=intersect_cld->inter_obj->subobj;
			intersect_cld->t=-Infinity;//initial
			if(curobj==NULL)
			{
				flag_stop=1;
			}
			else
			{
				flag_stop2=0;
				while(curobj!=NULL&&flag_stop2==0)
				//intersect_cld record closest intersection with child node
				{
					ObjIsect(ray,curobj,tmp_intersect,flag_outin_ray2obj);
					if(tmp_intersect->t > -DIS_EPSILON)
					{
						if(intersect_cld->t < -DIS_EPSILON||tmp_intersect->t < intersect_cld->t)
						{
							IsectCopy(tmp_intersect,intersect_cld);
						}
						if((intersect_cld->t-intersect->t)>-DIS_EPSILON&&(intersect_cld->t-intersect->t)<DIS_EPSILON)
						{
							flag_stop2=1;//stop search child only if get a contact child here
							IsectCopy(intersect_cld,intersect);
							if(flag_RTdebug_printf==1)
							{
								printf("    Find a neighbour's child in contact\n");
							}

							// update refr index stack and ray_level
							//tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
							//count_contact++;
							if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
							{
								//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
								tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
								count_contact++;
							}
							tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
							count_contact++;
						}
						else if((intersect_cld->t-intersect->t)<-DIS_EPSILON)
						{
							if(flag_errormsg==1)
							{
								printf("ErrorMsg-trace()-01: Ray is outside. detect intersect_cld is closer than intersect\n");
							}
							exit(1);
						}
					
					}
					else
					{
						if((tmp_intersect->t+Infinity) > DBL_EPSILON)
						{
							if(flag_errormsg==1)
							{
								printf("ErrorMsg-trace()-02: ObjIsect returns dis <= -DIS_EPSILON and it's not the intial \"-Infinity\"\n");
							}
							exit(1);
						}
					}
					curobj=curobj->nextobj;
				}
				if(flag_stop2!=1)
					flag_stop=1;//stop search cld's chl if no contact chd found
			}
		}while(flag_stop==0);
		if(flag_RTdebug_printf==1)
		{
			printf("    Search contact obj finish - child node\n");
		}
		
		//refr index stack done
		//calculate reflected and transmitted ray
		flag_outin_ray2obj=0;
		flag_outin_ray2obj_tr=1;
		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-MonteCarlo RefrIndex: Trace_dicto.c -- Line 739\n");
			printf("    refractive index stack:\n    ");
			for(loop_i=0;loop_i<count_contact;loop_i++)
			{
				printf("%e ",tmp_level_contact_refr[loop_i]);
			}
			printf("\n");
		}
		compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact_refr,count_contact,ray_next,flag_outin_ray2obj,flag_outin_ray2obj_tr);//flag_outin_ray2obj still=0 here
		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-ray_next: Trace_dicto.c -- Line 750\n");
			//printf("ray_next level:%d\n",ray_next->level);
			printf("    ray_next->P:[%e %e %e]\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
			printf("    ray_next->D:[%e %e %e]\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
			printf("    ray_next->I:%e\n",ray_next->I);
			printf("    ray_next->P_tri:%d\n",ray_next->P_tri);
			if(ray_next->trace_obj!=NULL)
			{
				printf("    ray_next->trace_obj->celltype:%d\n",ray_next->trace_obj->celltype);
				printf("    ray_next->trace_obj->cellname:%d\n",ray_next->trace_obj->cellname);
				printf("    ray_next->trace_obj->chlname:%d\n",ray_next->trace_obj->chlname);
				printf("    ray_next->trace_obj->vacuolename:%d\n",ray_next->trace_obj->vacuolename);
			}
			else
				printf("    ray_next->trace_obj->celltype:NULL\n");
		}
		if(flag_RT_file4plot==1)
		{
			fprintf(fout_file4plot,"%d %d %e %e %e\n",ray_i,ray_j,intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
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
			if(flag_RTdebug_printf==1)
			{
				printf("\n*** Trace() exit: ray discarded ***\n");
				printf("****************END****************\n");
			}
		}
		else
		{
			if(flag_RTdebug_printf==1)
			{
				printf("\n-------------------------------------------\n");
				printf("--- Trace(ray_next,ray_next->trace_obj) ---\n");
				printf("-------------------------------------------\n");
			}
			trace(ray_next,ray_next->trace_obj);
		}
	}
	else
	/***********************************************************************/
	/*1-2.[Intersection not found] intersect with mother node              */
	/*    and then search mother's neighbour, then neigbhour's child       */
	/*    until 1.NULL i.e. no children node; or 2.no contact with any objs*/
	/***********************************************************************/
	{
		if(obj==NULL)
			curobj=ray->P_obj;
		else
			curobj=obj->belongobj;

		flag_outin_ray2obj=1;//ray is inside relative to the mother node
		flag_outin_ray2obj_tr=0;
		ObjIsect(ray,curobj,intersect,flag_outin_ray2obj);//intersect with mother node
		if(intersect->t < -DIS_EPSILON)
		{
			if(flag_errormsg==1)
			{
				printf("[ErrorMsg] trace()-04: No intersection found when ObjIsect() with mother node\n");
				printf("    intersect->t=%e\n",intersect->t);
				printf("    intersect->inter_tri=%d\n",intersect->inter_tri);
				printf("    intersect->inter_pt:[%e %e %e]\n",intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
				printf("    intersect->inter_obj->celltype=%d\n",intersect->inter_obj->celltype);
				printf("    intersect->inter_obj->cellname=%d\n",intersect->inter_obj->cellname);
				printf("    intersect->inter_obj->chlname=%d\n",intersect->inter_obj->chlname);
				printf("    intersect->inter_obj->vacuolename=%d\n",intersect->inter_obj->vacuolename);				
			}
			exit(1);
		}
		else if(intersect->t < DIS_EPSILON)
		{
			if(flag_warningmsg==1)
			{
				printf("[WarningMsg] trace()-04: |intersect->t| < DIS_EPSILON found when ObjIsect() with mother node\n");
				printf("    intersect->t=%e\n",intersect->t);
				printf("    intersect->inter_tri=%d\n",intersect->inter_tri);
				printf("    intersect->inter_pt:[%e %e %e]\n",intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
				printf("    intersect->inter_obj->celltype=%d\n",intersect->inter_obj->celltype);
				printf("    intersect->inter_obj->cellname=%d\n",intersect->inter_obj->cellname);
				printf("    intersect->inter_obj->chlname=%d\n",intersect->inter_obj->chlname);
				printf("    intersect->inter_obj->vacuolename=%d\n",intersect->inter_obj->vacuolename);				
			}
		}

		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-Intersect: Trace_dicto.c -- Line 832\n");
			if(intersect==NULL)
				printf("    intersect=NULL - Impossible\n");
			else
			{
				printf("    intersect->t=%e\n",intersect->t);
				printf("    intersect->inter_tri=%d\n",intersect->inter_tri);
				printf("    intersect->inter_pt:[%e %e %e]\n",intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
				printf("    intersect->inter_obj->celltype=%d\n",intersect->inter_obj->celltype);
				printf("    intersect->inter_obj->cellname=%d\n",intersect->inter_obj->cellname);
				printf("    intersect->inter_obj->chlname=%d\n",intersect->inter_obj->chlname);
				printf("    intersect->inter_obj->vacuolename=%d\n",intersect->inter_obj->vacuolename);
			}
		}		

		if(flag_RT_file4plot==1)
		{
			fprintf(fout_file4plot,"%d %d %e %e %e\n",ray_i,ray_j,intersect->inter_pt[0],intersect->inter_pt[1],intersect->inter_pt[2]);
		}

		//calculate absorption - update 2021 Jun
		if(ray->outin_P==1)
		{
			////P of ray is on the inner side of surface, i.e. ray is from the mother node, and now intersect with mother node
			ray_next->I=ray->I*exp(-(ray->P_obj->sac_i)*(intersect->t));
			ray->P_obj->I_absorb_i+=(ray->I-ray_next->I);
			(ray->P_obj->S->tri+ray->P_tri)->I_absorb_tri_i+=(ray->I-ray_next->I);
			//NEW for RT training - added 2021 April
			if(flag_opt_absorbevents==1)
			{
				fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,intersect->t,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,0);
			}
		}
		else
		{
			////P of ray is on the outter side of surface, i.e. ray is reflected from a neighbour node, and now itnersect with mother node
			ray_next->I=ray->I*exp(-(ray->P_obj->sac_o)*(intersect->t));
			ray->P_obj->I_absorb_o+=(ray->I-ray_next->I);
			(ray->P_obj->S->tri+ray->P_tri)->I_absorb_tri_o+=(ray->I-ray_next->I);
			//NEW for RT training - added 2021 April
			if(flag_opt_absorbevents==1)
			{
				fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,intersect->t,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,0);
			}
		}

		// update refr index stack and ray_level
		count_contact=0;
		tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
		count_contact++;
		if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
		{
			//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
			tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
			count_contact++;
		}
		tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
		count_contact++;
		
		//search mother's mother etc
		IsectCopy(intersect,intersect_bk);//backup for boundary condition
		intersect_mth->inter_obj=curobj;//yes, obj can not be leaf here! but curobj can be
		if(flag_RTdebug_printf==1)
		{
			printf("\nDEBUGINFO-SearchContact: Trace_dicto.c -- Line 894\n");
			printf("    Search contact obj start - mother node\n");
		}
		flag_stop=0;
		flag_contact=0;
		do
		{
			curobj=intersect_mth->inter_obj->belongobj;
			if(curobj==NULL)//intersect->inter_obj=leaf
			{
				flag_stop=1;
			}
			else
			{
				ObjIsect(ray,curobj,intersect_mth,flag_outin_ray2obj);
				if((intersect_mth->t-intersect->t)>-DIS_EPSILON&&(intersect_mth->t-intersect->t)<DIS_EPSILON)
				{
					flag_stop=0;//find a mother's mother contact with mother...
					flag_contact=1;
					if(flag_RTdebug_printf==1)
					{
						printf("    Find intersect_tri contact with its mother\n");
					}
					IsectCopy(intersect_mth,intersect);

					// update refr index stack and ray_level
					//tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
					//count_contact++;
					if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
					{
						//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
						tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
						count_contact++;
					}
					tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
					count_contact++;
					
				}
				else if((intersect_mth->t-intersect->t)<-DIS_EPSILON)
				{
					if(flag_errormsg==1)
					{
						printf("[ErrorMsg] trace()-05: Ray is inside. detect intersect_mth is closer than intersect\n");
					}
					exit(1);
				}
				else
				{
					flag_stop=1;
				}
			}
		}while(flag_stop==0);
		
		if(flag_RTdebug_printf==1)
		{
			printf("    Search contact obj finish - mother node\n\n");
		}		

		//then search mother's neighbour
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
					tmp_level_contact_refr[1]=0.0001;//make sure total reflection; force ray back to the original obj no matter how many shared surfaces it has passed through
					flag_outin_ray2obj=1;
					flag_outin_ray2obj_tr=0;//but if found a contacted neigbour or even a neigbour's child, then flag_outin_ray2obj_tr becomes 1 again.
					if(flag_RTdebug_printf==1)
					{
						printf("\nDEBUGINFO-MonteCarlo RefrIndex: Trace_dicto.c -- Line 955\n");
						for(loop_i=0;loop_i<count_contact;loop_i++)
						{
							printf("%e ",tmp_level_contact_refr[loop_i]);
						}
						printf("\n");
					}
					compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact_refr,count_contact,ray_next,flag_outin_ray2obj,flag_outin_ray2obj_tr);//use intersect_bk to calculate ray_next
					if(flag_RTdebug_printf==1)
					{
						printf("\nDEBUGINFO-Mirror Boundary: Trace_dicto.c -- Line 965\n");
						printf("    ray_next->P:[%e %e %e]\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
						printf("    ray_next->D:[%e %e %e]\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
						printf("    ray_next->I:%e\n",ray_next->I);
						printf("    ray_next->P_tri:%d\n",ray_next->P_tri);
						if(ray_next->trace_obj!=NULL)
						{
							printf("    ray_next->trace_obj->celltype:%d\n",ray_next->trace_obj->celltype);
							printf("    ray_next->trace_obj->cellname:%d\n",ray_next->trace_obj->cellname);
							printf("    ray_next->trace_obj->chlname:%d\n",ray_next->trace_obj->chlname);
							printf("    ray_next->trace_obj->vacuolename:%d\n",ray_next->trace_obj->vacuolename);
						}
						else
							printf("    ray_next->trace_obj->celltype:NULL\n");
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
						if(flag_RTdebug_printf==1)
						{
							printf("\n*** Trace() exit: ray discarded ***\n");
							printf("****************END****************\n");
						}
					}
					else
					{
						if(flag_RTdebug_printf==1)
						{
							printf("\n-------------------------------------------\n");
							printf("--- Trace(ray_next,ray_next->trace_obj) ---\n");
							printf("-------------------------------------------\n");
						}
						//ray_next->trace_obj=ray->trace_obj;
						trace(ray_next,ray_next->trace_obj);
					}
					break;
				case 8:
				case 9:
				I_discard_Tr+=ray_next->I;
				//NEW for RT training - added 2021 July
				//although essentially no absorption, but this is for recalculate total reflectance and transmittance
				if(flag_opt_absorbevents==1)
				{
					// additional %d at last represents ray status, 0=absorb, 1=reflect, 2=transmit
					fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,0.0,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,2);
				}
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
				if(flag_RTdebug_printf==1)
				{
					printf("*** Trace() exit: transmitance ***\n");
					printf("****************END***************\n");
				}
				break;
				case 4:
				case 5:
				I_discard_Rf+=ray_next->I;
				//NEW for RT training - added 2021 July
				//although essentially no absorption, but this is for recalculate total reflectance and transmittance
				if(flag_opt_absorbevents==1)
				{
					// additional %d at last represents ray status, 0=absorb, 1=reflect, 2=transmit
					fprintf(fout_absorbevents,"%d %d %e %d %d %d %d %d %d %e %d\n",ray_i,ray_j,0.0,ray->P_obj->celltype,ray->P_obj->cellname,ray->P_obj->chlname,ray->P_obj->vacuolename,ray->P_tri,ray->outin_P,ray_next->I,1);
				}
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
				if(flag_RTdebug_printf==1)
				{
					printf("*** Trace() exit: reflectance ***\n");
					printf("****************END**************\n");
				}
				break;
			}
		}
		else//[1-2-2] and if mother node is not leaf boundary. then try search mother's neighbour
		{
			flag_outin_ray2obj=1;
			flag_outin_ray2obj_tr=0;
			//here if no more contact, then should finished with refr_index_o
			tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
			//otherwise, contact with neighbour or child of neighbour, then should finished with refr_index_i
			curobj=intersect->inter_obj->belongobj->subobj;
			flag_stop=0;
			flag_contact=0;
			if(flag_RTdebug_printf==1)
			{
				printf("\n    Search contact obj start - neighbour node\n");
			}
			while(curobj!=NULL&&flag_stop==0)
			{
				if(curobj!=intersect->inter_obj)//skip intersect->inter_obj
				{
					ObjIsect(ray,curobj,tmp_intersect,0);
					if((tmp_intersect->t-intersect->t)>-DIS_EPSILON&&(tmp_intersect->t-intersect->t)<DIS_EPSILON)
					{
						flag_stop=1;//impossible exist multi obj contact at the same position
						flag_contact=1;
						if(flag_RTdebug_printf==1)
						{
							printf("    Find intersect_tri contact with its neighbour\n");
						}
						IsectCopy(tmp_intersect,intersect);

						// update refr index stack and ray_level
						//tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
						//count_contact++;
						if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
						{
							//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
							tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
							count_contact++;
						}
						tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
						count_contact++;
							
						flag_outin_ray2obj_tr=1;
					}
				}
				curobj=curobj->nextobj;
			}
			if(flag_RTdebug_printf==1)
			{
				printf("    Search contact obj finish - neighbour node\n");
			}

			if(flag_contact==1)//find a contact neighbour, then search it's child and child's child etc
			{
				if(flag_RTdebug_printf==1)
				{
					printf("\n    Search contact obj start - child node\n");
				}
				intersect_cld->inter_obj=intersect->inter_obj;
				flag_stop=0;
				do
				{
					curobj=intersect_cld->inter_obj->subobj;
					intersect_cld->t=-Infinity;//initial
					if(curobj==NULL)
					{
						flag_stop=1;
					}
					else
					{
						flag_stop2=0;
						while(curobj!=NULL&&flag_stop==0)
						{
							ObjIsect(ray,curobj,tmp_intersect,0);
							if(tmp_intersect->t > -DIS_EPSILON)
							{
								if(intersect_cld->t < -DIS_EPSILON||tmp_intersect->t < intersect_cld->t)
								{
									IsectCopy(tmp_intersect,intersect_cld);
								}
								if((intersect_cld->t-intersect->t)>-DIS_EPSILON&&(intersect_cld->t-intersect->t)<DIS_EPSILON)
								{
									flag_stop2=1;//stop search child nodes only if find a contact child here
									IsectCopy(intersect_cld,intersect);
									if(flag_RTdebug_printf==1)
									{
										printf("    Find a neighbour's child in contact\n");
									}
									
									// update refr index stack and ray_level
									//tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_o;
									//count_contact++;
									if(intersect->inter_obj->refr_index_s>0.0)//i.e. refr_index_s!=-1, i.e. refr_index_s is not inactive
									{
										//in this way, avoid to specify additional refr index for cell wall and contacted chlo envelope
										tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_s;
										count_contact++;
									}
									tmp_level_contact_refr[count_contact]=intersect->inter_obj->refr_index_i;
									count_contact++;

									flag_outin_ray2obj_tr=1;
								}
								else if((intersect_cld->t-intersect->t)<-DIS_EPSILON)
								{
									if(flag_errormsg==1)
									{
										printf("ErrorMsg-trace()-06: Ray is outside. detect intersect_cld is closer than intersect\n");
									}
									exit(1);
								}
							}
							else
							{
								if((tmp_intersect->t+Infinity) > DBL_EPSILON)
								{
									if(flag_errormsg==1)
									{
										printf("ErrorMsg-trace()-07: ObjIsect returns dis <= -DIS_EPSILON and it's not the intial \"-Infinity\"\n");
									}
									exit(1);
								}
							}
							curobj=curobj->nextobj;
						}
						if(flag_stop2!=1)
							flag_stop=1;//stop search cld's chl if no contact chd found
					}
				}while(flag_stop==0);
				
				if(flag_RTdebug_printf==1)
				{
					printf("    Search contact obj finish - child node\n");
				}
			}
			
			//compute ray_next
			if(flag_RTdebug_printf==1)
			{
				printf("\nDEBUGINFO-MonteCarlo RefrIndex: Trace_dicto.c -- Line 1191\n");
				for(loop_i=0;loop_i<count_contact;loop_i++)
				{
					printf("%e ",tmp_level_contact_refr[loop_i]);
				}
				printf("\n");
			}
			compute_ray_next(ray,intersect_bk,intersect,tmp_level_contact_refr,count_contact,ray_next,flag_outin_ray2obj,flag_outin_ray2obj_tr);
			if(flag_RTdebug_printf==1)
			{
				printf("\nDEBUGINFO-ray_next: Trace_dicto.c -- Line 1201\n");
				//printf("    ray_next level:%d\n",ray_next->level);
				printf("    ray_next->P:[%e %e %e]\n",ray_next->P[0],ray_next->P[1],ray_next->P[2]);
				printf("    ray_next->D:[%e %e %e]\n",ray_next->D[0],ray_next->D[1],ray_next->D[2]);
				printf("    ray_next->I:%e\n",ray_next->I);
				printf("    ray_next->P_tri:%d\n",ray_next->P_tri);
				if(ray_next->trace_obj!=NULL)
				{
					printf("    ray_next->trace_obj->celltype:%d\n",ray_next->trace_obj->celltype);
					printf("    ray_next->trace_obj->cellname:%d\n",ray_next->trace_obj->cellname);
					printf("    ray_next->trace_obj->chlname:%d\n",ray_next->trace_obj->chlname);
					printf("    ray_next->trace_obj->vacuolename:%d\n",ray_next->trace_obj->vacuolename);
				}
				else
					printf("    ray_next->trace_obj->celltype:NULL\n");
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
				if(flag_RTdebug_printf==1)
				{
					printf("\n*** Trace() exit: ray discarded ***\n");
					printf("****************END****************\n");
				}
			}
			else
			{
				//if(flag_ray_next!=1)//should only change to transmitted ray
				//{
				//	ray_next->trace_obj=intersect->inter_obj->belongobj->subobj;
				//	ray_next->refr_index=intersect->inter_obj->refr_index_o;
				//}
				if(flag_RTdebug_printf==1)
				{
					printf("\n-------------------------------------------\n");
					printf("--- Trace(ray_next,ray_next->trace_obj) ---\n");
					printf("-------------------------------------------\n");
				}
				trace(ray_next,ray_next->trace_obj);
			}
		}	
	}
}
