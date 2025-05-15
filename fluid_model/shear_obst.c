/*
 *	Single-phase Lattice Boltzmann code
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <float.h>

#define i_3  0.33333333333333
#define t_3  0.66666666666667
#define i_6  0.16666666666667
#define i_9  0.11111111111111
#define f_9  0.44444444444444
#define i_12 0.083333333333333
#define i_36 0.027777777777778
#define t_i7 1.7142857142857
#define pi 3.14159265359

// Constants for random number generator (source: Numerical Recipes)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long seed;

// Random number generator

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	
	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

// Velocity distribution functions move along their respective velocity vectors (streaming step)
void bb_per_bc_ob(double *f,double *f_0,int *obst,int t_cs,int div,int y_s,int l_1,int l_2,int *fro,int *to,int bb_l){
	int j;
	
	// Horizontally pointing distributions simply move as a block
	memcpy(&f[t_cs+y_s],&f_0[t_cs],l_2);				// Move f_1s to the right
	memcpy(&f[3*t_cs],&f_0[3*t_cs+y_s],l_2);			// Move f_3s to the left
	
	// Now correct the positions of the exchanged distributions (those with diagonal unit vectors)
	memmove(&f[5*t_cs],&f[5*t_cs+1],l_1);				// Move leftmost column of f_5s up
	memmove(&f[7*t_cs-y_s],&f[7*t_cs-y_s+1],l_1);		// Move rightmost column of f_6s up
	memmove(&f[8*t_cs-y_s+1],&f[8*t_cs-y_s],l_1);		// Move rightmost column of f_7s down
	memmove(&f[8*t_cs+1],&f[8*t_cs],l_1);				// Move leftmost column of f_8s down
	
	// Other distributions must be moved one column at a time
	for(j=0;j<div;j++){
		memcpy(&f[2*t_cs+j*y_s],&f_0[2*t_cs+j*y_s+1],l_1);			// Move this column of f_2s up
		memcpy(&f[4*t_cs+j*y_s+1],&f_0[4*t_cs+j*y_s],l_1);			// Move this column of f_4s down
		if(j<div-1){
			memcpy(&f[5*t_cs+(j+1)*y_s],&f_0[5*t_cs+j*y_s+1],l_1);	// Move this column of f_5s up and to the right
			memcpy(&f[6*t_cs+j*y_s],&f_0[6*t_cs+(j+1)*y_s+1],l_1);	// Move this column of f_6s up and to the left
			memcpy(&f[7*t_cs+j*y_s+1],&f_0[7*t_cs+(j+1)*y_s],l_1);	// Move this column of f_7s down and to the left
			memcpy(&f[8*t_cs+(j+1)*y_s+1],&f_0[8*t_cs+j*y_s],l_1);	// Move this column of f_8s down and to the right
		}
		f[2*t_cs+(j+1)*y_s-1]=f_0[4*t_cs+(j+1)*y_s-1];		// f_4s bounce back to become f_2s
		f[4*t_cs+j*y_s]=f_0[2*t_cs+j*y_s];					// f_2s bounce back to become f_4s
		f[5*t_cs+(j+1)*y_s-1]=f_0[7*t_cs+(j+1)*y_s-1];		// f_7s bounce back to become f_5s
		f[6*t_cs+(j+1)*y_s-1]=f_0[8*t_cs+(j+1)*y_s-1];		// f_8s bounce back to become f_6s
		f[7*t_cs+j*y_s]=f_0[5*t_cs+j*y_s];					// f_5s bounce back to become f_7s
		f[8*t_cs+j*y_s]=f_0[6*t_cs+j*y_s];					// f_6s bounce back to become f_8s
	}
	
	/*// Bounce-back from obstacles
	for(j=0;j<bb_l;j++){
		f[to[j]]=f_0[fro[j]];
	}
	for(j=0;j<t_cs;j++){
		f[j]*=obst[y_s+j];
	}*/
}

// Calculate densities
void calc_rho(double *f,double *rho,int y_s,int div,int x_s,int t_cs){
	int i;
	for(i=0;i<t_cs;i++){
		rho[i]=f[i]+f[t_cs+i]+f[2*t_cs+i]+f[3*t_cs+i]+f[4*t_cs+i]+f[5*t_cs+i]+f[6*t_cs+i]+f[7*t_cs+i]+f[8*t_cs+i];
	}
}

void vel_force(double *f,double *rho,int div,int y_s,int t_cs,double lt,double U,int *e_is){
	int j,k,m=100;
	double e_du;
	
	for(j=0;j<div;j++){
		//f[j*y_s+m]=f_9*rho[j*y_s+m]*(1+lt);
		f[j*y_s+m]=f_9*(1+lt);
		
		for(k=1;k<5;k++){
			e_du=e_is[2*k]*U;
			//f[k*t_cs+j*y_s+m]=i_9*rho[j*y_s+m]*(1+3*e_du+4.5*e_du*e_du+lt);
			f[k*t_cs+j*y_s+m]=i_9*(1+3*e_du+4.5*e_du*e_du+lt);
			e_du=e_is[2*(k+4)]*U;
			//f[(k+4)*t_cs+j*y_s+m]=i_36*rho[j*y_s+m]*(1+3*e_du+4.5*e_du*e_du+lt);
			f[(k+4)*t_cs+j*y_s+m]=i_36*(1+3*e_du+4.5*e_du*e_du+lt);
		}
	}
}

// Now we calculate new momenta and equilibrium distributions at each lattice point
void eq_dists(double *f,double *f_0,double *rho,double *u,double *u_sq,int t_cs,int *e_is){
	int i,j,k;
	double l_t,e_du;
	// Calculate velocity components
	for(i=0;i<t_cs;i++){
		l_t=1/rho[i];
		u[i]=l_t*(f[t_cs+i]-f[3*t_cs+i]+f[5*t_cs+i]-f[6*t_cs+i]-f[7*t_cs+i]+f[8*t_cs+i]);
		u[t_cs+i]=l_t*(f[2*t_cs+i]-f[4*t_cs+i]+f[5*t_cs+i]+f[6*t_cs+i]-f[7*t_cs+i]-f[8*t_cs+i]);
		// Calculate new equilibrium distributions
		u_sq[i]=u[i]*u[i]+u[t_cs+i]*u[t_cs+i];
		l_t=-1.5*u_sq[i];
		f_0[i]=f_9*rho[i]*(1+l_t);
		for(k=1;k<5;k++){
			e_du=e_is[2*k]*u[i]+e_is[2*k+1]*u[t_cs+i];
			f_0[k*t_cs+i]=i_9*rho[i]*(1+3*e_du+4.5*e_du*e_du+l_t);
			e_du=e_is[2*(k+4)]*u[i]+e_is[2*(k+4)+1]*u[t_cs+i];
			f_0[(k+4)*t_cs+i]=i_36*rho[i]*(1+3*e_du+4.5*e_du*e_du+l_t);
		}
	}
}

// Apply collision step
void coll_step(double *f,double *f_0,double i_tau_v,int div,int y_s,int t_cs){	
	int i,k;
	
	 for(i=0;i<t_cs;i++){
         for(k=0;k<9;k++){
             f[k*t_cs+i]+=i_tau_v*(f_0[k*t_cs+i]-f[k*t_cs+i]);
         }
	 }
}

int main(int argc, char **argv){
	
	// Declare constants etc.
	int i,j,k,t,l=0,m;
	int num_procs,rank,div,t_cs,left,right,bb_l,*fro,*to;
	int x_s,y_s;
	int t_end=4e3,t_int=50;
	int e_is[18]={0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1};
	int l_1,l_2,*obst,o_t;
	double tau_v=0.6,Re=0.5e3,lt,C_s=pow(3,-0.5);
	double i_tau_v=1/tau_v,nu=C_s*C_s*(tau_v-0.5);
	double U=Re*nu/40.0;
	double l_t=-1.5*U*U;
	double *rho,*u,*u_sq,*u_x_all,*u_y_all,*f,*f_0;
	
	FILE *fout_1,*fin;
	MPI_Status status;
	MPI_Comm Comm_cart;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
	seed=-(signed)time(NULL)*rank*rank;
	
	left=rank-1;
	right=rank+1;
	if(rank==0){
		left=num_procs-1;
	}
	else if(rank==num_procs-1){
		right=0;
	}
	
	if(rank<3 && rank>1){
	U=0;
	}
	
	fin=fopen("obst.txt","r");
	fscanf(fin,"%d ",&y_s);
	fscanf(fin,"%d ",&x_s);
	div=x_s/num_procs;
	t_cs=div*y_s;
	l_1=(y_s-1)*sizeof(double); l_2=(div-1)*y_s*sizeof(double);
	
	// Assign memory for arrays
	rho=malloc(t_cs*sizeof(double));
	u=malloc(2*t_cs*sizeof(double)); u_sq=malloc(t_cs*sizeof(double));
	u_x_all=malloc(y_s*x_s*sizeof(double)); u_y_all=malloc(y_s*x_s*sizeof(double));
	f=malloc(9*t_cs*sizeof(double)); f_0=malloc(9*t_cs*sizeof(double));
	obst=malloc((t_cs+2*y_s)*sizeof(int));
	fro=malloc(9*t_cs*sizeof(int)); to=malloc(9*t_cs*sizeof(int));
	
	/*for(i=0;i<y_s*x_s;i++){
		fscanf(fin,"%d ",&o_t);
		if(i>=rank*t_cs && i<(rank+1)*t_cs){
			obst[i%t_cs+y_s]=o_t;
		}
	}*/
	fclose(fin);
	/*
	// Processes now exchange the obstacle columns on their edges
	MPI_Sendrecv(&obst[y_s],y_s,MPI_INT,left,1,&obst[t_cs+y_s],y_s,MPI_INT,right,1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv(&obst[t_cs],y_s,MPI_INT,right,2,&obst[0],y_s,MPI_INT,left,2,MPI_COMM_WORLD,&status);
	
	// Now generate the pointers for the bounce back from obstacles
	for(k=1;k<9;k++){
		for(i=1;i<y_s-1;i++){
			for(j=0;j<div;j++){
				if(obst[(j+1+e_is[2*k])*y_s+i+e_is[2*k+1]]==0){
					// To index
					m=3*(k==1)+4*(k==2)+1*(k==3)+2*(k==4)+7*(k==5)+8*(k==6)+5*(k==7)+6*(k==8);
					fro[l]=k*t_cs+j*y_s+i;
					to[l]=m*t_cs+j*y_s+i;
					l++;
				}
			}
		}
	}
	bb_l=l;
	fro=realloc(fro,l*sizeof(int)); to=realloc(to,l*sizeof(int));*/
	
	if(rank==0){
		fout_1=fopen("sh_out.txt","w");
		fprintf(fout_1,"%d %d %d %d %10.8f %10.8f ",y_s,x_s,t_end,t_int,tau_v,U);
	}
	
	// Set initial conditions
	for(i=0;i<t_cs;i++){
		u[i]=0.0;
		u[t_cs+i]=0.0;
		u_sq[i]=u[i]*u[i]+u[t_cs+i]*u[t_cs+i];
		//rho[i]=obst[y_s+i]*(1+0.05*2*(ran2(&seed)-0.5));
		rho[i]=1+0.02*2*(ran2(&seed)-0.5);
		f[i]=f_9*rho[i];
		for(k=1;k<5;k++){
			f[k*t_cs+i]=i_9*rho[i];
			f[(k+4)*t_cs+i]=i_36*rho[i];
		}
	}
	
	for(t=0;t<t_end;t++){
		memcpy(&f_0[0],&f[0],9*t_cs*sizeof(double));
		
		// Streaming step
		MPI_Sendrecv(&f_0[3*t_cs],y_s,MPI_DOUBLE,left,3,&f[4*t_cs-y_s],y_s,MPI_DOUBLE,right,3,MPI_COMM_WORLD,&status);
		MPI_Sendrecv(&f_0[6*t_cs],y_s,MPI_DOUBLE,left,6,&f[7*t_cs-y_s],y_s,MPI_DOUBLE,right,6,MPI_COMM_WORLD,&status);
		MPI_Sendrecv(&f_0[7*t_cs],y_s,MPI_DOUBLE,left,7,&f[8*t_cs-y_s],y_s,MPI_DOUBLE,right,7,MPI_COMM_WORLD,&status);
		MPI_Sendrecv(&f_0[2*t_cs-y_s],y_s,MPI_DOUBLE,right,1,&f[t_cs],y_s,MPI_DOUBLE,left,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv(&f_0[6*t_cs-y_s],y_s,MPI_DOUBLE,right,5,&f[5*t_cs],y_s,MPI_DOUBLE,left,5,MPI_COMM_WORLD,&status);
		MPI_Sendrecv(&f_0[9*t_cs-y_s],y_s,MPI_DOUBLE,right,8,&f[8*t_cs],y_s,MPI_DOUBLE,left,8,MPI_COMM_WORLD,&status);
		
		bb_per_bc_ob(f,f_0,obst,t_cs,div,y_s,l_1,l_2,fro,to,bb_l);
		calc_rho(f,rho,y_s,div,x_s,t_cs);
		vel_force(f,rho,div,y_s,t_cs,lt,U,e_is);
		eq_dists(f,f_0,rho,u,u_sq,t_cs,e_is);
		coll_step(f,f_0,i_tau_v,div,y_s,t_cs);
		
		// Gather and print whole array
		if(t%t_int==0){
			MPI_Gather(u,t_cs,MPI_DOUBLE,u_x_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(&u[t_cs],t_cs,MPI_DOUBLE,u_y_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if(rank==0){
			    for(i=0;i<y_s*x_s;i++){
					fprintf(fout_1,"%10.8f ",u_x_all[i]);
				}
				for(i=0;i<y_s*x_s;i++){
					fprintf(fout_1,"%10.8f ",u_y_all[i]);
				}
			}
		}
				
		/*if(t%t_int==0){
			MPI_Gather(u,t_cs,MPI_DOUBLE,u_x_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(&u[t_cs],t_cs,MPI_DOUBLE,u_y_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if(rank==0){
				for(i=0;i<y_s*x_s;i++){
					fprintf(fout_1,"%10.8f %10.8f ",u_x_all[i],u_y_all[i]);
				}
			}
		}*/
	}
	
	fclose(fout_1);
	MPI_Finalize();
}
