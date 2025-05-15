/*
 Multi-component Lattice Boltzmann simulation in which passive scalar chemical species and thermal energy
 are all advected by the fluid, which undergoes buoyancy-driven natural convection
 Created by Stuart Bartlett on 31/07/13.
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
#define rm_36 0.027777777777778

// Constants for random number generator (sourse: Numerical Recipes)
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

// All distribution functions move along their respective velocity vectors (streaming step)
void bb_per_bc(double *f,double *f_0,int *obst,int *fro,int *to,int t_cs,int div,int y_s,int l_1,int l_2,int bb_l){
    int i,j,k;
    
    for(k=0;k<18;k+=9){
        // Horizontally pointing distributions simply move as a block
        memcpy(&f[(k+1)*t_cs+y_s],&f_0[(k+1)*t_cs],l_2);	// Move f_1s to the right
        memcpy(&f[(k+3)*t_cs],&f_0[(k+3)*t_cs+y_s],l_2);	// Move f_3s to the left
        
        // Now correct the positions of the exchanged distributions (those with diagonal unit vectors)
        memmove(&f[(k+5)*t_cs],&f[(k+5)*t_cs+1],l_1);				// Move leftmost column of f_5s up
        memmove(&f[(k+7)*t_cs-y_s],&f[(k+7)*t_cs-y_s+1],l_1);		// Move rightmost column of f_6s up
        memmove(&f[(k+8)*t_cs-y_s+1],&f[(k+8)*t_cs-y_s],l_1);		// Move rightmost column of f_7s down
        memmove(&f[(k+8)*t_cs+1],&f[(k+8)*t_cs],l_1);				// Move leftmost column of f_8s down
        
        // Other distributions must be moved one column at a time
        for(j=0;j<div;j++){
            memcpy(&f[(k+2)*t_cs+j*y_s],&f_0[(k+2)*t_cs+j*y_s+1],l_1);			// Move this column of f_2s up
            memcpy(&f[(k+4)*t_cs+j*y_s+1],&f_0[(k+4)*t_cs+j*y_s],l_1);			// Move this column of f_4s down
            if(j<div-1){
                memcpy(&f[(k+5)*t_cs+(j+1)*y_s],&f_0[(k+5)*t_cs+j*y_s+1],l_1);	// Move this column of f_5s up and to the right
                memcpy(&f[(k+6)*t_cs+j*y_s],&f_0[(k+6)*t_cs+(j+1)*y_s+1],l_1);	// Move this column of f_6s up and to the left
                memcpy(&f[(k+7)*t_cs+j*y_s+1],&f_0[(k+7)*t_cs+(j+1)*y_s],l_1);	// Move this column of f_7s down and to the left
                memcpy(&f[(k+8)*t_cs+(j+1)*y_s+1],&f_0[(k+8)*t_cs+j*y_s],l_1);	// Move this column of f_8s down and to the right
            }
            /*if(k>8){
            f[(k+2)*t_cs+(j+1)*y_s-1]=f_0[(k+4)*t_cs+(j+1)*y_s-1];	// f_4s bounce back to become f_2s
            f[(k+4)*t_cs+j*y_s]=f_0[(k+2)*t_cs+j*y_s];				// f_2s bounce back to become f_4s
            f[(k+5)*t_cs+(j+1)*y_s-1]=f_0[(k+7)*t_cs+(j+1)*y_s-1];	// f_7s bounce back to become f_5s
            f[(k+6)*t_cs+(j+1)*y_s-1]=f_0[(k+8)*t_cs+(j+1)*y_s-1];	// f_8s bounce back to become f_6s
            f[(k+7)*t_cs+j*y_s]=f_0[(k+5)*t_cs+j*y_s];				// f_5s bounce back to become f_7s
            f[(k+8)*t_cs+j*y_s]=f_0[(k+6)*t_cs+j*y_s];				// f_6s bounce back to become f_8s
            }*/
        }
        // Bounce-back from obstacles
        for(j=0;j<bb_l;j++){
            f[k*t_cs+to[j]]=f_0[k*t_cs+fro[j]];
        }
    }
}

// Calculate densities and temperatures
void calc_dens(double *f,double *rho,double *eps,int *obst,double T_a,double T_b,double U_0,double W_0,int y_s,int div,int t_cs,int *e_is){
    int i,j,k;
    //double en_def,T_0,lt_u=-1.5*U_0*U_0,lt_w=lt_u*v_fac*v_fac,rho_0;
    double i_opw=1/(1+W_0),i_omw=1/(1-W_0),e_du,l_t_t=-1.5*(U_0*U_0+W_0*W_0),l_t_b=-1.5*W_0*W_0,en_def,rho_p,eps_p;
    
    // Velocity and temperature forcing
    for(j=0;j<div;j++){
        // Top row
        rho_p=(f[j*y_s]+f[t_cs+j*y_s]+f[3*t_cs+j*y_s]+2*(f[2*t_cs+j*y_s]+f[5*t_cs+j*y_s]+f[6*t_cs+j*y_s]))*i_opw;
        // e_4 direction
        e_du=-W_0;
        f[4*t_cs+j*y_s]=i_9*(1+3*e_du+4.5*e_du*e_du+l_t_t);             // Fluid
        f[13*t_cs+j*y_s]=i_9*T_b*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_t);    // Internal energy
        // e_7 direction
        e_du=-U_0-W_0;
        f[7*t_cs+j*y_s]=0.5*(f[t_cs+j*y_s]+f[2*t_cs+j*y_s]-f[3*t_cs+j*y_s]-f[4*t_cs+j*y_s]+2*f[5*t_cs+j*y_s]-rho_p*(U_0+W_0));
        f[16*t_cs+j*y_s]=i_36*T_b*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_t);   // Internal energy
        // e_8 direction
        e_du=U_0-W_0;
        f[8*t_cs+j*y_s]=rho_p*U_0-f[t_cs+j*y_s]+f[3*t_cs+j*y_s]-f[5*t_cs+j*y_s]+f[6*t_cs+j*y_s]+f[7*t_cs+j*y_s];
        f[17*t_cs+j*y_s]=i_36*T_b*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_t);  // Internal energy
        // Adjust boundary temperature
        eps_p=f[9*t_cs+j*y_s]+f[10*t_cs+j*y_s]+f[11*t_cs+j*y_s]+f[12*t_cs+j*y_s]+f[13*t_cs+j*y_s]+f[14*t_cs+j*y_s]+f[15*t_cs+j*y_s]+f[16*t_cs+j*y_s]+f[17*t_cs+j*y_s];
        en_def=rho_p*T_b-eps_p;
        f[13*t_cs+j*y_s]+=t_3*en_def;            // e_4 direction
        f[16*t_cs+j*y_s]+=i_6*en_def;            // e_7 direction
        f[17*t_cs+j*y_s]+=i_6*en_def;            // e_8 direction
        
        // Bottom row
        rho_p=(f[(j+1)*y_s-1]+f[t_cs+(j+1)*y_s-1]+f[3*t_cs+(j+1)*y_s-1]+2*(f[4*t_cs+(j+1)*y_s-1]+f[7*t_cs+(j+1)*y_s-1]+f[8*t_cs+(j+1)*y_s-1]))*i_omw;
        e_du=W_0;
        // e_2 direction
        f[2*t_cs+(j+1)*y_s-1]=i_9*(1+3*e_du+4.5*e_du*e_du+l_t_b);           // Fluid
        f[11*t_cs+(j+1)*y_s-1]=i_9*T_a*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_b);  // Internal energy
        // e_5 direction
        f[5*t_cs+(j+1)*y_s-1]=0.5*(-f[t_cs+(j+1)*y_s-1]-f[2*t_cs+(j+1)*y_s-1]+f[3*t_cs+(j+1)*y_s-1]+f[4*t_cs+(j+1)*y_s-1]+2*f[7*t_cs+(j+1)*y_s-1]+rho_p*W_0);
        f[14*t_cs+(j+1)*y_s-1]=i_36*T_a*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_b);  // Internal energy
        // e_6 direction
        f[6*t_cs+(j+1)*y_s-1]=f[t_cs+(j+1)*y_s-1]-f[3*t_cs+(j+1)*y_s-1]+f[5*t_cs+(j+1)*y_s-1]-f[7*t_cs+(j+1)*y_s-1]+f[8*t_cs+(j+1)*y_s-1];
        f[15*t_cs+(j+1)*y_s-1]=i_36*T_a*(1.5+1.5*e_du+4.5*e_du*e_du+l_t_b);  // Internal energy
        // Adjust boundary temperature
        eps_p=f[9*t_cs+(j+1)*y_s-1]+f[10*t_cs+(j+1)*y_s-1]+f[11*t_cs+(j+1)*y_s-1]+f[12*t_cs+(j+1)*y_s-1]+f[13*t_cs+(j+1)*y_s-1]+f[14*t_cs+(j+1)*y_s-1]+f[15*t_cs+(j+1)*y_s-1]+f[16*t_cs+(j+1)*y_s-1]+f[17*t_cs+(j+1)*y_s-1];
        en_def=rho_p*T_a-eps_p;
        f[11*t_cs+(j+1)*y_s-1]+=t_3*en_def;        // e_2 direction
        f[14*t_cs+(j+1)*y_s-1]+=i_6*en_def;        // e_5 direction
        f[15*t_cs+(j+1)*y_s-1]+=i_6*en_def;        // e_6 direction
    }
    
    // Calculate all densities and internal energies
    for(i=0;i<t_cs;i++){
        rho[i]=obst[y_s+i]*(f[i]+f[t_cs+i]+f[2*t_cs+i]+f[3*t_cs+i]+f[4*t_cs+i]+f[5*t_cs+i]+f[6*t_cs+i]+f[7*t_cs+i]+f[8*t_cs+i]);
        eps[i]=obst[y_s+i]*(f[9*t_cs+i]+f[10*t_cs+i]+f[11*t_cs+i]+f[12*t_cs+i]+f[13*t_cs+i]+f[14*t_cs+i]+f[15*t_cs+i]+f[16*t_cs+i]+f[17*t_cs+i])/rho[i];
    }
}

// Now we calculate new momenta and equilibrium distributions at each lattice point
void eq_dists(double *f,double *f_0,double *rho,double *eps,double *u,int t_cs,int *e_is){
    int i,j,k;
    double l_t,e_du,u_sq;
    
    // Calculate velocity components
    for(i=0;i<t_cs;i++){
        if(rho[i]){
            l_t=1/rho[i];
            u[i]=l_t*(f[t_cs+i]-f[3*t_cs+i]+f[5*t_cs+i]-f[6*t_cs+i]-f[7*t_cs+i]+f[8*t_cs+i]);
            u[t_cs+i]=l_t*(f[2*t_cs+i]-f[4*t_cs+i]+f[5*t_cs+i]+f[6*t_cs+i]-f[7*t_cs+i]-f[8*t_cs+i]);
            // Calculate new equilibrium distributions
            u_sq=u[i]*u[i]+u[t_cs+i]*u[t_cs+i];
            l_t=-1.5*u_sq;
            f_0[i]=f_9*rho[i]*(1+l_t);									// Fluid (rest particles)
            f_0[9*t_cs+i]=-t_3*rho[i]*eps[i]*u_sq;						// Internal energy (rest particles)
            for(k=1;k<5;k++){
                e_du=e_is[2*k]*u[i]+e_is[2*k+1]*u[t_cs+i];
                f_0[k*t_cs+i]=i_9*rho[i]*(1+3*e_du+4.5*e_du*e_du+l_t);						// Fluid (particle index k)
                f_0[(k+9)*t_cs+i]=i_9*rho[i]*eps[i]*(1.5+1.5*e_du+4.5*e_du*e_du+l_t);		// Internal energy (particle index k)
                e_du=e_is[2*(k+4)]*u[i]+e_is[2*(k+4)+1]*u[t_cs+i];
                f_0[(k+4)*t_cs+i]=i_36*rho[i]*(1+3*e_du+4.5*e_du*e_du+l_t);					// Fluid particles (particle index k+4)
                f_0[(k+13)*t_cs+i]=i_36*rho[i]*eps[i]*(3+6*e_du+4.5*e_du*e_du+l_t);			// Internal energy (particle index k+4)
            }
        }
    }
}

// Apply collision step
void coll_step(double *f,double *f_0,double *u,double *eps,int *obst,double T_0,double i_T_0,int *e_is,double b_g0,
               double i_tau_v,double i_tau_c,double c,int t_cs,int y_s,double *omega){
    int i,j,k;
    for(i=0;i<t_cs;i++){
        if(obst[y_s+i]){
            for(k=0;k<9;k++){
                // Fluid
                //f[k*t_cs+i]+=f_0[k*t_cs+i]*(i_tau_v+b_g0*i_T_0*(eps[i]-T_0)*c*(e_is[2*k+1]-u[t_cs+i]))-i_tau_v*f[k*t_cs+i];
                f[k*t_cs+i]+=i_tau_v*(f_0[k*t_cs+i]-f[k*t_cs+i]);             // Convection off
                // Internal energy
                f[(k+9)*t_cs+i]+=i_tau_c*(f_0[(k+9)*t_cs+i]-f[(k+9)*t_cs+i]);	// Collision
            }
        }
    }
}

int main(int argc, char **argv){
    // Declare constants etc.
    FILE *fout_1,*fin;
    int i,j,k,l=0,m,t,bb_l,num_procs,rank,div,t_cs,left,right;
    int t_end=1e6,t_int=t_end/20;
    int x_s,y_s,*obst,*fro,*to;//*emp,emp_l;
    // Read in system size
    fin=fopen("obst_lrg.txt","r");
    fscanf(fin,"%d %d ",&y_s,&x_s);
    int l_1=(y_s-1)*sizeof(double),l_2,l_3,e_is[18]={0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1};
    double i_nds=1.0/(y_s*x_s),Ra=1e6,Pr=7.2,Re=5e2;
    double T_a=1.5,T_b=0.5,T_0=1,i_T_0=1/T_0,c=sqrt(3*T_0);
    double H=y_s*c,tau_v=0.75,tau_c=0.5*(tau_v-0.5)/Pr+0.5,l_t,e_du;
    double nu=i_3*(tau_v-0.5)*c*c,chi=t_3*(tau_c-0.5)*c*c,i_tau_v=1/tau_v,i_tau_c=1/tau_c;
    double b_g0=Ra*chi*nu/(pow(H,3)*(T_a-T_b)),U_0=Re*nu/(1000*c*c),W_0=0.01;
    double *rho,*rho_all,*eps,*u,u_sq,*f,*f_0,*u_x_all,*u_y_all,*e_all,omega[9]={f_9,i_9,i_9,i_9,i_9,i_36,i_36,i_36,i_36};
    char d_out[16],run_st[3];
    
    MPI_Status status; MPI_Comm Comm_cart; MPI_Init(&argc,&argv); MPI_Comm_size(MPI_COMM_WORLD,&num_procs); MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    seed=-(signed)time(NULL)+rank;          // Seed random number generator
    div=x_s/num_procs;						// Calculate labour division
    t_cs=div*y_s;							// Total grid cells in block
    l_2=(div-1)*y_s*sizeof(double);         // Block memory size for horizontally moving velocities (e_1 and e_3)
    l_3=18*t_cs*sizeof(double);             // Block memory size for all distributions (mass and internal energy)
    
    left=rank-1;right=rank+1;
    if(rank==0){left=num_procs-1;}
    else if(rank==num_procs-1){right=0;}
    
    // Assign memory for arrays
    rho=malloc(t_cs*sizeof(double)); rho_all=malloc(y_s*x_s*sizeof(double));
    eps=malloc(t_cs*sizeof(double)); e_all=malloc(y_s*x_s*sizeof(double));
    u=malloc(2*t_cs*sizeof(double));
    u_x_all=malloc(y_s*x_s*sizeof(double)); u_y_all=malloc(y_s*x_s*sizeof(double));
    f=malloc(18*t_cs*sizeof(double)); f_0=malloc(18*t_cs*sizeof(double));
    obst=malloc((t_cs+2*y_s)*sizeof(int)); //emp=malloc(t_cs*sizeof(int));
    fro=malloc(9*t_cs*sizeof(int)); to=malloc(9*t_cs*sizeof(int));
    
    // Read in porous media array
    for(i=0;i<y_s*x_s;i++){
        fscanf(fin,"%d ",&k);
        //k=1;
        if(i>=rank*t_cs && i<(rank+1)*t_cs){
            obst[i%t_cs+y_s]=k;
        }
    }
    fclose(fin);
    
    // Processes now exchange the obstacle columns on their edges
    MPI_Sendrecv(&obst[y_s],y_s,MPI_INT,left,1,&obst[t_cs+y_s],y_s,MPI_INT,right,1,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(&obst[t_cs],y_s,MPI_INT,right,2,&obst[0],y_s,MPI_INT,left,2,MPI_COMM_WORLD,&status);
    
    // Now generate the pointers for the bounce back from obstacles
    for(k=1;k<9;k++){
        for(i=1;i<y_s-1;i++){
            for(j=0;j<div;j++){
                if(obst[(j+1+e_is[2*k])*y_s+i-e_is[2*k+1]]==0 && obst[(j+1)*y_s+i]==1){
                    // To index
                    m=3*(k==1)+4*(k==2)+1*(k==3)+2*(k==4)+7*(k==5)+8*(k==6)+5*(k==7)+6*(k==8);
                    fro[l]=k*t_cs+j*y_s+i;
                    to[l]=m*t_cs+j*y_s+i;
                    l++;
                }
            }
        }
        if(k!=4 && k<7){
            i=y_s-1;
            for(j=0;j<div;j++){
                if(obst[(j+1+e_is[2*k])*y_s+i-e_is[2*k+1]]==0 && obst[(j+1)*y_s+i]==1){
                    // To index
                    m=3*(k==1)+4*(k==2)+1*(k==3)+2*(k==4)+7*(k==5)+8*(k==6)+5*(k==7)+6*(k==8);
                    fro[l]=k*t_cs+j*y_s+i;
                    to[l]=m*t_cs+j*y_s+i;
                    l++;
                }
            }
        }
    }
    bb_l=l;         // Number of obstacle grid cells
    fro=realloc(fro,l*sizeof(int));
    to=realloc(to,l*sizeof(int));
    
    // Set uniform initial density, random temperature of mean T_0 and 0 velocity
    for(i=0;i<t_cs;i++){
        rho[i]=1;
        //eps[i]=obst[y_s+i]*T_0*(1+0.05*2*(ran2(&seed)-0.5));
        eps[i]=obst[y_s+i]*(T_b+(T_a-T_b)*((double)(i%y_s)/(y_s-1))+T_0*0.05*2*(ran2(&seed)-0.5));
        if(i%y_s<0.6*y_s){
            u[i]=obst[y_s+i]*U_0*(1-(double)(i%y_s)/(0.6*y_s-1))*(1+0.05*2*(ran2(&seed)-0.5));
        }
        else{
            u[i]=0;
        }
        // Set initial distributions
        u[t_cs+i]=obst[y_s+i]*W_0;
        u_sq=u[i]*u[i]+u[t_cs+i]*u[t_cs+i];
        l_t=-1.5*u_sq;
        f[i]=f_9*obst[y_s+i]*(1+l_t);									// Fluid (rest particles)
        f[9*t_cs+i]=-t_3*obst[y_s+i]*eps[i]*u_sq;						// Internal energy (rest particles)
        for(k=1;k<5;k++){
            e_du=e_is[2*k]*u[i]+e_is[2*k+1]*u[t_cs+i];
            f[k*t_cs+i]=i_9*obst[y_s+i]*(1+3*e_du+4.5*e_du*e_du+l_t);						// Fluid (particle index k)
            f[(k+9)*t_cs+i]=i_9*obst[y_s+i]*eps[i]*(1.5+1.5*e_du+4.5*e_du*e_du+l_t);		// Internal energy (particle index k)
            e_du=e_is[2*(k+4)]*u[i]+e_is[2*(k+4)+1]*u[t_cs+i];
            f[(k+4)*t_cs+i]=i_36*obst[y_s+i]*(1+3*e_du+4.5*e_du*e_du+l_t);					// Fluid particles (particle index k+4)
            f[(k+13)*t_cs+i]=i_36*obst[y_s+i]*eps[i]*(3+6*e_du+4.5*e_du*e_du+l_t);			// Internal energy (particle index k+4)
        }
    }
    
    // Run simulation
    for(t=0;t<=t_end;t++){
        // Gather and print whole array
        if(t%t_int==0){
            MPI_Gather(rho,t_cs,MPI_DOUBLE,rho_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gather(u,t_cs,MPI_DOUBLE,u_x_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gather(&u[t_cs],t_cs,MPI_DOUBLE,u_y_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gather(eps,t_cs,MPI_DOUBLE,e_all,t_cs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(rank==0){
                sprintf(run_st,"%d",t/t_int);
                strcpy(d_out,"opts_hyv_"); strcat(d_out,run_st); strcat(d_out,".txt");
                fout_1=fopen(d_out,"w");
                fprintf(fout_1,"%d %d %d %d %8.6f %8.6f %8.6f ",y_s,x_s,t_end,t_int,tau_v,tau_c,b_g0);
                for(i=0;i<y_s*x_s;i++){
                    fprintf(fout_1,"%8.6f ",rho_all[i]);
                }
                for(i=0;i<y_s*x_s;i++){
                    fprintf(fout_1,"%8.6f ",e_all[i]);
                }
                for(i=0;i<y_s*x_s;i++){
                    fprintf(fout_1,"%8.6f ",c*u_x_all[i]);
                }
                for(i=0;i<y_s*x_s;i++){
                    fprintf(fout_1,"%8.6f ",c*u_y_all[i]);
                }
                fclose(fout_1);
            }
        }
        
        // Make temporary copies of all fields
        memcpy(&f_0[0],&f[0],l_3);
        
        // Streaming step
        for(k=0;k<18;k+=9){
            // Exchange leftward moving distributions
            MPI_Sendrecv(&f_0[(k+3)*t_cs],y_s,MPI_DOUBLE,left,3+10*k,&f[(k+4)*t_cs-y_s],y_s,MPI_DOUBLE,right,3+10*k,MPI_COMM_WORLD,&status);
            MPI_Sendrecv(&f_0[(k+6)*t_cs],y_s,MPI_DOUBLE,left,6+20*k,&f[(k+7)*t_cs-y_s],y_s,MPI_DOUBLE,right,6+20*k,MPI_COMM_WORLD,&status);
            MPI_Sendrecv(&f_0[(k+7)*t_cs],y_s,MPI_DOUBLE,left,7+30*k,&f[(k+8)*t_cs-y_s],y_s,MPI_DOUBLE,right,7+30*k,MPI_COMM_WORLD,&status);
            // Exchange rightward moving distributions
            MPI_Sendrecv(&f_0[(k+2)*t_cs-y_s],y_s,MPI_DOUBLE,right,1+10*k,&f[(k+1)*t_cs],y_s,MPI_DOUBLE,left,1+10*k,MPI_COMM_WORLD,&status);
            MPI_Sendrecv(&f_0[(k+6)*t_cs-y_s],y_s,MPI_DOUBLE,right,5+20*k,&f[(k+5)*t_cs],y_s,MPI_DOUBLE,left,5+20*k,MPI_COMM_WORLD,&status);
            MPI_Sendrecv(&f_0[(k+9)*t_cs-y_s],y_s,MPI_DOUBLE,right,8+30*k,&f[(k+8)*t_cs],y_s,MPI_DOUBLE,left,8+30*k,MPI_COMM_WORLD,&status);
        }
        
        bb_per_bc(f,f_0,obst,fro,to,t_cs,div,y_s,l_1,l_2,bb_l);                             // Streaming and bounce back on boundaries
        calc_dens(f,rho,eps,obst,T_a,T_b,U_0,W_0,y_s,div,t_cs,e_is);                        // Calculate all new densities
        eq_dists(f,f_0,rho,eps,u,t_cs,e_is);                                                // Calculate new equilibrium distributions
        coll_step(f,f_0,u,eps,obst,T_0,i_T_0,e_is,b_g0,i_tau_v,i_tau_c,c,t_cs,y_s,omega);   // Collision step
    }
    MPI_Finalize();
}
