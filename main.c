#include <stdio.h>
#include <math.h>

/*
# Note
3-body simulation with (i) Euler and (ii) predictor-corrector algo
output: pos.xyz & data.csv

# Definitions
G    grav const (scaled)
n    max time steps
N    num of planets
T    total time (scaled) [physical]
dt   time step [physical]
t    iter step [code]
m    masses
x,vx,ax... are kinematic params for each body over time
ke,pe... are energies
note: require n > T/dt

# Scale
dist,mass,time are conveniently scaled
G        4*pi^2
dist     1 AU, R
mass     sun mass, M
time     earth period, 2*pi*R^1.5/sqrt(GM) = 1

# Setting
system:  sun-earth-mars
init:    aligned along x-axis, with suitable velocities
time:    in units of earth period
note:    energy, momentum should conserve

# Lookup
search " [##] " to look for changeable params
*/

#define G   4*M_PI*M_PI
#define n   10000
#define N   3
#define T   4 // [##]
#define dt  1e-3

int    t;
double m[N],
x[n][N],y[n][N],vx[n][N],vy[n][N],ax[n][N],ay[n][N],
ke[n][N],pe[n][N],ke_tot[n],pe_tot[n],KE,PE;

double d(int t, int i, int j){
	double dx=x[t][j]-x[t][i];
	double dy=y[t][j]-y[t][i];
	return sqrt(dx*dx+dy*dy);
}

double d3(int t, int i, int j){
	double dx=x[t][j]-x[t][i];
	double dy=y[t][j]-y[t][i];
	return pow(dx*dx+dy*dy,1.5);
}

void energy(){
	/* calc energies after all iter */
	int time=t;
	for(int t=0; t<time; t++){
		ke_tot[t]=0;
		pe_tot[t]=0;
		for(int i=0; i<N; i++){
			ke[t][i]=m[i]*(vx[t][i]*vx[t][i]+vy[t][i]*vy[t][i])/2;
			ke_tot[t]+=ke[t][i];
			pe[t][i]=0;
			for(int j=0; j<N; j++){
				if(i==j) continue;
				pe[t][i]-=m[i]*m[j]/d(t,i,j);
			}
			pe_tot[t]+=pe[t][i];
		}
		pe_tot[t]*=0.5;
	}
	KE=ke_tot[0]; // init KE; supposedly KE(t)=KE(0)
	PE=pe_tot[0]; // init PE; supposedly PE(t)=PE(0)
}

void init(){
	/* initialise */
	t=0;

	// system: sun-earth-mars [##]
	m[0]=1; // sun
	x[0][0]=0;
	y[0][0]=0;
	vx[0][0]=0;
	vy[0][0]=0;

	m[1]=3.0e-6; // earth
	x[0][1]=1;
	y[0][1]=0;
	vx[0][1]=0;
	vy[0][1]=2*M_PI;

	m[2]=3.3e-7; // mars
	x[0][2]=1.5;
	y[0][2]=0;
	vx[0][2]=0;
	vy[0][2]=2*M_PI/sqrt(1.5);
}

void update(){
	/* (i) Euler, (ii) predictor-corrector [##] */
	// predictor (Euler)
	for(int i=0; i<N; i++){
		// calc accel
		ax[t][i]=0;
		ay[t][i]=0;
		for(int j=0; j<N; j++){
			if(i==j) continue;
			ax[t][i]+=G*m[j]*(x[t][j]-x[t][i])/d3(t,i,j);
			ay[t][i]+=G*m[j]*(y[t][j]-y[t][i])/d3(t,i,j);
		}
		// update
		x[t+1][i]=x[t][i]+dt*vx[t][i];
		y[t+1][i]=y[t][i]+dt*vy[t][i];
		vx[t+1][i]=vx[t][i]+dt*ax[t][i];
		vy[t+1][i]=vy[t][i]+dt*ay[t][i];
	}

	// corrector
	// for simple Euler, comment whole section out
	double newx[N],newy[N];
	for(int i=0; i<N; i++){
		// calc accel
		ax[t+1][i]=0;
		ay[t+1][i]=0;
		for(int j=0; j<N; j++){
			if(i==j) continue;
			ax[t+1][i]+=G*m[j]*(x[t+1][j]-x[t+1][i])/d3(t+1,i,j);
			ay[t+1][i]+=G*m[j]*(y[t+1][j]-y[t+1][i])/d3(t+1,i,j);
		}
		// update
		newx[i]=x[t][i]+dt*(vx[t][i]+vx[t+1][i])/2;
		newy[i]=y[t][i]+dt*(vy[t][i]+vy[t+1][i])/2;
		vx[t+1][i]=vx[t][i]+dt*(ax[t][i]+ax[t+1][i])/2;
		vy[t+1][i]=vy[t][i]+dt*(ay[t][i]+ay[t+1][i])/2;
	}
	for(int i=0; i<N; i++){
		x[t+1][i]=newx[i];
		y[t+1][i]=newy[i];
	}

	t++;
}

void iter(){
	/* iterate */
	while(t*dt<T) update();
	energy();
}

void data(int num){
	/* print data */
	// num: num of snapshots to print
	int time=t,sep=time/num;
	FILE *f1=fopen("pos.xyz","w"),*f2=fopen("data.csv","w");
	fprintf(f2,"t,ke,pe,ke_err,pe_err\n");
	for(int t=0; t<time; t+=sep){
		fprintf(f1,"%d\n%d\n",N,t/sep);
		for(int i=0; i<N; i++)
			fprintf(f1,"H %lf %lf 0\n",x[t][i],y[t][i]);
		fprintf(f2,"%d,%.8e,%.8e,%.8e,%.8e\n",
			t/sep,ke_tot[t],pe_tot[t],fabs(ke_tot[t]/KE-1),fabs(pe_tot[t]/PE-1));
	}
	fclose(f1);
	fclose(f2);
}

int main(){
	init();
	iter();
	data(400);
	return 0;
}
