/* convergence study: ternary vs FP, undriven Gaussian pulse, fixed physical time */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#define NX 80
#define NY 80
#ifndef MC
#define MC 10000
#endif
#ifndef NSTEP
#define NSTEP 2000
#endif
#define M 100000
#define AMP 20.0
static const float CFL = (float)MC/(float)M;
static int8_t Ez[NX][NY],Hx[NX][NY],Hy[NX][NY];
static int32_t aEz[NX][NY],aHx[NX][NY],aHy[NX][NY];
static int32_t e1E[NX][NY],e2E[NX][NY],e1X[NX][NY],e2X[NX][NY],e1Y[NX][NY],e2Y[NX][NY];
static float fEz[NX][NY],fHx[NX][NY],fHy[NX][NY];
static long sat=0;
static inline int8_t enc(int32_t x,int32_t*e1,int32_t*e2){
  int32_t v=x+(3*(*e1)-(*e2))/2,s;
  if(v>=0)s=(v+M/2)/M; else s=-((-v+M/2)/M);
  if(s>1){s=1;sat++;} if(s<-1){s=-1;sat++;}
  *e2=*e1; *e1=v-s*M;
  if(*e1>2*M)*e1=2*M; if(*e1<-2*M)*e1=-2*M;
  return (int8_t)s;}
int main(void){int n,i,j;
 for(i=1;i<NX-1;i++)for(j=1;j<NY-1;j++){
   double g=10.0*exp(-(((i-40.0)*(i-40.0))+((j-40.0)*(j-40.0)))/50.0);
   fEz[i][j]=(float)g; aEz[i][j]=(int32_t)lrint(g*M/AMP);
   Ez[i][j]=enc(aEz[i][j],&e1E[i][j],&e2E[i][j]);}
 double E0=0; for(i=0;i<NX;i++)for(j=0;j<NY;j++)E0+=fEz[i][j]*fEz[i][j];
 for(n=0;n<NSTEP;n++){
  for(i=0;i<NX-1;i++)for(j=0;j<NY-1;j++){
    aHx[i][j]+=MC*(Ez[i][j]-Ez[i][j+1]); aHy[i][j]+=MC*(Ez[i+1][j]-Ez[i][j]);
    Hx[i][j]=enc(aHx[i][j],&e1X[i][j],&e2X[i][j]); Hy[i][j]=enc(aHy[i][j],&e1Y[i][j],&e2Y[i][j]);}
  for(i=1;i<NX-1;i++)for(j=1;j<NY-1;j++){
    aEz[i][j]+=MC*((Hy[i][j]-Hy[i-1][j])-(Hx[i][j]-Hx[i][j-1]));
    Ez[i][j]=enc(aEz[i][j],&e1E[i][j],&e2E[i][j]);}
  for(i=0;i<NX;i++){Ez[i][0]=0;Ez[i][NY-1]=0;}
  for(j=0;j<NY;j++){Ez[0][j]=0;Ez[NX-1][j]=0;}
  for(i=0;i<NX-1;i++)for(j=0;j<NY-1;j++){
    fHx[i][j]+=CFL*(fEz[i][j]-fEz[i][j+1]); fHy[i][j]+=CFL*(fEz[i+1][j]-fEz[i][j]);}
  for(i=1;i<NX-1;i++)for(j=1;j<NY-1;j++)
    fEz[i][j]+=CFL*((fHy[i][j]-fHy[i-1][j])-(fHx[i][j]-fHx[i][j-1]));
  for(i=0;i<NX;i++){fEz[i][0]=0;fEz[i][NY-1]=0;}
  for(j=0;j<NY;j++){fEz[0][j]=0;fEz[NX-1][j]=0;}}
 double num=0,den=0,Et=0,Ef=0;
 for(i=0;i<NX;i++)for(j=0;j<NY;j++){
   double t=aEz[i][j]*AMP/(double)M,f=fEz[i][j];
   num+=(t-f)*(t-f); den+=f*f; Et+=t*t; Ef+=f*f;}
 printf("CFL=%.5f steps=%6d  NRMSE(Ez field)=%8.4f  E_tern/E_fp=%8.4f  E_fp/E_fp(0)=%.5f  sat=%ld\n",
   CFL,NSTEP,sqrt(num/den),Et/Ef,Ef/E0,sat);

 /* Machine-readable row, appended so a sweep builds the table.
    Columns: CFL  steps  NRMSE  E_tern/E_fp  saturation_events        */
 { FILE *t = fopen("convergence.txt","a");
   if(t){ fprintf(t,"%.8f %d %.6f %.6f %ld\n",
                  CFL,NSTEP,sqrt(num/den),Et/Ef,sat); fclose(t); } }

 { char fn[64]; sprintf(fn,"cut_%d.txt",MC); FILE*g=fopen(fn,"w");
   for(j=0;j<NY;j++) fprintf(g,"%d %.6f %.6f\n",j,aEz[40][j]*AMP/(double)M,fEz[40][j]); fclose(g);}
 return 0;}
