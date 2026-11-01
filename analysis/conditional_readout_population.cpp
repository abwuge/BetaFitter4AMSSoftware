// Conditional timing density with readout-conditioned latent population.
// The event output still uses only the timing precision block, never a mass posterior.
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <omp.h>
namespace {
constexpr int K=13;
constexpr double NU=5.;
struct Node {double h,n,g[K],s[5];};
inline double log_sigmoid(double x){return x>=0 ? -std::log1p(std::exp(-x)) : x-std::log1p(std::exp(x));}
}
extern "C" double conditional_population_loss(
 const double *theta,int ne,int nq,int nf,const double *features,const double *y,const double *x,
 const double *moment0,const double *moment1,const double *r,const double *sigma,
 const double *q,const double *logw,const double *z,
 const double *u,const double *a,const double *b,
 double *gradient,double *event_nll,int threads,double rigidity_nu) {
    const int P=K+10*nf;const double *global_theta=theta;
    std::fill(gradient,gradient+P,0.);double total=0.;
    const double rscale=std::exp(theta[12]);
    #pragma omp parallel num_threads(threads)
    {
        std::vector<Node> nodes(2*nq);double local=0.;std::vector<double> grad(P,0.);
        #pragma omp for schedule(static)
        for(int e=0;e<ne;++e) {
            double local_theta[K];std::copy(global_theta,global_theta+K,local_theta);
            for(int h=0;h<10;++h)for(int j=0;j<nf;++j)local_theta[2+h]+=global_theta[K+h*nf+j]*features[e*nf+j];
            const double *theta=local_theta;
            const double d0=10*std::exp(theta[2]),d1=10*std::exp(theta[3]),d2=10*std::exp(theta[4]);
            const double l10=10*theta[5],l20=10*theta[6],l21=10*theta[7];
            const double logdet=3*std::log(10.)+theta[2]+theta[3]+theta[4];
            const double normal=std::lgamma((NU+3)/2)-std::lgamma(NU/2)-1.5*std::log(NU*std::acos(-1.))+logdet;
            double maxh=-std::numeric_limits<double>::infinity(),maxn=maxh;
            for(int species=0;species<2;++species)for(int k=0;k<nq;++k) {
                const int j=species*nq+k;Node &v=nodes[j];
                const double prior=theta[8]*z[k]+theta[9]*z[k]*z[k];
                const double hmix=theta[10]+theta[11]*z[k];
                const double pi=hmix>=0?1/(1+std::exp(-hmix)):std::exp(hmix)/(1+std::exp(hmix));
                const double rr=(r[e]-q[k])/(sigma[e]*rscale),rr2=rr*rr;
                double rlog,rscore;
                if(rigidity_nu>0){rlog=-.5*(rigidity_nu+1)*std::log1p(rr2/rigidity_nu);rscore=-1+(rigidity_nu+1)*rr2/(rigidity_nu+rr2);}
                else {rlog=-.5*rr2;rscore=-1+rr2;}
                v.h=logw[k]+prior+rlog+(species?log_sigmoid(hmix):log_sigmoid(-hmix));
                v.s[0]=z[k];v.s[1]=z[k]*z[k];v.s[2]=species-pi;v.s[3]=(species-pi)*z[k];v.s[4]=rscore;
                double delta[3],basis[3][2];
                for(int l=0;l<3;++l){
                    for(int g=0;g<2;++g)basis[l][g]=moment0[e*6+l*2+g]*a[j]+moment1[e*6+l*2+g]*b[j];
                    delta[l]=y[e*3+l]-x[e*3+l]*u[j]-theta[0]*basis[l][0]-theta[1]*basis[l][1];
                }
                const double v0=d0*delta[0]+l10*delta[1]+l20*delta[2];
                const double v1=d1*delta[1]+l21*delta[2],v2=d2*delta[2];
                const double radius=v0*v0+v1*v1+v2*v2;
                const double f=(NU+3)/(NU+radius);
                const double pe[3]={d0*v0,l10*v0+d1*v1,l20*v0+l21*v1+d2*v2};
                for(int g=0;g<2;++g)v.g[g]=f*(basis[0][g]*pe[0]+basis[1][g]*pe[1]+basis[2][g]*pe[2]);
                v.g[2]=1-f*delta[0]*v0*d0;v.g[3]=1-f*delta[1]*v1*d1;v.g[4]=1-f*delta[2]*v2*d2;
                v.g[5]=-10*f*delta[1]*v0;v.g[6]=-10*f*delta[2]*v0;v.g[7]=-10*f*delta[2]*v1;
                for(int h=0;h<5;++h)v.g[8+h]=v.s[h];
                v.n=v.h+normal-.5*(NU+3)*std::log1p(radius/NU);
                maxh=std::max(maxh,v.h);maxn=std::max(maxn,v.n);
            }
            double hn=0.,nn=0.,sg[5]={0},ng[K]={0};
            for(const Node &v:nodes){
                const double hw=std::exp(v.h-maxh),nw=std::exp(v.n-maxn);hn+=hw;nn+=nw;
                for(int h=0;h<5;++h)sg[h]+=hw*v.s[h];
                for(int h=0;h<K;++h)ng[h]+=nw*v.g[h];
            }
            const double loss=-(maxn+std::log(nn)-maxh-std::log(hn));
            local+=loss;if(event_nll)event_nll[e]=loss;
            for(int h=0;h<K;++h)grad[h]-=ng[h]/nn;
            for(int h=0;h<5;++h)grad[8+h]+=sg[h]/hn;
            for(int h=0;h<10;++h){
                const double score=-ng[2+h]/nn+(h>=6 ? sg[h-6]/hn : 0.);
                for(int j=0;j<nf;++j)grad[K+h*nf+j]+=score*features[e*nf+j];
            }
        }
        #pragma omp critical
        {total+=local;for(int h=0;h<P;++h)gradient[h]+=grad[h];}
    }
    for(int h=0;h<P;++h)gradient[h]/=ne;
    return total/ne;
}
