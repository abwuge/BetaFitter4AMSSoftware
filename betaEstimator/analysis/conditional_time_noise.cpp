// Conditional time-noise scales; fixed geometry/era features only.
// Scaling the first three precision-Cholesky columns preserves the marginal ADC
// covariance and the time-on-ADC regression; it changes conditional time noise.
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <omp.h>
namespace {
constexpr int D=7,K=39;constexpr double NU=5.,UREF=1.15;
struct Node{double h,n,g[K],s[5],mg[4],cg[3];};
inline double logsig(double x){return x>=0?-std::log1p(std::exp(-x)):x-std::log1p(std::exp(x));}
}
extern "C" double conditional_time_noise_loss(const double *theta,int ne,int nq,int nf,int nh,const double *features,const double *noise_features,
 const double *y,const double *x,const double *moment0,const double *moment1,const double *r,const double *sigma,
 const double *q,const double *logw,const double *z,const double *u,const double *a,const double *b,
 double *gradient,double *event_nll,int threads,double rigidity_nu){
    const int BASE=K+4*nf,P=BASE+3*nh;std::fill(gradient,gradient+P,0.);double total=0.,L[D][D]={{0}},logdet=0.;
    int ci[28],cj[28],kchol=0;
    for(int i=0;i<D;++i)for(int j=0;j<=i;++j){ci[kchol]=i;cj[kchol]=j;L[i][j]=10*(i==j?std::exp(theta[2+kchol]):theta[2+kchol]);if(i==j)logdet+=std::log(L[i][j]);++kchol;}
    const double constant=std::lgamma((NU+D)/2)-std::lgamma(NU/2)-.5*D*std::log(NU*std::acos(-1.))+logdet;
    const double rscale=std::exp(theta[34]);double slope[4];for(int l=0;l<4;++l)slope[l]=std::exp(theta[35+l]);
    #pragma omp parallel num_threads(threads)
    {
        std::vector<Node> nodes(2*nq);std::vector<double> grad(P,0.);double local=0.;
        #pragma omp for schedule(static)
        for(int e=0;e<ne;++e){
            double mean[4]={0};for(int l=0;l<4;++l)for(int f=0;f<nf;++f)mean[l]+=theta[K+l*nf+f]*features[(e*4+l)*nf+f];
            double scale[D]={1,1,1,1,1,1,1},event_logdet=0;
            for(int h=0;h<3;++h){double eta=0;for(int f=0;f<nh;++f)eta+=theta[BASE+h*nh+f]*noise_features[e*nh+f];scale[h]=std::exp(eta);event_logdet+=eta;}
            double maxh=-std::numeric_limits<double>::infinity(),maxn=maxh;
            for(int species=0;species<2;++species)for(int k=0;k<nq;++k){
                const int j=species*nq+k;Node &node=nodes[j];
                const double hmix=theta[32]+theta[33]*z[k];const double pi=hmix>=0?1/(1+std::exp(-hmix)):std::exp(hmix)/(1+std::exp(hmix));
                const double rr=(r[e]-q[k])/(sigma[e]*rscale),rr2=rr*rr;
                double rlog,rscore;
                if(rigidity_nu>0){rlog=-.5*(rigidity_nu+1)*std::log1p(rr2/rigidity_nu);rscore=-1+(rigidity_nu+1)*rr2/(rigidity_nu+rr2);}
                else{rlog=-.5*rr2;rscore=-1+rr2;}
                node.h=logw[k]+theta[30]*z[k]+theta[31]*z[k]*z[k]+rlog+(species?logsig(hmix):logsig(-hmix));
                node.s[0]=z[k];node.s[1]=z[k]*z[k];node.s[2]=species-pi;node.s[3]=(species-pi)*z[k];node.s[4]=rscore;
                double residual[D],basis[3][2];
                for(int l=0;l<3;++l){
                    for(int g=0;g<2;++g)basis[l][g]=moment0[e*6+l*2+g]*a[j]+moment1[e*6+l*2+g]*b[j];
                    residual[l]=y[e*D+l]-x[e*3+l]*u[j]-theta[0]*basis[l][0]-theta[1]*basis[l][1];
                }
                for(int l=0;l<4;++l)residual[3+l]=y[e*D+3+l]-mean[l]-slope[l]*(u[j]-UREF);
                double white[D]={0},pe[D]={0},radius=0.;
                for(int h=0;h<D;++h){for(int l=h;l<D;++l)white[h]+=L[l][h]*residual[l];white[h]*=scale[h];radius+=white[h]*white[h];}
                for(int h=0;h<D;++h)for(int l=0;l<=h;++l)pe[h]+=L[h][l]*scale[l]*white[l];
                const double f=(NU+D)/(NU+radius);
                for(int g=0;g<2;++g){node.g[g]=0.;for(int l=0;l<3;++l)node.g[g]+=f*basis[l][g]*pe[l];}
                for(int h=0;h<28;++h){const int i=ci[h],jj=cj[h];node.g[2+h]=(i==jj?1.:0.)-f*residual[i]*white[jj]*scale[jj]*(i==jj?L[i][jj]:10.);}
                for(int h=0;h<5;++h)node.g[30+h]=node.s[h];
                for(int l=0;l<4;++l){node.mg[l]=f*pe[3+l];node.g[35+l]=node.mg[l]*slope[l]*(u[j]-UREF);}
                for(int h=0;h<3;++h)node.cg[h]=1-f*white[h]*white[h];
                node.n=node.h+constant+event_logdet-.5*(NU+D)*std::log1p(radius/NU);
                maxh=std::max(maxh,node.h);maxn=std::max(maxn,node.n);
            }
            double hn=0.,nn=0.,sg[5]={0},ng[K]={0},mg[4]={0},cg[3]={0};
            for(const Node &v:nodes){const double hw=std::exp(v.h-maxh),nw=std::exp(v.n-maxn);hn+=hw;nn+=nw;
                for(int h=0;h<5;++h)sg[h]+=hw*v.s[h];for(int h=0;h<K;++h)ng[h]+=nw*v.g[h];for(int l=0;l<4;++l)mg[l]+=nw*v.mg[l];for(int h=0;h<3;++h)cg[h]+=nw*v.cg[h];}
            const double loss=-(maxn+std::log(nn)-maxh-std::log(hn));local+=loss;if(event_nll)event_nll[e]=loss;
            for(int h=0;h<K;++h)grad[h]-=ng[h]/nn;for(int h=0;h<5;++h)grad[30+h]+=sg[h]/hn;
            for(int l=0;l<4;++l)for(int f=0;f<nf;++f)grad[K+l*nf+f]-=mg[l]/nn*features[(e*4+l)*nf+f];
            for(int h=0;h<3;++h)for(int f=0;f<nh;++f)grad[BASE+h*nh+f]-=cg[h]/nn*noise_features[e*nh+f];
        }
        #pragma omp critical
        {total+=local;for(int h=0;h<P;++h)gradient[h]+=grad[h];}
    }
    for(int h=0;h<P;++h)gradient[h]/=ne;
    return total/ne;
}
