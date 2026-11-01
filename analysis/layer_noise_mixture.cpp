// Normalized mixture of a core seven-observable Student response and four
// single-layer timing-contamination responses, with physical latent states
// integrated using the already frozen conditional rigidity weights.
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <omp.h>

extern "C" double layer_noise_mixture_loss(const double *theta,int ne,int ns,const double *cache,const double *norm,
 const double *logprior,const double *h2,double logconstant,double *gradient,double *per_event,int threads){
    const int P=6;double logit[5]={0,theta[0],theta[1],theta[2],theta[3]},pi[5],maxlog=*std::max_element(logit,logit+5),sum=0;
    for(int k=0;k<5;++k){pi[k]=std::exp(logit[k]-maxlog);sum+=pi[k];}
    for(int k=0;k<5;++k){pi[k]/=sum;logit[k]=std::log(pi[k]);}
    const double tau2=std::exp(2*theta[4]),alpha2=std::exp(2*theta[5]);
    double b[5]={0},ratio[5]={0},det[5]={0};
    for(int k=1;k<5;++k){b[k]=tau2*alpha2*h2[k-1];ratio[k]=b[k]/(1+b[k]);det[k]=-.5*std::log1p(b[k]);}
    std::fill(gradient,gradient+P,0.);double total=0.;
    #pragma omp parallel num_threads(threads)
    {
        double local=0.;double grad[P]={0};std::vector<double> logs(ns*5),dt(ns*5),da(ns*5);
        #pragma omp for schedule(static)
        for(int e=0;e<ne;++e){
            double maximum=-std::numeric_limits<double>::infinity();
            for(int j=0;j<ns;++j){
                const double *c=cache+(e*ns+j)*6;
                const double first=alpha2*c[0],radius=first+c[1];
                for(int k=0;k<5;++k){
                    const int at=j*5+k;
                    const double removed=k?tau2*alpha2*alpha2*c[1+k]*c[1+k]/(1+b[k]):0;
                    const double rho=std::max(0.,radius-removed);
                    logs[at]=logprior[e*ns+j]-norm[e]+logit[k]+logconstant+3*theta[5]+det[k]-6*std::log1p(rho/5);
                    const double score=6/(5+rho);
                    dt[at]=-ratio[k]+score*2*removed/(1+b[k]);
                    da[at]=3-ratio[k]-score*(2*first-removed*(4-2*ratio[k]));
                    maximum=std::max(maximum,logs[at]);
                }
            }
            double denominator=0.,counts[5]={0},st=0.,sa=0.;
            for(int j=0;j<ns;++j)for(int k=0;k<5;++k){
                const int at=j*5+k;const double value=std::exp(logs[at]-maximum);
                denominator+=value;counts[k]+=value;st+=value*dt[at];sa+=value*da[at];
            }
            const double loss=-maximum-std::log(denominator);local+=loss;if(per_event)per_event[e]=loss;
            for(int k=1;k<5;++k)grad[k-1]+=pi[k]-counts[k]/denominator;
            grad[4]-=st/denominator;grad[5]-=sa/denominator;
        }
        #pragma omp critical
        {total+=local;for(int k=0;k<P;++k)gradient[k]+=grad[k];}
    }
    for(int k=0;k<P;++k)gradient[k]/=ne;
    return total/ne;
}
