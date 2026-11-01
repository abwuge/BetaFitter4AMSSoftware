#!/usr/bin/env python3
"""Build an event-dependent precision model without replacing the base library."""
import hashlib
import subprocess
from pathlib import Path
from research_flight_conditional import write,digest
from research_conditional_mc import ROOT,stamp


def replace_once(text,old,new):
    if text.count(old)!=1:raise ValueError('Unexpected base source structure: '+old[:80])
    return text.replace(old,new,1)


def main():
    source=ROOT/'analysis/flight_conditional_likelihood.cpp';out=ROOT/'results/flight_calibration_research/feature_noise'
    out.mkdir(exist_ok=True);target=ROOT/'analysis/feature_conditional_likelihood.cpp'
    if target.exists():raise FileExistsError(str(target))
    text=source.read_text();source_hash=digest(source)
    if source_hash!='81041530290fe376b4193222a0e4497ac3bbc3992c9c3c3a9087cf64cd205ed8':raise ValueError('Base likelihood source changed')
    text=replace_once(text,'conditional_loss(','conditional_feature_loss(')
    text=replace_once(text,'const double *theta,int ne,int nq,const double *y,const double *x,','const double *theta,int ne,int nq,int nf,const double *features,const double *y,const double *x,')
    text=replace_once(text,'    std::fill(gradient,gradient+K,0.);','    const int P=K+6*nf;const double *global_theta=theta;\n    std::fill(gradient,gradient+P,0.);')
    start=text.index('    const double d0=');end=text.index('    const double rscale=',start)
    local_precision=text[start:end];text=text[:start]+text[end:]
    text=replace_once(text,'std::vector<Node> nodes(2*nq);double local=0.,grad[K]={0};','std::vector<Node> nodes(2*nq);double local=0.;std::vector<double> grad(P,0.);')
    insert='''        for(int e=0;e<ne;++e) {
            double local_theta[K];std::copy(global_theta,global_theta+K,local_theta);
            for(int h=0;h<6;++h)for(int j=0;j<nf;++j)local_theta[2+h]+=global_theta[K+h*nf+j]*features[e*nf+j];
            const double *theta=local_theta;
'''+local_precision
    text=replace_once(text,'        for(int e=0;e<ne;++e) {\n',insert)
    text=replace_once(text,'            for(int h=0;h<5;++h)grad[8+h]+=sg[h]/hn;',
        '            for(int h=0;h<5;++h)grad[8+h]+=sg[h]/hn;\n            for(int h=0;h<6;++h)for(int j=0;j<nf;++j)grad[K+h*nf+j]-=ng[2+h]/nn*features[e*nf+j];')
    text=replace_once(text,'{total+=local;for(int h=0;h<K;++h)gradient[h]+=grad[h];}','{total+=local;for(int h=0;h<P;++h)gradient[h]+=grad[h];}')
    text=replace_once(text,'for(int h=0;h<K;++h)gradient[h]/=ne;','for(int h=0;h<P;++h)gradient[h]/=ne;')
    target.write_text(text)
    result=subprocess.run(['g++','-O3','-std=c++11','-fPIC','-shared','-fopenmp',str(target),'-o',str(out/'libfeature.so')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (out/'build.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError(result.stdout)
    write(out/'build.json',dict(complete=True,created_utc=stamp(),base_source_sha256=source_hash,feature_source_sha256=digest(target),library_sha256=digest(out/'libfeature.so'),base_library_replaced=False))
    print('FEATURE_MODEL_BUILT',flush=True)

if __name__=='__main__':main()
