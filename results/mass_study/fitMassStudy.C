#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TSystem.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

namespace MassStudy {
constexpr int P=8,B=4,NP=P+2*B;
const double EDGES[B+1]={2.5,3.,4.,5.,6.};
struct E {Long64_t run,event;int role,source,z,bin,richOK,tag;double R,Rref,u,ur,relR,mr;std::array<double,P>x;};
struct Model {std::array<double,P>center,scale;std::array<double,NP>p;int status;double nll,edm;};
double mass(int z,int k){return z==1?(k?1.875613:.938272):(k?3.727379:2.808392);}
double mu(const E&e,int k){double a=mass(e.z,k)/(e.z*e.R);return std::sqrt(1+a*a);}
double fraction(double a){return std::max(0.,std::min(1.,a));}
double quantile(std::vector<double> v,double p){if(v.empty())return NAN;std::sort(v.begin(),v.end());double j=p*(v.size()-1);size_t i=size_t(j);return v[i]+(j-i)*(v[std::min(i+1,v.size()-1)]-v[i]);}
double width(const std::vector<double>&v){return .5*(quantile(v,.84)-quantile(v,.16));}
std::array<double,P> features(const E&e,const Model&m){std::array<double,P>x;x[0]=1;for(int j=1;j<P;j++)x[j]=std::max(-4.,std::min(4.,(e.x[j]-m.center[j])/m.scale[j]));return x;}
double correction(const E&e,const Model&m){auto x=features(e,m);double c=0;for(int j=0;j<P;j++)c+=m.p[j]*x[j];return c;}
// Normalized t6 densities in measured inverse beta. Correction is independent
// of the time slope: du_corrected/du_raw = 1 at fixed geometry/contrasts.
// Tracker broadening uses a variance-matched t approximation; err(R^-1) is
// supplied at the first Kalman measurement, not propagated to z=0.
double logpdf(double u,double mean,double s,double relR){double d=(mean*mean-1)/mean;double v=s*s+(4./6)*d*d*relR*relR;double t=(u-mean)*(u-mean)/v;return std::lgamma(3.5)-std::lgamma(3.)-.5*std::log(6*M_PI*v)-3.5*std::log1p(t/6);}
double nllone(const E&e,double corrected,const double*p){double s=std::exp(p[P+e.bin]),f=fraction(p[P+B+e.bin]);double a=std::log1p(-f)+logpdf(corrected,mu(e,0),s,e.relR),b=std::log(f)+logpdf(corrected,mu(e,1),s,e.relR);double v=std::max(a,b);return -(v+std::log(std::exp(a-v)+std::exp(b-v)));}
double score(const E&e,double u,const Model&m){double s=std::exp(m.p[P+e.bin]);return logpdf(u,mu(e,1),s,e.relR)-logpdf(u,mu(e,0),s,e.relR);}
Model fit(const std::vector<E>&data,int mode,const Model*initial=nullptr,bool inject=false,bool optional=false){
  Model m{};if(data.size()<500)throw std::runtime_error("Too few fit events");
  if(initial)m=*initial;else{m.center[0]=0;m.scale[0]=1;for(int j=1;j<P;j++){std::vector<double>v;for(auto&e:data)v.push_back(e.x[j]);m.center[j]=quantile(v,.5);m.scale[j]=std::max(1e-6,.741301109*(quantile(v,.75)-quantile(v,.25)));}for(int j=0;j<P;j++)m.p[j]=0;for(int b=0;b<B;b++){m.p[P+b]=std::log(data[0].z==1?.04:.025);m.p[P+B+b]=data[0].z==1?.02:.85;}}
  std::vector<std::array<double,P>> xs;std::vector<E>fdata=data;std::vector<std::array<double,2>>means,tracker;
  for(auto&e:fdata){auto x=features(e,m);xs.push_back(x);if(inject)e.u+=.010*x[1]+.006*x[6];std::array<double,2>a,v;for(int k=0;k<2;k++){a[k]=mu(e,k);double sens=(a[k]*a[k]-1)/a[k];v[k]=(4./6)*sens*sens*e.relR*e.relR;}means.push_back(a);tracker.push_back(v);}
  const double norm=std::lgamma(3.5)-std::lgamma(3.)-.5*std::log(6*M_PI);
  auto evaluate=[&](const double*p,double*g){
    if(g)std::fill(g,g+NP,0.);double s2[B],ff[B];for(int b=0;b<B;b++){s2[b]=std::exp(2*p[P+b]);ff[b]=fraction(p[P+B+b]);}double value=0;
    for(size_t i=0;i<fdata.size();i++){
      int b=fdata[i].bin;double c=0;for(int j=0;j<P;j++)c+=p[j]*xs[i][j];double residual[2],variance[2],lp[2];for(int k=0;k<2;k++){residual[k]=fdata[i].u-c-means[i][k];variance[k]=s2[b]+tracker[i][k];lp[k]=norm-.5*std::log(variance[k])-3.5*std::log1p(residual[k]*residual[k]/(6*variance[k]));}
      double la=std::log1p(-ff[b])+lp[0],lb=std::log(ff[b])+lp[1],mx=std::max(la,lb),logmix=mx+std::log(std::exp(la-mx)+std::exp(lb-mx));value-=logmix;
      if(g){double densityRatio[2]={std::exp(lp[0]-logmix),std::exp(lp[1]-logmix)};double weight[2]={(1-ff[b])*densityRatio[0],ff[b]*densityRatio[1]},gc=0,gs=0;for(int k=0;k<2;k++){double r=residual[k],v=variance[k],den=6*v+r*r;gc-=weight[k]*7*r/den;gs+=weight[k]*s2[b]/v*(1-7*r*r/den);}for(int j=0;j<P;j++)g[j]+=gc*xs[i][j];g[P+b]+=gs;g[P+B+b]+=densityRatio[0]-densityRatio[1];}
    }
    for(int j=1;j<P;j++){value+=fdata.size()*.0005*std::pow(p[j]/.03,2);if(g)g[j]+=fdata.size()*.001*p[j]/(.03*.03);}return value;
  };
  auto objective=[&](const double*p){return evaluate(p,nullptr);};auto gradient=[&](const double*p,double*g){evaluate(p,g);};
  static bool checked[3]={false,false,false};if(!checked[data[0].z]){auto pp=m.p;for(int b=0;b<B;b++)pp[P+B+b]=.3;double grad[NP];gradient(pp.data(),grad);double maxerr=0;for(int j=0;j<NP;j++){double step=j<P?1e-6:1e-5;auto pa=pp,pb=pp;pa[j]+=step;pb[j]-=step;double fd=(objective(pa.data())-objective(pb.data()))/(2*step);maxerr=std::max(maxerr,std::abs(fd-grad[j])/(1+std::abs(fd)));}std::cout<<"gradient_check Z="<<data[0].z<<" max_relative_error="<<maxerr<<std::endl;if(maxerr>2e-4)throw std::runtime_error("Analytic likelihood gradient check failed");checked[data[0].z]=true;}
  std::unique_ptr<ROOT::Math::Minimizer> minim(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));ROOT::Math::GradFunctor fun(objective,NP,gradient);minim->SetFunction(fun);minim->SetMaxFunctionCalls(12000);minim->SetMaxIterations(3000);minim->SetTolerance(.1);minim->SetPrintLevel(0);minim->SetErrorDef(.5);
  for(int j=0;j<P;j++){std::string name="calibration_"+std::to_string(j);if(mode==0||(mode==1&&j>0))minim->SetFixedVariable(j,name,0);else if(mode==3)minim->SetFixedVariable(j,name,m.p[j]);else minim->SetLimitedVariable(j,name,m.p[j],.0003,-.08,.08);}
  for(int b=0;b<B;b++){minim->SetLimitedVariable(P+b,"log_scale_"+std::to_string(b),m.p[P+b],.02,std::log(.003),std::log(.15));minim->SetLimitedVariable(P+B+b,"heavy_fraction_"+std::to_string(b),m.p[P+B+b],.01,0.,1.);}
  bool ok=minim->Minimize();m.status=minim->Status();m.edm=minim->Edm()/data.size();m.nll=minim->MinValue()/data.size();std::copy(minim->X(),minim->X()+NP,m.p.begin());std::cout<<"fit Z="<<data[0].z<<" mode="<<mode<<" injected="<<inject<<" N="<<data.size()<<" ok="<<ok<<" status="<<m.status<<" nll="<<m.nll<<" edm="<<m.edm<<std::endl;
  bool valid=ok&&m.status==0&&std::isfinite(m.nll)&&std::isfinite(m.edm);for(double p:m.p)valid&=std::isfinite(p);if(!valid){if(optional){m.nll=INFINITY;std::cout<<"REJECTED optional restart"<<std::endl;return m;}throw std::runtime_error("Invalid/nonconverged mass fit");}
  for(int j=0;j<P;j++)if(std::abs(m.p[j])>.079)std::cout<<"BOUNDARY calibration "<<j<<std::endl;for(int b=0;b<B;b++)if(m.p[P+B+b]<1e-5||m.p[P+B+b]>1-1e-5)std::cout<<"BOUNDARY fraction bin "<<b<<std::endl;return m;
}
double auc(const std::vector<double>&scores,const std::vector<int>&labels){std::vector<std::pair<double,int>>a;size_t n0=0,n1=0;for(size_t i=0;i<scores.size();i++)if(labels[i]>=0){a.push_back({scores[i],labels[i]});n1+=labels[i]==1;n0+=labels[i]==0;}if(!n0||!n1)return NAN;std::sort(a.begin(),a.end());double rank=0;for(size_t i=0;i<a.size();){size_t j=i+1;while(j<a.size()&&a[j].first==a[i].first)j++;for(size_t k=i;k<j;k++)if(a[k].second==1)rank+=.5*(i+1+j);i=j;}return (rank-.5*n1*(n1+1))/(double(n0)*n1);}
double referenceU0(const E&e){return std::sqrt(1+std::pow(e.Rref/e.R,2)*(e.ur*e.ur-1));}
Model fitReference(const std::vector<E>&data,const Model&base){
  Model m=base;std::vector<std::array<double,P>>xx;std::vector<double>target;for(auto&e:data)if(e.richOK){xx.push_back(features(e,m));target.push_back(e.u-referenceU0(e));}if(target.size()<500)throw std::runtime_error("Too few RICH training references");double delta=3*width(target),variance=std::pow(width(target),2);
  auto obj=[&](const double*p){double value=0;for(size_t i=0;i<target.size();i++){double r=target[i];for(int j=0;j<P;j++)r-=p[j]*xx[i][j];double a=std::abs(r);value+=a<delta?.5*r*r:delta*(a-.5*delta);}for(int j=1;j<P;j++)value+=target.size()*.0005*p[j]*p[j];return value/variance;};
  std::unique_ptr<ROOT::Math::Minimizer>minim(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));ROOT::Math::Functor fun(obj,P);minim->SetFunction(fun);minim->SetMaxFunctionCalls(12000);minim->SetTolerance(.1);minim->SetErrorDef(.5);minim->SetPrintLevel(0);for(int j=0;j<P;j++)minim->SetLimitedVariable(j,"reference_"+std::to_string(j),0,.0003,-.08,.08);bool ok=minim->Minimize();if(!ok||minim->Status()!=0||!std::isfinite(minim->MinValue()))throw std::runtime_error("Failed RICH reference fit");std::copy(minim->X(),minim->X()+P,m.p.begin());std::cout<<"RICH reference fit Z="<<data[0].z<<" N="<<target.size()<<" status="<<minim->Status()<<" value="<<minim->MinValue()<<std::endl;return m;
}
std::string method(int k){return k==0?"raw":k==1?"offset":k==2?"mass_calibrated":"rich_reference";}
void saveModel(std::ofstream&out,int z,const char*name,const Model&m){for(int j=0;j<NP;j++)out<<z<<'\t'<<name<<'\t'<<j<<'\t'<<m.p[j]<<'\t'<<(j<P?m.center[j]:0)<<'\t'<<(j<P?m.scale[j]:1)<<'\t'<<m.status<<'\t'<<m.edm<<'\n';}
}

void fitMassStudy(const char*input="results/mass_study/flight.root",const char*dest="results/mass_study"){
  using namespace MassStudy;gSystem->mkdir(dest,true);std::string dir(dest);TFile in(input);auto t=dynamic_cast<TTree*>(in.Get("flight"));if(!t)throw std::runtime_error("No extracted flight tree");
  E e;int richOK,naf;double Rscalar,RK[28],err[21],chi[189],fx[8];
  t->SetBranchAddress("run",&e.run);t->SetBranchAddress("event",&e.event);t->SetBranchAddress("role",&e.role);t->SetBranchAddress("source",&e.source);t->SetBranchAddress("Z",&e.z);t->SetBranchAddress("u",&e.u);t->SetBranchAddress("richU",&e.ur);t->SetBranchAddress("richOK",&richOK);t->SetBranchAddress("NaF",&naf);t->SetBranchAddress("features",fx);t->SetBranchAddress("R",&Rscalar);t->SetBranchAddress("rigidityK",RK);t->SetBranchAddress("erigidity",err);t->SetBranchAddress("chis",chi);
  std::vector<E>data[3][2];Long64_t selected=0;std::set<Long64_t>runs[2];
  for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i);if(e.z!=1&&e.z!=2)continue;e.R=RK[8];e.Rref=RK[15];e.relR=std::abs(err[15]*e.R);if(!std::isfinite(e.R)||e.R<EDGES[0]||e.R>=EDGES[B]||!std::isfinite(e.relR)||e.relR<=0||e.relR>=.20||!std::isfinite(chi[171])||!std::isfinite(chi[172])||chi[171]<=0||chi[172]<=0||chi[171]>10||chi[172]>10)continue;
    e.bin=std::upper_bound(EDGES,EDGES+B+1,e.R)-EDGES-1;e.richOK=richOK&&naf&&std::isfinite(e.Rref)&&e.Rref>0&&std::isfinite(e.ur)&&e.ur>1/.97&&e.ur<1/.78;e.mr=e.richOK?e.z*e.Rref*std::sqrt(e.ur*e.ur-1):NAN;e.tag=-1;
    // RICH tags are only validation diagnostics, never passed to a fitter.
    // They share tracker information; they are not isotope truth labels.
    if(e.richOK){if(e.z==1){if(e.mr>.78&&e.mr<1.10)e.tag=0;if(e.mr>1.65&&e.mr<2.05)e.tag=1;}else{if(e.mr>2.60&&e.mr<3.02)e.tag=0;if(e.mr>3.52&&e.mr<3.94)e.tag=1;}}
    std::copy(fx,fx+P,e.x.begin());data[e.z][e.role].push_back(e);runs[e.role].insert(e.run);selected++;
  }
  for(auto r:runs[0])if(runs[1].count(r))throw std::runtime_error("Train validation run overlap");
  std::cout<<"Selected="<<selected<<" train_runs="<<runs[0].size()<<" validation_runs="<<runs[1].size()<<std::endl;
  std::ofstream models(dir+"/models.tsv");models<<std::setprecision(14)<<"Z\tmodel\tparameter\tvalue\tfeature_center\tfeature_scale\tfit_status\tedm\n";
  std::ofstream stats(dir+"/metrics.tsv");stats<<std::setprecision(10)<<"Z\tmethod\tN_train\tN_holdout\tN_rich\tN_tag_light\tN_tag_heavy\tmean_nll\trich_residual_median\trich_residual_W68\ttag_AUC\tlight_mass2_median\tlight_mass2_W68\theavy_mass2_median\theavy_mass2_W68\tseparation_W68\n";
  std::ofstream bins(dir+"/bins.tsv");bins<<std::setprecision(10)<<"Z\tmethod\tR_low\tR_high\tN_holdout\tmean_nll\ttrain_width_scale\ttrain_heavy_fraction\tN_tag_light\tN_tag_heavy\ttag_score_AUC\ttag_mass2_AUC\n";
  std::ofstream closure(dir+"/injection.tsv");closure<<std::setprecision(10)<<"Z\tN_holdout\tinjected_RMS\trecovered_difference_RMS\tresidual_W68\n";
  std::ofstream rclosure(dir+"/rigidity_injection.tsv");rclosure<<std::setprecision(10)<<"Z\tN_holdout\tR_relative_injection_RMS\tmass2_change_RMS_before\tmass2_change_RMS_after\tmass2_change_W68_before\tmass2_change_W68_after\n";
  std::ofstream boot(dir+"/bootstrap.tsv");boot<<std::setprecision(10)<<"Z\treplicate\trich_width_narrowing\tauc_difference\n";
  TFile result((dir+"/predictions.root").c_str(),"RECREATE");TTree pred("predictions","Frozen flight mass calibration; no MC truth");int z,source,tag,rich,bin;Long64_t run,event;double R,Rref,raw,offset,calibrated,richU,massRich,m2raw,m2offset,m2cal,scoreRaw,scoreCal,nearestMass,nearestU,injected,restored,richU0,referenceU,m2reference,scoreReference,Rinjected,m2Rbefore,m2Rafter;
  pred.Branch("Rinjected",&Rinjected);pred.Branch("mass2RigidityInjectionBefore",&m2Rbefore);pred.Branch("mass2RigidityInjectionAfter",&m2Rafter);
  pred.Branch("richU0",&richU0);pred.Branch("referenceU",&referenceU);pred.Branch("mass2Reference",&m2reference);pred.Branch("scoreReference",&scoreReference);
  pred.Branch("Z",&z);pred.Branch("source",&source);pred.Branch("run",&run);pred.Branch("event",&event);pred.Branch("R",&R);pred.Branch("Rref",&Rref);pred.Branch("rawU",&raw);pred.Branch("offsetU",&offset);pred.Branch("calibratedU",&calibrated);pred.Branch("richU",&richU);pred.Branch("massRich",&massRich);pred.Branch("richOK",&rich);pred.Branch("tag",&tag);pred.Branch("bin",&bin);pred.Branch("mass2Raw",&m2raw);pred.Branch("mass2Offset",&m2offset);pred.Branch("mass2Calibrated",&m2cal);pred.Branch("scoreRaw",&scoreRaw);pred.Branch("scoreCalibrated",&scoreCal);pred.Branch("nearestMass",&nearestMass);pred.Branch("nearestU",&nearestU);pred.Branch("injectedU",&injected);pred.Branch("restoredU",&restored);
  for(z=1;z<=2;z++){
    auto originalTrain=data[z][0];auto train=originalTrain;auto&val=data[z][1];std::mt19937 rng(2981+z);std::shuffle(train.begin(),train.end(),rng);if(train.size()>60000)train.resize(60000);
    // Strip validation reference fields from fitting inputs by construction.
    for(auto&v:train){v.richOK=0;v.tag=-1;v.ur=NAN;v.mr=NAN;v.Rref=NAN;}
    Model m0=fit(train,0);Model start=m0;for(int j=0;j<B;j++){start.p[P+j]+=std::log(1.2);start.p[P+B+j]=z==1?.15:.5;}Model m0b=fit(train,0,&start,false,true);if(m0b.nll<m0.nll)m0=m0b;Model m1=fit(train,1,&m0),m2=fit(train,2,&m1);start=m0;start.p[0]=-.015;for(int j=0;j<B;j++)start.p[P+B+j]=z==1?.15:.5;Model m2b=fit(train,2,&start,false,true);if(m2b.nll<m2.nll)m2=m2b;Model mi=fit(train,2,&m2,true);Model mr0=fitReference(originalTrain,m0),mr=fit(train,3,&mr0);std::array<Model,4>ms={m0,m1,m2,mr};
    auto rtrain=train;for(auto&v:rtrain){double factor=std::exp(.03*features(v,m2)[1]);v.R*=factor;v.relR*=factor;}Model mri=fit(rtrain,2,&m2);
    saveModel(models,z,"raw",m0);saveModel(models,z,"offset",m1);saveModel(models,z,"mass_calibrated",m2);saveModel(models,z,"injection",mi);saveModel(models,z,"rich_reference",mr);saveModel(models,z,"rigidity_injection",mri);models.flush();
    std::vector<double>ur[4],ll[4],a[4],b[4],sc[4],bs[4][B],bm[4][B];std::vector<int>labels,bl[B];std::vector<double>injections,closes,rinj,rclose0,rclose1;double bnll[4][B]={};int bn[B]={};std::set<Long64_t>validationRuns;for(auto&v:val)validationRuns.insert(v.run);
    for(auto&v:val){R=v.R;Rref=v.Rref;source=v.source;run=v.run;event=v.event;tag=v.tag;rich=v.richOK;bin=v.bin;raw=v.u;richU=v.ur;richU0=v.richOK?referenceU0(v):NAN;massRich=v.mr;offset=v.u-correction(v,m1);calibrated=v.u-correction(v,m2);referenceU=v.u-correction(v,mr);double uv[4]={raw,offset,calibrated,referenceU};bn[bin]++;labels.push_back(tag);
      bl[bin].push_back(tag);for(int k=0;k<4;k++){double n=nllone(v,uv[k],ms[k].p.data());ll[k].push_back(n);bnll[k][bin]+=n;if(v.richOK)ur[k].push_back(uv[k]-richU0);double mm=std::pow(v.z*v.R,2)*(uv[k]*uv[k]-1);if(v.tag==0)a[k].push_back(mm);if(v.tag==1)b[k].push_back(mm);sc[k].push_back(score(v,uv[k],ms[k]));bs[k][bin].push_back(sc[k].back());bm[k][bin].push_back(mm);}
      m2raw=std::pow(z*R,2)*(raw*raw-1);m2offset=std::pow(z*R,2)*(offset*offset-1);m2cal=std::pow(z*R,2)*(calibrated*calibrated-1);m2reference=std::pow(z*R,2)*(referenceU*referenceU-1);scoreRaw=sc[0].back();scoreCal=sc[2].back();scoreReference=sc[3].back();nearestMass=mass(z,std::abs(m2raw-mass(z,1)*mass(z,1))<std::abs(m2raw-mass(z,0)*mass(z,0)));nearestU=std::sqrt(1+std::pow(nearestMass/(z*R),2));
      auto x=features(v,m2);double delta=.010*x[1]+.006*x[6];injected=raw+delta;restored=injected-correction(v,mi);injections.push_back(delta);closes.push_back(restored-calibrated);
      double rf=std::exp(.03*x[1]);Rinjected=R*rf;E rv=v;rv.R=Rinjected;rv.relR*=rf;double ru=rv.u-correction(rv,mri);m2Rbefore=rf*rf*m2cal;m2Rafter=std::pow(z*Rinjected,2)*(ru*ru-1);rinj.push_back(rf-1);rclose0.push_back(m2Rbefore-m2cal);rclose1.push_back(m2Rafter-m2cal);pred.Fill();
    }
    for(int k=0;k<4;k++){double sep=(quantile(b[k],.5)-quantile(a[k],.5))/std::sqrt(std::pow(width(a[k]),2)+std::pow(width(b[k]),2));stats<<z<<'\t'<<method(k)<<'\t'<<train.size()<<'\t'<<val.size()<<'\t'<<ur[k].size()<<'\t'<<a[k].size()<<'\t'<<b[k].size()<<'\t'<<std::accumulate(ll[k].begin(),ll[k].end(),0.)/ll[k].size()<<'\t'<<quantile(ur[k],.5)<<'\t'<<width(ur[k])<<'\t'<<auc(sc[k],labels)<<'\t'<<quantile(a[k],.5)<<'\t'<<width(a[k])<<'\t'<<quantile(b[k],.5)<<'\t'<<width(b[k])<<'\t'<<sep<<'\n';for(int j=0;j<B;j++)bins<<z<<'\t'<<method(k)<<'\t'<<EDGES[j]<<'\t'<<EDGES[j+1]<<'\t'<<bn[j]<<'\t'<<bnll[k][j]/std::max(1,bn[j])<<'\t'<<std::exp(ms[k].p[P+j])<<'\t'<<fraction(ms[k].p[P+B+j])<<'\t'<<std::count(bl[j].begin(),bl[j].end(),0)<<'\t'<<std::count(bl[j].begin(),bl[j].end(),1)<<'\t'<<auc(bs[k][j],bl[j])<<'\t'<<auc(bm[k][j],bl[j])<<'\n';}
    auto rms=[](const std::vector<double>&v){double s=0;for(double x:v)s+=x*x;return std::sqrt(s/v.size());};closure<<z<<'\t'<<val.size()<<'\t'<<rms(injections)<<'\t'<<rms(closes)<<'\t'<<width(closes)<<'\n';rclosure<<z<<'\t'<<val.size()<<'\t'<<rms(rinj)<<'\t'<<rms(rclose0)<<'\t'<<rms(rclose1)<<'\t'<<width(rclose0)<<'\t'<<width(rclose1)<<'\n';
    // Paired whole-run bootstrap, fixed calibration and fixed tag windows.
    std::vector<Long64_t>src(validationRuns.begin(),validationRuns.end());std::mt19937 brng(8831+z);std::uniform_int_distribution<int>pick(0,src.size()-1);
    for(int bidx=0;bidx<200;bidx++){std::map<Long64_t,int>w;for(size_t i=0;i<src.size();i++)w[src[pick(brng)]]++;std::vector<double>r0,r2,s0,s2;std::vector<int>lab;for(auto&v:val){int n=w[v.run];for(int j=0;j<n;j++){if(v.richOK){double ru=referenceU0(v);r0.push_back(v.u-ru);r2.push_back(v.u-correction(v,m2)-ru);}if(v.tag>=0){s0.push_back(score(v,v.u,m0));s2.push_back(score(v,v.u-correction(v,m2),m2));lab.push_back(v.tag);}}}boot<<z<<'\t'<<bidx<<'\t'<<(1-width(r2)/width(r0))<<'\t'<<(auc(s2,lab)-auc(s0,lab))<<'\n';}
    std::cout<<"Z="<<z<<" holdout="<<val.size()<<" RICH="<<ur[0].size()<<" tags="<<a[0].size()<<","<<b[0].size()<<" width ratio="<<width(ur[2])/width(ur[0])<<" AUC="<<auc(sc[0],labels)<<" -> "<<auc(sc[2],labels)<<std::endl;
  }
  result.cd();pred.Write();TNamed info("protocol","Flight only; fixed eight affine-slope-invariant geometry/time-contrast features. No R, fitted beta, mass or RICH in correction features. t6 isotope mixture with four free fractions and scales. R=Kalman V6 inner z0; RICH diagnostic R at z=-70; radiator is lower. Tracker error propagation is approximate. Models frozen before holdout metrics. No MC opened.");info.Write();result.Close();std::cout<<"DONE mass study"<<std::endl;
}
