// Self-contained squared-error histogram gradient boosting. No external download.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
#include <exception>
struct Node {std::vector<int> rows; int feature=-1,cut=0,left=-1,right=-1,depth=0; double value=0,gain=-1; int bf=-1,bc=0;};
extern "C" int fit_hgb(const unsigned char *X,int n,int d,const double *y,const double *weight,const int *train,int nt,int trees,int maxleaves,int maxdepth,int minleaf,double lambda,double eta,const int *snap,int ns,double *predsnap,int *features,int *cuts,int *lefts,int *rights,double *values,double *initial){
 try {
  if(nt<2*minleaf||n<=0||d<=0||maxleaves<2)return 2;
  int maxnodes=2*maxleaves-1;
  double sy=0,sw=0;for(int j=0;j<nt;++j){int i=train[j];if(i<0||i>=n||!std::isfinite(y[i]))return 3;sy+=weight[i]*y[i];sw+=weight[i];}
  *initial=sy/sw;std::vector<double> pred(n,*initial),res(n,0);
  auto assess=[&](Node &node){
   double G=0,H=0;for(int i:node.rows){G+=weight[i]*res[i];H+=weight[i];}
   node.value=G/(H+lambda);node.gain=-1;node.bf=-1;
   if(node.depth>=maxdepth||(int)node.rows.size()<2*minleaf)return;
   double parent=G*G/(H+lambda);
   for(int j=0;j<d;++j){
    std::array<double,256> hs{},gs{};std::array<int,256> counts{};
    const unsigned char *col=X+(size_t)j*n;
    for(int i:node.rows){int b=col[i];counts[b]++;hs[b]+=weight[i];gs[b]+=weight[i]*res[i];}
    int count=0;double h=0,g=0;
    for(int b=0;b<255;++b){count+=counts[b];h+=hs[b];g+=gs[b];
     if(count<minleaf||(int)node.rows.size()-count<minleaf||h<=0||H-h<=0)continue;
     double gain=g*g/(h+lambda)+(G-g)*(G-g)/(H-h+lambda)-parent;
     if(gain>node.gain){node.gain=gain;node.bf=j;node.bc=b;}
    }
   }
  };
  for(int tree=0;tree<trees;++tree){
   for(int j=0;j<nt;++j){int i=train[j];res[i]=y[i]-pred[i];}
   std::vector<Node> nodes;nodes.reserve(maxnodes);nodes.emplace_back();nodes[0].rows.assign(train,train+nt);assess(nodes[0]);
   int leaves=1;
   while(leaves<maxleaves){
    int best=-1;double gain=1e-18;
    for(int j=0;j<(int)nodes.size();++j)if(nodes[j].left<0&&nodes[j].bf>=0&&nodes[j].gain>gain){gain=nodes[j].gain;best=j;}
    if(best<0)break;
    int f=nodes[best].bf,c=nodes[best].bc,depth=nodes[best].depth;
    Node l,r;l.depth=r.depth=depth+1;
    const unsigned char *col=X+(size_t)f*n;
    for(int i:nodes[best].rows){if(col[i]<=c)l.rows.push_back(i);else r.rows.push_back(i);}
    if((int)l.rows.size()<minleaf||(int)r.rows.size()<minleaf)return 4;
    int il=nodes.size(),ir=il+1;nodes[best].feature=f;nodes[best].cut=c;nodes[best].left=il;nodes[best].right=ir;std::vector<int>().swap(nodes[best].rows);
    assess(l);assess(r);nodes.push_back(std::move(l));nodes.push_back(std::move(r));++leaves;
   }
   for(int j=0;j<maxnodes;++j){size_t k=(size_t)tree*maxnodes+j;features[k]=-1;cuts[k]=0;lefts[k]=rights[k]=-1;values[k]=0;
    if(j<(int)nodes.size()){features[k]=nodes[j].feature;cuts[k]=nodes[j].cut;lefts[k]=nodes[j].left;rights[k]=nodes[j].right;values[k]=nodes[j].value;}
   }
   for(int i=0;i<n;++i){int j=0;while(nodes[j].left>=0)j=(X[(size_t)nodes[j].feature*n+i]<=nodes[j].cut?nodes[j].left:nodes[j].right);pred[i]+=eta*nodes[j].value;}
   for(int k=0;k<ns;++k)if(tree+1==snap[k])std::copy(pred.begin(),pred.end(),predsnap+(size_t)k*n);
   if((tree+1)%160==0){std::fprintf(stdout,"BOOST %d/%d n=%d d=%d\n",tree+1,trees,nt,d);std::fflush(stdout);}
  }
  return 0;
 }catch(const std::exception &e){std::fprintf(stderr,"HGB_ERROR %s\n",e.what());return 9;}
}
