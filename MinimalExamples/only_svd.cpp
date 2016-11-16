#include "itensor/all.h"
#include <iostream>

using namespace itensor;
using namespace std;

int
main()
{


  while(true)
  {
    //Make random ITensor with indices k,l,m
    auto k = Index("index k",5);
    auto l = Index("index l",15);
    auto m = Index("index m",5);
    auto T = randomTensor(k,l,m);
    ITensor U(k,m),S,V;
    svd(T,U,S,V,{"Cutoff",1E-12,"Maxm",10});
  }
}
