#include "itensor/all.h"
#include <iostream>

using namespace itensor;
using namespace std;

int
main()
{

  int N = 4;
  auto args = Args("Cutoff=",1E-12,"Maxm=",10);

  SpinHalf sites(N);

  auto state = InitState(sites);
  for (int i=1; i <= N; ++i)
    state.set(i, i % 2 != 0 ? "Up" : "Dn");
  auto psi = MPS(state);

  while(true)
  {
    for(int i=1; i < N; ++i)
    {
      ITensor AA = psi.A(i) * psi.A(i+1);

      //Make random ITensor with indices k,l,m
      auto k = Index("index k",5);
      auto l = Index("index l",15);
      auto m = Index("index m",5);
      auto T = randomTensor(k,l,m);
      ITensor U(k,m),S,V;
      svd(T,U,S,V,{"Cutoff",1E-12,"Maxm",10});
    }
  }
}
