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
  auto psi = IQMPS(state);

  int counter =0;
  while(true)
  {
    for(int i=1; i < N; ++i)
    {
      IQTensor AA = psi.A(i) * psi.A(i+1);
      IQTensor U(commonIndex(psi.A(i),AA,Site), commonIndex(psi.A(i),AA,Link) ), D, V;
      svd(AA, U, D, V, args);   
    }
  }
}
