#include "itensor/all.h"

using namespace itensor;
using namespace std;

int
main()
{

  int N = 4;
  double tstep = 0.1;
  auto args = Args("Cutoff=",1E-12,"Maxm=",10);

  SpinHalf sites(N);

  auto state = InitState(sites);
  for (int i=1; i <= N; ++i)
    state.set(i, i % 2 != 0 ? "Up" : "Dn");
  auto psi = IQMPS(state);


  vector<IQGate> gates;
  auto type = IQGate::tReal;
  for(int b = 1; b < N; ++b)
  {
    IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
    hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
    hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
    gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
  }
  for(int b = N-1; b >= 1; --b)
  {
    IQTensor hh = sites.op("Sz",b)*sites.op("Sz",b+1);
    hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
    hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
    gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
  }

  while(true)
  {
    for(auto& G : gates)
    {
      auto b = G.i1();
      // psi.position(b);
      IQTensor AA = psi.A(b) * psi.A(b+1) * G;
    }
  }
}
