#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/input.h"
#include "itensor/mps/dmrg.h"
#include "itensor/mps/tevol.h"
#include <iomanip>      // std::setprecision
using namespace std;
using namespace itensor;

int main(int argc, char* argv[])
{
  //
  //Parse the input file
  if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
  //
  //Get input parameters
  InputGroup basic(argv[1],"basic");
  auto N = basic.getInt("N"); //the 'M' stands for mandatory
  auto ttotal = basic.getReal("ttotal");
  auto tstep = basic.getReal("tstep");
  auto trotter = basic.getString("trotter");
  auto DoNormalize = basic.getString("DoNormalize");
  auto Normalize = basic.getString("Normalize");
  auto initialstate = basic.getString("initialstate");
  auto args = Args("Cutoff=",1E-9,"Maxm=",15);

  //
  //open the output file  
  string outputfilename1;
  string outputfilename2;
  outputfilename1.append(argv[1]);
  outputfilename1.append("_iqtmps_energy_output");
  outputfilename2.append(argv[1]);
  outputfilename2.append("_iqtmps_sz_output");    
  ofstream outputfile1;
  ofstream outputfile2;
  outputfile1.open (outputfilename1);
  outputfile2.open (outputfilename2);   
  //
  // Initialize the site degrees of freedom.
  SpinHalf sites(N);
  //
  // Generate the Hamiltonian
  // IQMPO H = Heisenberg(sites,opts); 
  //
  // Initialize the state
  InitState initState(sites);
  if (initialstate=="halfuphalfdown")
    {
      for(int i = 1; i <= N; ++i) 
        {
	  if(i<N/2+1)
	    //if(i%2 == 1)
            initState.set(i,"Up");
	  else
            initState.set(i,"Dn");
        }
    }
  if (initialstate=="'oneuponedown'")
    {
      for(int i = 1; i <= N; ++i) 
        {
	  //if(i<N/2+1)
	  if(i%2 == 1)
            initState.set(i,"Up");
	  else
            initState.set(i,"Dn");
        }
    }
  IQMPS psi(initState);

  //
  //create gates for time evolution
  vector<IQGate> gates;
  auto type = IQGate::tReal;
  //BondGate use the expansion of exp up to order 100 to calculate exp -iHt 
  if (trotter=="false")
    {
      for(int b = 1; b < N; ++b)
        {
	  IQTensor hh = 0.0 * sites.op("Sz",b)*sites.op("Sz",b+1);
	  hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
	  hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
	  gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
        }
      for(int b = N-1; b >= 1; --b)
        {
	  IQTensor hh = 0.0 * sites.op("Sz",b)*sites.op("Sz",b+1);
	  hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
	  hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
	  gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
        }
    }
  else
    {
      for(int b = 1; b < N; b+=2)
        {
	  IQTensor hh = 0.0 * sites.op("Sz",b)*sites.op("Sz",b+1);
	  hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
	  hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
	  gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
        }
      for(int b = N-2; b >= 2; b-=2)
        {
	  IQTensor hh = 0.0 * sites.op("Sz",b)*sites.op("Sz",b+1);
	  hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
	  hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
	  gates.push_back(IQGate(sites,b,b+1,type,tstep/1.,hh));
        }
      for(int b = 1; b < N; b+=2)
        {
	  IQTensor hh = 0.0 * sites.op("Sz",b)*sites.op("Sz",b+1);
	  hh += 0.5*sites.op("Sp",b)*sites.op("Sm",b+1);
	  hh += 0.5*sites.op("Sm",b)*sites.op("Sp",b+1);
	  gates.push_back(IQGate(sites,b,b+1,type,tstep/2.,hh));
        }
    }
  //
  //gateTEvol starts
  const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));    
  if(fabs(nt*tstep-ttotal) > 1E-9)
    {
      Error("Timestep not commensurate with total time");
    }
  Real tsofar = 0;
  Real tot_norm = psi.normalize();
  int step=0;
  for(int b=1; b<=N; ++b)
    {
      psi.position(b);
      Real obsSz= (dag(prime(psi.A(b),Site))*sites.op("Sz",b)*psi.A(b)).real();      
      //printfln("b=%4.0d %6d %23.20f",b,",     Sz =",obsSz);
      outputfile2 <<step<< "  " <<b<< "  " << setprecision(16)<< obsSz << endl;
    }
  for(int step = 1; step <= nt; ++step)
    {
      for(auto& G : gates)
	{
	  auto b = G.i1();
	  psi.position(b);
	  IQTensor AA = psi.A(b) * psi.A(b+1) * G;
	  AA.noprime();
	  IQTensor D;
	  Spectrum res = svd(AA,psi.Anc(b),D,psi.Anc(b+1),args);
	  if(DoNormalize=="false")
	    {
	      //  D *= 1./D.norm();
	    }  
	  psi.Anc(b+1) *= D;
	}
      if(Normalize=="false")
	{
	  tot_norm *= psi.normalize();
	}
      tsofar += tstep;
      //
      //prints energy and sz
      if(step % 1000 == 0)
        printfln("Step %d/%d",step,nt);

      for(int b=1; b<=N; ++b)
	{
	  psi.position(b);
	  Real obsSz= (dag(prime(psi.A(b),Site))*sites.op("Sz",b)*psi.A(b)).cplx().real();
	  //printfln("b=%4.0d %6d %23.20f",b,",     Sz =",obsSz);
	  outputfile2 <<step<< "  " <<b<< "  " << setprecision(16)<< obsSz << endl;
	}
    }
  outputfile1.close();
  outputfile2.close();
  return 0;
}
