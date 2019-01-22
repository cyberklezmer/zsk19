#include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
#include "mspp/msproblem.h"
//#include "mspp/cplex.h"
#include "zsk.h"
#include <mcheck.h>

using namespace mspp;


/*
template <typename O>
void alm1test(bool equivalent)
{
    almproblem mp(0.5,0.05,1);
    lvdistribution<double> d({{0.5500},{0.8500}, {1.0500},{1.3500}});
    using dist = fdprocessdistribution<lvdistribution<double>,noxi<vector<double>>>;
    dist pd({1},d,1);

    if(!equivalent)
    {
      std::cout <<"MALMtest DE direct..." << std::endl;

      desolution<almproblem,dist,lastxi<vector<double>>,O> x(mp,pd);
      if(fabs(x.obj()-1) > 0.001)
      {
          std::cerr << "alm1test: opt="
               << 1 << " expected, " << x.obj() << " achieved." << std::endl;
          throw;
      }

      sddpsolution<almproblem,dist,lastxi<vector<double>>,O> sx(mp,pd);
      if(fabs(sx.obj().lb()-x.obj()) > 0.05)
      {
          std::cerr << "alm1test: opt="
               << x.obj() << " expected, " <<
               sx.obj().lb() << "<" << sx.obj().ubm() <<
               " achieved." << std::endl;
          throw;
      }
    }
    else
    {
          std::cout <<"MALMtest DE equivalent..." << std::endl;

          mpmcvarequivalent<almproblem> ep(mp);

          desolution<mpmcvarequivalent<almproblem>,dist,
              lastxi<vector<double>>,O> x(ep,pd);
          if(fabs(x.obj()-1) > 0.001)
          {
              std::cerr << "alm1test: opt="
                   << 1 << " expected, " << x.obj() << " achieved." << std::endl;
              throw;
          }

          sddpsolution<mpmcvarequivalent<almproblem>,dist,lastxi<vector<double>>,O> sx(ep,pd);
          if(fabs(sx.obj().lb()-x.obj()) > 0.05)
          {
              std::cerr << "alm1test: opt="
                   << x.obj() << " expected, " <<
                   sx.obj().lb() << "<" << sx.obj().ubm() <<
                   " achieved." << std::endl;
              throw;
          }
    }
    std::cout << "passed" << std::endl;
}

template <typename O>
void almtest(unsigned int T=3, unsigned int nl=1)
{
    assert(T>=1);
    const double tol = 1e-5;
    double sol[]={
            0,0,0,1,1,0,1,0,1,0,1,0,1,
            0.727273,0.727273,0.272727,1,0.272727,
            1,0.272727,1,0.272727,1,0,0,1,1,1,1,1,
            1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,
            1,0,1,0.727273,0.727273,0.272727,1,0.272727,
            1,0.272727,1,0.272727,1};
    double obj = 0.906063;
    unsigned int k = sizeof(sol)/sizeof(sol[0]);


    using mydist = fhmcdistribution<mmcdistribution,almdistribution>;
    vector<mydist> d;

    mydist d1(mmcdistribution({{0.5,0.5}}),almdistribution(nl,{ 0.8, 1.1}));

    d.push_back(d1);

    if(T>1)
        d.push_back(mydist(mmcdistribution({{0.5,0.5,0},{0,0.5,0.5}}),
            almdistribution(nl,{ 0.8*0.8, 1.1*0.8, 1.1*1.1 })));

    if(T>2)
        d.push_back(mydist(mmcdistribution({{0.5,0.5,0,0},{0,0.5,0.5,0}
                                     ,{0,0,0.5,0.5}}),
          almdistribution(nl,{ 0.8*0.8*0.8, 0.8*0.8*1.1,
                                          0.8*1.1*1.1, 1.1*1.1*1.1 } )));

    using pdt = fdprocessdistribution<mydist,laststate<typename almdistribution::I_t>>;
    pdt pd({0,{1}},d);

    almproblem mp(0.5,0.05,T);

    std::cout <<"MALMtest DE..." << std::endl;

    desolution<almproblem,pdt,lastmdxi<vector<double>>,O> x(mp,pd);

    if(T==3 && nl==1 && fabs(x.obj()-obj) > tol)
    {
        std::cerr << "almtest: opt="
             << obj << " expected, " << x.obj() << " achieved." << std::endl;
        throw;
    }

    vector<double> xs;
    x.x()->exportlinear(xs);


    if(T==3 && nl==1)
        for(unsigned int i=0; i<k; i++)
        {
            if(fabs(xs[i]-sol[i]) > tol)
            {
                std::cerr << "x[" << i << "]="
                     << sol[i] << " expected, " << xs[i]
                     << " achieved." << std::endl;
                std::cerr << endl;
                throw;
            }
        }

    std::cout <<  "obj=." <<x.obj() << std::endl;

    std::cout <<"ALMtest Markov SDDP..." << std::endl;

    mpmcvarequivalent<almproblem> ep(mp);


//    msddpsolution<mpmcvarequivalent<almproblem>,pdt,lastmdxi<vector<double>>,O> sx(ep,pd);
    msddpsolution<almproblem,pdt,lastmdxi<vector<double>>,O> sx(mp,pd);

    if(fabs(sx.obj().lb()-x.obj()) > 0.1)
    {
        std::cerr << "malmtest: opt="
             << x.obj() << " expected, " <<
             sx.obj().lb() << "<" << sx.obj().ubm() <<
             " achieved." << std::endl;
        throw;
    }

    std::cout << "almtest: lb=" <<
         sx.obj().lb() << " < " << "opt=" << x.obj()
              << " < ub=" << sx.obj().ubm() << std::endl;
    std::cout <<  "passed." << std::endl;
}

*/

int main(int, char **)
{
    int res=mcheck(nullptr);
    if(res)
        cerr << "error " << res << " starting mcheck" << endl;

    std::ofstream log("zsk.log");
    sys::setlog(log);

    sys::seed(0);
//    using O=csvlpsolver<realvar>;
    using O=cplex<realvar>;


    assert(nstrikes<=maxstrikes);
    vector<vector<double>> etasqV(5);
    etasqV.push_back({epsvsq[0][0],epsvsq[0][1],0,0,0});
    etasqV.push_back({epsvsq[1][0],epsvsq[1][1],0,0,0});
    etasqV.push_back({0,0,varsigma,0,0});
    etasqV.push_back({0,0,0,varsigma,0});
    etasqV.push_back({0,0,0,0,varsigma});

    using stdd_t=ldistribution<double>;
    stdd_t stdd({-1,1});

    using etad_t = meanvardistribution<stdd_t>;

    etad_t etad({0,0,qcoef,qcoef,qcoef},etasqV,stdd);

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = iidprocessdistribution<etad_t>;

    etapd_t eta(dirac_t({0,0,v0[0],v0[1],v0[2]}),etad,T+1);

    double xim = -sigma*sigma/2.0;
    double xisd = sigma;

    arnormalprocessdistribution xipd(0,xim,xisd,1.0,T);

    vector<unsigned int> nx;
    for(unsigned int i=1;i<=T; i++)
        nx.push_back(i);

    using ha_t
       = chmcapproximation<arnormalprocessdistribution,onedcovering>;

    ha_t ha(xipd,nx);

    using xieta_t = xietaprocessdist<ha_t,etapd_t>;
    xieta_t xieta(ha,eta);

    using zeta_t = hmczeta<pair<double, vector<double>>,zskm>;


    return 0;
}
