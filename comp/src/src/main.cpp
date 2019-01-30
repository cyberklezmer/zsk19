#include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
#include "mspp/msproblem.h"
#include "mspp/mspptest.h"
//#include "mspp/cplex.h"
#include "zsk.h"
#include <mcheck.h>

using namespace mspp;



int main(int, char **)
{
//    int res=mcheck(nullptr);
//if(res)
//        cerr << "error " << res << " starting mcheck" << endl;

    std::ofstream logf("zsk.log");
    sys::setlog(logf);

    sys::seed(0);
    using O=csvlpsolver<realvar>;
//   using O=cplex<realvar>;


    assert(kappa<=maxkappa);
    vector<vector<double>> etasqV;
    etasqV.push_back({epsvsq[0][0],epsvsq[0][1],0,0,0});
    etasqV.push_back({epsvsq[1][0],epsvsq[1][1],0,0,0});
    etasqV.push_back({0,0,varsigma,0,0});
    etasqV.push_back({0,0,0,varsigma,0});
    etasqV.push_back({0,0,0,0,varsigma});

    unsigned int nl = 1;


    using stdd_t=ldistribution<double>;

    using etad_t = ldistribution<vector<double>>;

    etad_t etad({{0,0,qcoef,qcoef,qcoef}});

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = iidprocessdistribution<etad_t>;

    etapd_t etapd(dirac_t({0,0,v0[0],v0[1],v0[2]}),etad,T+1);


    vector<almdistribution> d;
    vector<mmcdistribution> m;

    double lspot0 = log(spot0);

//    m.push_back(mmcdistribution({{0.5,0.5}}));
vector<vector<double>> um(1);
um[0].push_back(1);
m.push_back(mmcdistribution(um));
    d.push_back(almdistribution(nl,{ lspot0-0.1, lspot0+0.1}));


    if(T>1)
    {
//        m.push_back(mmcdistribution({{0.5,0.5,0},{0,0.5,0.5}}));
m.push_back(mmcdistribution(um));
        d.push_back(almdistribution(nl,{  lspot0-0.2, lspot0, lspot0+0.2 }));
    }
    if(T>2)
    {
//        m.push_back(mmcdistribution({{0.5,0.5,0,0},{0,0.5,0.5,0}
//                                     ,{0,0,0.5,0.5}}));
m.push_back(mmcdistribution(um));

        d.push_back(almdistribution(nl,{ lspot0-0.3, lspot0-0.1,
                                          lspot0+0.1, lspot0+0.3 } ));
    }
    using hmcpd_t = hmcprocessdistribution<mmcdistribution,almdistribution>;

    hmcpd_t xipd(lspot0,m,d);


    using dxieta_t = dxietaprocessdist<hmcpd_t,etapd_t>;

    dxieta_t xieta(xipd,etapd);

    using zeta_t = hmczeta<pair<double, vector<double>>,zskm>;

    zskproblem p(0.5,0.05);

    static_assert(std::is_same<typename dxieta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);

    static_assert(std::is_same<typename zeta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);
    desolution<zskproblem,dxieta_t,zeta_t,O> sol(p,xieta);

//    vector<unsigned int> zd(T,zskm::zetasize);
//    msddpsolution<zskproblem,xieta_t,zeta_t,O> sx(p,xieta,zd);

    return 0;
}

/*

    using stdd_t=ldistribution<double>;
//    stdd_t stdd({-1,1});
    stdd_t stdd({0.0});

    using etad_t = meanvardistribution<stdd_t>;

    etad_t etad({0,0,qcoef,qcoef,qcoef},etasqV,stdd);

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = iidprocessdistribution<etad_t>;

    etapd_t etapd(dirac_t({0,0,v0[0],v0[1],v0[2]}),etad,T+1);




   double xim = -sigma*sigma/2.0;
    double xisd = sigma;

    arnormalprocessdistribution xipd(0,xim,xisd,1.0,T);
      logpd jako první

    vector<unsigned int> nx;
    for(unsigned int i=1;i<=T; i++)
        nx.push_back(i);

    using ha_t
       = chmcapproximation<arnormalprocessdistribution,onedcovering>;

    ha_t ha(xipd,nx);

    using xieta_t = xietaprocessdist<ha_t,etapd_t>;
    xieta_t xieta(ha,eta); */
