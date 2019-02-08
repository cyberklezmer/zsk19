//#define RISKNEUTRAL

#include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
#include "mspp/msproblem.h"
#include "mspp/mspptest.h"
//#include "mspp/cplex.h"
#include "zsk.h"
#include <mcheck.h>

using namespace mspp;

vectors<double> etasqV(unsigned int k)
{
    assert(k<=T);
    assert(k>0);

    vectors<double>r(5-k);
    for(unsigned int i=0; i<r.size(); i++)
        r[i]=vector<double>(5-k,0.0);
    r[0][0] = epsvsq[0][0];
    r[0][1] = epsvsq[0][1];
    r[1][0] = epsvsq[1][0];
    r[1][1] = epsvsq[1][1];

    for(unsigned int i=2; i<5-k; i++)
        r[i][i] = varsigma;
    return r;
}

struct dtestparams
{
    dtestparams() : T(3), trivialm(false),
        almleaves(1), etaleaves(1),delta(0.2),
        lambda(0.5), equivalent(false), sddp(false) {}

    unsigned int T;
    bool trivialm;
    unsigned int almleaves;
    unsigned int etaleaves;
    double delta;
    double lambda;
    bool equivalent;
    bool sddp;
};


void dtest(const dtestparams& p)
{
    unsigned int saveT = T;
    T = p.T;
    using stdd_t=ldistribution<double>;

    stdd_t stdd = p.etaleaves == 2 ? stdd_t(vector<double>({-1.0,1.0}))
                : (p.etaleaves == 3  ? stdd_t(vector<double>({-sqrt(1.5),0,sqrt(1.5)}))
                                       : stdd_t(vector<double>({0.0})));

    using etad_t = fmeanvardistribution<stdd_t>;
    vector<etad_t> dsts;
    for(unsigned int k=1; k<=T; k++)
    {
        vector<double> mean(5-k,0.0);
        for(unsigned int i=2; i<5-k; i++)
            mean[i]=qcoef;
        dsts.push_back(etad_t(mean,etasqV(k),stdd));
    }

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = processdistribution<etad_t,noxi<vector<double>>>;

    etapd_t etapd(dirac_t({0,0,v0[0],v0[1],v0[2]}),dsts);


    vector<almdistribution> d;
    vector<mmcdistribution> m;

    double lspot0 = log(spot0);

    if(p.trivialm)
    {
        for(unsigned int i=0; i<T; i++)
        {
            m.push_back(mmcdistribution(vector<vector<double>>({{1.0}})));
            d.push_back(almdistribution(p.almleaves,{ lspot0 },p.delta));
        }
    }
    else
    {
        double halfdelta = p.delta/2.0;
        m.push_back(mmcdistribution({{0.5,0.5}}));
        d.push_back(almdistribution(p.almleaves,{ lspot0-halfdelta, lspot0+halfdelta},p.delta));

        if(T>1)
        {
            m.push_back(mmcdistribution({{0.5,0.5,0},{0,0.5,0.5}}));
            d.push_back(almdistribution(p.almleaves,{  lspot0-p.delta, lspot0, lspot0+p.delta }, p.delta));
        }
        if(T>2)
        {
            m.push_back(mmcdistribution({{0.5,0.5,0,0},{0,0.5,0.5,0}
                                         ,{0,0,0.5,0.5}}));

            d.push_back(almdistribution(p.almleaves,{ lspot0-3*halfdelta, lspot0-halfdelta,
                                              lspot0+halfdelta, lspot0+3*halfdelta }, p.delta ));
        }
    }

    using hmcpd_t = hmcprocessdistribution<mmcdistribution,almdistribution>;

    hmcpd_t xipd(lspot0,m,d);

    using dxieta_t = dxietaprocessdist<hmcpd_t,etapd_t>;

    dxieta_t xieta(xipd,etapd);

    using zeta_t = hmczeta<pair<double, vector<double>>,zskm>;

    zskproblem pr(p.lambda,0.05);
    pr.fdebug = true;

    static_assert(std::is_same<typename dxieta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);

    static_assert(std::is_same<typename zeta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);
    desolution<zskproblem,dxieta_t,zeta_t,
cplex<realvar>
//csvlpsolver<realvar>
            > sol(pr,xieta);

    vectors<double> aves;
    vectors<unsigned int> cnts;

    sol.x()->print(sys::log());

    cout << "obj value=" << sol.obj() << endl;
    sol.x()->stats(cnts, aves);

#ifndef RISKNEUTRAL
    if(p.equivalent)
    {
        cout << "Solving an equivalent " << endl;

        mpmcvarequivalent<zskproblem> e(pr);

        desolution< mpmcvarequivalent<zskproblem>,dxieta_t,zeta_t,cplex<realvar>> esol(e,xieta);
        esol.x()->print(sys::log());
    }
#endif
    vector<double> sddpsol;
    if(p.sddp)
    {
        vector<unsigned int> dims(T,zskm::zetasize());

        msddpsolution<zskproblem,dxieta_t,zeta_t,cplex<realvar>> sx(pr,xieta,dims);
        double ratio = fabs(sx.obj().lb()/sol.obj());
        if(ratio < 0.98 && ratio > 1.02)
        {
            std::cerr << "dtest: opt="
                 << sol.obj() << " expected, " <<
                 sx.obj().lb() << "<" << sx.obj().ubm() <<
                 " achieved." << std::endl;
            throw;
        }
        sddpsol = *(sx.x());
    }

    vector<vector<std::string>> n;
    pr.varnames(n);

    for(unsigned int i=0; i< n.size(); i++)
    {
        for(unsigned int j=0; j< n[i].size(); j++)
        {
            cout << n[i][j] << "=" << aves[i][j] << " (" << cnts[i][j] << ")" ;
            if(i==0)
                cout << " sddp=" << sddpsol[j];
            cout << endl;
        }
    }

    double x = X0;
    double y = Y0;


    for(unsigned int i=0; i<T; i++)
    {
        double nx = xregconst+xregx*x + xregy*y;
        double ny = yregconst+yregx*x + yregy*y;
        cout << "(x,y)(" << i+1 << ")" << nx << "," << ny << "," << endl;
        x = nx;
        y = ny;
    }



//    desolution<zskproblem,dxieta_t,zeta_t,csvlpsolver<realvar>> solc(p,xieta);

//    vector<unsigned int> zd(T,zskm::zetasize);
//    msddpsolution<zskproblem,dxieta_t,zeta_t,O> sx(p,xieta,zd);
    T = saveT;

}


int main(int, char **)
{
//    int res=mcheck(nullptr);
//if(res)
//        cerr << "error " << res << " starting mcheck" << endl;

    try
    {

        std::ofstream logf("zsk.log");
        sys::setlog(logf);

        sys::seed(0);
    //    using O=csvlpsolver<realvar>;
       using O=cplex<realvar>;


        dtestparams p;

        p.T=3;
        p.trivialm=false;
        p.almleaves = 1;
        p.etaleaves = 1;
        p.lambda=0.5;
        p.delta = 0.2;
        p.sddp = false;
        dtest(p);
    }
    catch (mspp::exception& e)
    {
        cerr << "program throwed an exception: " << endl;
        cerr << e.msg() << " (#=" << e.erno() << ")" << endl;
        return 1;
    }

    return 0;
};

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
