//#define RISKNEUTRAL

#include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
#include "mspp/msproblem.h"
#include "mspp/mspptest.h"
#include "mspp/cplex.h"
#include "zsk.h"
#include <mcheck.h>

using namespace mspp;

#ifndef RISKNEUTRAL


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

template <typename O, typename S=cplex<realvar>>
void dtest(const dtestparams& p)
{
    constexpr bool mp = std::is_same<O,mpmcvar>::value;
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
        vector<double> mean(7-k,0.0);
        for(unsigned int i=4; i<7-k; i++)
            mean[i]=qcoef;
        dsts.push_back(etad_t(mean,etasqV(k,0 /*tbd*/),stdd));
    }

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = processdistribution<etad_t,noxi<vector<double>>>;

    etapd_t etapd(dirac_t({0,0,0,0,v0[0],v0[1],v0[2]}),dsts);

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

    zskproblem<O> pr(p.lambda,0.05);
    pr.fdebug = true;

    static_assert(std::is_same<typename dxieta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);

    static_assert(std::is_same<typename zeta_t::X_t,pair<unsigned int, pair<double, vector<double>>>>::value);
    desolution<zskproblem<mpmcvar>,dxieta_t,zeta_t,S>* sp = 0;
    vectors<double> aves;
    vectors<unsigned int> cnts;
    if constexpr(mp)
    {
        sp = new desolution<zskproblem<mpmcvar>,dxieta_t,zeta_t,S>(pr,xieta);


        sp->x()->print(sys::log());

        cout << "obj value=" << sp->obj() << endl;
        sp->x()->stats(cnts, aves);

    #ifndef RISKNEUTRAL
        if(p.equivalent)
        {
            cout << "Solving an equivalent " << endl;

            mpmcvarequivalent<zskproblem<mpmcvar>> e(pr);

            desolution<mpmcvarequivalent<zskproblem<mpmcvar>>,dxieta_t,zeta_t,cplex<realvar>> esol(e,xieta);
            esol.x()->print(sys::log());
        }
    #endif
    }
    vectors<double> sddpsol;
    vectors<double> sddperrs;
    if(p.sddp)
    {
        vector<unsigned int> dims(T,zskm::zetasize());

        msddpsolution<zskproblem<O>,dxieta_t,zeta_t,cplex<realvar>> sx(pr,xieta,dims);
        std::cout << "dtest:" <<
             sx.obj().lb() << "<" << sx.obj().ubm() <<
             " achieved." << std::endl;

        if constexpr(mp)
        {
            double ratio = fabs(sx.obj().lb()/sp->obj());
            if(ratio < 0.98 && ratio > 1.02)
            {
                std::cerr << "dtest: opt="
                     << sp->obj() << " expected, " <<
                     sx.obj().lb() << "<" << sx.obj().ubm() <<
                     " achieved." << std::endl;
                throw;
            }
        }
        sddpsol = sx.x()->means;
        sddperrs = sx.x()->sterrs;
    }

    vector<vector<std::string>> n;
    pr.varnames(n);

    for(unsigned int i=0; i< n.size(); i++)
    {
        for(unsigned int j=0; j< n[i].size(); j++)
        {
            cout << n[i][j] << ": ";
            if(mp)
                 cout << aves[i][j] << " (" << cnts[i][j] << ")" ;
            if(p.sddp)
                cout << " sddp=" << sddpsol[i][j] << " (" << sddperrs[i][j] << ") ";
            if(mp && p.sddp)
                cout << "e=" << fabs(aves[i][j] - sddpsol[i][j]);
            cout << endl;
        }
    }

/*    double x = X0;
    double y = Y0;

    for(unsigned int i=0; i<T; i++)
    {
        double nx = xregconst+xregx*x + xregy*y;
        double ny = yregconst+yregx*x + yregy*y;
        cout << "(x,y)(" << i+1 << ")" << nx << "," << ny << "," << endl;
        x = nx;
        y = ny;
    }*/



//    desolution<zskproblem,dxieta_t,zeta_t,csvlpsolver<realvar>> solc(p,xieta);

//    vector<unsigned int> zd(T,zskm::zetasize);
//    msddpsolution<zskproblem,dxieta_t,zeta_t,O> sx(p,xieta,zd);
    T = saveT;

}

struct compparams
{
    string id = string("nan");
    string comment = string("");
    unsigned int T = ::T;
    double varrho = ::varrho;
    unsigned int patoms = 3;
    double lambda = 0.5;
    double omega = 0;
};

class expmapping : public mapping<double,double>
{
public:
   virtual double operator()(const double& x) const
   {
        return exp(x);
   }
};

using cha_t = chmcapproximation<arnormalprocessdistribution,onedcovering,true>;
using dha_t=dhmcapproximation<arnormalprocessdistribution,onedcovering,1,true,expmapping>;

template <typename O,typename ha_t, bool nox, bool nod>
void cont(const compparams& pars, std::ofstream& res,
          bool headers = false  )
{
    constexpr bool contin = std::is_same<ha_t,cha_t>::value;
    unsigned int saveT = T;
    T=pars.T;
    double savevarrho = varrho;
    varrho = pars.varrho;
    using stdd_t=stdnormaldistribution;

    ostringstream rdn;
    rdn << pars.id << ".csv";
    std::ofstream rdet(rdn.str());

    using etad_t = meanvardistribution<stdd_t>;

    vector<etad_t> dsts;
    for(unsigned int k=1; k<=T; k++)
    {
        vector<double> mean(7-k,0.0);
        for(unsigned int i=4; i<7-k; i++)
            mean[i]=qcoef;
        dsts.push_back(etad_t(mean,etasqV(k,pars.omega)));
    }

    using dirac_t = diracdistribution<vector<double>>;
    using etapd_t = processdistribution<etad_t,noxi<vector<double>>>;

    etapd_t eta(dirac_t({0,0,0,0,v0[0],v0[1],v0[2]}),dsts);

    double xim = -sigma*sigma/2.0;
    double xisd = sigma;

    arnormalprocessdistribution xipd(log(spot0),xim,xisd,1.0,T);

    vector<unsigned int> nx;
    unsigned int na = pars.patoms;
    for(unsigned int i=1;i<T; i++)
    {
        nx.push_back(na);
        na = static_cast<unsigned int>(round(static_cast<double>(na) * sqrt(2.0)));
    }
    nx.push_back(na);

    vector<double> trends;
    if constexpr(contin)
    {
        trends = vector<double>(T,xim);
    }
    ha_t ha(xipd,nx,trends);

//cout << "prs" << endl;
//    cout << ha.d(1).first().m() << endl;
//cout << endl;
    if constexpr(std::is_same<ha_t,dha_t>::value)
    {
        printtreed(cout,ha);
    }

/*
    for(unsigned int i=1; i<=T; i++)
    {
        const mmcdistribution& m = ha.md(i);
        vectors<probability> p = m.m();
        cout << "stage=" << i << endl;
        for(unsigned int j=0; j < p.size(); j++)
        {
            for(unsigned int k=0; k < p[j].size(); k++)
                cout << p[j][k] << " ";
            cout << endl;
        }
    }*/

    using xieta_t = xietaprocessdist<ha_t,etapd_t>;
    xieta_t xieta(ha,eta);

    vector<unsigned int> dims(T,zskm::zetasize());

    using zeta_t = hmczeta<pair<double, vector<double>>,zskm>;

    zskproblem<O,nox,nod> p(pars.lambda,0.05);

    vector<vector<std::string>> n;
    p.varnames(n);

    for(unsigned int i=0; i< n.size();i++)
    {
        for(unsigned int j=0; j< n[i].size(); j++)
        {
            rdet << n[i][j] << ",";
        }
    }
    rdet << endl;


    if(headers)
    {
        res << "id,lambda,omega";
        for(unsigned int i=0; i< n.size();i++)
        {
            for(unsigned int j=0; j< n[i].size(); j++)
            {
                res << "," << n[i][j];
            }
        }
        res << ",time,lb,ubm,ubb,states..." << endl;
     }

    msddpsolution<zskproblem<O,nox,nod>,xieta_t,zeta_t,cplex<realvar>> sx(p,xieta,dims);

    res << pars.id << "," << pars.lambda  << "," << pars.omega;

    for(unsigned int i=0; i< n.size();i++)
    {       
        for(unsigned int j=0; j< n[i].size(); j++)
        {
            res << "," << sx.x()->means[i][j];
        }
    }

    res << "," << sx.obj().time()
        << "," << sx.obj().lb() << "," << sx.obj().ubm()
        << "," << sx.obj().ubb();

    for(unsigned int i=0; i <T; i++)
        res << "," << nx[i];
    res << "," << pars.comment;
    res << endl;

    bool any = true;
    for(unsigned int k=0; any ; k++)
    {
        any = false;
        for(unsigned int i=0; i< n.size();i++)
        {

            for(unsigned int j=0; j< n[i].size(); j++)
                if(k < sx.x()->x[i].size())
                {
                    rdet  << sx.x()->x[i][k][j] << ",";
                    any = true;
                }
                else
                    rdet << ",";
        }
        rdet << endl;
    }

    varrho = savevarrho;
    T = saveT;
}

#endif // RISKNEUTRAL

//continuous 5 stages
extern vector<double> x01;
extern vector<probability> aprs1;
extern vector<vector<double>> adata1;
extern vector<double> x02;
extern vector<probability> aprs2;
extern vector<vector<double>> adata2;
extern vector<double> x03;
extern vector<probability> aprs3;
extern vector<vector<double>> adata3;


void atest(const vector<double>& x0,
    const vector<probability>& aprs,
    const vector<vector<double>>& adata,
    unsigned int npp
)
{
    unsigned int saveT = T;
    T=1;
    vector<atom<vector<double>>> a;
    double ave = 0;
    for(unsigned int k=0; k<adata.size()/100; k++)
        for(unsigned int i=0; i<npp; i++)
        {
            vector<double> x = adata[100*k+i];
            probability p = aprs[k]/(double) npp;
            ave += p * x[0];
            a.push_back({x,p});
        }
    cout << "ave = " << ave << endl;
    ldistribution<vector<double>> d(a);

    using pd_t = fdprocessdistribution<ldistribution<vector<double>>, noxi<vector<double>>>;

    pd_t pd(x0,d,1);


#ifdef RISKNEUTRAL
    zskproblem<expectation> p(0.1,0.05);
    desolution<zskproblem<expectation>,pd_t,
            lastxi<vector<double>>,
//csvlpsolver<realvar>
            cplex<realvar>
            > x(p,pd);
#else
    zskproblem<mpmcvar> p(0.1,0.05);
    desolution<zskproblem<mpmcvar>,pd_t,
            lastxi<vector<double>>,
//csvlpsolver<realvar>
            cplex<realvar>
            > x(p,pd);
#endif



    x.x()->print(cout);
    cout << x.obj() << endl;
    T = saveT;
}

template <bool nod>
void mainanal( std::ofstream& res)
{
    compparams p;
    p.T = 3;
    p.patoms = 15;

    bool headers = true;
    for(double lambda = 0.1; lambda < 1.001; lambda+=0.45) //0.45)
        for(double omega = 0.0; omega < 1.001; omega+=0.5 ) //0.5)
        {
            ostringstream s;
            s << "lambda" << lambda << "omega"
              << omega << "n" << p.patoms;
            p.lambda = lambda;
            p.omega = omega;
            p.id = s.str();
            cont<nestedmcvar,dha_t,false,nod>(p, res, headers);
            headers = false;
            s << "nox";
            p.id = s.str();
            cont<nestedmcvar,dha_t,true,nod>(p, res, headers);
        }
}



void prelim( std::ofstream& res)
{
    compparams p;
    p.T = 3;
    p.patoms = 10;
    p.varrho = 0.96;

    bool headers = true;
    for(double lambda = 0.1; lambda < 1.001; lambda+=0.9) //0.45)
        for(double omega = 0.0; omega < 1.001; omega+=1 ) //0.5)
        {
            ostringstream s;
            s << "lambda" << lambda << "omega"
              << omega << "n" << p.patoms;
            p.lambda = lambda;
            p.omega = omega;
            p.id = s.str();
            cont<nestedmcvar,dha_t,false,false>(p, res, headers);
            headers = false;
            s << "nox";
            p.id = s.str();
            cont<nestedmcvar,dha_t,true,false>(p, res, headers);
        }
}


int main(int, char **)
{
//    int res=mcheck(nullptr);
//if(res)
//        cerr << "error " << res << " starting mcheck" << endl;

    std::ofstream res("results.csv");

    try
    {
        std::ofstream logf("zsk.log");
        sys::setlog(logf);

        sys::seed(0);
    //    using O=csvlpsolver<realvar>;
       using O=cplex<realvar>;
#ifndef RISKNEUTRAL
        if constexpr(0) // reproducing high upper bound
        {
            dtestparams p;

            p.T=1;
            p.trivialm=true;//false
            p.almleaves = 1;//2;
            p.etaleaves = 1;//2
            p.lambda=1;
            p.delta = 0.2;
            p.sddp = true;
            dtest<mpmcvar /*nestedmcvar*/,
//csvlpsolver<realvar>
                    cplex<realvar>
                    >(p);
        }

        if constexpr(0) // simple test
        {
            dtestparams p;

            p.T=1;
            p.trivialm=true;
            p.almleaves = 1;
            p.etaleaves = 1;
            p.lambda=0.1;
            p.delta = 0.2;
            p.sddp = true;
            dtest<mpmcvar/*,csvlpsolver<realvar>*/>(p);
        }

        compparams p;
        if constexpr(0)
            prelim(res);
        if constexpr(1)
            mainanal<true>(res);
        if  constexpr(0)
        {
            p.T = 1;
            p.patoms = 15;

            p.lambda = 0.1;
            p.omega = 0;
            cont<nestedmcvar,dha_t,false,false>(p, res, true);
         }
#endif // RISKNEUTRAL
        if  constexpr(0)
            atest(x02,aprs2,adata2,100);
    }
    catch (mspp::exception& e)
    {
        cerr << "program throwed an exception: " << endl;
        cerr << e.msg() << " (#=" << e.erno() << ")" << endl;
        res << "program throwed an exception: " << endl;
        res << e.msg() << " (#=" << e.erno() << ")" << endl;
        return 1;
    }

    return 0;
};



