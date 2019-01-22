#ifndef ZSK_H
#define ZSK_H

#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/msproblem.h"
#include "mspp/hmcapprox.h"

using namespace mspp;


// paste from params.cpp

extern	const double varrho	;
extern	const double iota	;

extern	const vector<vector<double>> epsvsq
                ;
extern	const double X0	;
extern	const double Y0	;

extern	const double xregconst	;
extern	const double yregconst	;
extern	const double xregx	;
extern	const double yregx	;
extern	const double xregy	;
extern	const double yregy	;

extern	const double spot0	;
extern	const double sigma 	;



extern	const double qcoef	;
extern	const double varsigma	;
extern	const vector<double> v0

                ;

extern	const double rf	;
extern	const double bsc	;
extern	const double bsl	;
extern	const double bsq	;

extern	const double strikes[]



                ;


// end of paste

const unsigned int T = 2;
const unsigned int nstrikes=1;

extern const double maxstrikes;

probability Phi(double x)
{
    boost::math::normal norm;
    return boost::math::cdf(norm,x);
}

double bsf(double spot, double strike, double tau)
{
    double x = strike-spot;
    double vol = bsc+bsl*x+bsq*x*x;
    double d1 = (log(spot/strike)+(rf+vol*vol/2.0)*tau)/vol/sqrt(tau);
    double d2 = d1-vol*sqrt(tau);
    double res = spot*Phi(d1)-strike*exp(-rf*tau)*Phi(d2);
}

class zskproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,vector<double>,realvar,lastx>
{
    static std::vector<unsigned int> makeps(unsigned int t)
    {
        std::vector<unsigned int> res = {1};
        for(unsigned int i=1; i<=t; i++)
            res.push_back(2);
        return res;
    }
public:
    zskproblem(double lambda, double alpha, unsigned int T) :
        msproblem<mpmcvar, linearfunction,
                linearmsconstraint,vector<double>,realvar,lastx>
        (makeps(T),mpmcvar(lambda,alpha))
    {}


    virtual void f_is(unsigned int k,
                      const vector<double>& zeta,
                      linearfunction& f
                      ) const
    {
        if(k)
            f.setc(0,zeta[0]);
        else
            f.setc(0,1);
    }

    virtual  void x_is(
            unsigned int k,
            const vector<double>& zeta,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g
            ) const
    {
        if(k)
            xs[1].setlimits(0,1);
        else
            xs[0].setlimits(0,1);
        if(!k)
            return;
        if(k==1)
            g.add(linearmsconstraint({1.0,1.0,-1.0},constraint::eq, 0.0));
        else
            g.add(linearmsconstraint({0.0,1.0,1.0,-1.0},constraint::eq, 0.0));
        if(k==this->T())
            xs[1].setlimits(1,1);
    }
    double minf_is(unsigned int) const
    {
        return -1000;
    }
    double maxf_is(unsigned int) const
    {
        return 1000;
    }
};


class zskm: public mapping<pair<double,vector<double>>,vector<double>>
{
public:

  enum { spot, xeps, yeps, firstfuture, f1=firstfuture, f2, f3, firstoption };

  virtual vector<double> operator() (const pair<double,vector<double>>& a) const
  {
    double p=exp(a.first);
    vector<double> r({ p,a.second[0],a.second[1],
                    p*exp( a.second[2]), p*exp( 2* a.second[3]),
                    p*exp( 3* a.second[4])});
    for(unsigned int i=0; i<nstrikes; i++)
    {
        for(unsigned int tau=1; tau<=T; tau++)
        {
            double b = bsf(p,strikes[i],tau);
            r.push_back(b);
        }
    }
  }
};



#endif // ZSK_H
