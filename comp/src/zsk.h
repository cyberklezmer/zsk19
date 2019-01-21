#ifndef ZSK_H
#define ZSK_H

#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/msproblem.h"
#include "mspp/hmcapprox.h"

using namespace mspp;

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


extern	const vector<vector<double>> etavsq
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




extern	const double qcoef	;
extern	const double varsigma	;
extern	const vector<double> f0

                ;

extern	const double rf	;
extern	const double bsc	;
extern	const double bsl	;
extern	const double bsq	;




#endif // ZSK_H