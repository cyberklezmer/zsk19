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

extern	const double K[]



                ;

extern	const double maxkappa	;

extern	const double r[]

                ;


// end of paste


const unsigned int T = 2;
const unsigned int kappa=1;
const unsigned int nfprices = 3;

class almdistribution:
        public fdistribution<double, unsigned int>
{
public:
    using I_t = double;
    using C_t = unsigned int ;
    almdistribution(unsigned int n, const vector<double>& sc) :
         fn(n), fsc(sc) {}
private:
    virtual atom<double> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c < fsc.size());
        return { fsc[c] - 0.5 +  (double) (2*i+1) / (double) (2*fn), 1.0/ (double) fn };
    }

    virtual unsigned int natoms_is(const unsigned int& ) const
    { return fn; }

    unsigned int fn;
    vector<double> fsc;
};



probability inline Phi(double x)
{
    boost::math::normal norm;
    return boost::math::cdf(norm,x);
}

double inline bsf(double spot, double strike, double tau)
{
    double x = strike-spot;
    double vol = bsc+bsl*x+bsq*x*x;
    double d1 = (log(spot/strike)+(rf+vol*vol/2.0)*tau)/vol/sqrt(tau);
    double d2 = d1-vol*sqrt(tau);
    double res = spot*Phi(d1)-strike*exp(-rf*tau)*Phi(d2);
    return res;
}

class zskm: public mapping<pair<double,vector<double>>,vector<double>>
{
public:

  enum { espot, exeps, eyeps, efirstfuture, ef1=efirstfuture, ef2, ef3, efirstoption };

  static constexpr unsigned int zetasize = efirstfuture + nfprices + kappa * T;

  static double P(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize);
      return zeta[espot];
  }

  static double Q(const vector<double> zeta, unsigned int ttm)
  {
      assert(zeta.size()==zetasize);
      assert(efirstfuture+ttm-1<zeta.size());
      return zeta[efirstfuture+ttm-1];
  }

  static double Xeps(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize);
      return zeta[exeps];
  }

  static double Yeps(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize);
      return zeta[eyeps];
  }


  static double B(const vector<double> zeta,unsigned int ttm,
                                                  unsigned int i)
  {
      assert(zeta.size()==zetasize);
      assert(ttm>0);
      assert(efirstoption+kappa*(ttm-1)+i<zeta.size());
      return zeta[efirstoption+kappa*(ttm-1)+i];
  }


  virtual vector<double> operator() (const pair<double,vector<double>>& a) const
  {
    double p=exp(a.first);
    vector<double> r({ p,a.second[0],a.second[1],
                    p*exp( a.second[2]), p*exp( 2* a.second[3]),
                    p*exp( 3* a.second[4])});
    for(unsigned int tau=1; tau<=T; tau++)
    {
        for(unsigned int i=0; i<kappa; i++)
        {
            double b = bsf(p,K[i],tau);
            r.push_back(b);
        }
    }
    assert(r.size()==zetasize);
    return r;
  }
};

using dc=expectation; //mpmcvar;

class zskproblem: public msproblem<dc, linearfunction,
        linearmsconstraint,vector<double>,realvar,lastx>
{
public:
    /// concerns all the stages except for the last one
    enum vars {zt, xval, yval, et, firstpositive = et, firstft,
                             nfixed = firstft };


    static unsigned int nft(unsigned int k)
    {
        assert(k<=::T);
        return ::T-k;
    }
    static unsigned int nphit(unsigned int k)
    {
        assert(k<=::T);
        return (::T-k) * kappa;
    }
    static unsigned int phirelpos(unsigned int k,
                        unsigned int tau, unsigned int strike)
    {
        assert(k<tau);
        return (tau-k-1)*kappa + strike;
    }
    static unsigned int nvars(unsigned int k)
    {
        assert(k<=::T);
        return nfixed + nft(k)+nphit(k);
    }

    static double gf(unsigned int k)
    {
        assert(k);
        assert(k<=::T);
        return r[k-1];
    }

    unsigned int findex(unsigned int k,unsigned int tau) const
    {
        assert(tau <= this->T());
        assert(k <= this->T());
        assert(tau);
        assert(tau>k);
        return firstft + (tau-1-k);
    }
    unsigned int phiindex(unsigned int k,unsigned int tau, unsigned int i)
      const
    {
        assert(tau <= this->T());
        assert(k <= this->T());
        assert(tau);
        assert(tau>k);
        assert(i<kappa);
        return firstft + nft(k) + phirelpos(k,tau,i);
    }
private:
    static std::vector<unsigned int> makeps()
    {
        std::vector<unsigned int> res;
        for(unsigned int i=0; i<=::T; i++)
            res.push_back(nvars(i));
        return res;
    }
public:
    zskproblem(double lambda, double alpha) :
        msproblem<dc, linearfunction,
                linearmsconstraint,vector<double>,realvar,lastx>
        (makeps(),
expectation()
//      dc(lambda,alpha)
         )
    {
    }

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        std::ostringstream s;
        switch(i)
        {
           case zt:
              s << "z";
              break;
           case xval:
             s << "X";
             break;
           case yval:
             s << "Y";
             break;
           case et:
             s << "e";
             break;
           default:
           {
              unsigned int off = i-firstft;
              if(off<nft(stage))
                 s << "f" << stage + off + 1;
              else
              {
                 unsigned int offf = off - nft(stage);
                 unsigned int tau = offf / kappa + stage + 1;
                 unsigned int i = offf % kappa + 1;
                 s << "phi" << tau << "(" << i << ")";
              }
            }
        }
        return s.str();
    }


    virtual void f_is(unsigned int k,
                      const vector<double>&,
                      linearfunction& f
                      ) const
    {
        if(k==0)
            f.setc(zt,-1);
        else
            f.setc(zt,-pow(varrho,k));
    }

    virtual  void x_is(
            unsigned int k,
            const vector<double>& zeta,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g
            ) const
    {
cout << "P=" << zskm::P(zeta)
     << " Xeps=" << zskm::Xeps(zeta)
     << " Yeps=" << zskm::Yeps(zeta) ;
if(!k)
    cout << " F=" << zskm::Q(zeta,1)
         << " B=" << zskm::B(zeta,1,0);
cout << endl;

        unsigned int T = this->T();

        unsigned int i=0;
        for(; i<firstpositive; i++)
            xs[i].setlimits();
        for(; i<xs.size(); i++)
            xs[i].setlimits(0);

        unsigned int ls = k ? this->xdim(k)+this->xdim(k-1) : this->xdim(k);
        unsigned int toff = k ? this->xdim(k-1) : 0;

        if(k==T)
            xs[et].setlimits(0,0);

        // zt

        vector<double> zl(ls,0.0);
        double zr = 0.0;

        zl[toff+zt]=-1;

        double P=zskm::P(zeta);
        if(k==0)
            zl[toff+et]=-P;
        else
        {
//enum vars {zt, xval, yval, et, firstpositive = et, firstft, nfixed = firstft };
            zl[toff+xval]=1;

            zl[toff+et]=-P;
            zl[et] = P;
            zr-= gf(k)*P;
            zl[findex(k-1,k)]=P;
            for(unsigned int i=0; i<kappa; i++)
                zl[phiindex(k-1,k,i)]=P;
            zl[toff+yval] = -P;

            for(unsigned int i=0; i<kappa; i++)
                zl[phiindex(k-1,k,i)] = -min(P,K[i]);
        }

        double df = varrho;
        for(unsigned int tau=k+1; tau<=T; tau++)
        {
            double Q = zskm::Q(zeta,tau-k);
            zl[toff+findex(k,tau)] = -df*Q;
            if(k)
                zl[findex(k-1,tau)] = df*Q;
            df *= varrho;
        }
        for(unsigned int tau=k+1; tau<=T; tau++)
            for(unsigned int i=0; i<kappa; i++)
            {
                double B = zskm::B(zeta,tau-k,i);
                zl[toff+phiindex(k,tau,i)] = -B;
                if(k)
                    zl[phiindex(k-1,tau,i)] = B;
            }

        g.add(linearmsconstraint(zl,constraint::eq, zr));

        // x,y

        vector<double> xl(ls,0.0);
        double xr = 0.0;
        vector<double> yl(ls,0.0);
        double yr = 0.0;

        if(k==0)
        {
            xl[xval] = 1;
            xr = X0;
            yl[yval] = 1;
            yr = Y0;
        }
        else
        {
            xl[toff+xval] = 1;
            xl[xval]=-xregx;
            xl[yval]=-xregy;
            xr = xregconst+zskm::Xeps(zeta);


            yl[toff+yval] = 1;
            yl[xval]=-yregx;
            yl[yval]=-yregy;
            yr = yregconst+zskm::Yeps(zeta);;
        }

        g.add(linearmsconstraint(xl,constraint::eq, xr));
        g.add(linearmsconstraint(yl,constraint::eq, yr));

    }
    double minf_is(unsigned int) const
    {
        return -1e7;
    }
    double maxf_is(unsigned int) const
    {
        return 1e7;
    }
};




#endif // ZSK_H
