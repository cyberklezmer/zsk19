#ifndef ZSK_H
#define ZSK_H

#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/msproblem.h"
#include "mspp/hmcapprox.h"

using namespace mspp;


// paste from params.ods

extern	unsigned int T	;
extern	double varrho	;
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

extern	const double L[]



                ;




extern	const double r[]

                ;




const unsigned int nfprices=	3	;
const unsigned int kappa=	2	;


// end of paste



class almdistribution:
        public fdistribution<double, unsigned int>
{
public:
    using I_t = double;
    using C_t = unsigned int ;
    almdistribution(unsigned int n, const vector<double>& sc, double delta) :
         fn(n), fsc(sc), fdelta(delta) {}
private:
    virtual atom<double> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c < fsc.size());
//cout << "c = " << c << " D=" << fsc[c]- 0.5 +  (double) (2*i+1) / (double) (2*fn), 1.0/ (double) fn;
//cout << endl;
        return { fsc[c] - fdelta * (0.5 +  (double) (2*i+1) / (double) (2*fn)), 1.0/ (double) fn };
    }

    virtual unsigned int natoms_is(const unsigned int& ) const
    { return fn; }

    unsigned int fn;
    vector<double> fsc;
    double fdelta;
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
//cout << "bsf(" << spot << "," << strike << ","
//       << vol << "," << tau << ")=" << res << endl;

    return res;
}

double inline csf(double spot, double strike, double tau)
{
    double pv = strike * exp(-rf*tau);
    double res = pv - spot + bsf(spot,strike,tau);
    return res;
}


inline vectors<double> etasqV(unsigned int k, double omega)
{
    assert(k<=T);
    assert(k>0);

    vectors<double>r(7-k);
    for(unsigned int i=0; i<r.size(); i++)
        r[i]=vector<double>(7-k,0.0);

    r[0][0] = sqrt(1-omega)*epsvsq[0][0];
    r[0][1] = sqrt(1-omega)*epsvsq[0][1];
    r[1][0] = sqrt(1-omega)*epsvsq[1][0];
    r[1][1] = sqrt(1-omega)*epsvsq[1][1];

    r[2][2] = sqrt(omega)*epsvsq[0][0];
    r[2][3] = sqrt(omega)*epsvsq[0][1];
    r[3][2] = sqrt(omega)*epsvsq[1][0];
    r[3][3] = sqrt(omega)*epsvsq[1][1];

    for(unsigned int i=4; i<7-k; i++)
        r[i][i] = varsigma;
    return r;
}


class zskm: public mapping<pair<double,vector<double>>,vector<double>>
{
public:

  enum { espot,  exg, eyg, exnextf, eynextf, efirstfuture,
         ef1=efirstfuture, ef2, ef3, efirstoption };

  static unsigned int zetasize()
  {
      return efirstfuture + nfprices + 2 * kappa * T;
  }

  static double P(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize());
      return zeta[espot];
  }

  static double Q(const vector<double> zeta, unsigned int ttm)
  {
      assert(zeta.size()==zetasize());
      assert(efirstfuture+ttm-1<zeta.size());
//cout << "ttm=" << ttm << " Q=" << zeta[efirstfuture+ttm-1] << endl;
      return zeta[efirstfuture+ttm-1];
  }

  static double Xnextf(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize());
      return zeta[exnextf];
  }

    static double Ynextf(const vector<double> zeta)
  {
      assert(zeta.size()==zetasize());
      return zeta[eynextf];
  }

    static double Xg(const vector<double> zeta)
    {
        assert(zeta.size()==zetasize());
        return zeta[exg];
    }

    static double Yg(const vector<double> zeta)
    {
        assert(zeta.size()==zetasize());
        return zeta[eyg];
    }


  static double B(const vector<double> zeta,unsigned int ttm,
                                                  unsigned int i)
  {
      assert(zeta.size()==zetasize());
      assert(ttm>0);
      assert(efirstoption+kappa*(ttm-1)+i<zeta.size());
      return zeta[efirstoption+kappa*(ttm-1)+i];
  }

  static double C(const vector<double> zeta,unsigned int ttm,
                                                  unsigned int i)
  {
      assert(zeta.size()==zetasize());
      assert(ttm>0);
      assert(efirstoption+kappa*T+kappa*(ttm-1)+i<zeta.size());
      return zeta[efirstoption+kappa*T+kappa*(ttm-1)+i];
  }


  virtual vector<double> operator() (const pair<double,vector<double>>& a) const
  {
    double p=exp(a.first);
//cout << "a.first = " << a.first << endl;
//cout << "Px = " << p << endl;
    vector<double> r({ p,a.second[0],a.second[1],a.second[2],a.second[3] });
    unsigned int is=4;
    for(unsigned int i=0; i<nfprices; i++)
    {
        if(is+i < a.second.size())
            r.push_back(p*exp( (i+1)*a.second[is+i]));
        else
            r.push_back(0);
    }
    for(unsigned int tau=1; tau<=T; tau++)
    {
        for(unsigned int i=0; i<kappa; i++)
        {
            double b = bsf(p,K[i],tau);
            r.push_back(b);
        }
    }
    for(unsigned int tau=1; tau<=T; tau++)
    {
        for(unsigned int i=0; i<kappa; i++)
        {
            double c = csf(p,L[i],tau);
            r.push_back(c);
        }
    }
    assert(r.size()==zetasize());
    return r;
  }
};


template <typename O, bool nox = true, bool nod = false>
class zskproblem: public msproblem<O, linearfunction,
        linearmsconstraint,vector<double>,realvar,lastx>
{
public:
    double F() const
    {  // tbd do params, pokud se bude pouivat
        double res = T * 32639 - sqrt(T) * 8811 * 2;
        if(!nox)
            res += T * 227286 - sqrt(T) * 26339*2;
    }

    bool fdebug;
    static constexpr double instnbound = 1e12;
    /// concerns all the stages except for the last one
    enum vars {zt, xval, yval, det, et, vt, wt=vt,  fx, nlast=fx, fy, firstft,
                             nfixed = firstft };


    static unsigned int nft(unsigned int k)
    {
        assert(k<=::T);
        return ::T-k;
    }

    static unsigned int ndft(unsigned int k)
    {
        assert(k<=::T);
        return k ? nft(k) : 0;
    }

    static unsigned int nphit(unsigned int k)
    {
        assert(k<=::T);
        return (::T-k) * kappa;
    }

    static unsigned int npsit(unsigned int k)
    {
        return nphit(k);
    }


    static unsigned int phirelpos(unsigned int k,
                        unsigned int tau, unsigned int strike)
    {
        assert(k<tau);
        return (tau-k-1)*kappa + strike;
    }

    static unsigned int psirelpos(unsigned int k,
                        unsigned int tau, unsigned int strike)
    {
        return phirelpos(k,tau,strike);
    }


    static unsigned int nvars(unsigned int k)
    {
        assert(k<=::T);
        if(k==T)
            return nlast;
        else
            return nfixed + nft(k) + ndft(k) + nphit(k) + npsit(k);
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

    unsigned int dfindex(unsigned int k,unsigned int tau) const
    {
        assert(tau <= this->T());
        assert(k <= this->T());
        assert(tau);
        assert(tau>k);
        assert(k);
        return firstft + nft(k)+(tau-1-k);
    }

    unsigned int phiindex(unsigned int k,unsigned int tau, unsigned int i)
      const
    {
        assert(tau <= this->T());
        assert(k <= this->T());
        assert(tau);
        assert(tau>k);
        assert(i<kappa);
        return firstft + nft(k) + ndft(k) + phirelpos(k,tau,i);
    }
    unsigned int psiindex(unsigned int k,unsigned int tau, unsigned int i)
      const
    {
        assert(tau <= this->T());
        assert(k <= this->T());
        assert(tau);
        assert(tau>k);
        assert(i<kappa);
        return firstft + nft(k) + ndft(k) + nphit(k) + psirelpos(k,tau,i);
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
        msproblem<O, linearfunction,
                linearmsconstraint,vector<double>,realvar,lastx>
        (makeps(),
#ifdef RISKNEUTRAL
         expectation()
#else
      O(lambda,alpha)
#endif
         ), fdebug(false)
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
           case fx:
             s << "xf";
             break;
           case yval:
             s << "Y";
             break;
           case fy:
              s << "yf";
              break;
           case et:
              s << "e";
              break;
           case det:
             s << "ds";
             break;
           case wt:
              s << (stage==this->T() ?  "w" : "v");
              break;
          break;
           default:
           {
              unsigned int off = i-firstft;
              if(off<nft(stage))
                 s << "f" << stage + off + 1;
              else if(off<nft(stage)+ndft(stage))
                  s << "df" << stage + off - nft(stage) + 1;
              else if(off<nft(stage)+ndft(stage)+nphit(stage))
              {
                 unsigned int offf = off - nft(stage)-ndft(stage);
                 unsigned int tau = offf / kappa + stage + 1;
                 unsigned int i = offf % kappa + 1;
                 s << "phi" << tau << "(" << i << ")";
              }
              else
              {
                 unsigned int offf = off - nft(stage)-ndft(stage)-nphit(stage);
                 unsigned int tau = offf / kappa + stage + 1;
                 unsigned int i = offf % kappa + 1;
                 s << "psi" << tau << "(" << i << ")";
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
        {
            f.setc(zt,-pow(varrho,k));
            if(k==this->T())
                f.setc(wt,pow(varrho,k));
        }
    }

    virtual void x_is(
            unsigned int k,
            const vector<double>& zeta,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g
            ) const
    {
        if(fdebug)
        {
            sys::log() << k;
            for(unsigned int i=0; i<zeta.size(); i++)
                sys::log() << "," << zeta[i];
            sys::log() << endl;
        }
/*                              zskm::P(zeta)
     << " Xeps=" << zskm::Xeps(zeta)
     << " Yeps=" << zskm::Yeps(zeta) ;
if(!k)
    cout << " F=" << zskm::Q(zeta,1)
         << " B=" << zskm::B(zeta,1,0);
cout << endl;*/

        unsigned int T = this->T();

        unsigned int i=0;
        xs[zt].setlimits();
        xs[xval].setlimits();
        xs[yval].setlimits();
        xs[det].setlimits();
        xs[et].setpositive();
        xs[vt].setlimits();
        for(; i<=det; i++)
            xs[i].setlimits();
        if(k==0)
        {
            xs[fx].setlimits(0,0);
            xs[fy].setlimits(0,0);
        }
        else if(k<T)
        {
            xs[fx].setlimits(zskm::Xnextf(zeta),zskm::Xnextf(zeta));
            xs[fy].setlimits(zskm::Ynextf(zeta),zskm::Ynextf(zeta));
        }
        xs[et].setlimits(0,instnbound);
        for(i=firstft; i<firstft+nft(k); i++)
            xs[i].setlimits(0,nod ? 0 : instnbound);
        for( ;i<firstft+nft(k)+ndft(k); i++)
            xs[i].setlimits(0,nod ? 0 : instnbound);
        for(; i<xs.size(); i++)
            xs[i].setlimits(0,nod ? 0 : instnbound);

        unsigned int ls = k ? this->xdim(k)+this->xdim(k-1) : this->xdim(k);
        unsigned int toff = k ? this->xdim(k-1) : 0;

        if(k==T)
            xs[et].setlimits(0,0);

        if(k==0)
            xs[vt].setlimits(0,0);
        // zt

        vector<double> zl(ls,0.0);
        double zr = 0.0;

        vector<double> del(ls,0.0);
        double der = 0.0;

        del[toff+det]=1;
        if(k==0)
            del[toff+et] = -1;
        else
        {
            del[toff+et]=-1;
            del[et] = 1;
            der-= gf(k);
            del[findex(k-1,k)]=1;
            for(unsigned int i=0; i<kappa; i++)
                del[phiindex(k-1,k,i)]=1;
            for(unsigned int i=0; i<kappa; i++)
                del[psiindex(k-1,k,i)]=-1;
            del[toff+yval] = -1;
        }

        double P=zskm::P(zeta);
        zl[toff+zt]=-1;

        if(k==0)
            zl[toff+et]=-P;
        else
        {
            zl[toff+det]=-P;
            if(!nox)
                zl[toff+xval]=1000;
            for(unsigned int i=0; i<kappa; i++)
                zl[phiindex(k-1,k,i)] += -min(P,K[i]);
            for(unsigned int i=0; i<kappa; i++)
                zl[psiindex(k-1,k,i)] += max(P,L[i]);
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
        {
            for(unsigned int i=0; i<kappa; i++)
            {
                double B = zskm::B(zeta,tau-k,i);
                zl[toff+phiindex(k,tau,i)] = -B;
                if(k)
                    zl[phiindex(k-1,tau,i)] = B;
            }
            for(unsigned int i=0; i<kappa; i++)
            {
                double C = zskm::C(zeta,tau-k,i);
                zl[toff+psiindex(k,tau,i)] = -C;
                if(k)
                    zl[psiindex(k-1,tau,i)] = C;
            }
        }
        g.add(linearmsconstraint(del,constraint::eq, der));
        g.add(linearmsconstraint(zl,constraint::eq, zr));

        if(k)
            for(unsigned int tau=k+1; tau<=T; tau++)
            {
                vector<double> dfl(ls,0.0);
                dfl[toff+dfindex(k,tau)]=1;
                dfl[toff+findex(k,tau)]=-1;
                dfl[findex(k-1,tau)]=1;
                g.add(linearmsconstraint(dfl,constraint::eq, 0.0));
            }
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
            xl[fx]=-1;
            xr = xregconst+zskm::Xg(zeta);

            yl[toff+yval] = 1;
            yl[xval]=-yregx;
            yl[yval]=-yregy;
            yl[fy]=-1;
            yr = yregconst+zskm::Yg(zeta);
        }

        g.add(linearmsconstraint(xl,constraint::eq, xr));
        g.add(linearmsconstraint(yl,constraint::eq, yr));

        if(k>0)
        {
            if(k<T)
            {
                vector<double> vl(ls,0.0);
                vl[toff+vt] = -1;
                vl[vt] = 1.0 / varrho;
                vl[toff+zt] = 1;
                g.add(linearmsconstraint(vl,constraint::eq, 0));
            }
            else // k==T
            {
                if(::iota > 0)
                {
                    vector<double> ul(ls,0.0);
                    ul[toff+wt] = 1.0 / ::iota;
                    ul[vt] = 1.0 / varrho;
                    ul[toff+zt] = 1;
                    g.add(linearmsconstraint(ul,constraint::geq, F()));
                }
                xs[wt].setpositive();
//xs[wt].setlimits(P,P);
            }
        }

    }
    double minf_is() const
    {
        return  -1e15;
    }
    double maxf_is() const
    {
        return 1e15;
    }
};




#endif // ZSK_H
