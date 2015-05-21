// $Id: wilson_flow_w.cc,v 1.12 2011/12/11 17:23:22 cmcneile Exp cmcneile $
/*! \file
 *  \brief Code for Wilson flow
 *   
 *   A collection of routines to compute the
 *   Wilson flow.
 *
 *   Essentially the code implements appendix C
 *   of
 *
 *    Properties and uses of the Wilson flow in lattice QCD.
 *    Martin Luscher 
 *    Published in JHEP 1008 (2010) 071 , arXiv:1006.4518 
 *
 *    The code calls the existing stout smearing routines
 *    -- which already exist in chroma.
 *
 *   The Wilson flow can be used to determine the lattice spacing.
 *

 *    Adaptive step size code follows arXiv:1301.4388

 *
 * See below for additional information about the Wilson flow.
 * 
 * Continuous smearing of Wilson Loops.
 * Robert Lohmayer, Herbert Neuberger. arXiv:1110.3522. 
 * Published in PoS LATTICE2011 (2011) 249 
 *
 *
 *
 *
 */

#include "meas/glue/wilson_flow_w.h"
#include "meas/glue/mesfield.h"
#include "util/gauge/stout_utils.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

//using namespace Chroma;
namespace Chroma
{



  /**

   **/


  void measure_wilson_gauge(multi1d<LatticeColorMatrix> & u,
			    Real & gspace, Real & gtime,
			    int jomit)
  {

    multi1d<LatticeColorMatrix> field_st(10) ;
    LatticeColorMatrix tmp ; 

    mesField( field_st,u); 

    int offset = 0;

    gtime  = 0.0 ;
    gspace = 0.0 ;

    Double tr ; 

    //  for(int mu=0; mu < Nd-1; ++mu)

    for(int mu=0; mu < Nd; ++mu)
    {
      for(int nu=mu+1; nu < Nd; ++nu)
      {
	//	  tr = real(sum(trace(field_st[offset]))) ;
	//	  tr = imag(sum(trace(field_st[offset]))) ;
	tmp = field_st[offset] * field_st[offset] ;
	tr = real(sum(trace(tmp))) ; 

	// Real tt = 2.0 * tr ; 
	//	  cout << "DEBUG " << mu << " " << nu  << " " << tt << endl ;

	if (nu==jomit)
	{
	  gtime += 2.0*(tr);
	}
	else
	{
	  gspace += 2.0*(tr);
	}

	++offset  ; 
      }
    }

    gspace /= -2.0*Layout::vol() ;
    gtime /= -2.0*Layout::vol() ;

  }


  void wilson_flow_one_step(multi1d<LatticeColorMatrix> & u, Real rho)
  {
    int mu, dir;
    multi1d<LatticeColorMatrix> dest(Nd);
    multi1d<LatticeColorMatrix> next(Nd);


    // -------------------------------------

    multi1d<bool> smear_in_this_dirP(4) ;
    multi2d<Real> rho_a(4,4) ;
    multi2d<Real> rho_b1(4,4) ;
    multi2d<Real> rho_b2(4,4) ;
    multi2d<Real> rho_c(4,4) ;

    for (mu = 0; mu <= Nd-1; mu++)
    {
      smear_in_this_dirP(mu) = true ;
      for (dir = 0; dir <= Nd-1; dir++)
      {
	rho_a[mu][dir] = rho * 0.25 ;

	rho_b1[mu][dir] = rho * 8.0/9.0 ;
	rho_b2[mu][dir] = rho * 17.0/36.0 ;

	rho_c[mu][dir] = rho * 3.0/4.0 ;

      }
    }


    Stouting::smear_links(u, dest,smear_in_this_dirP, rho_a);

    LatticeColorMatrix  Q, QQ, C ;
    LatticeColorMatrix  Q2, QQ2  ;

    multi1d<LatticeColorMatrix> Q0(Nd);
    multi1d<LatticeColorMatrix> Q1(Nd);


    multi1d<LatticeComplex> f;   // routine will resize these


    for (mu = 0; mu <= Nd-1; mu++)
    {
      Stouting::getQsandCs(dest, Q1[mu],QQ,C,mu,smear_in_this_dirP,rho_b1) ;
      Stouting::getQsandCs(u   , Q0[mu],QQ,C,mu,smear_in_this_dirP,rho_b2) ;

      Q = Q1[mu] - Q0[mu] ;
      QQ = Q * Q ;
      Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
          
      // Assemble the stout links exp(iQ)U_{mu} 
      next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      

    }

    for (mu = 0; mu <= Nd-1; mu++)
    {
      u[mu]    =  next[mu] ;
      dest[mu] =  next[mu] ;      
    }


    for (mu = 0; mu <= Nd-1; mu++)
    {
      Stouting::getQsandCs(dest, Q2,QQ,C,mu,smear_in_this_dirP,rho_c) ;

      Q = Q2 - Q1[mu] + Q0[mu] ;
      QQ = Q * Q ;
      Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
          
      // Assemble the stout links exp(iQ)U_{mu} 
      next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      

    }


    for (mu = 0; mu <= Nd-1; mu++)
    {
      u[mu]    =  next[mu] ;
    }



  }


  void wilson_flow(XMLWriter& xml,
		   multi1d<LatticeColorMatrix> & u, int nstep, 
		   Real  wflow_eps, int jomit)
  {
    Real gact4i, gactij;
    int dim = nstep + 1 ;
    multi1d<Real> gact4i_vec(dim);
    multi1d<Real> gactij_vec(dim);
    multi1d<Real> step_vec(dim);



    measure_wilson_gauge(u,gactij,gact4i,jomit) ;
    gact4i_vec[0] = gact4i ;
    gactij_vec[0] = gactij ;
    step_vec[0] = 0.0 ;

    //  QDPIO::cout << "WFLOW " << 0.0 << " " << gact4i << " " << gactij <<  endl ; 

    QDPIO::cout << "START_ANALYZE_wflow" << endl ; 
    QDPIO::cout << "WFLOW time gact4i gactij" << endl ; 

    for(int i=0 ; i < nstep ; ++i)
    {
      wilson_flow_one_step(u,wflow_eps) ;

      measure_wilson_gauge(u,gactij,gact4i,jomit) ;
      gact4i_vec[i+1] = gact4i ;
      gactij_vec[i+1] = gactij ;


      Real xx = (i + 1) * wflow_eps ;
      QDPIO::cout << "WFLOW " << xx << " " << gact4i << " " << gactij <<  endl ; 

      step_vec[i+1] = xx ;

    }
    QDPIO::cout << "END_ANALYZE_wflow" << endl ; 

    push(xml, "wilson_flow_results");
    write(xml,"wflow_step",step_vec) ; 
    write(xml,"wflow_gact4i",gact4i_vec) ; 
    write(xml,"wflow_gactij",gactij_vec) ; 
    pop(xml);  // elem

  }




/////////////////


  struct Eps
  {
    Real prev_;
    Real this_;
    Real next_;
    Real max_;
    Real cut_;
  }; 

  void wilson_flow_one_step_adapt(multi1d<LatticeColorMatrix> & u, Eps & rho, Real tol)
  {
    int mu, dir;
    multi1d<LatticeColorMatrix> dest(Nd);
    multi1d<LatticeColorMatrix> next(Nd);

    rho.prev_ = rho.this_;

    // -------------------------------------

    multi1d<bool> smear_in_this_dirP(4) ;
    multi2d<Real> rho_a(4,4) ;
    multi2d<Real> rho_b1(4,4) ;
    multi2d<Real> rho_b2(4,4) ;
    multi2d<Real> rho_c(4,4) ;
    multi2d<Real> rho_d0(4,4) ;
    multi2d<Real> rho_d1(4,4) ;

    multi1d<LatticeColorMatrix> u0(u);
    multi1d<LatticeColorMatrix> du(Nd);

    LatticeColorMatrix  Q, QQ, C ;
    LatticeColorMatrix  Q2, QQ2  ;

    multi1d<LatticeColorMatrix> Q0(Nd);
    multi1d<LatticeColorMatrix> Q1(Nd);

    multi1d<LatticeComplex> f;   // routine will resize these
 
    bool go=true;

    while (go) {

      for (mu = 0; mu <= Nd-1; mu++)
      {
        smear_in_this_dirP(mu) = true ;
        for (dir = 0; dir <= Nd-1; dir++)
        {
          rho_a[mu][dir] = rho.next_ * 0.25 ;
  
          rho_b1[mu][dir] = rho.next_ * 8.0/9.0 ;
          rho_b2[mu][dir] = rho.next_ * 17.0/36.0 ;
  
          rho_c[mu][dir] = rho.next_  * 3.0/4.0 ;
  
          rho_d0[mu][dir] = rho.next_ ;
          rho_d1[mu][dir] = rho.next_ * 2.0 ;
        }
      }
  
  
      // This is W1
      Stouting::smear_links(u0, dest,smear_in_this_dirP, rho_a);
  
      for (mu = 0; mu <= Nd-1; mu++)
      {
  
        // this is the second order integrator
        Stouting::getQsandCs(dest, Q1[mu],QQ,C,mu,smear_in_this_dirP,rho_d1) ;
        Stouting::getQsandCs(u0   , Q0[mu],QQ,C,mu,smear_in_this_dirP,rho_d0) ;
  
        Q = Q1[mu] - Q0[mu] ;
        QQ = Q * Q ;
        Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
            
        du=(f[0] + f[1]*Q + f[2]*QQ)*u[mu];      
        // this is the second order integrator
  
        Stouting::getQsandCs(dest, Q1[mu],QQ,C,mu,smear_in_this_dirP,rho_b1) ;
        Stouting::getQsandCs(u0   , Q0[mu],QQ,C,mu,smear_in_this_dirP,rho_b2) ;
  
        Q = Q1[mu] - Q0[mu] ;
        QQ = Q * Q ;
        Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
            
        // Assemble the stout links exp(iQ)U_{mu} 
        // This is W2
        next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      
  
      }
  
      for (mu = 0; mu <= Nd-1; mu++)
      {
        u[mu]    =  next[mu] ;
        dest[mu] =  next[mu] ;      
      }
  
  
      for (mu = 0; mu <= Nd-1; mu++)
      {
        Stouting::getQsandCs(dest, Q2,QQ,C,mu,smear_in_this_dirP,rho_c) ;
  
        Q = Q2 - Q1[mu] + Q0[mu] ;
        QQ = Q * Q ;
        Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
            
        // Assemble the stout links exp(iQ)U_{mu} 
        // This is W3
        next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      
  
      }
  
      multi1d<Real> dist(Nd);
      for (mu = 0; mu <= Nd-1; mu++)
      {
        u[mu]    =  next[mu] ;
        du[mu] = u[mu] - u0[mu];
        dist[mu] = globalMax( sqrt( localNorm2( du[mu]) ) /(Nc*Nc) );
      }
  
      Real max_dist = dist[0];
      if ( toDouble(dist[1]) > toDouble(dist[0]) ) max_dist = dist[1];
      if ( toDouble(dist[2]) > toDouble(dist[1]) ) max_dist = dist[2];
      if ( toDouble(dist[3]) > toDouble(dist[2]) ) max_dist = dist[3];
 
//      QDPIO::cout << "eps : " << rho.next_ << " , max_dist : " << max_dist  << endl ; 

      if ( toDouble(max_dist) < toDouble(tol) ) {
        go=false;
        rho.this_ = rho.next_;
      } else {
        u = u0;
      }
  
      rho.next_ = rho.next_ * 0.95*pow(toDouble(tol/max_dist), 1./3.);
      if ( toDouble(rho.next_) > toDouble(rho.max_) ) { rho.next_ = rho.max_; }
      if ( toDouble(rho.next_) > toDouble(rho.cut_-rho.this_) ) { rho.next_ = rho.cut_ - rho.this_; }

    }

  }





  Real wilson_flow(XMLWriter& xml,
		   multi1d<LatticeColorMatrix> & u, Real wtime, Real eps_init, Real tol, 
		   int jomit)
  {
    Real gact4i, gactij;
    int dim = 10000 ;
    multi1d<Real> gact4i_vec(dim);
    multi1d<Real> gactij_vec(dim);
    multi1d<Real> step_vec(dim);
    multi1d<Real> gact4i_vec_tmp(dim);
    multi1d<Real> gactij_vec_tmp(dim);
    multi1d<Real> step_vec_tmp(dim);

    measure_wilson_gauge(u,gactij,gact4i,jomit) ;
    gact4i_vec[0] = gact4i ;
    gactij_vec[0] = gactij ;
    step_vec[0] = 0.0 ;


    QDPIO::cout << "START_ANALYZE_wflow" << endl ; 
    QDPIO::cout << "WFLOW time gact4i gactij" << endl ; 

    Eps eps;
    eps.prev_ = 0.0;
    eps.this_ = 0.0;
    eps.next_ = eps_init;
    eps.max_ = 0.2;
    eps.cut_ = wtime; 
    Real t(0.0);
    int counter(1);
    QDPIO::cout << "WFLOW " << 0.0 << " " << gact4i << " " << gactij <<  endl ; 
    while (toFloat(wtime-t) > 1e-7)
    {
      wilson_flow_one_step_adapt(u, eps, tol) ;
      t += eps.this_ ;
      eps.cut_ = wtime - t;

      measure_wilson_gauge(u,gactij,gact4i,jomit) ;
      gact4i_vec[counter] = gact4i ;
      gactij_vec[counter] = gactij ;

      QDPIO::cout << "WFLOW " << t << " " << gact4i << " " << gactij <<  endl ; 

      step_vec[counter] = t ;

      // If we hit the end of the multid, the resize it
      if ( (counter+1)%dim==0) {
        // This nonsense is needed since resize() kills the data within the multi1d

        for (int k=0; k<counter+1; ++k) {
          gact4i_vec_tmp[k] = gact4i_vec[k];
          gactij_vec_tmp[k] = gactij_vec[k];
          step_vec_tmp[k] = step_vec[k];
        }

        gact4i_vec.resize(counter+1+dim);
        gactij_vec.resize(counter+1+dim);
        step_vec.resize(counter+1+dim);

        for (int k=0; k<counter+1; ++k) {
          gact4i_vec[k] = gact4i_vec_tmp[k];
          gactij_vec[k] = gactij_vec_tmp[k];
          step_vec[k] = step_vec_tmp[k];
        }

        gact4i_vec_tmp.resize(counter+1+dim);
        gactij_vec_tmp.resize(counter+1+dim);
        step_vec_tmp.resize(counter+1+dim);

      }
      counter++;
    }
    {
      // This nonsense is needed since resize() kills the data within the multi1d
      for (int k=0; k<counter; ++k) {
        gact4i_vec_tmp[k] = gact4i_vec[k];
        gactij_vec_tmp[k] = gactij_vec[k];
        step_vec_tmp[k] = step_vec[k];
      }
      gact4i_vec.resize(counter);
      gactij_vec.resize(counter);
      step_vec.resize(counter);
      for (int k=0; k<counter; ++k) {
        gact4i_vec[k] = gact4i_vec_tmp[k];
        gactij_vec[k] = gactij_vec_tmp[k];
        step_vec[k] = step_vec_tmp[k];
      }
    }

    QDPIO::cout << "END_ANALYZE_wflow" << endl ; 

    push(xml, "wilson_flow_results");
    write(xml,"wflow_step",step_vec) ; 
    write(xml,"wflow_gact4i",gact4i_vec) ; 
    write(xml,"wflow_gactij",gactij_vec) ; 
    pop(xml);  // elem

    // returns the last full step size (should be eps_next_, but that one gets chopped);
    // eps.this_ is close enough
    return eps.prev_;

  }




}  // end namespace Chroma


