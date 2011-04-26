// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ELASTICITY_WEAK_FORMS_H
#define __H2D_ELASTICITY_WEAK_FORMS_H

#include "../integrals/integrals_h1.h"

/* Default weak form for linear elasticity (Lame equations)
   with Dirichlet and/or zero Neumann BC

   Nonzero Neumann and Newton boundary conditions can be enabled
   by creating a descendant and adding surface forms to it.
*/

namespace WeakFormsElasticity {

namespace VolumetricMatrixForms {

class DefaultLinearXX : public WeakForm::MatrixFormVol
{
public:
  DefaultLinearXX(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
  DefaultLinearXX(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM, area), lambda(lambda), mu(mu), gt(gt) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    if (gt == HERMES_PLANAR)
        return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    else if (gt == HERMES_AXISYM_X)
    {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
            result += wt[i] * (e->x[i] * lambda * (u->dx[i] * v->dx[i] +
                                         u->val[i]/e->x[i] * v->dx[i] +
                                         u->dx[i] * v->val[i]/e->x[i] +
                                         1/sqr(e->x[i]) * u->val[i] * v->val[i]) +
                               e->x[i] * mu * (2 * u->dx[i] * v->dx[i] +
                                         2 * 1/sqr(e->x[i]) * u->val[i] * v->val[i] +
                                         u->dy[i] * v->dy[i]));
        return result;
    }
    else
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
          result += wt[i] * (e->x[i] * lambda * (u->dx[i] * v->dx[i] +
                                       u->val[i]/e->x[i] * v->dx[i] +
                                       u->dx[i] * v->val[i]/e->x[i] +
                                       1/sqr(e->x[i]) * u->val[i] * v->val[i]) +
                             e->x[i] * mu * (2 * u->dx[i] * v->dx[i] +
                                       2 * 1/sqr(e->x[i]) * u->val[i] * v->val[i] +
                                       u->dy[i] * v->dy[i]));
      return result;
    }
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  virtual WeakForm::MatrixFormVol* clone() {
    return new DefaultLinearXX(*this);
  }

private:
  double lambda, mu;
  GeomType gt;
};

class DefaultLinearXY : public WeakForm::MatrixFormVol
{
public:
  DefaultLinearXY(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
  DefaultLinearXY(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), lambda(lambda), mu(mu), gt(gt) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    if (gt == HERMES_PLANAR)
        return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
                mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
    else if (gt == HERMES_AXISYM_X)
    {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
        result += wt[i] * (e->y[i] * lambda * (u->dx[i] * v->dy[i] +
                           u->dx[i] * v->val[i]/e->y[i]) +
                           e->y[i] * mu * u->dy[i] * v->dx[i]);
        return result;
    }
    else
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
          result += wt[i] * (e->x[i] * lambda * (u->dy[i] * v->dx[i] +
                             u->dy[i] * v->val[i]/e->x[i]) +
                             e->x[i] * mu * u->dx[i] * v->dy[i]);
      return result;
    }
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  virtual WeakForm::MatrixFormVol* clone() {
    return new DefaultLinearXY(*this);
  }

private:
  double lambda, mu;
  GeomType gt;
};

class DefaultLinearYY : public WeakForm::MatrixFormVol
{
public:
  DefaultLinearYY(unsigned int i, unsigned int j, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu), gt(gt) { }
  DefaultLinearYY(unsigned int i, unsigned int j, std::string area, double lambda, double mu, GeomType gt = HERMES_PLANAR)
    : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM, area), lambda(lambda), mu(mu), gt(gt) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    if (gt == HERMES_PLANAR)
        return  mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    else if (gt == HERMES_AXISYM_X)
    {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
        result += wt[i] * (e->y[i] * lambda * (u->dx[i] * v->dx[i]) +
                           e->y[i] * mu * (u->dy[i] * v->dy[i] +
                           2 * u->dx[i] * v->dx[i]));
        return result;
    }
    else
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
          result += wt[i] * (e->x[i] * lambda * (u->dy[i] * v->dy[i]) +
                             e->x[i] * mu * (u->dx[i] * v->dx[i] +
                             2 * u->dy[i] * v->dy[i]));
      return result;
    }
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  virtual WeakForm::MatrixFormVol* clone() {
    return new DefaultLinearYY(*this);
  }

private:
  double lambda, mu;
  GeomType gt;
};


/*
class MultiComponentDefaultVolumetricMatrixFormLinearSym : public WeakForm::MultiComponentMatrixFormVol
{
public:
  MultiComponentDefaultVolumetricMatrixFormLinearSym(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double lambda, double mu)
    : WeakForm::MultiComponentMatrixFormVol(coordinates, HERMES_SYM), lambda(lambda), mu(mu) { }
  MultiComponentDefaultVolumetricMatrixFormLinearSym(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area, double lambda, double mu)
    : WeakForm::MultiComponentMatrixFormVol(coordinates, HERMES_SYM, area), lambda(lambda), mu(mu) { }

  template<typename Real, typename Scalar>
  void matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                   Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const {
    Scalar result_0_0 = 0;
    Scalar result_1_1 = 0;
    for (int i = 0; i < n; i++) {
      result_0_0 += wt[i] * (lambda + 2*mu) * u->dx[i] * v->dx[i];
      result_0_0 += wt[i] * mu * u->dy[i] * v->dy[i];

      result_1_1 += wt[i] * mu * u->dx[i] * v->dx[i];
      result_1_1 += wt[i] * (lambda + 2*mu) * u->dy[i] * v->dy[i];
    }
    result.push_back(result_0_0);
    result.push_back(result_1_1);
  }

  virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                     Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const {
    matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext, result);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const {
    Hermes::vector<Ord> result;
    matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext, result);

    // Choose the maximum order.
    Ord to_return = result[0];
    if(result[1] > to_return)
      to_return = result[1];
    return to_return;
  }

private:
  double lambda, mu;
};
*/
};
};

#endif
