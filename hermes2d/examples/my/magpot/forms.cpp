template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  double perm=4*M_PI*1e-7; // permeabilita vakua
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += 1/perm * wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] 
                                 - u->dx[i]/e->x[i] * v->val[i] 
                                + u->val[i]/sqr(e->x[i]) * v->val[i]
    ); 
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += -wt[i] * ( F<Real, Scalar>(e->x[i], e->y[i]) * v->val[i]);          
  }

  return result;
}
