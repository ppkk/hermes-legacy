#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::SurfaceMatrixForms;
using namespace WeakFormsH1::SurfaceVectorForms;
  
/* Weak forms */
class CustomWeakFormAcoustics : public WeakForm
{ 
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho, double sound_speed, double omega)
  : WeakForm(1) {
    scalar ii =  cplx(0.0, 1.0);
    // Volumetric terms.
    add_matrix_form(new DefaultLinearDiffusion(0, 0, 1.0/rho));
    add_matrix_form(new DefaultLinearMass(0, 0, - sqr(omega) / rho / sqr(sound_speed)));
    // Term generated by the Newton condition at outlet.
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0,  bdy_newton, -ii * omega / rho / sound_speed));
  };
};
