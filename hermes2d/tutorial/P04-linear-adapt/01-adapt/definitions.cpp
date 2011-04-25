#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_motor, double eps_motor, 
                        std::string mat_air, double eps_air) : WeakForm(1)
  {
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, mat_motor, eps_motor));
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, mat_air, eps_air));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, mat_motor, eps_motor));
    add_vector_form(new DefaultResidualDiffusion(0, mat_air, eps_air));
  };
};
