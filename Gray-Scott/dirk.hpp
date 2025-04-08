#pragma once

#include "integrator.hpp"
#include <sundials/sundials_context.hpp>
#include <arkode/arkode_arkstep.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_spbcgs.h>

class DIRK : public Integrator
{
private:
  void *arkode{};
  SUNMatrix mat{};
  SUNLinearSolver ls{};

  static int jacobian(sunrealtype, N_Vector y, N_Vector, SUNMatrix jac, void *const data, N_Vector, N_Vector, N_Vector) noexcept
  {
    const auto &model = *static_cast<Model*>(data);
    model.jacobian(N_VGetArrayPointer(y), SUNSparseMatrix_IndexValues(jac), SUNSparseMatrix_IndexPointers(jac), SUNSparseMatrix_Data(jac));
    return 0;
  }

public:
  DIRK(const Model &model, const N_Vector y, const sundials::Context &ctx, const Opts &opts) noexcept
  {
    arkode = ARKStepCreate(nullptr, f, model.t0(), y, ctx);
    ARKodeSetFixedStep(arkode, opts.dt);
    ARKodeSStolerances(arkode, opts.rel_tol, opts.abs_tol);
    ARKodeSetOrder(arkode, opts.order);
    ARKodeSetStopTime(arkode, model.tf());
    ARKodeSetMaxNumSteps(arkode, -1);
    ARKodeSetUserData(arkode, (void *)&model);
    ARKodeSetPredictorMethod(arkode, 1);
    ARKodeSetDeduceImplicitRhs(arkode, true);

    mat = SUNSparseMatrix(model.dim(), model.dim(), model.nnz(), CSR_MAT, ctx);
    ls = SUNLinSol_SPBCGS(y, SUN_PREC_NONE, -1, ctx);
    ARKodeSetLinearSolver(arkode, ls, mat);
    ARKodeSetJacFn(arkode, jacobian);
  }

  ~DIRK() noexcept
  {
    SUNMatDestroy(mat);
    SUNLinSolFree(ls);
    ARKodeFree(&arkode);
  }

  void solve(const N_Vector y) const noexcept
  {
    sunrealtype tret;
    ARKodeEvolve(arkode, std::numeric_limits<sunrealtype>::max(), y, &tret, ARK_NORMAL);
  }

  void print_stats() const noexcept
  {
    ARKodePrintAllStats(arkode, stdout, SUN_OUTPUTFORMAT_TABLE);
  }
};
