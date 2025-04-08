#pragma once

#include "integrator.hpp"
#include <sundials/sundials_context.hpp>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_lsrkstep.h>

class ERK : public Integrator
{
private:
  void *arkode{};

public:
  ERK(const Model &model, const N_Vector y, const sundials::Context &ctx, const Opts &opts) noexcept
  {
    if (opts.low_storage) {
      arkode = LSRKStepCreateSTS(f, model.t0(), y, ctx);
      LSRKStepSetDomEigFn(arkode, dom_eig);
    } else {
      arkode = ERKStepCreate(f, model.t0(), y, ctx);
      ARKodeSetOrder(arkode, opts.order);
    }
    
    ARKodeSetFixedStep(arkode, opts.dt);
    ARKodeSStolerances(arkode, opts.rel_tol, opts.abs_tol);
    ARKodeSetStopTime(arkode, model.tf());
    ARKodeSetMaxNumSteps(arkode, -1);
    ARKodeSetUserData(arkode, (void *)&model);
  }

  ~ERK() noexcept
  {
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
