#pragma once

#include "integrator.hpp"
#include <sundials/sundials_context.hpp>
#include <arkode/arkode_splittingstep.h>
#include <arkode/arkode_erkstep.h>

class Splitting : public Integrator
{
private:
  void *arkode_diff{};
  void *arkode_splitting{};
  std::array<SUNStepper, 3> steppers{};

  template <typename F>
  struct Content
  {
    sunrealtype t{};
    N_Vector v;
    const F f;
    const Model &model;

    Content(const N_Vector tmpl, const F f, const Model &model) : v(N_VClone(tmpl)), f(f), model(model) {}

    ~Content() { N_VDestroy(v); }

    static Content& from_stepper(SUNStepper s)
    {
      void* content = nullptr;
      SUNStepper_GetContent(s, &content);
      return *static_cast<Content*>(content);
    }
  };

  template <typename F>
  static SUNStepper create_stepper(const N_Vector y, const F f, const Model &model, const sundials::Context& ctx)
  {
    SUNStepper stepper{};
    SUNStepper_Create(ctx, &stepper);
    SUNStepper_SetContent(stepper, new Content<F>(y, f, model));

    auto reset = [](const SUNStepper s, const sunrealtype tR, const N_Vector vR)
    {
      auto& content = Content<F>::from_stepper(s);
      content.t     = tR;
      N_VScale(SUN_RCONST(1.0), vR, content.v);
      return 0;
    };
    SUNStepper_SetResetFn(stepper, reset);

    auto empty_func = [](auto...) { return 0; };
    SUNStepper_SetStopTimeFn(stepper, empty_func);
    SUNStepper_SetStepDirectionFn(stepper, empty_func);
    SUNStepper_SetFullRhsFn(stepper, empty_func);

    auto evolve = [](const SUNStepper s, const sunrealtype tout, const N_Vector vret, sunrealtype* const tret)
    {
      auto& content = Content<F>::from_stepper(s);
      content.f(content.model, tout - content.t, N_VGetArrayPointer(content.v), N_VGetArrayPointer(vret));
      *tret = tout;
      return 0;
    };
    SUNStepper_SetEvolveFn(stepper, evolve);

    auto destroy = [](SUNStepper s)
    {
      delete &Content<F>::from_stepper(s);
      return 0;
    };
    SUNStepper_SetDestroyFn(stepper, destroy);
    return stepper;
  }

public:
  Splitting(const Model &model, const N_Vector y, const sundials::Context &ctx, const Opts &opts)
  {
    steppers[0] = create_stepper(y, [](const Model &model, const sunrealtype dt, const sunrealtype *const vin, sunrealtype *const vout) {
      model.evolve_reaction_u(dt, vin, vout);
    }, model, ctx);

    steppers[1] = create_stepper(y, [](const Model &model, const sunrealtype dt, const sunrealtype *const vin, sunrealtype *const vout) {
      model.evolve_reaction_v(dt, vin, vout);
    }, model, ctx);

    if (opts.low_storage) {
      arkode_diff = LSRKStepCreateSTS(f_diffusion, model.t0(), y, ctx);
      LSRKStepSetDomEigFn(arkode_diff, dom_eig);
    } else {
      arkode_diff = ERKStepCreate(f_diffusion, model.t0(), y, ctx);
      ARKodeSetOrder(arkode_diff, opts.order);
    }
    ARKodeSetUserData(arkode_diff, (void*) &model);
    ARKodeSetFixedStep(arkode_diff, opts.dt);
    ARKodeCreateSUNStepper(arkode_diff, &steppers[2]);

    arkode_splitting = SplittingStepCreate(steppers.data(), steppers.size(), model.t0(), y, ctx);
    ARKodeSetOrder(arkode_splitting, opts.order);
    ARKodeSetFixedStep(arkode_splitting, opts.dt);
    ARKodeSetStopTime(arkode_splitting, model.tf());
    ARKodeSetMaxNumSteps(arkode_splitting, -1);
    ARKodeSetUserData(arkode_splitting, (void *)&model);
  }

  ~Splitting() noexcept
  {
    ARKodeFree(&arkode_diff);
    ARKodeFree(&arkode_splitting);
    for (auto stepper : steppers) {
      SUNStepper_Destroy(&stepper);
    }
  }

  void solve(const N_Vector y) const noexcept
  {
    sunrealtype tret;
    ARKodeEvolve(arkode_splitting, std::numeric_limits<sunrealtype>::max(), y, &tret, ARK_NORMAL);
  }

  void print_stats() const noexcept
  {
    ARKodePrintAllStats(arkode_splitting, stdout, SUN_OUTPUTFORMAT_TABLE);
    ARKodePrintAllStats(arkode_diff, stdout, SUN_OUTPUTFORMAT_TABLE);
  }
};