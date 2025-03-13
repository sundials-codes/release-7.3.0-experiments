#pragma once
#include "sundials/sundials_nvector.h"
#include "gray_scott_model.hpp"
#include "options.hpp"

class Integrator {
protected:
  using Model = GrayScottModel<sunrealtype, sunindextype>;
  using Opts = Options<sunrealtype, sunindextype>;

  static int f(sunrealtype, N_Vector y, N_Vector f, void *const data) noexcept
  {
    const auto &model = *static_cast<Model*>(data);
    model.f(N_VGetArrayPointer(y), N_VGetArrayPointer(f));
    return 0;
  }

  static int f_diffusion(sunrealtype, N_Vector y, N_Vector f, void *const data) noexcept
  {
    const auto &model = *static_cast<Model*>(data);
    model.f_diffusion(N_VGetArrayPointer(y), N_VGetArrayPointer(f));
    return 0;
  }

  static int f_reaction(sunrealtype, N_Vector y, N_Vector f, void *const data) noexcept
  {
    const auto &model = *static_cast<Model*>(data);
    model.f_reaction(N_VGetArrayPointer(y), N_VGetArrayPointer(f));
    return 0;
  }

  static int dom_eig(sunrealtype, N_Vector, N_Vector, sunrealtype *lambdaR, sunrealtype *lambdaI, void *data, N_Vector, N_Vector, N_Vector) noexcept
  {
    const auto& model = *static_cast<Model*>(data);
    *lambdaR = -model.spectral_radius_estimate();
    *lambdaI = 0;
    return 0;
  }

public:
  virtual ~Integrator() noexcept = default;
  virtual void solve(N_Vector y) const noexcept = 0;
  virtual void print_stats() const noexcept = 0;
};
