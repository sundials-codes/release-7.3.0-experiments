#include "gray_scott_model.hpp"
#include "options.hpp"
#include "erk.hpp"
#include "dirk.hpp"
#include "splitting.hpp"
#include <stdexcept>
#include <omp.h>
#include <nvector/nvector_openmp.h>
#include <chrono>

int main(int argc, char *argv[]) {
  const Options<sunrealtype, sunindextype> opts(argc, argv);
  const sundials::Context ctx;

  const GrayScottModel<sunrealtype, sunindextype> model(opts.grid_pts_1d, opts.threads);

  std::printf("Solving with %d threads\n", opts.threads);
  const auto y = N_VNew_OpenMP(model.dim(), opts.threads, ctx);
  model.initial_condition(N_VGetArrayPointer(y));

  const auto integrator = [&]() -> std::unique_ptr<Integrator> {
    if (opts.method == "ERK") {
      return std::make_unique<ERK>(model, y, ctx, opts);
    } else if (opts.method == "DIRK") {
      return std::make_unique<DIRK>(model, y, ctx, opts);
    } else if (opts.method == "Splitting") {
      return std::make_unique<Splitting>(model, y, ctx, opts);
    } else {
      throw std::invalid_argument("Invalid method name");
    }
  }();

  const auto start = std::chrono::high_resolution_clock::now();
  integrator->solve(y);
  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> duration = end - start;

  std::printf("Runtime: %g\n\n", duration.count());
  integrator->print_stats();

  FILE *f = fopen(opts.out_file.c_str(), "w");
  std::fprintf(f, "%g\n", duration.count());
  N_VPrintFile(y, f);
  fclose(f);

  N_VDestroy(y);

  return 0;
}
