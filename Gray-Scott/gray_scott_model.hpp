#pragma once

#include <array>
#include <cmath>
#include <cstddef>

template <typename Real = double, typename Index = std::size_t>
class GrayScottModel
{
private:
  static constexpr Index components = 2;
  static constexpr Index nnz_per_row = 6;
  static constexpr std::array<Real, 2> diff{2.0e-5L, 1.0e-5L};
  static constexpr Real a = 0.04L;
  static constexpr Real b = 0.06L;

  static constexpr std::array<Real, 2> tspan{0, 3500};

  const Index grid_pts_1d;
  const Index grid_pts;
  const Real dx;
  const std::array<Real, 2> scaled_diff;
  const int threads;

  constexpr Index index(const Index i, const Index j, const Index comp) const noexcept
  {
    return components * (grid_pts_1d * i + j) + comp;
  }

  template <typename F>
  void stencil_iterate(const F f) const noexcept {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts_1d; i++)
    {
      const auto i_prev = i == 0 ? grid_pts_1d - 1 : i - 1;
      const auto i_next = i == grid_pts_1d - 1 ? 0 : i + 1;

      f(i_prev, i, i_next, grid_pts_1d - 1, 0, 1);
      for (Index j = 1; j < grid_pts_1d - 1; j++) {
        f(i_prev, i, i_next, j - 1, j, j + 1);
      }
      f(i_prev, i, i_next, grid_pts_1d - 2, grid_pts_1d - 1, 0);
    }
  }

public:
  constexpr GrayScottModel(const Index grid_pts_1d, const int threads = 1) noexcept
      : grid_pts_1d(grid_pts_1d),
        grid_pts(grid_pts_1d * grid_pts_1d),
        dx(Real(2) / grid_pts_1d),
        scaled_diff{diff[0] / std::pow(dx, 2), diff[1] / std::pow(dx, 2)},
        threads(threads)
  {
  }

  [[nodiscard]] constexpr Index dim() const noexcept
  {
    return components * grid_pts;
  }

  void initial_condition(Real *const y0) const noexcept
  {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts_1d; i++)
    {
      for (Index j = 0; j < grid_pts_1d; j++)
      {
        const auto x = 2 * static_cast<Real>(i) / grid_pts_1d - 1;
        const auto y = 2 * static_cast<Real>(j) / grid_pts_1d - 1;
        constexpr Real xp = 0.05L;
        constexpr Real yp = 0.02L;
        y0[index(i, j, 0)] = 1 - std::exp(-80 * (std::pow(x + xp, 2) + std::pow(y + yp, 2)));
        y0[index(i, j, 1)] = std::exp(-80 * (std::pow(x - xp, 2) + std::pow(y - yp, 2)));
      }
    }
  }

  [[nodiscard]] constexpr Real t0() const noexcept
  {
    return tspan[0];
  }

  [[nodiscard]] constexpr Real tf() const noexcept
  {
    return tspan[1];
  }

  void f(const Real *const y, Real *const f) const noexcept
  {
    f_diffusion(y, f);
    f_reaction<true>(y, f);
  }

  constexpr Index nnz() const noexcept
  {
    return dim() * nnz_per_row;
  }

  void jacobian(const Real *const y, Index *const indexvals, Index *const indexptrs, Real *const jac_data) const noexcept
  {
    // In CSR format
    stencil_iterate([=](const auto i_prev, const auto i, const auto i_next, const auto j_prev, const auto j, const auto j_next){
      for (Index comp = 0; comp < components; comp++)
      {
        const auto row = index(i, j, 0);
        const auto u = y[row];
        const auto v = y[row + 1];

        const auto idx_val = nnz_per_row * row;
        indexptrs[row] = idx_val;
        indexptrs[row + 1] = idx_val + nnz_per_row;

        indexvals[idx_val] = index(i_prev, j, 0);
        indexvals[idx_val + 1] = index(i, j_prev, 0);
        indexvals[idx_val + 2] = row;
        indexvals[idx_val + 3] = index(i, j, 1);
        indexvals[idx_val + 4] = index(i, j_next, 0);
        indexvals[idx_val + 5] = index(i_next, j, 0);

        indexvals[idx_val + 6] = index(i_prev, j, 1);
        indexvals[idx_val + 7] = index(i, j_prev, 1);
        indexvals[idx_val + 8] = row;
        indexvals[idx_val + 9] = index(i, j, 1);
        indexvals[idx_val + 10] = index(i, j_next, 1);
        indexvals[idx_val + 11] = index(i_next, j, 1);

        jac_data[idx_val] = scaled_diff[0];
        jac_data[idx_val + 1] = scaled_diff[0];
        jac_data[idx_val + 2] = -4 * scaled_diff[0] - std::pow(v, 2) - a;
        jac_data[idx_val + 3] = -2 * u * v;
        jac_data[idx_val + 4] = scaled_diff[0];
        jac_data[idx_val + 5] = scaled_diff[0];

        jac_data[idx_val + 6] = scaled_diff[1];
        jac_data[idx_val + 7] = scaled_diff[1];
        jac_data[idx_val + 8] = std::pow(v, 2);
        jac_data[idx_val + 9] = -4 * scaled_diff[1] + 2 * u * v - (a + b);
        jac_data[idx_val + 10] = scaled_diff[1];
        jac_data[idx_val + 11] = scaled_diff[1];
      }
    });
    indexptrs[dim()] = nnz();
  }

  void f_diffusion(const Real *const y, Real *const f) const noexcept
  {
    stencil_iterate([=](const auto i_prev, const auto i, const auto i_next, const auto j_prev, const auto j, const auto j_next){
      for (Index comp = 0; comp < components; comp++)
      {
        f[index(i, j, comp)] = scaled_diff[comp] * (y[index(i_prev, j, comp)] + y[index(i, j_prev, comp)] - 4 * y[index(i, j, comp)] + y[index(i, j_next, comp)] + y[index(i_next, j, comp)]);
      }
    });
  }

  [[nodiscard]] constexpr Real spectral_radius_estimate() const noexcept {
    return std::max(diff[0], diff[1]) * 2 * grid_pts;
  }

  template <bool Increment = false>
  void f_reaction(const Real *const y, Real *const f) const noexcept
  {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts; i++)
    {
      const auto u = y[2 * i];
      const auto v = y[2 * i + 1];
      const auto uv = u * std::pow(v, 2);
      const auto fu = -uv + a * (1 - u);
      const auto fv = uv - (a + b) * v;
      if constexpr (Increment)
      {
        f[2 * i] += fu;
        f[2 * i + 1] += fv;
      }
      else
      {
        f[2 * i] = fu;
        f[2 * i + 1] = fv;
      }
    }
  }

  void evolve_reaction_u(const Real dt, const Real *const y, Real *const yNext) const noexcept {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts; i++)
    {
      const auto u = y[2 * i];
      const auto v = y[2 * i + 1];
      const auto v2 = std::pow(v, 2);
      const auto av2 = a + v2;
      const auto expExpr = std::exp(dt * av2);
      yNext[2 * i] = (a * (expExpr - 1 + u) + u * v2) / (expExpr * av2);
      yNext[2 * i + 1] = v;
    }
  }

  void evolve_reaction_v(const Real dt, const Real *const y, Real *const yNext) const noexcept {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts; i++)
    {
      const auto u = y[2 * i];
      const auto v = y[2 * i + 1];
      const auto ab = a + b;
      const auto uv = u * v;
      yNext[2 * i] = u;
      yNext[2 * i + 1] = ab * v / (uv + std::exp(ab * dt) * (ab - uv));
    }
  }

  template <typename Matrix>
  void jacobian_reaction(const Real *const y, const Matrix mat) const noexcept {
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (Index i = 0; i < grid_pts_1d; i++)
    {
      for (Index j = 0; j < grid_pts_1d; j++)
      {
        const auto uIdx = index(i, j, 0);
        const auto vIdx = index(i, j, 1);
        const auto u = y[uIdx];
        const auto v = y[vIdx];
        mat(uIdx, uIdx) = -std::pow(v, 2) - a;
        mat(uIdx, vIdx) = -2 * u * v;
        mat(vIdx, uIdx) = std::pow(v, 2);
        mat(vIdx, vIdx) = 2 * u * v - (a + b);
      }
    }
  }
};
