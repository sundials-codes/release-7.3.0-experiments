#pragma once

#include <array>
#include <cmath>

template <typename Real = double, typename Index = std::int64_t>
class GrayScottModel
{
private:
  static constexpr Index components = 2;
  static constexpr Index nnz_per_row = 6;
  static constexpr std::array<Real, 2> diff{2.0e-5L, 1.0e-5L};
  static constexpr Real a = 0.04L;
  static constexpr Real b = 0.06L;

  static constexpr std::array<Real, 2> tspan{0, 3500};

  const Index grid_pts_1D;
  const Index grid_pts;
  const Real dx;
  const Real dx2;
  const int threads;

  constexpr Index index(const Index i, const Index j, const Index comp) const noexcept
  {
    return components * (grid_pts_1D * i + j) + comp;
  }

public:
  constexpr GrayScottModel(const Index grid_pts_1D, const int threads = 1) noexcept
      : grid_pts_1D(grid_pts_1D),
        grid_pts(grid_pts_1D * grid_pts_1D),
        dx(Real(2) / grid_pts_1D),
        dx2(std::pow(dx, -2)),
        threads(threads)
  {
  }

  [[nodiscard]] constexpr Index dim() const noexcept
  {
    return components * grid_pts;
  }

  constexpr void initial_condition(Real *const y0) const noexcept
  {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
    for (Index i = 0; i < grid_pts_1D; i++)
    {
      for (Index j = 0; j < grid_pts_1D; j++)
      {
        const auto x = 2 * static_cast<Real>(i) / grid_pts_1D - 1;
        const auto y = 2 * static_cast<Real>(j) / grid_pts_1D - 1;
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

  constexpr void f(const Real *const y, Real *const f) const noexcept
  {
    f_diffusion(y, f);
    f_reaction<true>(y, f);
  }

  constexpr Index nnz() const noexcept
  {
    return dim() * nnz_per_row;
  }

  constexpr void jacobian(const Real *const y, Index *const indexvals, Index *const indexptrs, Real *const jac_data) const noexcept
  {
    // In CSR format
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
    for (Index i = 0; i < grid_pts_1D; i++)
    {
      const auto i_prev = i == 0 ? grid_pts_1D - 1 : i - 1;
      const auto i_next = i == grid_pts_1D - 1 ? 0 : i + 1;
      for (Index j = 0; j < grid_pts_1D; j++)
      {
        const auto j_prev = j == 0 ? grid_pts_1D - 1 : j - 1;
        const auto j_next = j == grid_pts_1D - 1 ? 0 : j + 1;

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

        jac_data[idx_val] = diff[0] * dx2;
        jac_data[idx_val + 1] = diff[0] * dx2;
        jac_data[idx_val + 2] = -4 * diff[0] * dx2 - std::pow(v, 2) - a;
        jac_data[idx_val + 3] = -2 * u * v;
        jac_data[idx_val + 4] = diff[0] * dx2;
        jac_data[idx_val + 5] = diff[0] * dx2;

        jac_data[idx_val + 6] = diff[1] * dx2;
        jac_data[idx_val + 7] = diff[1] * dx2;
        jac_data[idx_val + 8] = std::pow(v, 2);
        jac_data[idx_val + 9] = -4 * diff[1] * dx2 + 2 * u * v - (a + b);
        jac_data[idx_val + 10] = diff[1] * dx2;
        jac_data[idx_val + 11] = diff[1] * dx2;
      }
    }
    indexptrs[dim()] = nnz();
  }

  constexpr void f_diffusion(const Real *const y, Real *const f) const noexcept
  {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
    for (Index i = 0; i < grid_pts_1D; i++)
    {
      const auto i_prev = i == 0 ? grid_pts_1D - 1 : i - 1;
      const auto i_next = i == grid_pts_1D - 1 ? 0 : i + 1;
      for (Index j = 0; j < grid_pts_1D; j++)
      {
        const auto j_prev = j == 0 ? grid_pts_1D - 1 : j - 1;
        const auto j_next = j == grid_pts_1D - 1 ? 0 : j + 1;

        for (Index comp = 0; comp < components; comp++)
        {
          f[index(i, j, comp)] = diff[comp] * dx2 * (-4 * y[index(i, j, comp)] + y[index(i_prev, j, comp)] + y[index(i_next, j, comp)] + y[index(i, j_prev, comp)] + y[index(i, j_next, comp)]);
        }
      }
    }
  }

  [[nodiscard]] constexpr Real spectral_radius_estimate() const noexcept {
    return std::max(diff[0], diff[1]) * 2 * grid_pts;
  }

  template <bool Increment = false>
  constexpr void f_reaction(const Real *const y, Real *const f) const noexcept
  {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
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

  constexpr void evolve_reaction_u(const Real dt, const Real *const y, Real *const yNext) const noexcept {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
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

  constexpr void evolve_reaction_v(const Real dt, const Real *const y, Real *const yNext) const noexcept {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
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
  constexpr void jacobian_reaction(const Real *const y, const Matrix mat) const noexcept {
    #pragma omp parallel for schedule(static) default(none) num_threads(threads)
    for (Index i = 0; i < grid_pts_1D; i++)
    {
      for (Index j = 0; j < grid_pts_1D; j++)
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
