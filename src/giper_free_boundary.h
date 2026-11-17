#pragma once

#include "free_boundary.h"

namespace mc3d {
class HyperFreeBoundary : public FreeBoundary {
 private:
  bool stochastic_rounding_ = false;
 public:
  void SetStochasticRounding(bool enabled) { stochastic_rounding_ = enabled; }
  HyperFreeBoundary();
  HyperFreeBoundary(Point pstn, Point nrml, unsigned int np, double S, double T,
                    double alpha = 0);
  ~HyperFreeBoundary();
  BoundaryInflow GenerateInflow(double dt) override;
};
}  // namespace mc3d
