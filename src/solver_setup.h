#pragma once

#include <memory>

#include "boundary.h"
#include "config.h"
#include "geometry.h"

namespace mc3d {

double ReferenceParticleDensity(const SimulationConfig& config);

struct BoundaryDescriptor {
  Point position;
  Point normal;
  Point mixing;
};

BoundaryDescriptor MakeBoundaryDescriptor(BoundaryFace face,
                                          const SimulationConfig& cfg);

std::unique_ptr<Boundary> MakeBoundary(const BoundaryDescriptor& descriptor,
                                       BoundaryType type,
                                       const SimulationConfig& cfg);

void ApplyGeometryPostProcessing(const SimulationConfig& cfg,
                                 Geometry& geometry);

std::unique_ptr<Geometry> BuildGeometry(const SimulationConfig& cfg);

}  // namespace mc3d
