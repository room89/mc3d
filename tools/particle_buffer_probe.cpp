// Build separately against baseline/current mc3d_core. Capacities in bytes, not RSS.
#include "cell_cluster.h"
#include <iostream>
#include <type_traits>
namespace mc3d {
class CellClusterTestAccess {
 public:
  template<class Cluster> static void Probe(Cluster& c, unsigned outside) {
    c.SetApex({0,0,0}); c.SetSize(64,1,1);
    c.Initialize(64,1,1,0,.04,.5,std::make_unique<Geometry>(),0,0,1);
    constexpr std::size_t count=1048576;
    for (std::size_t i=0;i<count;++i)
      c.cells[i%64].GetBuffer().emplace_back(Point(i%100<outside ? -1. : 63.5,.5,.5),Point(double(i),0,0));
    std::size_t local_size=0,local_capacity=0;
    for (auto& cell:c.cells) {
      local_size+=cell.GetBuffer().size()*sizeof(Particle);
      local_capacity+=cell.GetBuffer().capacity()*sizeof(Particle);
    }
    std::size_t gathered_size=0,gathered_capacity=0,metadata=0,temp_size=0,temp_capacity=0;
    if constexpr (std::is_void_v<decltype(c.GatherParticles())>) {
      c.GatherParticles();
      gathered_size=c.partile_buffer.size()*sizeof(Particle);
      gathered_capacity=c.partile_buffer.capacity()*sizeof(Particle);
      // Old implementation reserves every input chunk before filtering.
      temp_size=(count/100*outside+std::min<std::size_t>(count%100,outside))*sizeof(Particle);
      temp_capacity=count*sizeof(Particle);
      c.DistributeParticles(c.partile_buffer);
    } else {
      auto sources=c.GatherParticles();
      gathered_size=c.partile_buffer.size()*sizeof(Particle);
      gathered_capacity=c.partile_buffer.capacity()*sizeof(Particle);
      metadata=sources.capacity()*sizeof(typename decltype(sources)::value_type);
      auto storage=c.DistributeSources(sources,c.partile_buffer);
      temp_size=storage.size*sizeof(Particle);
      temp_capacity=storage.capacity*sizeof(Particle);
    }
    std::cout << "{\"outside_percent\":" << outside << ",\"source_size_bytes\":" << local_size
              << ",\"source_capacity_bytes\":" << local_capacity
              << ",\"gathered_size_bytes\":" << gathered_size << ",\"gathered_capacity_bytes\":" << gathered_capacity
              << ",\"metadata_capacity_bytes\":" << metadata << ",\"temporary_size_bytes\":" << temp_size
              << ",\"temporary_capacity_bytes\":" << temp_capacity
              << ",\"final_size_bytes\":" << c.partile_buffer.size()*sizeof(Particle)
              << ",\"final_capacity_bytes\":" << c.partile_buffer.capacity()*sizeof(Particle) << "}\n";
  }
};
}
int main(int argc,char** argv) {
  const auto outside=argc>1?std::stoul(argv[1]):0;
  if(outside>100) return 1;
  mc3d::CellCluster c(4,123);
  mc3d::CellClusterTestAccess::Probe(c,outside);
}
