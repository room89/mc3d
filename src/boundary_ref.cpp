#include "boundary_ref.h"


boundary_ref::boundary_ref(void)
{
}


boundary_ref::~boundary_ref(void)
{
}

int boundary_ref::bondary_condition(deque<particle> *cluster_particle)
{
	return ref->bondary_condition(cluster_particle);
}
