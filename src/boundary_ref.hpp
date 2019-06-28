#pragma once
#include "boundary.h"

struct boundary_ref
{
public:
	boundary *ref;
	int bondary_condition(deque<particle> *cluster_particle);
	boundary_ref(void);
	~boundary_ref(void);
};

