#include "timer.h"

timer::timer(void)
{
}

timer::~timer(void)
{
}

void timer::checkpoint()
{
	time a;
	a.data = MPI_Wtime();
	time_deque.push_back(a);
}

double timer::calc_av()
{
	double av_t = 0;
	deque<time>::iterator timer_iter = time_deque.begin();
	while(timer_iter != time_deque.end())
	{
		av_t += timer_iter->data;
		timer_iter++;
	}
	av_t -= timer_iter->data * time_deque.size();
	av_t /= time_deque.size();
	return av_t;
}