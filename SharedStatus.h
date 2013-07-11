#ifndef SHARED_STATUS_H
#define SHARED_STATUS_H

#include <vector>
#include <numeric>
#include <assert.h>
#include <math.h>

class SharedStatus
{
public:
	std::vector<double> norm_reqInfo_reqUpdate_time, norm_reqUpdate_reqData_time,
		vtkdax_norm_time, vtkdax_mc_time, vtkdax_total_time,
		dax_mc_nores_time, dax_mc_res_time,
		vtk_mc_time;
	std::vector<double> norm_copyPointsToDev_time, norm_copyCellsToDev_time, norm_copyToMem_time;

	void print() ;

private:
	template <typename T>
	T average(std::vector<T> v) { return v.size()? std::accumulate(v.begin(), v.end(), (T)0) / (double)v.size() : (T)NAN; }

    SharedStatus()
    {

    }

    ///////////// static //////////////
private:
    static SharedStatus *self;
public:
    static void init( )
    {
        if ( NULL == SharedStatus::self)
        	SharedStatus::self = new SharedStatus();
    }
    static inline SharedStatus *getInstance() {
    	assert(SharedStatus::self);
    	return SharedStatus::self;
    }

};

#endif
