#include <stdio.h>
#include "SharedStatus.h"


SharedStatus *SharedStatus::self = NULL;

void SharedStatus::print()
{
	printf("Average vtkdax mc filter: %lf, normal filter: %lf, total: %lf\n",
			average(vtkdax_mc_time), average(vtkdax_norm_time), average(vtkdax_total_time) );
	printf("Average dax mc no resolution: %lf\n", average(dax_mc_nores_time));
	printf("Average dax mc with resolution: %lf\n", average(dax_mc_res_time));
	printf("Average vtk mc: %lf\n", average(vtk_mc_time));

	printf("Average normal filter tx time: %lf\n", 
		average(norm_copyPointsToDev_time)+
		average(norm_copyCellsToDev_time)+
		average(norm_copyToMem_time)
	);
	printf("Average host->dev time: copy points: %lf, copy cells: %lf\n",
			average(norm_copyPointsToDev_time), average(norm_copyCellsToDev_time));
	printf("Average dev->host time: %lf\n", average(norm_copyToMem_time));

}
