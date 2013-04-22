
#ifndef piston_scalar_image3d_h
#define piston_scalar_image3d_h


#ifdef PISTON_ENABLED
#include <piston/image3d.h>
#endif PISTON_ENABLED

struct piston_scalar_image3d : piston::image3d<thrust::device_system_tag>
{
  typedef thrust::device_vector<dax::Scalar> PointDataContainer;
  PointDataContainer point_data_vector;
  typedef PointDataContainer::iterator PointDataIterator;

  piston_scalar_image3d(dax::Id xsize, dax::Id ysize, dax::Id zsize,
                        const std::vector<dax::Scalar> &data)
    : piston::image3d< thrust::device_system_tag >(xsize, ysize, zsize),
      point_data_vector(data)
  {
    assert(this->NPoints == this->point_data_vector.size());
  }

  PointDataIterator point_data_begin() {
    return this->point_data_vector.begin();
  }
  PointDataIterator point_data_end() {
    return this->point_data_vector.end();
  }
};

#endif
