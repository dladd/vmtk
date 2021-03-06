/*=========================================================================

  Program:   VMTK
  Module:    $RCSfile: vtkvmtkPolyDataCylinderHarmonicMappingFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:44 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// .NAME vtkvmtkPolyDataCylinderHarmonicMappingFilter - ..
// .SECTION Description
// ..

#ifndef __vtkvmtkPolyDataCylinderHarmonicMappingFilter_h
#define __vtkvmtkPolyDataCylinderHarmonicMappingFilter_h

#include "vtkvmtkPolyDataHarmonicMappingFilter.h"
#include "vtkvmtkWin32Header.h"

class VTK_VMTK_DIFFERENTIAL_GEOMETRY_EXPORT vtkvmtkPolyDataCylinderHarmonicMappingFilter : public vtkvmtkPolyDataHarmonicMappingFilter
{
public:
  static vtkvmtkPolyDataCylinderHarmonicMappingFilter* New();
  vtkTypeRevisionMacro(vtkvmtkPolyDataCylinderHarmonicMappingFilter,vtkvmtkPolyDataHarmonicMappingFilter);

protected:
  vtkvmtkPolyDataCylinderHarmonicMappingFilter();
  ~vtkvmtkPolyDataCylinderHarmonicMappingFilter();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkvmtkPolyDataCylinderHarmonicMappingFilter(const vtkvmtkPolyDataCylinderHarmonicMappingFilter&);  // Not implemented.
  void operator=(const vtkvmtkPolyDataCylinderHarmonicMappingFilter&);  // Not implemented.
};

#endif

