// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
// This file is part of RBDyn.
//
// RBDyn is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RBDyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with RBDyn.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// includes
// std
#include <vector>

// RBDyn
#include "Joint.h"

#include <rbdyn/config.hh>

namespace rbd
{
class MultiBody;
struct MultiBodyConfig;

/**
	* Integrate joint configuration using a second-order "leapfrog" method.
	* @param type Joint type.
	* @param alpha Joint velocity vector.
	* @param step Integration step.
	* @param q Joint configuration vector.
	*/
class RBDYN_DLLAPI LeapFrogIntegration
{
public:
	LeapFrogIntegration(const MultiBody& mb);
/**
	* Use the leapfrog method to integrate.
	* @param mb MultiBody used has model.
	* @param mbc Use alphaD, alpha and q. Fill alpha and q.
	* @param step Integration step.
	*/
RBDYN_DLLAPI void leapfrogIntegration(const MultiBody& mb, MultiBodyConfig& mbc, double step);

/// safe version of @see eulerIntegration.
RBDYN_DLLAPI void sLeapfrogIntegration(const MultiBody& mb, MultiBodyConfig& mbc, double step);

private:
	void jointIntegration(Joint::Type type, std::vector<double>& alpha,
	const std::vector<double>& alphaD,
	std::vector<double>& alphaDi,
	double step, std::vector<double>& q);

	std::vector<std::vector<double>> alphaDi_;
};

} // namepsace rbd
