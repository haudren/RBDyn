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

// associated header
#include "LeapFrogIntegration.h"

// includes
// RBDyn
#include "MultiBody.h"
#include "MultiBodyConfig.h"

namespace rbd
{

namespace{

Eigen::Matrix<double, 4, 3> angularVelToQuatVel(const std::vector<double>& q)
{
	return 0.5*(Eigen::Matrix<double, 4, 3>() <<
		-q[1], -q[2], -q[3],
		q[0], -q[3], q[2],
		q[3], q[0], -q[1],
		-q[2], q[1], q[0]).finished();
}

}

LeapFrogIntegration::LeapFrogIntegration(const MultiBody& mb)
{
	alphaDi_.clear();
	alphaDi_.resize(mb.nrJoints());
	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		alphaDi_[i] = mb.joint(i).zeroDof();
	}
}

void LeapFrogIntegration::jointIntegration(Joint::Type type, std::vector<double>& alpha,
	const std::vector<double>& alphaD, std::vector<double>& alphaDi,
	double step, std::vector<double>& q)
{
	switch(type)
	{
		case Joint::Rev:
		case Joint::Prism:
		{
			// alpha_{i+1/2}
			alpha[0] += alphaDi[0]*step/2;
			q[0] += alpha[0]*step;
			alpha[0] += alphaD[0]*step/2;
			break;
		}

		/// @todo manage reverse joint
		case Joint::Planar:
		{
			throw std::runtime_error("Planar joint not implemented yet!");
			double q1Step = q[2]*alpha[0] + alpha[1];
			double q2Step = -q[1]*alpha[0] + alpha[2];
			q[0] += alpha[0]*step;
			q[1] += q1Step*step;
			q[2] += q2Step*step;
			break;
		}

		case Joint::Cylindrical:
		{
			throw std::runtime_error("Cylindrical joint not implemented yet!");
			q[0] += alpha[0]*step;
			q[1] += alpha[1]*step;
			break;
		}

		/// @todo manage reverse joint
		case Joint::Free:
		{
			Eigen::Vector3d v;
			alpha[3] += alphaDi[3]*step/2;
			alpha[4] += alphaDi[4]*step/2;
			alpha[5] += alphaDi[5]*step/2;

			v << alpha[3], alpha[4], alpha[5];
			// v is in FS coordinate. We have to put it back in FP coordinate.
			v = QuatToE(q).transpose()*v;

			q[4] += v[0]*step;
			q[5] += v[1]*step;
			q[6] += v[2]*step;

			alpha[3] += alphaD[3]*step/2;
			alpha[4] += alphaD[4]*step/2;
			alpha[5] += alphaD[5]*step/2;
			// don't break, we go in spherical
		}
		/// @todo manage reverse joint
		case Joint::Spherical:
		{
			Eigen::Vector4d qi;
			Eigen::Vector3d w;
			qi << q[0], q[1], q[2], q[3];

			alpha[0] += alphaDi[0]*step/2;
			alpha[1] += alphaDi[1]*step/2;
			alpha[2] += alphaDi[2]*step/2;

			w << alpha[0], alpha[1], alpha[2];

			qi.noalias() += angularVelToQuatVel(q)*w*step;
			qi.normalize();

			alpha[0] += alphaD[0]*step/2;
			alpha[1] += alphaD[1]*step/2;
			alpha[2] += alphaD[2]*step/2;

			q[0] = qi(0);
			q[1] = qi(1);
			q[2] = qi(2);
			q[3] = qi(3);
			break;
		}

		case Joint::Fixed:
		default:
		;
	}
}

void LeapFrogIntegration::leapfrogIntegration(const MultiBody& mb, MultiBodyConfig& mbc, double step)
{
	const std::vector<Joint>& joints = mb.joints();

	// integrate
	for(std::size_t i = 0; i < joints.size(); ++i)
	{
		jointIntegration(joints[i].type(), mbc.alpha[i], mbc.alphaD[i], alphaDi_[i], step, mbc.q[i]);
		alphaDi_[i] = mbc.alphaD[i];
	}
}

void LeapFrogIntegration::sLeapfrogIntegration(const MultiBody& mb, MultiBodyConfig& mbc, double step)
{
	checkMatchQ(mb, mbc);
	checkMatchAlpha(mb, mbc);
	checkMatchAlphaD(mb, mbc);

	leapfrogIntegration(mb, mbc, step);
}

} // namespace rbd

