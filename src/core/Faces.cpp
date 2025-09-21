//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-04 22:09:56 gtaubin>
//------------------------------------------------------------------------
//
// Faces.cpp
//
// Written by: <Your Name>
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <algorithm>
#include <math.h>
#include "Faces.hpp"

Faces::Faces(const int nV, const vector<int>& coordIndex) {
	_nV = nV;

	_coordIndexReference = &coordIndex;

	_faceIndex.push_back(0);
	for (int idx = 0; idx < coordIndex.size(); idx++) {
		if (coordIndex[idx] == -1)
			_faceIndex.push_back(idx+1);
	}
}

int Faces::getNumberOfVertices() const {
	return _nV;
}

int Faces::getNumberOfFaces() const {
	return _faceIndex.size()-1;
}

int Faces::getNumberOfCorners() const {
	return _coordIndexReference->size();
}

int Faces::getFaceSize(const int iF) const {
	if (iF == _faceIndex.size()-1) return 0;
	return _faceIndex[iF+1] - _faceIndex[iF] - 1;
}

int Faces::getFaceFirstCorner(const int iF) const {
	if (iF >= _faceIndex.size()-1) return -1;
	return _faceIndex[iF];
}

int Faces::getFaceVertex(const int iF, const int j) const {
	return (*_coordIndexReference)[_faceIndex[iF]+j];
}

int Faces::getCornerFace(const int iC) const {
	if ((*_coordIndexReference)[iC] == -1) return -1;
	return _faceIndex.end() - lower_bound(_faceIndex.begin(), _faceIndex.end(), iC);
}

int Faces::getNextCorner(const int iC) const {
	if ((*_coordIndexReference)[iC] == -1) return -1;
	if ((*_coordIndexReference)[iC+1] != -1) return (*_coordIndexReference)[iC+1];
	return getFaceFirstCorner(getCornerFace(iC));
}

