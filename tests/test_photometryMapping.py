# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

import abc
import unittest
import lsst.utils.tests

import lsst.afw.geom
import lsst.jointcal.photometryMappings
import lsst.jointcal.photometryTransform
import lsst.jointcal.star


CHEBYSHEV_T = [
    lambda x: 1,
    lambda x: x,
    lambda x: 2*x**2 - 1,
    lambda x: (4*x**2 - 3)*x,
    lambda x: (8*x**2 - 8)*x**2 + 1,
    lambda x: ((16*x**2 - 20)*x**2 + 5)*x,
]


class PhotometryMappingTestBase:
    def setUp(self):
        self.value = 5.0
        self.valueErr = 2.0

        baseStar0 = lsst.jointcal.star.BaseStar(0, 0, 1, 2)
        self.star0 = lsst.jointcal.star.MeasuredStar(baseStar0)
        baseStar1 = lsst.jointcal.star.BaseStar(1, 2, 3, 4)
        self.star1 = lsst.jointcal.star.MeasuredStar(baseStar1)
        self.star1.setXFocal(2)
        self.star1.setYFocal(3)


class PhotometryMappingTestCase(PhotometryMappingTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(PhotometryMappingTestCase, self).setUp()
        self.scale = 3
        transform = lsst.jointcal.photometryTransform.FluxTransformSpatiallyInvariant(self.scale)
        self.mapping = lsst.jointcal.photometryMappings.PhotometryMapping(transform)

    def test_getNpar(self):
        result = self.mapping.getNpar()
        self.assertEqual(result, 1)

    def _test_offsetParams(self, delta, expect):
        self.mapping.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.mapping.getTransform().getParameters())

    def test_transform(self):
        result = self.mapping.transform(self.star0, self.value)
        self.assertEqual(result, self.value*self.scale)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.array([0.0])
        self._test_offsetParams(delta, np.array([self.scale]))
        delta -= 1
        self._test_offsetParams(delta, self.scale-delta)

    def test_computeParameterDerivatives(self):
        """Test that the derivative of a spatially invariant transform is always the same."""
        result = self.mapping.computeParameterDerivatives(self.star0, self.value)
        self.assertEqual(self.value, result)
        result = self.mapping.computeParameterDerivatives(self.star1, self.value)
        self.assertEqual(self.value, result)
        transform = lsst.jointcal.FluxTransformSpatiallyInvariant(1000.0)
        mapping = lsst.jointcal.PhotometryMapping(transform)
        result = mapping.computeParameterDerivatives(self.star0, self.value)
        self.assertEqual(self.value, result)

    def test_getMappingIndices(self):
        """A mapping with one invariant transform has one index"""
        self.mapping.setIndex(5)
        result = self.mapping.getMappingIndices()
        self.assertEqual(result, [5])


class ChipVisitPhotometryMappingTestCase(PhotometryMappingTestBase, abc.ABC):
    def setUp(self):
        super().setUp()
        self.bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(-5, -6), lsst.afw.geom.Point2D(7, 8))
        self.order = 1
        self.coefficients = np.array([[5, 2], [3, 0]], dtype=float)
        self.chipScale = 2
        self.visitScale = 3
        self.chipIndex = 5
        self.visitIndex = 1000

    def _initMappings(self, InvariantTransform, ChebyTransform, ChipVisitMapping):
        """Initialize self.mappingInvariants and self.mappingCheby.
        Call after setUp().

        Parameters
        ----------
        InvariantTransform : `PhotometryTransformSpatiallyInvariant`-type
            The PhotometryTransformSpatiallyInvariant-derived class to construct
            invariant transforms for.
        ChebyTransform : `PhotometryTransform`-type
            The PhotometryTransformChebyshev-derived class to construct
            2d transforms for.
        ChipVisitMapping : `PhotometryMapping`-type
            The PhotometryMapping-derived class to construct for both mappings.
        """
        # self.mappingInvariants has two trivial transforms in it, to serve
        # as a simpler test of functionality.
        chipTransform = InvariantTransform(self.chipScale)
        chipMapping = lsst.jointcal.PhotometryMapping(chipTransform)
        chipMapping.setIndex(self.chipIndex)
        visitTransform = InvariantTransform(self.visitScale)
        visitMapping = lsst.jointcal.PhotometryMapping(visitTransform)
        visitMapping.setIndex(self.visitIndex)
        self.mappingInvariants = ChipVisitMapping(chipMapping, visitMapping)
        self.mappingInvariants.setWhatToFit(True, True)  # default to fitting both

        # self.mappingCheby is a more realistic mapping, with two components:
        # spatially-invariant per chip and a chebyshev per visit.
        # Need a new chipMapping, as it stores shared_ptr to the transform.
        chipTransform = InvariantTransform(self.chipScale)
        chipMapping = lsst.jointcal.PhotometryMapping(chipTransform)
        chipMapping.setIndex(self.chipIndex)
        visitTransform2 = ChebyTransform(self.coefficients, self.bbox)
        visitMapping2 = lsst.jointcal.PhotometryMapping(visitTransform2)
        visitMapping2.setIndex(self.visitIndex)
        self.mappingCheby = ChipVisitMapping(chipMapping, visitMapping2)
        self.mappingCheby.setWhatToFit(True, True)  # default to fitting both

    def test_getNpar(self):
        result = self.mappingInvariants.getNpar()
        self.assertEqual(result, 2)
        # order 1 implies 3 parameters, plus one for the chip mapping
        result = self.mappingCheby.getNpar()
        self.assertEqual(result, 4)

    def _evaluate_chebyshev(self, x, y):
        """Evaluate the chebyshev defined by self.coefficients at (x,y)"""
        # sx, sy: transform from self.bbox range to [-1, -1]
        cx = (self.bbox.getMinX() + self.bbox.getMaxX())/2.0
        cy = (self.bbox.getMinY() + self.bbox.getMaxY())/2.0
        sx = 2.0 / self.bbox.getWidth()
        sy = 2.0 / self.bbox.getHeight()
        result = 0
        for j in range(self.order+1):
            Ty = CHEBYSHEV_T[j](sy*(y - cy))
            for i in range(0, self.order-j+1):
                Tx = CHEBYSHEV_T[i](sx*(x - cx))
                result += self.coefficients[j, i]*Tx*Ty
        return result

    def _computeChebyshevDerivative(self, star):
        """Return the derivatives w.r.t. the Chebyshev components."""
        cx = (self.bbox.getMinX() + self.bbox.getMaxX())/2.0
        cy = (self.bbox.getMinY() + self.bbox.getMaxY())/2.0
        sx = 2.0 / self.bbox.getWidth()
        sy = 2.0 / self.bbox.getHeight()
        Tx = np.array([CHEBYSHEV_T[i](sx*(star.getXFocal() - cx))
                      for i in range(self.order+1)], dtype=float)
        Ty = np.array([CHEBYSHEV_T[i](sy*(star.getYFocal() - cy))
                      for i in range(self.order+1)], dtype=float)
        expect = []
        for j in range(len(Ty)):
            for i in range(0, self.order-j+1):
                expect.append(Ty[j]*Tx[i])
        return np.array(expect)

    @abc.abstractmethod
    def _computeVisitDerivative(self, star):
        """Return the derivative w.r.t. the chebyshev visit component."""
        pass

    @abc.abstractmethod
    def _computeChipDerivative(self, star):
        """Return the derivative w.r.t. the chip component."""
        pass

    def test_getMappingIndices(self):
        """There are npar indices in a constrained mapping."""
        expect = [self.chipIndex, self.visitIndex]
        result = self.mappingInvariants.getMappingIndices()
        self.assertEqual(result, expect)

        # npar - 1 because the chip mapping has the 1st parameter
        expect = [self.chipIndex, ] + list(range(self.visitIndex,
                                                 self.visitIndex + self.mappingCheby.getNpar() - 1))
        result = self.mappingCheby.getMappingIndices()
        self.assertEqual(result, expect)

    def _test_transform_mappingInvariants(self, star, expect):
        result = self.mappingInvariants.transform(star, self.value)
        self.assertEqual(result, expect)

    def _test_transform_mappingCheby(self, star, expect):
        result = self.mappingCheby.transform(star, self.value)
        self.assertEqual(result, expect)

    def _test_computeParameterDerivatives(self, star, expectInvariant):
        """Test self.mappingInvariants and self.mappingCheby transforming star.
        expectCheby is calculated from _computeChipDerivative and
        _computeChebyshevDerivative.
        """
        result = self.mappingInvariants.computeParameterDerivatives(star, self.value)
        self.assertFloatsAlmostEqual(result, expectInvariant)

        # the chip derivative is a single number
        expectCheby = [self._computeChipDerivative(self.star1)]
        # the Chebyshev Derivatives are a list, so we have to use extend
        expectCheby.extend(self._computeVisitDerivative(self.star1))
        expectCheby = np.array(expectCheby)
        result = self.mappingCheby.computeParameterDerivatives(star, self.value)
        self.assertFloatsAlmostEqual(result, expectCheby)

    def _test_setWhatToFit(self, fittingChips, fittingVisits, nPar, indices, derivatives):
        """
        Parameters
        ----------
        fittingChips : `bool`
            Are we fitting the chip component?
            Passed to ``self.mappingCheby.setWhatToFit()``.
        fittingVisits : `bool`
            Are we fitting the visit component?
            Passed to ``self.mappingCheby.setWhatToFit()``.
        nPar : `int`
            Expected result from ``self.mappingCheby.getNpar()``.
        indices : `list`
            Expected result from ``self.mappingCheby.getMappingIndices()``.
        derivatives : `list`
            Expected result from ``self.mappingCheby.computeParameterDerivatives()``.
        """
        self.mappingCheby.setWhatToFit(fittingChips, fittingVisits)
        self.assertEqual(self.mappingCheby.getNpar(), nPar)
        self.assertEqual(self.mappingCheby.getMappingIndices(), indices)
        result = self.mappingCheby.computeParameterDerivatives(self.star1, self.value)
        self.assertFloatsAlmostEqual(result, derivatives)

    def test_setWhatToFit(self):
        """Test that mapping methods behave correctly when chip and/or visit
        fitting is disabled.

        The "fit both" case (True, True) is tested by all of the above tests.
        """
        # Using mappingCheby so getNpar() will distinguish chips (1 param) from visits (3 params).

        # fit nothing means 0 parameters and no indices
        self._test_setWhatToFit(False, False, 0, [], [])

        # fit just chips means 1 parameter and one index [self.chipIndex]
        self._test_setWhatToFit(True, False, 1, [self.chipIndex],
                                np.array([self._computeChipDerivative(self.star1)]))

        # fit just visits means 3 parameters (order 1) and 3 indices starting at self.visitIndex
        self._test_setWhatToFit(False, True, 3, list(range(self.visitIndex, self.visitIndex+3)),
                                np.array([self._computeVisitDerivative(self.star1)]))


class ChipVisitFluxMappingTestCase(ChipVisitPhotometryMappingTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self._initMappings(lsst.jointcal.FluxTransformSpatiallyInvariant,
                           lsst.jointcal.FluxTransformChebyshev,
                           lsst.jointcal.ChipVisitFluxMapping)

    def _computeVisitDerivative(self, star):
        return self._computeChebyshevDerivative(star) * self.value * self.chipScale

    def _computeChipDerivative(self, star):
        return self.value * self._evaluate_chebyshev(star.getXFocal(), star.getYFocal())

    def test_transform(self):
        expect = self.value * self.chipScale * self.visitScale
        self._test_transform_mappingInvariants(self.star0, expect)
        # The doubly-spatially invariant mapping should be independent of star position.
        self._test_transform_mappingInvariants(self.star1, expect)

        expect = self.value * self.chipScale * self._evaluate_chebyshev(self.star0.getXFocal(),
                                                                        self.star0.getYFocal())
        self._test_transform_mappingCheby(self.star0, expect)
        expect = self.value * self.chipScale * self._evaluate_chebyshev(self.star1.getXFocal(),
                                                                        self.star1.getYFocal())
        self._test_transform_mappingCheby(self.star1, expect)

    def test_computeParameterDerivatives(self):
        expectInvariant = np.array([self.value*self.visitScale, self.value*self.chipScale])
        self._test_computeParameterDerivatives(self.star1, expectInvariant)


class ChipVisitMagnitudeMappingTestCase(ChipVisitPhotometryMappingTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self._initMappings(lsst.jointcal.MagnitudeTransformSpatiallyInvariant,
                           lsst.jointcal.MagnitudeTransformChebyshev,
                           lsst.jointcal.ChipVisitMagnitudeMapping)

    def _computeVisitDerivative(self, star):
        return self._computeChebyshevDerivative(star)

    def _computeChipDerivative(self, star):
        # Magnitude chip derivative is always identically 1:
        #     d(M(m))/d(m0)=1 where M(m) = m + m0
        return 1.0

    def test_transform(self):
        expect = self.value + self.chipScale + self.visitScale
        self._test_transform_mappingInvariants(self.star0, expect)
        # The doubly-spatially invariant mapping should be independent of star position.
        self._test_transform_mappingInvariants(self.star1, expect)

        expect = self.value + self.chipScale + self._evaluate_chebyshev(self.star0.getXFocal(),
                                                                        self.star0.getYFocal())
        self._test_transform_mappingCheby(self.star0, expect)

        expect = self.value + self.chipScale + self._evaluate_chebyshev(self.star1.getXFocal(),
                                                                        self.star1.getYFocal())
        self._test_transform_mappingCheby(self.star1, expect)

    def test_computeParameterDerivatives(self):
        # the parameter derivative of a spatially invariant magnitude transform is always 1.
        expectInvariant = np.array([1.0, 1.0])
        self._test_computeParameterDerivatives(self.star1, expectInvariant)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
