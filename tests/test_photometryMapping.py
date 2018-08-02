import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.geom
import lsst.jointcal.photometryMappings
import lsst.jointcal.photometryTransfo
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
        self.instFlux = 5.0
        self.instFluxErr = 2.0

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
        transfo = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(self.scale)
        self.mapping = lsst.jointcal.photometryMappings.PhotometryMapping(transfo)

    def test_getNpar(self):
        result = self.mapping.getNpar()
        self.assertEqual(result, 1)

    def _test_offsetParams(self, delta, expect):
        self.mapping.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.mapping.getTransfo().getParameters())

    def test_transformFlux(self):
        result = self.mapping.transform(self.star0, self.instFlux)
        self.assertEqual(result, self.instFlux*self.scale)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.array([0.0])
        self._test_offsetParams(delta, np.array([self.scale]))
        delta -= 1
        self._test_offsetParams(delta, self.scale-delta)

    def test_computeParameterDerivatives(self):
        """Test that the derivative of a spatially invariant transform is always the same."""
        result = self.mapping.computeParameterDerivatives(self.star0, self.instFlux)
        self.assertEqual(self.instFlux, result)
        result = self.mapping.computeParameterDerivatives(self.star1, self.instFlux)
        self.assertEqual(self.instFlux, result)
        transfo = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(1000.0)
        mapping = lsst.jointcal.photometryMappings.PhotometryMapping(transfo)
        result = mapping.computeParameterDerivatives(self.star0, self.instFlux)
        self.assertEqual(self.instFlux, result)

    def test_getMappingIndices(self):
        """A mapping with one invariant transfo has one index"""
        self.mapping.setIndex(5)
        result = self.mapping.getMappingIndices()
        self.assertEqual(result, [5])


class ChipVisitPhotometryMappingTestCase(PhotometryMappingTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(ChipVisitPhotometryMappingTestCase, self).setUp()
        self.bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(-5, -6), lsst.afw.geom.Point2D(7, 8),
                                        invert=False)
        self.order = 1
        self.coefficients = np.array([[5, 2], [3, 0]], dtype=float)
        self.chipScale = 2
        self.visitScale = 3
        self.chipIndex = 5
        self.visitIndex = 1000
        chipTransfo = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(self.chipScale)
        chipMapping = lsst.jointcal.photometryMappings.PhotometryMapping(chipTransfo)
        chipMapping.setIndex(self.chipIndex)
        visitTransfo = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(self.visitScale)
        visitMapping = lsst.jointcal.photometryMappings.PhotometryMapping(visitTransfo)
        visitMapping.setIndex(self.visitIndex)
        self.mappingInvariants = lsst.jointcal.photometryMappings.ChipVisitPhotometryMapping(chipMapping,
                                                                                             visitMapping)
        self.mappingInvariants.setWhatToFit(True, True)  # default to fitting both

        # Have to make a new chipMapping, as it stores shared_ptr to the transfo.
        chipTransfo = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(self.chipScale)
        chipMapping = lsst.jointcal.photometryMappings.PhotometryMapping(chipTransfo)
        chipMapping.setIndex(self.chipIndex)
        visitTransfo2 = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.coefficients,
                                                                                   self.bbox)
        visitMapping2 = lsst.jointcal.photometryMappings.PhotometryMapping(visitTransfo2)
        visitMapping2.setIndex(self.visitIndex)
        self.mappingCheby = lsst.jointcal.photometryMappings.ChipVisitPhotometryMapping(chipMapping,
                                                                                        visitMapping2)
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

    def test_transformFlux(self):
        result = self.mappingInvariants.transform(self.star0, self.instFlux)
        self.assertEqual(result, self.instFlux*self.chipScale*self.visitScale)

        result = self.mappingCheby.transform(self.star0, self.instFlux)
        expect = self.instFlux*self.chipScale*self._evaluate_chebyshev(self.star0.getXFocal(),
                                                                       self.star0.getYFocal())
        self.assertEqual(result, expect)

        result = self.mappingCheby.transform(self.star1, self.instFlux)
        expect = self.instFlux*self.chipScale*self._evaluate_chebyshev(self.star1.getXFocal(),
                                                                       self.star1.getYFocal())
        self.assertEqual(result, expect)

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
                expect.append(Ty[j]*Tx[i]*self.instFlux*self.chipScale)
        return expect

    def _computeChipDerivative(self, star):
        """Return the derivative w.r.t. the chip component."""
        return self.instFlux*self._evaluate_chebyshev(star.getXFocal(), star.getYFocal())

    def test_computeParameterDerivatives(self):
        result = self.mappingInvariants.computeParameterDerivatives(self.star0, self.instFlux)
        expect = np.array([self.instFlux*self.visitScale, self.instFlux*self.chipScale])
        self.assertFloatsAlmostEqual(result, expect)

        expect = [self._computeChipDerivative(self.star1), ]
        expect.extend(self._computeChebyshevDerivative(self.star1))
        result = self.mappingCheby.computeParameterDerivatives(self.star1, self.instFlux)
        self.assertFloatsAlmostEqual(np.array(expect), result)

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
        result = self.mappingCheby.computeParameterDerivatives(self.star1, self.instFlux)
        self.assertFloatsEqual(result, derivatives)

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
                                [self._computeChipDerivative(self.star1)])

        # fit just visits means 3 parameters (order 1) and 3 indices starting at self.visitIndex
        self._test_setWhatToFit(False, True, 3, list(range(self.visitIndex, self.visitIndex+3)),
                                [self._computeChebyshevDerivative(self.star1)])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
