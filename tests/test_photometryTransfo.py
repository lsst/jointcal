import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.geom
import lsst.jointcal.photometryTransfo


CHEBYSHEV_T = [
    lambda x: 1,
    lambda x: x,
    lambda x: 2*x**2 - 1,
    lambda x: (4*x**2 - 3)*x,
    lambda x: (8*x**2 - 8)*x**2 + 1,
    lambda x: ((16*x**2 - 20)*x**2 + 5)*x,
]


class PhotometryTransfoTestBase:
    def setUp(self):
        self.value = 5.0
        self.point = [1., 5.]


class SpatiallyInvariantTestBase(PhotometryTransfoTestBase):
    """Tests for PhotometryTransfoSpatiallyInvariant.
     Subclasses need to call setUp to define:
         self.transfo1 == a default initalized PhotometryTransfoSpatiallyInvariant.
         self.transfo2 == a transfo initialized with self.t2init.
    """
    def setUp(self):
        super().setUp()
        self.t2init = 1000.0

    def _test_transform(self, transfo, expect):
        result = transfo.transform(self.point[0], self.point[1], self.value)
        self.assertEqual(result, expect)

    def _offsetParams(self, delta, value, expect):
        self.transfo1.offsetParams(delta)
        result = self.transfo1.transform(self.point[0], self.point[1], value)
        self.assertFloatsAlmostEqual(result, expect)

    def _test_offsetParams(self, expect):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        # check that offset by 0 doesn't change anything.
        delta = np.zeros(1, dtype=float)
        self._offsetParams(delta, self.value, self.value)

        # offset by +1 should result in `expect`
        delta -= 1
        self._offsetParams(delta, self.value, expect)

    def test_clone(self):
        clone1 = self.transfo1.clone()
        self.assertEqual(self.transfo1.getParameters(), clone1.getParameters())
        clone2 = self.transfo2.clone()
        self.assertEqual(self.transfo2.getParameters(), clone2.getParameters())
        self.assertNotEqual(clone1.getParameters(), clone2.getParameters())

    def _test_computeParameterDerivatives(self, expect):
        """The derivative of a spatially invariant transform is always the same.
        Should be indepdendent of position
        """
        result = self.transfo1.computeParameterDerivatives(1, 2, self.value)
        self.assertEqual(expect, result)
        result = self.transfo1.computeParameterDerivatives(-5, -100, self.value)
        self.assertEqual(expect, result)
        result = self.transfo2.computeParameterDerivatives(-1000, 150, self.value)
        self.assertEqual(expect, result)


class FluxTransfoSpatiallyInvariantTestCase(SpatiallyInvariantTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.transfo1 = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant()
        self.transfo2 = lsst.jointcal.photometryTransfo.FluxTransfoSpatiallyInvariant(self.t2init)

    def test_transform(self):
        self._test_transform(self.transfo1, self.value)
        self._test_transform(self.transfo2, self.value*self.t2init)

    def test_offsetParams(self):
        self._test_offsetParams(self.value*2)

    def test_computeParameterDerivatives(self):
        """Should be indepdendent of position, and equal to the flux."""
        self._test_computeParameterDerivatives(self.value)


class MagnitudeTransfoSpatiallyInvariantTestCase(SpatiallyInvariantTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.transfo1 = lsst.jointcal.photometryTransfo.MagnitudeTransfoSpatiallyInvariant()
        self.transfo2 = lsst.jointcal.photometryTransfo.MagnitudeTransfoSpatiallyInvariant(self.t2init)

    def test_transform(self):
        self._test_transform(self.transfo1, self.value)
        self._test_transform(self.transfo2, self.value + self.t2init)

    def test_offsetParams(self):
        self._test_offsetParams(self.value + 1)

    def test_computeParameterDerivatives(self):
        """Should always be identically 1."""
        self._test_computeParameterDerivatives(1.0)


class PhotometryTransfoChebyshevTestCase(PhotometryTransfoTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(-5, -6), lsst.afw.geom.Point2D(7, 8))
        self.order1 = 2
        self.transfo1 = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.order1, self.bbox)
        self.order2 = 1
        self.coefficients = np.array([[5, 3], [4, 0]], dtype=float)
        self.transfo2 = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.coefficients,
                                                                                   self.bbox)

    def test_getNpar(self):
        self.assertEqual(self.transfo1.getNpar(), 6)
        self.assertEqual(self.transfo2.getNpar(), 3)

    def _evaluate_chebyshev(self, x, y):
        """Evaluate the chebyshev defined by self.coefficients at (x,y)"""
        # sx, sy: transform from self.bbox range to [-1, -1]
        cx = (self.bbox.getMinX() + self.bbox.getMaxX())/2.0
        cy = (self.bbox.getMinY() + self.bbox.getMaxY())/2.0
        sx = 2.0 / self.bbox.getWidth()
        sy = 2.0 / self.bbox.getHeight()
        result = 0
        order = len(self.coefficients)
        for j in range(order):
            for i in range(0, order-j):
                Tx = CHEBYSHEV_T[i](sx*(x - cx))
                Ty = CHEBYSHEV_T[j](sy*(y - cy))
                result += self.coefficients[j, i]*Tx*Ty
        return result

    def test_transform(self):
        result = self.transfo1.transform(self.point[0], self.point[1], self.value)
        self.assertEqual(result, self.value)  # transfo1 is the identity

        result = self.transfo2.transform(self.point[0], self.point[1], self.value)
        expect = self.value*self._evaluate_chebyshev(self.point[0], self.point[1])
        self.assertEqual(result, expect)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(self.transfo1.getNpar(), dtype=float)
        expect = np.zeros((self.order1+1, self.order1+1), dtype=float)
        expect[0, 0] = 1
        self.transfo1.offsetParams(delta)
        # nothing should have changed if we transform a delta of 0
        self.assertFloatsAlmostEqual(expect, self.transfo1.getCoefficients())

        delta[0] = 1
        delta[1] = -2
        delta[2] = -3
        delta[3] = -4
        delta[4] = -5
        delta[5] = -6
        expect[0, 0] = 0
        expect[0, 1] = 2
        expect[0, 2] = 3
        expect[1, 0] = 4
        expect[1, 1] = 5
        expect[2, 0] = 6
        self.transfo1.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.transfo1.getCoefficients())

    def test_clone(self):
        clone1 = self.transfo1.clone()
        self.assertFloatsEqual(self.transfo1.getParameters(), clone1.getParameters())
        self.assertEqual(self.transfo1.getOrder(), clone1.getOrder())
        self.assertEqual(self.transfo1.getBBox(), clone1.getBBox())
        clone2 = self.transfo2.clone()
        self.assertFloatsEqual(self.transfo2.getParameters(), clone2.getParameters())
        self.assertEqual(self.transfo2.getOrder(), clone2.getOrder())
        self.assertEqual(self.transfo2.getBBox(), clone2.getBBox())

    def test_computeParameterDerivatives(self):
        cx = (self.bbox.getMinX() + self.bbox.getMaxX())/2.0
        cy = (self.bbox.getMinY() + self.bbox.getMaxY())/2.0
        sx = 2.0 / self.bbox.getWidth()
        sy = 2.0 / self.bbox.getHeight()
        result = self.transfo1.computeParameterDerivatives(self.point[0], self.point[1], self.value)
        Tx = np.array([CHEBYSHEV_T[i](sx*(self.point[0] - cx)) for i in range(self.order1+1)], dtype=float)
        Ty = np.array([CHEBYSHEV_T[i](sy*(self.point[1] - cy)) for i in range(self.order1+1)], dtype=float)
        expect = []
        for j in range(len(Ty)):
            for i in range(0, self.order1-j+1):
                expect.append(Ty[j]*Tx[i]*self.value)
        self.assertFloatsAlmostEqual(np.array(expect), result)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
