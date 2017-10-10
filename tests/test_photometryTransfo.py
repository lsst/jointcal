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


class PhotometryTransfoTestBase():
    def setUp(self):
        self.instFlux = 1.0
        self.point = [1., 5.]


class PhotometryTransfoSpatiallyInvariantTestCase(PhotometryTransfoTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(PhotometryTransfoSpatiallyInvariantTestCase, self).setUp()
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant()

    def test_transform(self):
        result = self.transfo.transform(self.point[0], self.point[1], self.instFlux)
        self.assertEqual(result, self.instFlux)

    def _test_offsetParams(self, delta, expect):
        self.transfo.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.transfo.transform(1, 2, 1))

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(1, dtype=float)
        self._test_offsetParams(delta, np.array(1.0))
        delta -= 1
        self._test_offsetParams(delta, np.array(2.0))

    def test_clone(self):
        clone1 = self.transfo.clone()
        self.assertEqual(self.transfo.getParameters(), clone1.getParameters())
        transfo2 = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant(1234)
        clone2 = transfo2.clone()
        self.assertEqual(transfo2.getParameters(), clone2.getParameters())
        self.assertNotEqual(clone1.getParameters(), clone2.getParameters())

    def test_computeParameterDerivatives(self):
        """The derivative of a spatially invariant transform is always the same.
        Should be indepdendent of position, flux, and value.
        """
        result = self.transfo.computeParameterDerivatives(1, 2, self.instFlux)
        self.assertEqual(self.instFlux, result)
        result = self.transfo.computeParameterDerivatives(-5, -100, self.instFlux)
        self.assertEqual(self.instFlux, result)
        transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant(1000.0)
        result = transfo.computeParameterDerivatives(1, 2, self.instFlux)
        self.assertEqual(self.instFlux, result)


class PhotometryTransfoChebyshevTestCase(PhotometryTransfoTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(PhotometryTransfoChebyshevTestCase, self).setUp()
        self.bbox = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(-5, -6), lsst.afw.geom.Point2D(7, 8))
        self.degree1 = 2
        self.transfo1 = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.degree1, self.bbox)
        self.degree2 = 1
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
        degree = len(self.coefficients)
        for j in range(degree):
            for i in range(0, degree-j):
                Tx = CHEBYSHEV_T[i](sx*(x - cx))
                Ty = CHEBYSHEV_T[j](sy*(y - cy))
                result += self.coefficients[j, i]*Tx*Ty
        return result

    def test_transform(self):
        result = self.transfo1.transform(self.point[0], self.point[1], self.instFlux)
        self.assertEqual(result, self.instFlux)  # transfo1 is the identity

        result = self.transfo2.transform(self.point[0], self.point[1], self.instFlux)
        expect = self.instFlux*self._evaluate_chebyshev(self.point[0], self.point[1])
        self.assertEqual(result, expect)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(self.transfo1.getNpar(), dtype=float)
        expect = np.zeros((self.degree1+1, self.degree1+1), dtype=float)
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
        self.assertEqual(self.transfo1.getDegree(), clone1.getDegree())
        self.assertEqual(self.transfo1.getBBox(), clone1.getBBox())
        clone2 = self.transfo2.clone()
        self.assertFloatsEqual(self.transfo2.getParameters(), clone2.getParameters())
        self.assertEqual(self.transfo2.getDegree(), clone2.getDegree())
        self.assertEqual(self.transfo2.getBBox(), clone2.getBBox())

    def test_computeParameterDerivatives(self):
        cx = (self.bbox.getMinX() + self.bbox.getMaxX())/2.0
        cy = (self.bbox.getMinY() + self.bbox.getMaxY())/2.0
        sx = 2.0 / self.bbox.getWidth()
        sy = 2.0 / self.bbox.getHeight()
        result = self.transfo1.computeParameterDerivatives(self.point[0], self.point[1], self.instFlux)
        Tx = np.array([CHEBYSHEV_T[i](sx*(self.point[0] - cx)) for i in range(self.degree1+1)], dtype=float)
        Ty = np.array([CHEBYSHEV_T[i](sy*(self.point[1] - cy)) for i in range(self.degree1+1)], dtype=float)
        expect = []
        for j in range(len(Ty)):
            for i in range(0, self.degree1-j+1):
                expect.append(Ty[j]*Tx[i]*self.instFlux)
        self.assertFloatsAlmostEqual(np.array(expect), result)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
