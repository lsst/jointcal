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

"""Tests of the AstrometryTransform objects and its helpers."""
import numpy as np

import unittest
import lsst.utils.tests

import lsst.geom
import lsst.log
import lsst.jointcal
from lsst.jointcal.astrometryTransform import AstrometryTransformPolynomial, inversePolyTransform


class AstrometryTransformPolynomialBase:
    def setUp(self):
        self.longMessage = True
        np.random.seed(100)

        # TRACE level so we can see all the inverse iteration steps.
        lsst.log.setLevel('', lsst.log.TRACE)

        # default initialized gtansfoPoly should be the identity
        self.polyIdentity = AstrometryTransformPolynomial(9)

        self.poly2 = AstrometryTransformPolynomial(2)
        # A "reasonable" polynomial from fitting some data.
        poly2Str = "1\norder 2\n-0.167892516757 -5.19165870491e-05 2.0273225968e-07 2.75150063324e-11 -1.41377300894e-11 7.74131250823e-12 -0.355492046755 3.64496426108e-07 5.18028912201e-05 -1.38429955144e-11 -2.23437887691e-11 2.23415522405e-11 \n"  # noqa: E501
        self.poly2.read(poly2Str)

        # NOTE: taken from the "visit" component of a model, so needs to be
        # evaluated on a different bbox.
        self.poly3 = AstrometryTransformPolynomial(3)
        poly3Str = "1\norder 3\n-0.110057982995 -0.00384939195727 2.22344576374e-05 1.01023479109e-09 -5.2469672592e-09 1.50668618207e-09 1.46346365118e-09 -1.63797654003e-10 1.27721649626e-09 -2.05428888876e-12 -0.115317494814 2.21413899258e-05 0.00385039699775 -5.72570678996e-09 -9.59649028554e-09 1.41401351683e-08 -1.97886824954e-10 -1.30841763583e-09 -1.95789465504e-11 -1.26379172569e-09 "  # noqa: E501
        self.poly3.read(poly3Str)

        self.poly9 = AstrometryTransformPolynomial(9)
        # A "reasonable" polynomial from fitting some data.
        poly9Str = "1\norder 9\n-0.16793632503 -5.19113500529e-05 4.6484637807e-07 1.41209689574e-09 -8.13063885535e-10 -5.52862561337e-10 -7.20917258572e-12 1.44268609502e-12 1.01158414497e-12 5.81014911789e-13 1.72417265293e-14 -1.39482997394e-15 -1.34837756422e-15 -5.5315543147e-16 -3.26263945732e-16 -2.27835189163e-17 1.54838068036e-19 1.36520025882e-18 4.39218993176e-19 1.09711054987e-19 1.07909937222e-19 1.75422369904e-20 9.46438635874e-22 -7.47653768123e-22 -4.21000295132e-22 -1.63261283932e-23 2.18164463104e-23 -2.29125892915e-23 -7.79605618485e-24 -8.96899627417e-25 1.8328850322e-25 2.429872115e-25 1.57383988209e-26 -1.25289729437e-26 -1.60798488421e-26 3.50802217715e-27 1.84635563103e-27 3.27562810981e-28 4.22037315014e-30 -5.54664361872e-29 -2.8639429519e-29 1.1383384424e-29 1.04059046024e-30 3.18765755999e-30 -3.9564850064e-31 -1.8028544084e-31 -3.82041218173e-32 -1.70229355716e-32 1.32310477371e-32 9.28363194124e-35 1.5990001923e-33 -1.31379494936e-33 6.88708291597e-35 -2.23831688762e-34 2.33466866348e-35 -0.355489287839 1.33735110546e-06 5.17833395824e-05 -6.95722643803e-09 -1.04122329759e-09 6.40573984039e-12 2.20559605231e-11 5.41569680032e-12 4.5540873349e-13 1.26533773481e-13 -3.78843184326e-14 -1.40950401053e-14 -1.79770054411e-15 -8.32768571948e-18 -2.33287234831e-16 3.78385302886e-17 2.09977398841e-17 2.56668227438e-18 6.94279232614e-19 -1.62606627953e-19 2.04795207418e-19 -2.24588539576e-20 -1.79824251677e-20 -2.83870287858e-21 -3.26896753863e-22 -2.59328106552e-22 1.11975449899e-22 -9.70654488549e-23 7.70881954106e-24 8.70008772668e-24 2.08267357514e-24 -3.206329176e-26 1.25458989726e-25 5.34275512133e-26 -3.49011241433e-26 2.55423596751e-26 -1.37582742645e-27 -2.21887269529e-27 -7.87749245556e-28 4.74161700401e-29 -2.40765221773e-29 -1.32783315748e-29 -7.20736097924e-30 5.46646062406e-30 -3.51762994693e-30 9.24376571391e-32 2.36641447777e-31 1.05246663621e-31 7.65324655394e-33 -1.07479816221e-32 6.75168688604e-33 -7.73208895042e-34 6.25567958067e-34 -3.54400370636e-34 1.97908370676e-34 \n"  # noqa: E501
        self.poly9.read(poly9Str)

        # make a grid of points to evaluate on
        self._makePoints(0, 1000, 0, 1000, 200)

    def _makePoints(self, minX, maxX, minY, maxY, num):
        """Sets self.points to a 2d grid of num points from min->max."""
        self.frame = lsst.jointcal.frame.Frame(lsst.jointcal.star.Point(minX, minY),
                                               lsst.jointcal.star.Point(maxX, maxY))
        num = 200
        xx = np.linspace(minX, maxX, num)
        yy = np.linspace(minY, maxY, num)
        self.points = []
        for x in xx:
            for y in yy:
                self.points.append(lsst.geom.Point2D(x, y))


class InversePolyTransformTestCase(AstrometryTransformPolynomialBase, lsst.utils.tests.TestCase):
    def checkInverse(self, poly, inverse, maxDiff):
        """Test that ``astrometryTransformPolynomial(inverse(point))==point`` to within maxDiff."""
        results = []
        for point in self.points:
            # TODO: Fix these "Point"s once DM-4044 is done.
            tempPoint = lsst.jointcal.star.Point(point[0], point[1])
            result = inverse.apply(poly.apply(tempPoint))
            results.append(lsst.geom.Point2D(result.x, result.y))

        self.assertPairListsAlmostEqual(results, self.points, maxDiff=maxDiff)

    def testInversePolyIdentity(self):
        precision = 1e-8
        inverse = inversePolyTransform(self.polyIdentity, self.frame, precision)
        self.checkInverse(self.polyIdentity, inverse, precision)

    def testInversePoly2(self):
        precision = 1e-6
        inverse = inversePolyTransform(self.poly2, self.frame, precision)
        self.checkInverse(self.poly2, inverse, 1e-7)

    def testInversePoly3(self):
        # Different bbox for this one, because it is a focal plane to tangent plane transform.
        minX = 14.6927
        maxX = 42.38
        minY = -62.5486
        maxY = -0.323848
        self._makePoints(minX, maxX, minY, maxY, 200)

        precision = 1e-7
        inverse = inversePolyTransform(self.poly3, self.frame, precision, maxOrder=5)
        self.checkInverse(self.poly3, inverse, 3e-8)

    def testInversePoly9(self):
        precision = 1e-6
        inverse = inversePolyTransform(self.poly9, self.frame, precision,
                                       maxOrder=11, nSteps=100)
        self.checkInverse(self.poly9, inverse, 4e-5)

    def testNotEnoughPoints(self):
        with self.assertRaises(RuntimeError):
            inversePolyTransform(self.poly2, self.frame, 1e-4, nSteps=2)


class AstrometryTransformPolynomialTestCase(AstrometryTransformPolynomialBase, lsst.utils.tests.TestCase):
    def checkToAstMap(self, poly, inverseMaxDiff=1e-6):
        """Test that AstrometryTransformPolynomial.toAstMap() gives accurate results.

        Parameters
        ----------
        poly : `lsst.jointcal.astrometryTransform.AstrometryTransformPolynomial`
            The polynomial to be tested.
        inverseMaxDiff : `float`
            Required accuracy on inverse polynomial.
            See `lsst.afw.geom.utils.assertPairsAlmostEqual`.
        """
        # maxDiff should be small, because this should be a near-exact conversion,
        # modulo implementation details.
        maxDiff = 1e-10

        astMap = poly.toAstMap(self.frame)
        expects = []
        forwards = []
        inverses = []
        for point in self.points:
            # TODO: Fix these "Point"s once DM-4044 is done.
            tempPoint = lsst.jointcal.star.Point(point[0], point[1])
            expect = poly.apply(tempPoint)
            expects.append(lsst.geom.Point2D(expect.x, expect.y))
            result = astMap.applyForward(point)
            forwards.append(result)
            inverses.append(astMap.applyInverse(result))

        self.assertPairListsAlmostEqual(forwards, expects, maxDiff=maxDiff)
        self.assertPairListsAlmostEqual(inverses, self.points, maxDiff=inverseMaxDiff)

    def testToAstMapIdentity(self):
        self.checkToAstMap(self.polyIdentity)

    def testToAstMapOrder2(self):
        self.checkToAstMap(self.poly2)

    def testToAstMapOrder3(self):
        # Different bbox for this one, because it was fit on a focal plane.
        minX = 14.6927
        maxX = 42.38
        minY = -62.5486
        maxY = -0.323848
        self._makePoints(minX, maxX, minY, maxY, 200)

        self.checkToAstMap(self.poly3, inverseMaxDiff=5e-6)

    def testToAstMapOrder9(self):
        # looser tolerance: 9th order polynomials are harder to get a good inverse for.
        self.checkToAstMap(self.poly9, inverseMaxDiff=4e-5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
