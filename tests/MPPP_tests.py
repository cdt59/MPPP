import unittest
from src2.ImProcessor import ImProcessor
from src2.Image import Image


class TestImageProcessing(unittest.TestCase):

    def setUp(self):
        """
        To instantiate the ImageProcessor and the Image objects
        """
        self.processor = ImProcessor()
        self.zcam_l = Image(
            ".\\data\\Sol_0709\\ZL0_0709_0729889008_069RAD_N0332864ZCAM07114_1100LMA01.IMG")
        self.zcam_r = Image(
            ".\\data\\Sol_0709\\ZR0_0709_0729888971_069RAD_N0332864ZCAM07114_0340LMA01.IMG")
        self.flf = Image(
            ".\\data\\Sol_0709\\FLF_0711_0730074611_304RAD_N0332864FHAZ00206_0A0295J01.IMG")
        self.frf = Image(
            ".\\data\\Sol_0709\\FRF_0711_0730074611_304RAD_N0332864FHAZ00206_0A0295J01.IMG")
        self.nlf = Image(
            ".\\data\\Sol_0709\\NLF_0709_0729883381_848RAD_N0332864SAPP00601_0A00LLJ01.IMG")
        self.nrf = Image(
            ".\\data\\Sol_0709\\NRF_0709_0729889968_787RAD_N0332864NCAM03709_0A0195J01.IMG")
        self.nlmv = Image(
            ".\\data\\Sol_0709\\NLMV0709_0729882743_168RAD_N0332841TRAV00262_0A02LLJ01.IMG")
        self.nrmv = Image(
            ".\\data\\Sol_0709\\NRMV0709_0729882743_168RAD_N0332841TRAV00262_0A02LLJ01.IMG")

    def test_ftau(self):
        self.assertAlmostEqual(self.processor.find_ftau(
            self.zcam_l, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(self.processor.find_ftau(
            self.zcam_r, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(self.processor.find_ftau(
            self.flf, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(self.processor.find_ftau(
            self.frf, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(self.processor.find_ftau(
            self.nlf, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(self.processor.find_ftau(
            self.nrf, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(
            self.processor.find_ftau(self.nlmv, 0.5), 0.42324086, 5)
        self.assertAlmostEqual(
            self.processor.find_ftau(self.nrmv, 0.5), 0.42324086, 5)

    # write a test suite for calculate_padding in ImProcessor
    def test_calculate_padding(self):
        self.assertEqual(self.processor.calculate_padding(self.zcam_l),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.zcam_r),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.flf),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.frf),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.nlf),
                         (1928, 1912, 1448, 1432))
        self.assertEqual(self.processor.calculate_padding(self.nrf),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.nlmv),
                         (0, 0, 0, 0))
        self.assertEqual(self.processor.calculate_padding(self.nrmv),
                         (0, 0, 0, 0))


if __name__ == '__main__':
    unittest.main()
