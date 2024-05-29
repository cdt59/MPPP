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
        self.processor.single_to_triple_channels(self.zcam_l)
        self.zcam_l.pad_im = True

        self.zcam_r = Image(
            ".\\data\\Sol_0709\\ZR0_0709_0729888971_069RAD_N0332864ZCAM07114_0340LMA01.IMG")
        self.processor.single_to_triple_channels(self.zcam_r)
        self.zcam_r.pad_im = True

        self.flf = Image(
            ".\\data\\Sol_0709\\FLF_0711_0730074611_304RAD_N0332864FHAZ00206_0A0295J01.IMG")
        self.processor.single_to_triple_channels(self.flf)
        self.flf.pad_im = True

        self.frf = Image(
            ".\\data\\Sol_0709\\FRF_0711_0730074611_304RAD_N0332864FHAZ00206_0A0295J01.IMG")
        self.processor.single_to_triple_channels(self.frf)
        self.frf.pad_im = True

        self.nlf = Image(
            ".\\data\\Sol_0709\\NLF_0709_0729883381_848RAD_N0332864SAPP00601_0A00LLJ01.IMG")
        self.processor.single_to_triple_channels(self.nlf)
        self.nlf.pad_im = True

        self.nrf = Image(
            ".\\data\\Sol_0709\\NRF_0709_0729889968_787RAD_N0332864NCAM03709_0A0195J01.IMG")
        self.processor.single_to_triple_channels(self.nrf)
        self.nrf.pad_im = True

        self.nlmv = Image(
            ".\\data\\Sol_0709\\NLMV0709_0729882743_168RAD_N0332841TRAV00262_0A02LLJ01.IMG")
        self.processor.single_to_triple_channels(self.nlmv)
        self.nlmv.pad_im = True

        self.nrmv = Image(
            ".\\data\\Sol_0709\\NRMV0709_0729882743_168RAD_N0332841TRAV00262_0A02LLJ01.IMG")
        self.processor.single_to_triple_channels(self.nrmv)
        self.nrmv.pad_im = True

    # def test_ftau(self):
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.zcam_l, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.zcam_r, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.flf, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.frf, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.nlf, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(self.processor.find_ftau(
    #         self.nrf, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(
    #         self.processor.find_ftau(self.nlmv, 0.5), 0.42324086, 5)
    #     self.assertAlmostEqual(
    #         self.processor.find_ftau(self.nrmv, 0.5), 0.42324086, 5)

    # def test_calculate_padding(self):
    #     self.assertEqual(self.processor.calculate_padding(self.zcam_l),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.zcam_r),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.flf),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.frf),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.nlf),
    #                      (1928, 1912, 1448, 1432))
    #     self.assertEqual(self.processor.calculate_padding(self.nrf),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.nlmv),
    #                      (0, 0, 0, 0))
    #     self.assertEqual(self.processor.calculate_padding(self.nrmv),
    #                      (0, 0, 0, 0))

    # def test_make_image_mask(self):
    #     zcam_l_padded_image = self.processor.pad_image(
    #         self.zcam_l.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.zcam_l, (0, 0, 0, 0)).shape,
    #                      zcam_l_padded_image.shape[:2])

    #     zcam_r_padded_image = self.processor.pad_image(
    #         self.zcam_r.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.zcam_r, (0, 0, 0, 0)).shape,
    #                      zcam_l_padded_image.shape[:2])

    #     flf_padded_image = self.processor.pad_image(
    #         self.flf.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.flf, (0, 0, 0, 0)).shape,
    #                      flf_padded_image.shape[:2])

    #     frf_padded_image = self.processor.pad_image(
    #         self.frf.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.frf, (0, 0, 0, 0)).shape,
    #                      frf_padded_image.shape[:2])

    #     nlf_padded_image = self.processor.pad_image(
    #         self.nlf.image, (1928, 1912, 1448, 1432))
    #     self.assertEqual(self.processor.make_image_mask(
    #         self.nlf, (1928, 1912, 1448, 1432)).shape, nlf_padded_image.shape[:2])

    #     nrf_padded_image = self.processor.pad_image(
    #         self.nrf.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.nrf, (0, 0, 0, 0)).shape,
    #                      nrf_padded_image.shape[:2])

    #     nlmv_padded_image = self.processor.pad_image(
    #         self.nlmv.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.nlmv, (0, 0, 0, 0)).shape,
    #                      nlmv_padded_image.shape[:2])

    #     nrmv_padded_image = self.processor.pad_image(
    #         self.nrmv.image, (0, 0, 0, 0))
    #     self.assertEqual(self.processor.make_image_mask(self.nrmv, (0, 0, 0, 0)).shape,
    #                      nrmv_padded_image.shape[:2])

    # tests for color_brightness_correction
    # def test_color_brightness_correction(self):
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.zcam_l, 0.5, 12, 1, 1).shape, self.zcam_l.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.zcam_r, 0.5, 12, 1, 1).shape, self.zcam_r.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.flf, 0.5, 12, 1, 1).shape, self.flf.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.frf, 0.5, 12, 1, 1).shape, self.frf.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.nlf, 0.5, 12, 1, 1).shape, self.nlf.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.nrf, 0.5, 12, 1, 1).shape, self.nrf.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.nlmv, 0.5, 12, 1, 1).shape, self.nlmv.image.shape)
    #     self.assertEqual(self.processor.color_brightness_correction(
    #         self.nrmv, 0.5, 12, 1, 1).shape, self.nrmv.image.shape)

    # test process_image()
    def test_process_image_pipeline(self):
        self.setUp()

        self.processor.process_image(self.zcam_l)
        self.processor.process_image(self.zcam_r)
        self.processor.process_image(self.flf)
        self.processor.process_image(self.frf)
        self.processor.process_image(self.nlf)
        self.processor.process_image(self.nrf)
        self.processor.process_image(self.nlmv)
        self.processor.process_image(self.nrmv)


if __name__ == '__main__':
    unittest.main()
