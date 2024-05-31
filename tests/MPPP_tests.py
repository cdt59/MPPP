import unittest
import os
from src2.ImProcessor import ImProcessor
from src2.Image import Image
from src.MPPP import image
import matplotlib.image
import numpy as np


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
    # def test_process_image_pipeline(self):
    #     self.setUp()

    #     self.processor.process_image(self.zcam_l)
    #     self.processor.process_image(self.zcam_r)
    #     self.processor.process_image(self.flf)
    #     self.processor.process_image(self.frf)
    #     self.processor.process_image(self.nlf)
    #     self.processor.process_image(self.nrf)
    #     self.processor.process_image(self.nlmv)
    #     self.processor.process_image(self.nrmv)

    def test_image_process(self):
        # for each image, compare the processing of ImProcessor and MPPP.Image
        self.setUp()

        scale = 12
        scale_red = 1.0
        scale_blue = 1.0
        clip_low = 0.01
        gamma = 2
        pad_im = True
        save_im = False
        save_mask = False
        find_offsets_mode = 0

        # ZCAM_LEFT
        zcam_l_old = image(os.path.join(
            "data/sol_0709", self.zcam_l.filename))
        zcam_l_old.scale = scale
        zcam_l_old.scale_red = scale_red
        zcam_l_old.scale_blue = scale_blue
        zcam_l_old.clip_low = clip_low
        zcam_l_old.gamma = gamma
        zcam_l_old.pad_im = pad_im
        zcam_l_old.save_im = save_im
        zcam_l_old.save_mask = save_mask
        zcam_l_old.find_offsets_mode = find_offsets_mode

        zcam_l_old.image_process()
        old_proc_zcam_l = zcam_l_old.im8

        self.processor.process_image(self.zcam_l)
        new_proc_zcam_l = self.zcam_l.im8

        self.assertEqual(old_proc_zcam_l.shape, new_proc_zcam_l.shape)
        self.assertTrue((old_proc_zcam_l == new_proc_zcam_l).all())

        old_mask_zcam_l = zcam_l_old.mask_im
        new_mask_zcam_l = self.zcam_l.mask_image

        # TODO: remove this block
        # # mask verbose information for dubugging
        # mask_diff = old_mask_zcam_l - new_mask_zcam_l
        # print(f"Mask total pixel difference: {mask_diff.sum().sum()}")
        # print(f"Mask max difference: {mask_diff.max().max()}")
        # print(f"Mask min difference: {mask_diff.min().min()}")
        # print(f"mask shape: {old_mask_zcam_l.shape}")

        # # save mask and differences
        # matplotlib.image.imsave("tests/old_mask_zcam_l.png", old_mask_zcam_l)
        # matplotlib.image.imsave("tests/new_mask_zcam_l.png", new_mask_zcam_l)
        # matplotlib.image.imsave("tests/mask_diff.png", -mask_diff)

        self.assertEqual(old_mask_zcam_l.shape, new_mask_zcam_l.shape)
        self.assertTrue((old_mask_zcam_l == new_mask_zcam_l).all())

        # ZCAM_right
        zcam_r_old = image(os.path.join(
            "data/sol_0709", self.zcam_r.filename))
        zcam_r_old.scale = scale
        zcam_r_old.scale_red = scale_red
        zcam_r_old.scale_blue = scale_blue
        zcam_r_old.clip_low = clip_low
        zcam_r_old.gamma = gamma
        zcam_r_old.pad_im = pad_im
        zcam_r_old.save_im = save_im
        zcam_r_old.save_mask = save_mask
        zcam_r_old.find_offsets_mode = find_offsets_mode

        zcam_r_old.image_process()
        old_proc_zcam_r = zcam_r_old.im8

        self.processor.process_image(self.zcam_r)
        new_proc_zcam_r = self.zcam_r.im8

        self.assertEqual(old_proc_zcam_r.shape, new_proc_zcam_r.shape)
        self.assertTrue((old_proc_zcam_r == new_proc_zcam_r).all())

        # FLF
        flf_old = image(os.path.join(
            "data/sol_0709", self.flf.filename))
        flf_old.scale = scale
        flf_old.scale_red = scale_red
        flf_old.scale_blue = scale_blue
        flf_old.clip_low = clip_low
        flf_old.gamma = gamma
        flf_old.pad_im = pad_im
        flf_old.save_im = save_im
        flf_old.save_mask = save_mask
        flf_old.find_offsets_mode = find_offsets_mode

        flf_old.image_process()
        old_proc_flf = flf_old.im8

        self.processor.process_image(self.flf)
        new_proc_flf = self.flf.im8

        self.assertEqual(old_proc_flf.shape, new_proc_flf.shape)
        self.assertTrue((old_proc_flf == new_proc_flf).all())

        old_mask_zcam_r = zcam_r_old.mask_im
        new_mask_zcam_r = self.zcam_r.mask_image
        self.assertEqual(old_mask_zcam_r.shape, new_mask_zcam_r.shape)
        self.assertTrue((old_mask_zcam_r == new_mask_zcam_r).all())

        old_mask_flf = flf_old.mask_im
        new_mask_flf = self.flf.mask_image
        self.assertEqual(old_mask_flf.shape, new_mask_flf.shape)
        self.assertTrue((old_mask_flf == new_mask_flf).all())

        # FRF
        frf_old = image(os.path.join(
            "data/sol_0709", self.frf.filename))
        frf_old.scale = scale
        frf_old.scale_red = scale_red
        frf_old.scale_blue = scale_blue
        frf_old.clip_low = clip_low
        frf_old.gamma = gamma
        frf_old.pad_im = pad_im
        frf_old.save_im = save_im
        frf_old.save_mask = save_mask
        frf_old.find_offsets_mode = find_offsets_mode

        frf_old.image_process()
        old_proc_frf = frf_old.im8

        self.processor.process_image(self.frf)
        new_proc_frf = self.frf.im8

        self.assertEqual(old_proc_frf.shape, new_proc_frf.shape)
        self.assertTrue((old_proc_frf == new_proc_frf).all())

        old_mask_frf = frf_old.mask_im
        new_mask_frf = self.frf.mask_image
        self.assertEqual(old_mask_frf.shape, new_mask_frf.shape)
        self.assertTrue((old_mask_frf == new_mask_frf).all())

        # NLF
        nlf_old = image(os.path.join(
            "data/sol_0709", self.nlf.filename))
        nlf_old.scale = scale
        nlf_old.scale_red = scale_red
        nlf_old.scale_blue = scale_blue
        nlf_old.clip_low = clip_low
        nlf_old.gamma = gamma
        nlf_old.pad_im = pad_im
        nlf_old.save_im = save_im
        nlf_old.save_mask = save_mask
        nlf_old.find_offsets_mode = find_offsets_mode

        nlf_old.image_process()
        old_proc_nlf = nlf_old.im8

        self.processor.process_image(self.nlf)
        new_proc_nlf = self.nlf.im8

        self.assertEqual(old_proc_nlf.shape, new_proc_nlf.shape)
        self.assertTrue((old_proc_nlf == new_proc_nlf).all())

        old_mask_nlf = nlf_old.mask_im
        new_mask_nlf = self.nlf.mask_image
        self.assertEqual(old_mask_nlf.shape, new_mask_nlf.shape)
        self.assertTrue((old_mask_nlf == new_mask_nlf).all())

        # NRF
        nrf_old = image(os.path.join(
            "data/sol_0709", self.nrf.filename))
        nrf_old.scale = scale
        nrf_old.scale_red = scale_red
        nrf_old.scale_blue = scale_blue
        nrf_old.clip_low = clip_low
        nrf_old.gamma = gamma
        nrf_old.pad_im = pad_im
        nrf_old.save_im = save_im
        nrf_old.save_mask = save_mask
        nrf_old.find_offsets_mode = find_offsets_mode

        nrf_old.image_process()
        old_proc_nrf = nrf_old.im8

        self.processor.process_image(self.nrf)
        new_proc_nrf = self.nrf.im8

        self.assertEqual(old_proc_nrf.shape, new_proc_nrf.shape)
        self.assertTrue((old_proc_nrf == new_proc_nrf).all())

        old_mask_nrf = nrf_old.mask_im
        new_mask_nrf = self.nrf.mask_image
        self.assertEqual(old_mask_nrf.shape, new_mask_nrf.shape)
        self.assertTrue((old_mask_nrf == new_mask_nrf).all())

        # NLMV
        nlmv_old = image(self.nlmv.IMG_path)
        nlmv_old.scale = scale
        nlmv_old.scale_red = scale_red
        nlmv_old.scale_blue = scale_blue
        nlmv_old.clip_low = clip_low
        nlmv_old.gamma = gamma
        nlmv_old.pad_im = pad_im
        nlmv_old.save_im = save_im
        nlmv_old.save_mask = save_mask
        nlmv_old.find_offsets_mode = find_offsets_mode

        nlmv_old.image_process()
        old_proc_nlmv = nlmv_old.im

        self.processor.process_image(self.nlmv)
        new_proc_nlmv = self.nlmv.proc_image

        self.assertEqual(old_proc_nlmv.shape, new_proc_nlmv.shape)
        self.assertTrue((old_proc_nlmv == new_proc_nlmv).all())

        old_mask_nlmv = nlmv_old.mask_im
        new_mask_nlmv = self.nlmv.mask_image
        self.assertEqual(old_mask_nlmv.shape, new_mask_nlmv.shape)
        self.assertTrue((old_mask_nlmv == new_mask_nlmv).all())

        # NRMV
        nrmv_old = image(self.nrmv.IMG_path)
        nrmv_old.scale = scale
        nrmv_old.scale_red = scale_red
        nrmv_old.scale_blue = scale_blue
        nrmv_old.clip_low = clip_low
        nrmv_old.gamma = gamma
        nrmv_old.pad_im = pad_im
        nrmv_old.save_im = save_im
        nrmv_old.save_mask = save_mask
        nrmv_old.find_offsets_mode = find_offsets_mode

        nrmv_old.image_process()
        old_proc_nrmv = nrmv_old.im

        self.processor.process_image(self.nrmv)
        new_proc_nrmv = self.nrmv.proc_image

        self.assertEqual(old_proc_nrmv.shape, new_proc_nrmv.shape)
        self.assertTrue((old_proc_nrmv == new_proc_nrmv).all())

        old_mask_nrmv = nrmv_old.mask_im
        new_mask_nrmv = self.nrmv.mask_image
        self.assertEqual(old_mask_nrmv.shape, new_mask_nrmv.shape)
        self.assertTrue((old_mask_nrmv == new_mask_nrmv).all())


if __name__ == '__main__':
    unittest.main()
