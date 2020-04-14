from uimfpy.UIMFReader import *
from uimfpy.utils import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def find_mz_pixel_idx(mz_arr, min_mz, max_mz, npixels):
    pixel_idx = npixels*(mz_arr-min_mz)/(max_mz-min_mz)
    if np.sum(pixel_idx < 0)>0: print("[ERR] mz_arr has some mass values less than min_mz")
    if np.sum(pixel_idx >= npixels)>0: pixel_idx[pixel_idx >= npixels] = npixels-1
    return pixel_idx.astype(np.int)

class UIMFViz(UIMFReader):
    """UIMFViz"""
    def __init__(self, uimf_file, TIC_threshold=None):
        super().__init__(uimf_file, TIC_threshold)

    def draw_2d_xic_in_target_mz(self, peak_region, target_mz, ppm_window=1000, isotopes = 4, npixels_scans=200,
        npixels_mz=100, padding_frames=10, padding_scans=10, verbose=False, fout=None):

        min_scan = peak_region.scan_start-padding_frames
        max_scan = peak_region.scan_end+padding_frames
        min_frame = peak_region.frame_start-padding_scans
        max_frame = peak_region.frame_end+padding_scans

        frame_nums=list(range(min_frame,max_frame))
        scan_nums=list(range(min_scan,max_scan))
        
        bins = self.get_mz_peaks(frame_nums=frame_nums,scan_nums=scan_nums)

        if isinstance(target_mz, list): __target_mz = target_mz
        else: __target_mz = [target_mz]
        
        for _mz in __target_mz:
            self.__plot_mz_by_scan(bins, _mz, ppm_window, max_scan, min_scan, npixels_scans, npixels_mz, verbose=verbose, fout=fout)
            self.__plot_2dxic_frame_by_scan(bins, _mz, 50, max_frame, min_frame, max_scan, min_scan, verbose=verbose, fout=fout)
        
        # heavy_mz = heavy_info['mono_mass']
        # light_mz = heavy_info['light_mono_mz']
        
        # plot_2d_map(bins, heavy_mz, 1000, max_scan, min_scan, num_pixels_scans, num_pixels_scans)
        # plot_2d_map(bins, light_mz, 1000, max_scan, min_scan, num_pixels_scans, num_pixels_scans)
        
        # heavy_pepseq = heavy_info['pepseq']
        # light_pepseq = heavy_info['light_pepseq']
        
        # mzbins_list = []
        # for iso in range(isotopes+1):
        #     mzbins_list.append(find_mzbin_infos(mzbins_by_mz, heavy_pepseq, heavy_mz, iso, heavy_info['z']))
        #     mzbins_list.append(find_mzbin_infos(mzbins_by_mz, light_pepseq, light_mz, iso, heavy_info['z']))
        
        # fig, ax = plt.subplots(2, 2, sharey=False)
        # cmap=plt.get_cmap("tab20")
        
        # fig_iso, axs_iso = plt.subplots(2, isotopes+1)
        # for i,info in enumerate(mzbins_list):
        #     if info not in mzbins_by_mz: continue
        #     print("info in mzbins_by_mz", info)
            
        #     _mzbins = mzbins_by_mz[info]
        #     if len(_mzbins) > 0:
        #         _xic2d = xic_matrix(xic_by_mzbin, _mzbins, reader.num_frames, reader.num_scans, normalize=False)
        #         _vec = _xic2d[min_frame:max_frame, min_scan:max_scan]
        #         axs_iso[i%2, info[3]].imshow(_vec)
        #         axs_iso[i%2, info[3]].set_axis_off()
        #         if info[0]==light_pepseq:
        #             color = 'b'
        #             ax[0, 1].plot(_vec.sum(axis=0), color='b')
        #             ax[1, 1].plot(_vec.sum(axis=1), color='b')
        #         else:
        #             ax[0, 0].plot(_vec.sum(axis=0), color='r')
        #             ax[1, 0].plot(_vec.sum(axis=1), color='r')
        # plt.show()
    
    def __plot_2dxic_frame_by_scan(self, bins, mz, ppm, max_frame, min_frame, max_scan, min_scan, verbose=False, fout=None):
        min_mz = mz*(1-1e-6*ppm)
        max_mz = mz*(1+1e-6*ppm)

        if verbose: print('[__plot_2dxic_frame_by_scan] min_mz:{}, max_mz:{}'.format(min_mz, max_mz))
        mat_2d = np.zeros((max_frame-min_frame, max_scan-min_scan))

        for f, s in bins:
            mz_arr = bins[f,s]['mz']
            _error = np.abs(ppm_error(mz, mz_arr))
            idx_good = (_error<= ppm)
            if np.sum(idx_good) > 0:
                int_sum = np.sum(bins[f,s]['int'][idx_good])
                mat_2d[max_frame-1-f, s-min_scan] += int_sum
        
        if verbose: print('[__plot_2dxic_frame_by_scan] mat_2d.min:{}, mat_2d.max:{}'.format(mat_2d.min(), mat_2d.max()))
        plt.imshow(mat_2d)
        plt.axis('off')
        im_ratio = (max_frame-min_frame)/(max_scan-min_scan)
        plt.colorbar(fraction=0.046*im_ratio, pad=0.04)
        if fout: plt.savefig('{}_{:.3f}_2dxic_frame_by_scan.pdf'.format(fout, mz))
        plt.show()

    def __plot_mz_by_scan(self, bins, mz, ppm, max_scan, min_scan, npixels_scans, npixels_mz, verbose=False, fout=None):
        min_mz = mz*(1-1e-6*ppm)
        max_mz = mz*(1+1e-6*ppm)

        if verbose: print('[__plot_mz_by_scan] min_mz:{}, max_mz:{}'.format(min_mz, max_mz))
        
        scan_width = max_scan-min_scan
        mat_2d = np.zeros((npixels_mz, npixels_scans))

        for f, s in bins:
            mz_arr = bins[f,s]['mz']
            _error = np.abs(ppm_error(mz, mz_arr))
            idx_good = (_error<= ppm)
            if np.sum(idx_good) > 0:
                f_mz = mz_arr[idx_good]
                f_int = bins[f,s]['int'][idx_good]

                mz_pixel_idx = find_mz_pixel_idx(f_mz, min_mz, max_mz, npixels_mz)
                
                if (npixels_scans > scan_width):
                    scan_idx1 = max(0, int(npixels_scans*(s-0.5-min_scan)/scan_width))
                    scan_idx2 = min(npixels_scans-1, int(npixels_scans*(s+0.5-min_scan)/scan_width))
                    for i in range(scan_idx1, scan_idx2):
                        mat_2d[npixels_mz-1-mz_pixel_idx, i] += f_int
                else:
                    scan_idx = int(npixels_scans*(s-min_scan)/scan_width)
                    mat_2d[npixels_mz-1-mz_pixel_idx, scan_idx] += f_int
        
        if verbose: print('[__plot_mz_by_scan] mat_2d.min:{}, mat_2d.max:{}'.format(mat_2d.min(), mat_2d.max()))
        plt.imshow(mat_2d, norm=colors.LogNorm(vmin=1, vmax=mat_2d.max()))
        plt.axis('off')
        im_ratio = npixels_mz/npixels_scans
        plt.colorbar(fraction=0.046*im_ratio, pad=0.04)
        if fout: plt.savefig('{}_{:.3f}_mz_by_scan.pdf'.format(fout, mz))
        plt.show()

        plt.plot(np.flip(mat_2d.sum(axis=1)))
        if fout: plt.savefig('{}_{:.3f}_mz_profile.pdf'.format(fout, mz))
        plt.show()
        
        # plt.figure()
        # ax = plt.gca()
        # im = ax.imshow(mat_2d)

        # # create an axes on the right side of ax. The width of cax will be 5%
        # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="5%", pad=0.05)

        # plt.colorbar(im, cax=cax)
