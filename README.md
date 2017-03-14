# Catalog Match

Given input `ra, dec` coordinates, uses the [astroquery][1] package to download
the same frame from a selected catalog.

Both catalogs are then matched using astropy's [match_to_catalog_sky][2].

Available catalogs in `astroquery`:

```
allwise_p3as_psd                AllWISE Source Catalog
allwise_p3as_mep                AllWISE Multiepoch Photometry Table
allwise_p3as_psr                AllWISE Reject Table
allwise_p3as_cdd                AllWISE Atlas Metadata Table
allwise_p3am_xrf                AllWISE Frame Cross-Reference Table
allwise_p3al_lod                AllWISE Atlas Inventory Table
allwise_p3am_cdd                AllWISE Atlas Image Inventory Table
neowiser_p1bs_psd               NEOWISE-R Single Exposure (L1b) Source Table
neowiser_p1ba_mch               NEOWISE-R Known Solar System Object Possible Association List ( Caution )
neowiser_p1bs_frm               NEOWISE-R Single Exposure (L1b) Frame Metadata Table
neowiser_p1bl_lod               NEOWISE-R Single Exposure (L1b) Scan Inventory Table
neowiser_p1bm_frm               NEOWISE-R Single Exposure (L1b) Image Inventory Table
allsky_4band_p3as_psd           WISE All-Sky Source Catalog
allsky_4band_p1bs_psd           WISE All-Sky Single Exposure (L1b) Source Table
allsky_4band_p1ba_mch           WISE All-Sky Known Solar System Object Possible Association List ( Caution )
allsky_4band_p3as_psr           WISE All-Sky Reject Table
allsky_4band_p3as_cdd           WISE All-Sky Atlas Metadata Table
allsky_4band_p3am_xrf           WISE All-Sky Frame Cross-Reference Table
allsky_4band_p1bs_frm           WISE All-Sky Single Exposure (L1b) Frame Metadata Table
allsky_4band_p1bl_lod           WISE All-Sky Single Exposure (L1b) Scan Inventory Table
allsky_4band_p3al_lod           WISE All-Sky Atlas Inventory Table
allsky_4band_p1bm_frm           WISE All-Sky Single Exposure (L1b) Image Inventory Table
allsky_4band_p3am_cdd           WISE All-Sky Atlas Image Inventory Table
allsky_3band_p3as_psd           WISE 3-Band Cryo Source Working Database ( Readme)
allsky_3band_p1bs_psd           WISE 3-Band Cryo Single Exposure (L1b) Source Table
allsky_3band_p1ba_mch           WISE 3-Band Cryo Known Solar System Object Possible Association List ( Caution )
allsky_3band_p3as_cdd           WISE 3-Band Cryo Atlas Metadata Table
allsky_3band_p3am_xrf           WISE 3-Band Cryo Frame Cross-Reference Table
allsky_3band_p1bs_frm           WISE 3-Band Cryo Single Exposure (L1b) Frame Metadata Table
allsky_3band_p1bl_lod           WISE 3-Band Cryo Single Exposure (L1b) Scan Inventory Table
allsky_3band_p3al_lod           WISE 3-Band Cryo Atlas Inventory Table
allsky_3band_p1bm_frm           WISE 3-Band Cryo Single Exposure (L1b) Image Inventory Table
allsky_3band_p3am_cdd           WISE 3-Band Cryo Atlas Image Inventory Table
allsky_2band_p1bs_psd           WISE Post-Cryo Single Exposure (L1b) Source Table
allsky_2band_p1bs_frm           WISE Post-Cryo Single Exposure (L1b) Frame Metadata Table
allsky_2band_p1ba_mch           WISE Post-Cryo Single Exposure (L1b) Known SSO Possible Association List ( Caution )
allsky_2band_p1bl_lod           WISE Post-Cryo Single Exposure (L1b) Scan Inventory Table
allsky_2band_p1bm_frm           WISE Post-Cryo Single Exposure (L1b) Image Inventory Table
prelim_2band_p1bs_psd           WISE Preliminary Post-Cryo Single Exposure (L1b)  Source Table (Superseded)
prelim_2band_p1ba_mch           WISE Preliminary Post-Cryo Solar System Object Possible Association List ( Caution , Superseded)
prelim_2band_p1bs_frm           WISE Preliminary Post-Cryo Single Exposure (L1b) Frame Metadata Table (Superseded)
prelim_2band_p1bl_lod           WISE Preliminary Post-Cryo Single Exposure (L1b) Scan Inventory Table (Superseded)
prelim_2band_p1bm_frm           WISE Preliminary Post-Cryo Single Exposure (L1b) Image Inventory Table (Superseded)
prelim_p3as_psd                 WISE Preliminary Release Source Catalog (Superseded)
prelim_p1bs_psd                 WISE Preliminary Release Single Exposure (L1b) Source Table (Superseded)
prelim_p1ba_mch                 WISE Preliminary Release Known Solar System Object Possible Association List ( Caution , Superseded)
prelim_p3as_cdd                 WISE Preliminary Release Atlas Metadata Table (Superseded)
prelim_p3am_xrf                 WISE Preliminary Release Frame Cross-Reference Table (Superseded)
prelim_p1bs_frm                 WISE Preliminary Release Single Exposure (L1b) Frame Metadata Table (Superseded)
prelim_p1bl_lod                 WISE Preliminary Release Single Exposure (L1b) Scan Inventory Table (Superseded)
prelim_p3al_lod                 WISE Preliminary Release Atlas Inventory Table (Superseded)
prelim_p1bm_frm                 WISE Preliminary Release Single Exposure (L1b) Image Inventory Table (Superseded)
prelim_p3am_cdd                 WISE Preliminary Release Atlas Image Inventory Table (Superseded)
fp_psc                          2MASS All-Sky Point Source Catalog (PSC)
fp_xsc                          2MASS All-Sky Extended Source Catalog (XSC)
lga_v2                          The 2MASS Large Galaxy Atlas
fp_scan_dat                     2MASS All-Sky Survey Scan Info (Read Me!)
fp_coadd_dat                    2MASS All-Sky Survey Atlas Image Info
pt_src_rej                      2MASS Survey Point Source Reject Table
wsdb_info                       2MASS Survey Merged Point Source Information Table
wsdb_link                       2MASS Survey Merged Point Source Link Table
ext_src_rej                     2MASS Survey Extended Source Reject Table
ewsdbf_info                     2MASS Survey Merged Extended Source Information Table
ewsdbf_link                     2MASS Survey Merged Extended Source Link Table
pscan_dat                       2MASS Survey Scan Info
coadd_dat                       2MASS Survey Atlas Image Info
pt_src_6x2                      2MASS 6X w/LMC/SMC Point Source Working Database / Catalog ( Read Me! )
sixxf_info                      2MASS 6X w/LMC/SMC Merged Point Source Information Table
sixxf_link                      2MASS 6X w/LMC/SMC Merged Point Source Link Table
ext_src_6x2                     2MASS 6X w/LMC/SMC Extended Source Working Database / Catalog ( Read Me! )
esixxf_info                     2MASS 6X w/LMC/SMC Merged Extended Source Information Table
esixxf_link                     2MASS 6X w/LMC/SMC Merged Extended Source Link Table
pscan_dat_6x2                   2MASS 6X w/LMC/SMC Scan Info
coadd_dat_6x2                   2MASS 6X w/LMC/SMC Atlas Image Info
pt_src_c                        2MASS Calibration Point Source Working Database
cf_info                         2MASS Calibration Merged Point Source Information Table
cf_link                         2MASS Calibration Merged Point Source Link Table
ext_src_c                       2MASS Calibration Extended Source Working Database
ecf_info                        2MASS Calibration Merged Extended Source Information Table
ecf_link                        2MASS Calibration Merged Extended Source Link Table
pscan_dat_c                     2MASS Calibration Scan Info
coadd_dat_c                     2MASS Calibration Atlas Image Info
deepcal_src                     2MASS Combined Calibration Field Source Table
pt_src_sc                       2MASS LMC/SMC Calibration Point Source Working Database
scf_info                        2MASS LMC/SMC Calibration Merged Point Source Information Table
scf_link                        2MASS LMC/SMC Calibration Merged Point Source Link Table
ext_src_sc                      2MASS LMC/SMC Calibration Extended Source Working Database
escf_info                       2MASS LMC/SMC Calibration Merged Extended Source Information Table
escf_link                       2MASS LMC/SMC Calibration Merged Extended Source Link Table
pscan_dat_sc                    2MASS LMC/SMC Calibration Scan Info
coadd_dat_sc                    2MASS LMC/SMC Calibration Atlas Image Info
pt_src_cat                      2MASS Second Incremental Release Point Source Catalog (PSC)
ext_src_cat                     2MASS Second Incremental Release Extended Source Catalog (XSC)
scan_dat_2                      2MASS Second Incremental Release Survey Scan Info
pt_src_cat1                     2MASS First Incremental Release Point Source Catalog (PSC)
ext_src_cat1                    2MASS First Incremental Release Extended Source Catalog (XSC)
scan_dat                        2MASS First Incremental Release Survey Scan Info
pts_samp_cat                    2MASS Sampler Point Source Catalog (PSC)
exts_samp_cat                   2MASS Sampler Extended Source Catalog (XSC)
irasfsc                         IRAS Faint Source Catalog v2.0 (FSC)
iraspsc                         IRAS Point Source Catalog v2.1 (PSC)
irasgal                         IRAS Cataloged Galaxies and Quasars
irasssc                         IRAS Serendipitous Survey Catalog
irassss                         IRAS Small Scale Structure Catalog
iras_ao                         IRAS Additional Observations (AO) Catalog
irasfscr                        IRAS Faint Source Catalog Rejects
iraspscr                        IRAS Point Source Catalog Rejects
iraspscw                        IRAS PSC joined with WSDB
iraspsch                        IRAS PSC joined with HCON and WSDB
irascatalog                     IRAS 1.2-Jy Redshift Survey
comsight                        IRAS Asteroid and Comet Survey
astsight                        IRAS Minor Planet Survey
summary                         IRAS Large Galaxies Catalog
ercscbm                         Planck ERCSC Bandmerged Catalog
planckphz                       Planck List of High Redshift Source Candidates
com_pccs2_030                   Planck PCCS2 30GHz Catalog
com_pccs2_044                   Planck PCCS2 44GHz Catalog
com_pccs2_070                   Planck PCCS2 70GHz Catalog
com_pccs2_100                   Planck PCCS2 100GHz Catalog
com_pccs2_143                   Planck PCCS2 143GHz Catalog
com_pccs2_217                   Planck PCCS2 217GHz Catalog
com_pccs2_353                   Planck PCCS2 353GHz Catalog
com_pccs2_545                   Planck PCCS2 545GHz Catalog
com_pccs2_857                   Planck PCCS2 857GHz Catalog
com_pccs2e_100                  Planck PCCS2E 100GHz Catalog (lower reliability)
com_pccs2e_143                  Planck PCCS2E 143GHz Catalog (lower reliability)
com_pccs2e_217                  Planck PCCS2E 217GHz Catalog (lower reliability)
com_pccs2e_353                  Planck PCCS2E 353GHz Catalog (lower reliability)
com_pccs2e_545                  Planck PCCS2E 545GHz Catalog (lower reliability)
com_pccs2e_857                  Planck PCCS2E 857GHz Catalog (lower reliability)
com_pccs2_sz_mmf1               Planck PR2 Sunyaev-Zeldovich Cluster MMF1 List
com_pccs2_sz_mmf3               Planck PR2 Sunyaev-Zeldovich Cluster MMF3 List
com_pccs2_sz_pws                Planck PR2 Sunyaev-Zeldovich Cluster PwS List
com_pccs2_sz_union              Planck PR2 Sunyaev-Zeldovich Cluster UNION List
com_pccs2_gcc                   Planck Catalog of Galactic Cold Clumps
com_pccs1_030                   Planck PCCS 30GHz Catalog
com_pccs1_044                   Planck PCCS 44GHz Catalog
com_pccs1_070                   Planck PCCS 70GHz Catalog
com_pccs1_100                   Planck PCCS 100GHz Catalog
com_pccs1_143                   Planck PCCS 143GHz Catalog
com_pccs1_217                   Planck PCCS 217GHz Catalog
com_pccs1_353                   Planck PCCS 353GHz Catalog
com_pccs1_545                   Planck PCCS 545GHz Catalog
com_pccs1_857                   Planck PCCS 857GHz Catalog
com_pccs1_sz_mmf1               Planck Sunyaev-Zeldovich Cluster MMF1 List
com_pccs1_sz_mmf3               Planck Sunyaev-Zeldovich Cluster MMF3 List
com_pccs1_sz_pws                Planck Sunyaev-Zeldovich Cluster PwS List
com_pccs1_sz_union2             Planck Sunyaev-Zeldovich Cluster UNION List v2.1
ercsc_f030_e                    Planck ERCSC 30GHz Catalog
ercsc_f044_e                    Planck ERCSC 44GHz Catalog
ercsc_f070_e                    Planck ERCSC 70GHz Catalog
ercsc_f100_e                    Planck ERCSC 100GHz Catalog
ercsc_f143_e                    Planck ERCSC 143GHz Catalog
ercsc_f217_e                    Planck ERCSC 217GHz Catalog
ercsc_f353_e                    Planck ERCSC 353GHz Catalog
ercsc_f545_e                    Planck ERCSC 545GHz Catalog
ercsc_f857_e                    Planck ERCSC 857GHz Catalog
ecc                             Planck Early Cold Core Source List (ECC)
esz                             Planck Early Sunyaev-Zeldovich Cluster List (ESZ)
slphotdr4                       SEIP Source List
slicovv2                        SEIP IRAC Coverage Table
slmcovv2                        SEIP MIPS Coverage Table
sltracev2                       SEIP Traceback Table
irs_enhv211                     IRS Enhanced Products
a1763t2                         Abell 1763 Source Catalog
a1763t3                         Abell 1763 MIPS 70 micron Catalog
dr4_clouds_full                 C2D Fall '07 Full CLOUDS Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_clouds_hrel                 C2D Fall '07 High Reliability (HREL) CLOUDS Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_clouds_ysoc                 C2D Fall '07 candidate Young Stellar Objects (YSO) CLOUDS Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_off_cloud_full              C2D Fall '07 Full OFF-CLOUD Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_off_cloud_hrel              C2D Fall '07 High Reliability (HREL) OFF-CLOUD Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_off_cloud_ysoc              C2D Fall '07 candidate Young Stellar Objects (YSO) OFF-CLOUD Catalog (CHA_II, LUP, OPH, PER, SER)
dr4_cores_full                  C2D Fall '07 Full CORES Catalog
dr4_cores_hrel                  C2D Fall '07 High Reliability (HREL) CORES Catalog
dr4_cores_ysoc                  C2D Fall '07 candidate Young Stellar Objects (YSO) CORES Catalog
dr4_stars_full                  C2D Fall '07 Full STARS Catalog
dr4_stars_hrel                  C2D Fall '07 High Reliability (HREL) STARS Catalog
dr4_stars_ysoc                  C2D Fall '07 candidate Young Stellar Objects (YSO) STARS Catalog
dr4_MM                          C2D Fall '07 Millimeter (MM) Sources Catalog (OPH, PER, SER Clouds)
dr4_trans1_full                 C2D Fall '07 Perseus Epoch 1 Transient Sources FULL Catalog
dr4_trans2_full                 C2D Fall '07 Perseus Epoch 2 Transient Sources FULL Catalog
clash36_v2                      CLASH 3.6 micron Catalog
clash45_v2                      CLASH 4.5 micron Catalog
clash58_v2                      CLASH 5.8 micron Catalog
clash80_v2                      CLASH 8.0 micron Catalog
csi2264t1                       CSI 2264 Object Table
csi2264t2                       CSI 2264 CoRoT Light Curves
csi2264t3                       CSI 2264 Spitzer Light Curves
cygx_arch                       Cygnus-X Archive
cygx_cat                        Cygnus-X Catalog
dustingsfull                    DUSTiNGS Full Catalog
dustingsgsc                     DUSTiNGS Good Source Catalog
feps_phot                       FEPS Photometry Catalog
xfls_i1m                        Extragalactic FLS IRAC Channel 1 Main Field Catalog
xfls_i2m                        Extragalactic FLS IRAC Channel 2 Main Field Catalog
xfls_i3m                        Extragalactic FLS IRAC Channel 3 Main Field Catalog
xfls_i4m                        Extragalactic FLS IRAC Channel 4 Main Field Catalog
xfls_i1v                        Extragalactic FLS IRAC Channel 1 Verification Field Catalog
xfls_i2v                        Extragalactic FLS IRAC Channel 2 Verification Field Catalog
xfls_i3v                        Extragalactic FLS IRAC Channel 3 Verification Field Catalog
xfls_i4v                        Extragalactic FLS IRAC Channel 4 Verification Field Catalog
xfls_iallm                      Extragalactic FLS IRAC Bandmerged Main Field Catalog
xfls_iallv                      Extragalactic FLS IRAC Bandmerged Verification Field Catalog
xfls_m1t2                       Extragalactic FLS MIPS 24 micron Extended Source Catalog
xfls_m1t4                       Extragalactic FLS MIPS 24 micron Calibration Star Catalog
xfls_m1t5                       Extragalactic FLS MIPS 24 micron Point Source Catalog
xfls_m2                         Extragalactic FLS MIPS 70 micron Catalog
xfls_m3                         Extragalactic FLS MIPS 160 micron Catalog
xfls_kpno                       Extragalactic FLS KPNO R-band Source List
fls_release_v2_mmt_spectra      FLS MMT/Hectospec Spectroscopic Catalog (V2)
fls_release_v2_photom           FLS SDSS and MIPS Astrometric and Photometric Catalog (V2)
fls_release_v2_sdss_spectra     FLS SDSS Spectroscopic Catalog (V2)
xfls_w2                         Extragalactic FLS WIYN/Hydra Spectroscopic Catalog
xfls_w3                         Extragalactic FLS WIYN/Hydra Line Strength and Equivalent Width Catalog
xfls_w4                         Extragalactic FLS WIYN/Hydra Line Ratios and Extinction Catalog
galcen_psc                      Point Source in a Spitzer/IRAC Survey of the Galactic Center (Ramirez et al. 2008)
glimpse_s07_ar                  GLIMPSE I Spring '07 Archive (more complete, less reliable)
glimpse_s07                     GLIMPSE I Spring '07 Catalog (highly reliable)
glimpse2_v2arc                  GLIMPSE II Spring '08 Archive (more complete, less reliable)
glimpse2_v2cat                  GLIMPSE II Spring '08 Catalog (highly reliable)
glimpse2ep1a08                  GLIMPSE II Epoch 1 December '08 Archive (more complete, less reliable)
glimpse2ep1c08                  GLIMPSE II Epoch 1 December '08 Catalog (highly reliable)
glimpse2ep2a09                  GLIMPSE II Epoch 2 November '09 Archive (more complete, less reliable)
glimpse2sub                     GLIMPSEII Subarray Source List
glimpse2ep2mra09                GLIMPSE II Epoch 2 November '09 More Reliable Archive (more reliable)
glimpse3d_v2arc                 GLIMPSE 3D, 2007-2009 Archive (more complete, less reliable),(Erratum)
glimpse3d_v1cat_tbl             GLIMPSE 3D, 2007-2009 Catalog (highly reliable)
glimpse3dep1a                   GLIMPSE 3D Epoch 1 Archive (more complete, less reliable)
glimpse3dep1c                   GLIMPSE 3D Epoch 1 Catalog (highly reliable)
glimpse3dep2a                   GLIMPSE 3D Epoch 2 Archive (more complete, less reliable)
glimpse3dep2mra                 GLIMPSE 3D Epoch 2 More Reliable Archive (more complete, less reliable)
glimpse360a                     GLIMPSE360 Archive (more complete, less reliable)
glimpse360c                     GLIMPSE360 Catalog (highly reliable)
velcara                         Vela-Carina Archive (more complete, less reliable)
velcarc                         Vela-Carina Catalog (highly reliable)
deepglimpsea                    Deep GLIMPSE Archive (more complete, less reliable)
deepglimpsec                    Deep GLIMPSE Catalog (highly reliable)
glimpsesmoga                    SMOG Archive (more complete, less reliable)
glimpsesmogc                    SMOG Catalog (highly reliable)
glimpsecygxa                    GLIMPSE Cygnus-X Archive (more complete, less reliable)
glimpsecygxc                    GLIMPSE Cygnus-X Catalog (highly reliable)
goodsn_mips24                   GOODS-N MIPS 24 micron Photometry Catalog
goods_mips24                    GOODS-S MIPS 24 micron Photometry Catalog
goodsnirs16                     GOODS-N IRS 16 micron Photometry Catalog
goodssirs16                     GOODS-S IRS 16 micron Photometry Catalog
m31irac                         M31 IRAC Catalog
mipslg                          MIPS Local Galaxies Catalog
mipsgala                        MIPSGAL Archive
mipsgalc                        MIPSGAL Catalog
s4gcat                          Spitzer Survey of Stellar Structure in Galaxies (S4G)
safires70                       Spitzer Archival Far-Infrared Extragalactic Survey (SAFIRES) MIPS 70 micron Catalog
safires160                      Spitzer Archival Far-Infrared Extragalactic Survey (SAFIRES) MIPS 160 micron Catalog
sagearciracv2                   SAGE Winter '08 IRAC Epoch 1 and Epoch 2 Archive (more complete, less reliable)
sagecatiracv2                   SAGE Winter '08 IRAC Epoch 1 and Epoch 2 Catalog (more reliable)
sage_ar_irac                    SAGE IRAC Single Frame + Mosaic Photometry Archive (more complete, less reliable)
sage_cat_irac                   SAGE IRAC Single Frame + Mosaic Photometry Catalog (more reliable)
sage_ar_irac_e1e2               SAGE IRAC Epoch 1 and Epoch 2 Archive (more complete, less reliable)
sage_cat_irac_e1e2              SAGE IRAC Epoch 1 and Epoch 2 Catalog (more reliable)
sage_ar_irac_match              SAGE IRAC Matched Epoch Catalog (more complete, less reliable)
sage_cat_irac_match             SAGE IRAC Matched Epoch Archive (more reliable)
sage_ar_irac_off                SAGE IRAC Offset Position Epoch 1 and Epoch 2 Archive (more complete, less reliable)
sage_cat_irac_off               SAGE IRAC Offset Position Epoch 1 and Epoch 2 Catalog (more reliable)
sagecatmips24v2                 SAGE Winter '08 MIPS 24 um Epoch 1 and Epoch 2 Catalog (more reliable)
sage_cat_m24                    SAGE MIPS 24 um Epoch 1 and Epoch 2 Catalog (more reliable)
sage_full_m24                   SAGE MIPS 24 um Epoch 1 and Epoch 2 Full List (more complete, less reliable)
sage_cat_m24_match              SAGE MIPS 24 um Matched Epoch Catalog (more reliable)
sage_full_m24_match             SAGE MIPS 24 um Matched Epoch Full List (more complete, less reliable)
sage_cat_m70                    SAGE MIPS 70 um Combined Epoch Catalog (more reliable)
sage_full_m70                   SAGE MIPS 70 um Combined Epoch Full List (more complete, less reliable)
sage_cat_m160                   SAGE MIPS 160 um Combined Epoch Catalog (more reliable)
sage_full_m160                  SAGE MIPS 160 um Combined Epoch Catalog (more complete, less reliable)
sagefull                        SAGE-Var LMC Full Catalog
sagevar                         SAGE-Var LMC Variable Catalog
sagesmc_iracadr3                SAGE-SMC IRAC Single Frame + Mosaic Photometry Archive v1.5
sagesmc_iraccdr3                SAGE-SMC IRAC Single Frame + Mosaic Photometry Catalog v1.5
sagesmc_iracep1a                SAGE-SMC IRAC Epoch 1 Archive (less reliable)
sagesmc_iracep1c                SAGE-SMC IRAC Epoch 1 Catalog (more reliable)
sagesmc_iraca                   SAGE-SMC IRAC Epoch 0, Epoch 1, and Epoch 2 Archive (less reliable)
sagesmc_iracc                   SAGE-SMC IRAC Epoch 0, Epoch 1, and Epoch 2 Catalog (more reliable)
sagesmc_mips24ep1c              SAGE-SMC MIPS 24um Epoch 1 Catalog (more reliable)
sagesmc_mips24ep1f              SAGE-SMC MIPS 24um Epoch 1 Full List (less reliable)
sagesmc_mips24f                 SAGE-SMC MIPS 24 um Epoch 0, Epoch 1, and Epoch 2 Full List (more complete, less reliable)
sagesmc_mips24c                 SAGE-SMC MIPS 24 um Epoch 0, Epoch 1, and Epoch 2 Catalog (more reliable)
sagesmc_mips70f                 SAGE-SMC MIPS 70um Combined Epoch Full List (more complete, less reliable)
sagesmc_mips70c                 SAGE-SMC MIPS 70um Combined Epoch Catalog (more reliable)
sagesmc_mips160f                SAGE-SMC MIPS 160um Combined Epoch Full List (more complete, less reliable)
sagesmc_mips160c                SAGE-SMC MIPS 160um Combined Epoch Catalog (more reliable)
sagesmcfull                     SAGE-Var SMC Full Catalog
sagesmcvar                      SAGE-Var SMC Variable Catalog
ssid2                           SAGE-Spec ID Search
sass_v3                         SASS October 2011 Catalog
scosmos_irac_0407               S-COSMOS IRAC 4-channel Photometry Catalog June 2007 (README)
scosmos_mips_24_go2             S-COSMOS MIPS 24um MAIN Photometry Catalog June 2007 ((Aug 2008: Important Flux-correction Note))
scosmos_mips_24_go2_deep        S-COSMOS MIPS 24um DEEP Photometry Catalog June 2007 ((Aug 2008: Important Flux-correction Note))
scosmos_mips_70_v3              S-COSMOS MIPS 70um Photometry Catalog v3 Jan 2009
scosmos_mips_160_v3             S-COSMOS MIPS 160um Photometry Catalog v3 Jan 2009
sdwfs_ch1_stack                 SDWFS Aug'09 DR1.1 IRAC 3.6um-Selected Total Coadd Stack
sdwfs_ch2_stack                 SDWFS Aug '09 DR1.1 IRAC 4.5um-Selected Total Coadd Stack
sdwfs_ch3_stack                 SDWFS Aug '09 DR1.1 IRAC 5.8um-Selected Total Coadd Stack
sdwfs_ch4_stack                 SDWFS Aug '09 DR1.1 IRAC 8.0um-Selected Total Coadd Stack
sdwfs_ch1_epoch1                SDWFS Aug '09 DR1.1 IRAC 3.6um-Selected 3x30sec Coadd, epoch 1 (Jan '04)
sdwfs_ch2_epoch1                SDWFS Aug '09 DR1.1 IRAC 4.5um-Selected 3x30sec Coadd, epoch 1 (Jan '04)
sdwfs_ch3_epoch1                SDWFS Aug '09 DR1.1 IRAC 5.8um-Selected 3x30sec Coadd, epoch 1 (Jan '04)
sdwfs_ch4_epoch1                SDWFS Aug '09 DR1.1 IRAC 8.0um-Selected 3x30sec Coadd, epoch 1 (Jan '04)
sdwfs_ch1_epoch2                SDWFS Aug '09 DR1.1 IRAC 3.6um-Selected 3x30sec Coadd, epoch 2 (Aug '07)
sdwfs_ch2_epoch2                SDWFS Aug '09 DR1.1 IRAC 4.5um-Selected 3x30sec Coadd, epoch 2 (Aug '07)
sdwfs_ch3_epoch2                SDWFS Aug '09 DR1.1 IRAC 5.8um-Selected 3x30sec Coadd, epoch 2 (Aug '07)
sdwfs_ch4_epoch2                SDWFS Aug '09 DR1.1 IRAC 8.0um-Selected 3x30sec Coadd, epoch 2 (Aug '07)
sdwfs_ch1_epoch3                SDWFS Aug '09 DR1.1 IRAC 3.6um-Selected 3x30sec Coadd, epoch 3 (Feb '08)
sdwfs_ch2_epoch3                SDWFS Aug '09 DR1.1 IRAC 4.5um-Selected 3x30sec Coadd, epoch 3 (Feb '08)
sdwfs_ch3_epoch3                SDWFS Aug '09 DR1.1 IRAC 5.8um-Selected 3x30sec Coadd, epoch 3 (Feb '08)
sdwfs_ch4_epoch3                SDWFS Aug '09 DR1.1 IRAC 8.0um-Selected 3x30sec Coadd, epoch 3 (Feb '08)
sdwfs_ch1_epoch4                SDWFS Aug '09 DR1.1 IRAC 3.6um-Selected 3x30sec Coadd, epoch 4 (Mar '08)
sdwfs_ch2_epoch4                SDWFS Aug '09 DR1.1 IRAC 4.5um-Selected 3x30sec Coadd, epoch 4 (Mar '08)
sdwfs_ch3_epoch4                SDWFS Aug '09 DR1.1 IRAC 5.8um-Selected 3x30sec Coadd, epoch 4 (Mar '08)
sdwfs_ch4_epoch4                SDWFS Aug '09 DR1.1 IRAC 8.0um-Selected 3x30sec Coadd, epoch 4 (Mar '08)
sdwfs_lcurve                    SDWFS Light Curve Catalog
sdwfs_var                       SDWFS Variability Catalog
sepm24                          SEP MIPS 24 micron Point Source Catalog
sepm70                          SEP MIPS 70 micron Point Source Catalog
sepmext                         SEP MIPS Extended Source Catalog
sepirac                         SEP IRAC-based Multiwavelength Photometric Catalog
servscdfsi1                     SERVS CDFS 3.6 micron Catalog
servscdfsi2                     SERVS CDFS 4.5 micron Catalog
servscdfsi12                    SERVS CDFS 2-band Catalog (highly reliable)
servseni1                       SERVS ELAIS N1 3.6 micron Catalog
servseni2                       SERVS ELAIS N1 4.5 micron Catalog
servseni12                      SERVS ELAIS N1 2-band Catalog (highly reliable)
servsesi1                       SERVS ELAIS S1 3.6 micron Catalog
servsesi2                       SERVS ELAIS S1 4.5 micron Catalog
servsesi12                      SERVS ELAIS S1 2-band Catalog (highly reliable)
servslhi1                       SERVS Lockman Hole 3.6 micron Catalog
servslhi2                       SERVS Lockman Hole 4.5 micron Catalog
servslhi12                      SERVS Lockman Hole 2-band Catalog (highly reliable)
servsxmmi1                      SERVS XMM-LSS 3.6 micron Catalog
servsxmmi2                      SERVS XMM-LSS 4.5 micron Catalog
servsxmmi12                     SERVS XMM-LSS 2-band Catalog (highly reliable)
shelacomb                       SHELA Combined Epoch IRAC Catalog
shelaep1                        SHELA Epoch 1 IRAC Catalog
shelaep2                        SHELA Epoch 2 IRAC Catalog
shelaep3                        SHELA Epoch 3 IRAC Catalog
shelasdss                       SHELA-SDSS Stripe 82 Catalog
simple                          SIMPLE Photometry Catalog
spiesch1                        SpIES 3.6 micron-only Catalog
spiesch2                        SpIES 4.5 micron-only Catalog
spiesch12                       SpIES Dual-band Catalog
spuds_irac                      SpUDS IRAC Catalog
spuds_mips                      SpUDS MIPS Catalog
ssdf1                           SSDF IRAC Ch1 Catalog
ssdf2                           SSDF IRAC Ch2 Catalog
chandra_cat_f05                 SWIRE CDFS Region Fall '05  Spitzer Catalog
chandra_24_cat_f05              SWIRE CDFS Region 24um Fall '05 Spitzer Catalog
chandra_70_cat_f05              SWIRE CDFS Region 70um Fall '05 Spitzer Catalog
chandra_160_cat_f05             SWIRE CDFS Region 160um Fall '05 SWIRE Spitzer Catalog
elaisn1_cat_s05                 SWIRE ELAIS N1 Region Spring '05 Spitzer Catalog
elaisn1_24_cat_s05              SWIRE ELAIS N1 Region 24um Spring '05 Spitzer Catalog
elaisn1_70_cat_s05              SWIRE ELAIS N1 Region 70um Spring '05 Spitzer Catalog
elaisn1_160_cat_s05             SWIRE ELAIS N1 Region 160um Spring '05 Spitzer Catalog
elaisn2_cat_s05                 SWIRE ELAIS N2 Region Spring '05 Spitzer Catalog
elaisn2_24_cat_s05              SWIRE ELAIS N2 Region 24um Spring '05 Spitzer Catalog
elaisn2_70_cat_s05              SWIRE ELAIS N2 Region 70um Spring '05 Spitzer Catalog
elaisn2_160_cat_s05             SWIRE ELAIS N2 Region 160um Spring '05 Spitzer Catalog
elaiss1_cat_f05                 SWIRE ELAIS S1 Region Fall '05 SWIRE Spitzer Catalog
elaiss1_24_cat_f05              SWIRE ELAIS S1 Region 24um Fall '05 Spitzer Catalog
elaiss1_70_cat_f05              SWIRE ELAIS S1 Region 70um Fall '05 Spitzer Catalog
elaiss1_160_cat_f05             SWIRE ELAIS S1 Region 160um Fall '05 Spitzer Catalog
lockman_cat_s05                 SWIRE Lockman Region Spring '05 SWIRE Spitzer Catalog
swire_lhisod                    SWIRE Lockman Hole ISOCAM Deep Field Catalog
swire_lhisos                    SWIRE Lockman Hole ISOCAM Shallow Field Catalog
lockman_24_cat_s05              SWIRE Lockman Region 24um Spring '05 Spitzer Catalog
lockman_70_cat_s05              SWIRE Lockman Region 70um Spring '05 Spitzer Catalog
lockman_160_cat_s05             SWIRE Lockman Region 160um Spring '05 Spitzer Catalog
xmm_cat_s05                     SWIRE XMM_LSS Region Spring '05 Spitzer Catalog
xmm_24_cat_s05                  SWIRE XMM_LSS Region 24um Spring '05 Spitzer Catalog
xmm_70_cat_s05                  SWIRE XMM_LSS Region 70um Spring '05 Spitzer Catalog
xmm_160_cat_s05                 SWIRE XMM_LSS Region 160um Spring '05 Spitzer Catalog
taurus_2008_2_1                 Taurus Catalog October 2008 v2.1
ysoggd1215obj                   YSOVAR GGD 12-15 Object Table
ysoggd1215lc                    YSOVAR GGD 12-15 Light Curve Table
ysoi20050obj                    YSOVAR IRAS 20050+2720 Object Table
ysoi20050lc                     YSOVAR IRAS 20050+2720 Light Curve Table
ysol1688obj                     YSOVAR L1688 Object Table
ysol1688lc                      YSOVAR L1688 Light Curve Table
yson1333obj                     YSOVAR NGC1333 Object Table
yson1333lc                      YSOVAR NGC1333 Light Curve Table
msxc6                           The Midcourse Space Experiment (MSXC6)
msxc6_rej                       The Midcourse Space Experiment (MSXC6) Rejects
cosmos_chandra_bsc21            Chandra-COSMOS Bright Source Catalog v2.1
galex_emphot_v3                 GALEX/COSMOS Prior-based Photometry Catalog June 2008
acs_iphot_sep07                 COSMOS ACS I-band photometry catalog September 2007
cosmos_ib_phot                  COSMOS Intermediate and Broad Band Photometry Catalog 2008
cosmos_phot                     COSMOS Photometry Catalog January 2006
cosmos_zphot_mag25              COSMOS Photometric Redshift Catalog Fall 2008 (README - mag 25 limited)
cosmos_morph_cassata_1_1        COSMOS Cassata Morphology Catalog v1.1
cosmos_morph_tasca_1_1          COSMOS Tasca Morphology Catalog v1.1
morphology_2005                 COSMOS Morphology Catalog 2005
cosmos_morph_col_1              COSMOS Zamojski Morphology Catalog v1.0
cosmos_morph_zurich_1           COSMOS Zurich Structure and Morphology Catalog v1.0
scosmos_mips_24_go3             S-COSMOS MIPS 24 Photometry Catalog October 2008
cosmos_vla_deep_may2010         COSMOS VLA Deep Catalog May 2010
cosmos327                       COSMOS VLA 327 MHz Catalog
cosmos_xmm_2                    COSMOS XMM Point-like Source Catalog v2.0
cosmos_xgroups                  COSMOS X-ray Group Catalog
cosmos_xgal                     COSMOS X-ray Group Member Catalog
bolocamv21                      BOLOCAM Galactic Plane Survey Catalog v2.1
bolocam_gps_v1_0_1              BOLOCAM Galactic Plane Survey Catalog
bolocamdist                     BOLOCAM Galactic Plane Survey Distance Catalog
akari_fis                       Akari/FIS Bright Source Catalogue
akari_irc                       Akari/IRC Point Source Catalogue
usno_b1                         USNO-B1 (United States Naval Observatory B1.0 Catalog)
ucac4_sources                   USNO CCD Astrograph Catalog (UCAC4)
urat1                           The First USNO Robotic Astrometric Telescope Catalog (URAT1)
denis3                          DENIS 3rd Release (Sep. 2005)
spsc250                         SPIRE Point Source Catalog: 250 microns
spsc350                         SPIRE Point Source Catalog: 350 microns
spsc500                         SPIRE Point Source Catalog: 500 microns
spscxid                         SPIRE Point Source Catalog Cross-Reference Matrix
acmccat                         ACMC Catalog
dunes                           DUNES Catalog
hgoodsn                         GOODS North Catalog
hgoodss                         GOODS South Catalog
heritagelphot                   HERITAGE LMC Band-Matched Catalog
heritagelclass                  HERITAGE LMC Band-Matched Classification Table
heritagel100                    HERITAGE LMC PACS 100 micron Catalog
heritagel160                    HERITAGE LMC PACS 160 micron Catalog
heritagel250                    HERITAGE LMC SPIRE 250 micron Catalog
heritagel350                    HERITAGE LMC SPIRE 350 micron Catalog
heritagel500                    HERITAGE LMC SPIRE 500 micron Catalog
heritagesphot                   HERITAGE SMC Band-Matched Catalog
heritagesclass                  HERITAGE SMC Band-Matched Classification Table
heritages100                    HERITAGE SMC PACS 100 micron Catalog
heritages160                    HERITAGE SMC PACS 160 micron Catalog
heritages250                    HERITAGE SMC SPIRE 250 micron Catalog
heritages350                    HERITAGE SMC SPIRE 350 micron Catalog
heritages500                    HERITAGE SMC SPIRE 500 micron Catalog
hermesscat250                   HerMES 250 micron StarFinder Catalog
hermesscat250sxt                HerMES 250 micron SUSSEXtractor Catalog
hermesscat350                   HerMES 350 micron StarFinder Catalog
hermesscat350sxt                HerMES 350 micron SUSSEXtractor Catalog
hermesscat500                   HerMES 500 micron StarFinder Catalog
hermesscat500sxt                HerMES 500 micron SUSSEXtractor Catalog
hermesxid250                    HerMES Band-merged Catalog (250 micron positions)
hermesxid24                     HerMES Band-merged Catalog (24 micron positions)
hopsdata                        Herschel Orion Protostars Survey SED Data Catalog
hopsfits                        Herschel Orion Protostars Survey SED Fits Catalog
hpdpsso                         HPDP SSO Table
pep100                          PEP PACS 100 micron Catalog
pep160                          PEP PACS 160 micron Catalog
pepprior                        PEP PACS Extractions Using MIPS 24 micron Priors
pepgoodss70                     PEP PACS 70 micron GOODS-S Catalog
pepxid                          PEP PACS and MIPS Cross-IDs Catalog
pep250                          PEP SPIRE 250 micron Catalog
pep350                          PEP SPIRE 350 micron Catalog
pep500                          PEP SPIRE 500 micron Catalog
peplh24                         PEP Lockman Hole MIPS 24 micron Comparison Catalog
ppmxl                           PPMXL: A Proper Motion Catalog Combining USNO-B and 2MASS
brava                           BRAVA Catalog
musyc_photz                     MUSYC Photometric Redshift Catalog
musyc_phot                      MUSYC Photometry Catalog
ptf_lightcurves                 PTF Lightcurve Table
ptf_sources                     PTF Sources Catalog
ptf_objects                     PTF Objects
ptfphotcalcat                   PTF Photometric Calibrator Catalog
gaia_source_dr1                 Gaia Catalogue
schemas                         Available schemas at NASA/IPAC Infrared Science Archive
tables                          Available tables at NASA/IPAC Infrared Science Archive
columns                         Available columns at NASA/IPAC Infrared Science Archive
keys                            Available Keys at NASA/IPAC Infrared Science Archive
key_columns                     Available Key Columns at NASA/IPAC Infrared Science Archive
```



[1]: https://astroquery.readthedocs.io
[2]: http://docs.astropy.org/en/stable/coordinates/matchsep.html#matching-catalogs