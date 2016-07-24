"""
Run integration tests from pixsim through redshifts

python -m desispec.test.integration_test
"""
from __future__ import absolute_import, print_function
import sys
import os
import random
import time
import subprocess as sp
import glob

import numpy as np
from astropy.io import fits

import desispec.pipeline as pipe
import desispec.io as io
import desispec.log as logging

import desispec.scripts.makebricks as makebricks


#- prevent nose from trying to run this test since it takes too long
__test__ = False


def check_env():
    """
    Check required environment variables; raise RuntimeException if missing
    """
    log = logging.get_logger()
    #- template locations
    missing_env = False
    if 'DESI_BASIS_TEMPLATES' not in os.environ:
        log.warning('missing $DESI_BASIS_TEMPLATES needed for simulating spectra'.format(name))
        missing_env = True

    if not os.path.isdir(os.getenv('DESI_BASIS_TEMPLATES')):
        log.warning('missing $DESI_BASIS_TEMPLATES directory')
        log.warning('e.g. see NERSC:/project/projectdirs/desi/spectro/templates/basis_templates/v1.0')
        missing_env = True

    for name in (
        'DESI_SPECTRO_SIM', 'DESI_SPECTRO_REDUX', 'PIXPROD', 'SPECPROD', 'DESIMODEL'):
        if name not in os.environ:
            log.warning("missing ${0}".format(name))
            missing_env = True

    if missing_env:
        log.warning("Why are these needed?")
        log.warning("    Simulations written to $DESI_SPECTRO_SIM/$PIXPROD/")
        log.warning("    Raw data read from $DESI_SPECTRO_DATA/")
        log.warning("    Spectro pipeline output written to $DESI_SPECTRO_REDUX/$SPECPROD/")
        log.warning("    Templates are read from $DESI_BASIS_TEMPLATES")

    #- Wait until end to raise exception so that we report everything that
    #- is missing before actually failing
    if missing_env:
        log.critical("missing env vars; exiting without running pipeline")
        sys.exit(1)

    #- Override $DESI_SPECTRO_DATA to match $DESI_SPECTRO_SIM/$PIXPROD
    os.environ['DESI_SPECTRO_DATA'] = os.path.join(os.getenv('DESI_SPECTRO_SIM'), os.getenv('PIXPROD'))


# Simulate raw data

def sim(night, nspec=5, clobber=False):
    """
    Simulate data as part of the integration test.

    Args:
        night (str): YEARMMDD
        nspec (int, optional): number of spectra to include
        clobber (bool, optional): rerun steps even if outputs already exist
        
    Raises:
        RuntimeError if any script fails
    """
    log = logging.get_logger()

    # Create input fibermaps, spectra, and pixel-level raw data

    for expid, flavor in zip([0,1,2], ['flat', 'arc', 'dark']):
        cmd = "newexp-desi --flavor {flavor} --nspec {nspec} --night {night} --expid {expid}".format(
            expid=expid, flavor=flavor, nspec=nspec, night=night)
        fibermap = io.findfile('fibermap', night, expid)
        simspec = '{}/simspec-{:08d}.fits'.format(os.path.dirname(fibermap), expid)
        inputs = []
        outputs = [fibermap, simspec]
        if pipe.runcmd(cmd, inputs=inputs, outputs=outputs, clobber=clobber) != 0:
            raise RuntimeError('pixsim newexp failed for {} exposure {}'.format(flavor, expid))

        cmd = "pixsim-desi --preproc --nspec {nspec} --night {night} --expid {expid}".format(expid=expid, nspec=nspec, night=night)
        inputs = [fibermap, simspec]
        outputs = list()
        outputs.append(fibermap.replace('fibermap-', 'simpix-'))
        for camera in ['b0', 'r0', 'z0']:
            pixfile = io.findfile('pix', night, expid, camera)
            outputs.append(pixfile)
            #outputs.append(os.path.join(os.path.dirname(pixfile), os.path.basename(pixfile).replace('pix-', 'simpix-')))
        if pipe.runcmd(cmd, inputs=inputs, outputs=outputs, clobber=clobber) != 0:
            raise RuntimeError('pixsim failed for {} exposure {}'.format(flavor, expid))

<<<<<<< HEAD
    #-----
    #- Extract

    waverange = dict(
        b = "3570,5940,1.0",
        r = "5630,7740,1.0",
        z = "7440,9830,1.0",
        )
    for expid in [0,1,2]:
        for channel in ['b', 'r', 'z']:
            camera = channel+'0'
            pixfile = io.findfile('pix', night, expid, camera)
            fiberfile = io.findfile('fibermap', night, expid)
            psffile = '{}/data/specpsf/psf-{}.fits'.format(os.getenv('DESIMODEL'), channel)
            framefile = io.findfile('frame', night, expid, camera)
            # cmd = "exspec -i {pix} -p {psf} --specmin 0 --nspec {nspec} -w {wave} -o {frame}".format(
            #     pix=pixfile, psf=psffile, wave=waverange[channel], frame=framefile, **params)
            cmd = "desi_extract_spectra -i {pix} -p {psf} -f {fibermap} --specmin 0 --nspec {nspec} -o {frame}".format(
                pix=pixfile, psf=psffile, frame=framefile, fibermap=fiberfile, **params)

            inputs = [pixfile, psffile, fiberfile]
            outputs = [framefile,]
            if runcmd(cmd, inputs, outputs, clobber) != 0:
                raise RuntimeError('extraction failed for {} expid {}'.format(camera, expid))

    #-----
    #- Fiber flat
    expid = 0
    for channel in ['b', 'r', 'z']:
        camera = channel+"0"
        framefile = io.findfile('frame', night, expid, camera)
        fiberflat = io.findfile('fiberflat', night, expid, camera)
        fibermap = io.findfile('fibermap', night, expid)  # for QA
        qafile = io.findfile('qa_calib', night, expid, camera)
        qafig = io.findfile('qa_flat_fig', night, expid, camera)
        cmd = "desi_compute_fiberflat --infile {frame} --fibermap {fibermap} --outfile {fiberflat} --qafile {qafile} --qafig {qafig}".format(
            frame=framefile, fibermap=fibermap, fiberflat=fiberflat, qafile=qafile, qafig=qafig, **params)
        inputs = [framefile,fibermap,]
        outputs = [fiberflat,qafile,qafig,]
        if runcmd(cmd, inputs, outputs, clobber) != 0:
            raise RuntimeError('fiberflat failed for '+camera)

    #-----
    #- Sky model
    flat_expid = 0
    expid = 2
    for channel in ['b', 'r', 'z']:
        camera = channel+"0"
        framefile = io.findfile('frame', night, expid, camera)
        fibermap = io.findfile('fibermap', night, expid)
        fiberflat = io.findfile('fiberflat', night, flat_expid, camera)
        skyfile = io.findfile('sky', night, expid, camera)
        qafile = io.findfile('qa_data', night, expid, camera)
        qafig = io.findfile('qa_sky_fig', night, expid, camera)
        cmd="desi_compute_sky --infile {frame} --fibermap {fibermap} --fiberflat {fiberflat} --outfile {sky} --qafile {qafile} --qafig {qafig}".format(
            frame=framefile, fibermap=fibermap, fiberflat=fiberflat, sky=skyfile, qafile=qafile, qafig=qafig, **params)
        inputs = [framefile, fibermap, fiberflat]
        outputs = [skyfile, qafile, qafig,]
        if runcmd(cmd, inputs, outputs, clobber) != 0:
            raise RuntimeError('sky model failed for '+camera)


    #-----
    #- Fit standard stars
    if 'STD_TEMPLATES' in os.environ:
        std_templates = os.getenv('STD_TEMPLATES')
    else:
        std_templates = os.getenv('DESI_ROOT')+'/spectro/templates/star_templates/v1.0/stdstar_templates_v1.0.fits'

    stdstarfile = io.findfile('stdstars', night, expid, spectrograph=0)
    cmd = """desi_fit_stdstars --spectrograph 0 \
      --fibermap {fibermap} \
      --fiberflatexpid {flat_expid} \
      --models {std_templates} --outfile {stdstars}""".format(
        flat_expid=flat_expid, fibermap=fibermap, std_templates=std_templates,
        stdstars=stdstarfile)

    inputs = [fibermap, std_templates]
    outputs = [stdstarfile,]
    if runcmd(cmd, inputs, outputs, clobber) != 0:
        raise RuntimeError('fitting stdstars failed')


    #-----
    #- Flux calibration
    for channel in ['b', 'r', 'z']:
        camera = channel+"0"
        framefile = io.findfile('frame', night, expid, camera)
        fibermap  = io.findfile('fibermap', night, expid)
        fiberflat = io.findfile('fiberflat', night, flat_expid, camera)
        skyfile   = io.findfile('sky', night, expid, camera)
        calibfile = io.findfile('calib', night, expid, camera)
        qafile = io.findfile('qa_data', night, expid, camera)
        qafig = io.findfile('qa_flux_fig', night, expid, camera)

        #- Compute flux calibration vector
        cmd = """desi_compute_fluxcalibration \
          --infile {frame} --fibermap {fibermap} --fiberflat {fiberflat} --sky {sky} \
          --models {stdstars} --outfile {calib} --qafile {qafile} --qafig {qafig}""".format(
            frame=framefile, fibermap=fibermap, fiberflat=fiberflat, sky=skyfile,
            stdstars=stdstarfile, calib=calibfile, qafile=qafile, qafig=qafig
            )
        inputs = [framefile, fibermap, fiberflat, skyfile, stdstarfile]
        outputs = [calibfile, qafile, qafig]
        if runcmd(cmd, inputs, outputs, clobber) != 0:
            raise RuntimeError('flux calibration failed for '+camera)

        #- Apply the flux calibration to write a cframe file
        cframefile = io.findfile('cframe', night, expid, camera)
        cmd = """desi_process_exposure \
          --infile {frame} --fiberflat {fiberflat} --sky {sky} --calib {calib} \
          --outfile {cframe}""".format(frame=framefile, fibermap=fibermap,
            fiberflat=fiberflat, sky=skyfile, calib=calibfile, cframe=cframefile)
        inputs = [framefile, fiberflat, skyfile, calibfile]
        outputs = [cframefile, ]
        if runcmd(cmd, inputs, outputs, clobber) != 0:
            raise RuntimeError('combining calibration steps failed for '+camera)

    #-----
    #- Bricks
    inputs = list()
    for camera in ['b0', 'r0', 'z0']:
        inputs.append( io.findfile('cframe', night, expid, camera) )

    outputs = list()
    fibermap = io.read_fibermap(io.findfile('fibermap', night, expid))
    bricks = set(fibermap['BRICKNAME'])
    for b in bricks:
        for channel in ['b', 'r', 'z']:
            outputs.append( io.findfile('brick', brickname=b, band=channel))
=======
    return
>>>>>>> 9f7a6f2c35a3dba59c29aee81ec75a224615bbb0


<<<<<<< HEAD
    #-----
    #- Redshifts!
    for b in bricks:
        inputs = [io.findfile('brick', brickname=b, band=channel) for channel in ['b', 'r', 'z']]
        zbestfile = io.findfile('zbest', brickname=b)
        outputs = [zbestfile, ]
        cmd = "desi_zfind --brick {} -o {}".format(b, zbestfile)
        if runcmd(cmd, inputs, outputs, clobber) != 0:
            raise RuntimeError('redshifts failed for brick '+b)
    # ztruth QA
    qafile = io.findfile('qa_ztruth', night)
    qafig = io.findfile('qa_ztruth_fig', night)
    cmd = "desi_qa_zfind --night {night} --qafile {qafile} --qafig {qafig}".format(
        night=night, qafile=qafile, qafig=qafig)
    inputs = []
    outputs = [qafile, qafig]
    if runcmd(cmd, inputs, outputs, clobber) != 0:
        raise RuntimeError('redshift QA failed for night '+night)

    #-----
    #- Did it work?
    #- (this combination of fibermap, simspec, and zbest is a pain)
    simdir = os.path.dirname(io.findfile('fibermap', night=night, expid=expid))
=======
def integration_test(night=None, nspec=5, clobber=False):
    """Run an integration test from raw data simulations through redshifts
    
    Args:
        night (str, optional): YEARMMDD, defaults to current night
        nspec (int, optional): number of spectra to include
        clobber (bool, optional): rerun steps even if outputs already exist
        
    Raises:
        RuntimeError if any script fails
      
    """
    log = logging.get_logger()
    log.setLevel(logging.DEBUG)

    # YEARMMDD string, rolls over at noon not midnight
    # TODO: fix usage of night to be something other than today
    if night is None:
        #night = time.strftime('%Y%m%d', time.localtime(time.time()-12*3600))
        night = "20160726"

    # check for required environment variables
    check_env()

    # simulate inputs
    sim(night, nspec=nspec, clobber=clobber)

    # create production

    # FIXME:  someday run PSF estimation too...
    ### com = "desi_pipe --env env.txt --spectrographs 0 --fakeboot --fakepsf"
    rawdir = os.path.join(os.getenv('DESI_SPECTRO_SIM'), os.getenv('PIXPROD'))    
    com = "desi_pipe --spectrographs 0 --fakeboot --fakepsf --raw {}".format(rawdir)
    sp.check_call(com, shell=True)

    # raw and production locations

    rawdir = os.path.abspath(io.rawdata_root())
    proddir = os.path.abspath(io.specprod_root())

    # Modify options file to restrict the spectral range

    optpath = os.path.join(proddir, "run", "options.yaml")
    opts = pipe.read_options(optpath)
    opts['extract']['specmin'] = 0
    opts['extract']['nspec'] = nspec
    pipe.write_options(optpath, opts)

    # run the generated shell scripts

    # FIXME:  someday run PSF estimation too...

    # print("Running bootcalib script...")
    # com = os.path.join(proddir, "run", "scripts", "bootcalib_all.sh")
    # sp.check_call(["bash", com])

    # print("Running specex script...")
    # com = os.path.join(proddir, "run", "scripts", "specex_all.sh")
    # sp.check_call(["bash", com])

    # print("Running psfcombine script...")
    # com = os.path.join(proddir, "run", "scripts", "psfcombine_all.sh")
    # sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "extract_all.sh")
    print("Running extraction script "+com)
    sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "fiberflat-procexp_all.sh")
    print("Running calibration script "+com)
    sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "bricks.sh")
    print("Running makebricks script "+com)
    sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "zfind_all.sh")
    print("Running zfind script "+com)
    sp.check_call(["bash", com])

    # #-----
    # #- Did it work?
    # #- (this combination of fibermap, simspec, and zbest is a pain)
    expid = 2
    fmfile = io.findfile('fibermap', night=night, expid=expid)
    fibermap = io.read_fibermap(fmfile)
    simdir = os.path.dirname(fmfile)
>>>>>>> 9f7a6f2c35a3dba59c29aee81ec75a224615bbb0
    simspec = '{}/simspec-{:08d}.fits'.format(simdir, expid)
    siminfo = fits.getdata(simspec, 'METADATA')

    brickdirs = glob.glob(os.path.join(proddir, "bricks", "*"))
    bricks = [ os.path.basename(x) for x in brickdirs ]

    print()
    print("--------------------------------------------------")
    print("Brick     True  z        ->  Class  z        zwarn")
    # print("3338p190  SKY   0.00000  ->  QSO    1.60853   12   - ok")
    for b in bricks:
        zbest = io.read_zbest(io.findfile('zbest', brickname=b))        
        for i in range(len(zbest.z)):
            if zbest.spectype[i] == 'ssp_em_galaxy':
                objtype = 'GAL'
            elif zbest.spectype[i] == 'spEigenStar':
                objtype = 'STAR'
            else:
                objtype = zbest.spectype[i]

            z, zwarn = zbest.z[i], zbest.zwarn[i]

            j = np.where(fibermap['TARGETID'] == zbest.targetid[i])[0][0]
            truetype = siminfo['OBJTYPE'][j]
            truez = siminfo['REDSHIFT'][j]
            dv = 3e5*(z-truez)/(1+truez)
            if truetype == 'SKY' and zwarn > 0:
                status = 'ok'
            elif zwarn == 0:
                if truetype == 'LRG' and objtype == 'GAL' and abs(dv) < 150:
                    status = 'ok'
                elif truetype == 'ELG' and objtype == 'GAL' and abs(dv) < 150:
                    status = 'ok'
                elif truetype == 'QSO' and objtype == 'QSO' and abs(dv) < 750:
                    status = 'ok'
                elif truetype == 'STD' and objtype == 'STAR':
                    status = 'ok'
                else:
                    status = 'OOPS'
            else:
                status = 'OOPS'
            print('{0}  {1:4s} {2:8.5f}  -> {3:5s} {4:8.5f} {5:4d}  - {6}'.format(
                b, truetype, truez, objtype, z, zwarn, status))

    print("--------------------------------------------------")




if __name__ == '__main__':
    integration_test()
