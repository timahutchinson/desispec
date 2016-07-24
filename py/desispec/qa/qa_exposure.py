"""
Classes to organize and execute QA for a DESI exposure
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from desispec.log import get_logger

log=get_logger()

<<<<<<< HEAD
        Notes:

        Attributes:
            All input args become object attributes.
        """
        # Parse from frame.meta
        if frame is not None:
            # Parse from meta, if possible
            try:
                flavor = frame.meta['FLAVOR']
            except:
                pass
            else:
                try:
                    camera = frame.meta['CAMERA']
                except KeyError:
                    pass

        assert flavor in ['none', 'flat', 'arc', 'dark', 'bright', 'bgs', 'mws', 'lrg', 'elg', 'qso']
        self.flavor = flavor
        self.camera = camera
        
        # Initialize data
        if in_data is None:
            self.data = dict(flavor=self.flavor, camera=self.camera)
        else:
            assert isinstance(in_data,dict)
            self.data = in_data

    def init_qatype(self, qatype, param, re_init=False):
        """Initialize parameters for a given qatype
        qatype: str  
          Type of QA to be performed (e.g. SKYSUB)
        param: dict
          Dict of parameters to guide QA
        re_init: bool, (optional)
          Re-initialize parameter dict
          Code will always add new parameters if any exist
        """
        # Fill and return if not set previously or if re_init=True
        if (qatype not in self.data.keys()) or re_init: 
            self.data[qatype] = {}
            self.data[qatype]['PARAM'] = param
            return

        # Update the new parameters only
        for key in param.keys():
            if key not in self.data[qatype]['PARAM'].keys():
                self.data[qatype]['PARAM'][key] = param[key]

    def init_fiberflat(self, re_init=False):
        """Initialize parameters for FIBERFLAT QA
        QA method is desispec.fiberflat.qa_fiberflat

        Parameters:
        ------------
        re_init: bool, (optional)
          Re-initialize FIBERFLAT parameter dict
        """
        #
        assert self.flavor in ['flat']

        # Standard FIBERFLAT input parameters
        fflat_dict = dict(MAX_N_MASK=20000,  # Maximum number of pixels to mask
                          MAX_SCALE_OFF=0.05,  # Maximum offset in counts (fraction)
                          MAX_OFF=0.15,       # Maximum offset from unity
                          MAX_MEAN_OFF=0.05,  # Maximum offset in mean of fiberflat
                          MAX_RMS=0.02,      # Maximum RMS in fiberflat
                          )
        # Init
        self.init_qatype('FIBERFLAT', fflat_dict, re_init=re_init)

    def init_fluxcalib(self, re_init=False):
        """ Initialize parameters for FLUXCALIB QA
        Args:
            re_init: bool, (optional)
              Re-initialize  parameter dict

        Returns:

        """
        assert self.flavor in ['dark','bright','bgs','mws','lrg','elg','qso']

        # Standard FLUXCALIB input parameters
        flux_dict = dict(ZP_WAVE=0.,        # Wavelength for ZP evaluation (camera dependent)
                         MAX_ZP_OFF=0.2,    # Max offset in ZP for individual star
                         )

        if self.camera[0] == 'b':
            flux_dict['ZP_WAVE'] = 4800.  # Ang
        elif self.camera[0] == 'r':
            flux_dict['ZP_WAVE'] = 6500.  # Ang
        elif self.camera[0] == 'z':
            flux_dict['ZP_WAVE'] = 8250.  # Ang
        else:
            raise ValueError("Not ready for this camera!")

        # Init
        self.init_qatype('FLUXCALIB', flux_dict, re_init=re_init)

    def init_skysub(self, re_init=False):
        """Initialize parameters for SkySub QA 
        QA method is desispec.sky.qa_skysub

        Parameters:
        ------------
        re_init: bool, (optional)
          Re-initialize SKYSUB parameter dict
        """
        #
        assert self.flavor in ['dark','bright','bgs','mws','lrg','elg','qso']

        # Standard SKYSUB input parameters
        sky_dict = dict(
            PCHI_RESID=0.05, # P(Chi^2) limit for bad skyfiber model residuals
            PER_RESID=95.,   # Percentile for residual distribution
            )
        # Init
        self.init_qatype('SKYSUB', sky_dict, re_init=re_init)

    def run_qa(self, qatype, inputs, clobber=True):
        """Run QA tests of a given type
        Over-writes previous QA of this type, unless otherwise specified

        qatype: str  
          Type of QA to be performed (e.g. SKYSUB)
        inputs: tuple
          Set of inputs for the tests
        clobber: bool, optional [True]
          Over-write previous QA 
        """
        from desispec.sky import qa_skysub
        from desispec.fiberflat import qa_fiberflat
        from desispec.fluxcalibration import qa_fluxcalib

        # Check for previous QA if clobber==False
        if not clobber:
            # QA previously performed?
            if 'QA' in self.data[qatype].keys():
                return
        # Run
        if qatype == 'SKYSUB':
            # Expecting: frame, skymodel
            assert len(inputs) == 2
            # Init parameters (as necessary)
            self.init_skysub()
            # Run
            qadict = qa_skysub(self.data[qatype]['PARAM'],
                inputs[0], inputs[1])
        elif qatype == 'FIBERFLAT':
            # Expecting: frame, fiberflat
            assert len(inputs) == 2
            # Init parameters (as necessary)
            self.init_fiberflat()
            # Run
            qadict = qa_fiberflat(self.data[qatype]['PARAM'], inputs[0], inputs[1])
        elif qatype == 'FLUXCALIB':
            # Expecting: frame, fluxcalib, individual_outputs (star by star)
            assert len(inputs) == 3
            # Init parameters (as necessary)
            self.init_fluxcalib()
            # Run
            qadict = qa_fluxcalib(self.data[qatype]['PARAM'],
                                  inputs[0], inputs[1], inputs[2])
        else:
            raise ValueError('Not ready to perform {:s} QA'.format(qatype))
        # Update
        self.data[qatype]['QA'] = qadict

    def __repr__(self):
        """
        Print formatting
        """
        return ('{:s}: camera={:s}, flavor={:s}'.format(
                self.__class__.__name__, self.camera, self.flavor))
=======
>>>>>>> 9f7a6f2c35a3dba59c29aee81ec75a224615bbb0

class QA_Exposure(object):
    def __init__(self, expid, night, flavor, specprod_dir=None, in_data=None, **kwargs):
        """
        Class to organize and execute QA for a DESI Exposure

        x.flavor, x.data
        
        Args:
            expid: int -- Exposure ID
            night: str -- YYYYMMDD
            flavor: str
              exposure type (e.g. flat, arc, science)
            specprod_dir(str): Path containing the exposures/ directory to use. If the value
                is None, then the value of :func:`specprod_root` is used instead.
            in_data: dict, optional -- Input data
              Mainly for reading from disk

        Notes:

        Attributes:
            All input args become object attributes.
        """
        flavors = ['none', 'flat', 'arc', 'gray', 'dark', 'bright', 'bgs', 'mws', 'lrg', 'elg', 'qso']
        assert flavor in flavors
        if flavor in ['dark', 'bright', 'bgs', 'mws', 'lrg', 'elg', 'qso']:
            self.type = 'data'
        else:
            self.type = 'calib'

        self.expid = expid
        self.night = night
        self.specprod_dir = specprod_dir
        self.flavor = flavor

        if in_data is None:
            self.data = dict(flavor=self.flavor, expid=self.expid,
                             night=self.night, frames={})
            self.load_qa_data(**kwargs)
        else:
            assert isinstance(in_data,dict)
            self.data = in_data

    def fluxcalib(self, outfil):
        """ Perform QA on fluxcalib results for an Exposure

        Args:
            outfil: str -- Filename for PDF  (may automate)

        Independent results for each channel
        """
        from . import qa_plots
        # Init
        if 'FLUXCALIB' not in self.data.keys():
            self.data['FLUXCALIB'] = {}
        # Loop on channel
        cameras = self.data['frames'].keys()
        for channel in ['b','r','z']:
            # Init
            if channel not in self.data['FLUXCALIB'].keys():
                self.data['FLUXCALIB'][channel] = {}
            # Load
            ZPval = []
            for camera in cameras:
                if camera[0] == channel:
                    ZPval.append(self.data['frames'][camera]['FLUXCALIB']['QA']['ZP'])
            # Measure RMS
            if len(ZPval) > 0:
                self.data['FLUXCALIB'][channel]['ZP_RMS'] = np.std(ZPval)

        # Figure
        qa_plots.exposure_fluxcalib(outfil, self.data)

    def load_qa_data(self, remove=False):
        """ Load the QA data files for a given exposure (currently yaml)
        Args:
            remove: bool, optional
              Remove QA frame files
        """

        from desispec import io as desiio
        qafiles = desiio.get_files(filetype='qa_'+self.type, night=self.night,
                                  expid=self.expid,
                                  specprod_dir=self.specprod_dir)
        #import pdb; pdb.set_trace()
        # Load into frames
        for camera,qadata_path in qafiles.iteritems():
            qa_frame = desiio.load_qa_frame(qadata_path)
            # Remove?
            if remove:
                #import pdb; pdb.set_trace()
                os.remove(qadata_path)
            # Test
            for key in ['expid','night']:
                assert getattr(qa_frame,key) == getattr(self, key)
            # Save
            self.data['frames'][camera] = qa_frame.qa_data

    def __repr__(self):
        """ Print formatting
        """
        return ('{:s}: night={:s}, expid={:d}, type={:s}, flavor={:s}'.format(
                self.__class__.__name__, self.night, self.expid, self.type, self.flavor))
