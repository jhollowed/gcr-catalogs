"""
Alpha Q galaxy catalog class.
"""
from __future__ import division
import os
import numpy as np
import h5py
from astropy.cosmology import FlatLambdaCDM
from GCR import BaseGenericCatalog
import pdb
from .alphaq import AlphaQGalaxyCatalog

__all__ = ['AlphaQClusterCatalog']

class AlphaQClusterCatalog(AlphaQGalaxyCatalog):
    """
    The galaxy cluster catalog. Inherits AlphaQGalaxyCatalog, overloading select methods.

    The AlphaQ cluster catalog is structured in the following way: under the root hdf group, there
    is a group per each halo with SO mass above 1e14 M_sun/h. Each of these groups contains the same
    datasets as the original AlphaQ galaxy catalog, but with only as many rows as member galaxies for
    the halo in question. Each group has attributes which contain halo-wide quantities, such as mass,
    position, etc.

    This class offers filtering on any halo quantity (group attribute), as seen in all three of the
    methods of this class (all the group attributes are iterated over in contexts concerning the
    pre-filtering). The valid filtering quantities are:
    {'host_halo_mass', 'sod_halo_cdelta', 'sod_halo_cdelta_error', 'sod_halo_c_acc_mass',
     'fof_halo_tag', 'halo_index', 'halo_step', 'halo_ra', 'halo_dec', 'halo_z',
     'halo_z_err', 'sod_halo_radius', 'sod_halo_mass', 'sod_halo_ke', 'sod_halo_vel_disp'}
    """


    def _subclass_init(self, filename, **kwargs):
        super(AlphaQClusterCatalog, self)._subclass_init(filename, **kwargs)
        #with h5py.File(self._file, 'r') as fh:
            #pdb.set_trace()
            #self._native_filter_quantities = set(fh[next(fh.keys().__iter__())].attrs)


    def _iter_native_dataset(self, native_filters=None):
        with h5py.File(self._file, 'r') as fh:
            for key in fh['clusters']:
                halo = fh['clusters'][key]

                if native_filters and not all(f[0](*(halo.attrs[k] for k in f[1:])) for f in native_filters):
                    continue

                def native_quantity_getter(native_quantity):
                    raise NotImplementedError

                yield native_quantity_getter
    

    def _generate_native_quantity_list(self):
        with h5py.File(self._file, 'r') as fh:
            # quantities are the same for all halos - we use the first here
            hgroup = fh['clusters']['halo_1']
            hattrs = list(hgroup.attrs.keys())
            hobjects = []
            #get all the names of objects in this tree
            hgroup.visit(hobjects.append)
            #filter out the group objects and keep the dataset objects
            hdatasets = [hobject for hobject in hobjects if type(hgroup[hobject]) != h5py.Group]
            hdatasets = hdatasets + hattrs
            pdb.set_trace()
            native_quantities = set(hdatasets)
        return native_quantities


    def _iter_native_dataset(self, native_filters=None):
        assert not native_filters, '*native_filters* is not supported'
        with h5py.File(self._file, 'r') as fh:
            def native_quantity_getter(native_quantity):
                return fh['galaxyProperties/{}'.format(native_quantity)].value
            yield native_quantity_getter


    def _get_native_quantity_info_dict(self, quantity, default=None):
        with h5py.File(self._file,'r') as fh:
            if 'halo1/'+quantity not in fh['clusters']:
                return default
            else:
                info_dict = dict()
                for key in fh['clusters']['halo_1/'+quantity].attrs:
                    info_dict[key] = fh['clusters']['halo_1/'+quantity].attrs[key]
                return info_dict
