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

        # add halo property quantity modifiers
        for key in self._quantity_modifiers.keys():
            self._quantity_modifiers[key] = 'galaxyProperties/{}'.format(self._quantity_modifiers[key])
        self._quantity_modifiers.update({
            'GAL_VELOCITY_DISPERSION':  'haloProperties/gal_vel_disp',
            'GAL_VELOCITY_DISPERSION_OBS': 'haloProperties/gal_vel_disp_obs',
            'PARTICLE_VELOCITY_DISPERSION': 'haloProperties/sod_halo_vel_disp',
            'HALO_REDSHIFT': 'haloProperties/halo_z',
            'HALO_REDSHIFT_ERROR': 'haloProperties/halo_z_err',
            'HALO_DEC': 'haloProperties/halo_dec',
            'HALO_RA': 'haloProperties_halo_ra',
            'HALO_ID': 'haloProperties/halo_index',
            'HALO_MASS': 'haloProperties/halo_mass',
            'HALO_SOD_TAG': 'haloProperties/fof_halo_tag',
            'HALO_SOD_MASS_M200': 'haloProperties/sod_halo_mass',
            'HALO_SOD_RADIUS': 'haloProperties/sod_halo_radius',
            'HALO_SOD_CDELTA': 'haloProperties/sod_halo_cdelta',
            'HALO_SOD_CDELTA_ERROR': 'haloProperties/sod_halo_cdelta_error',
            'HALO_SOD_KE': 'haloProperties/sod_halo_ke'
        })

    
    def _iter_native_dataset(self, native_filters=None):
        assert not native_filters, '*native_filters* is not supported'
        with h5py.File(self._file, 'r') as fh:
            def native_quantity_getter(native_quantity):
                halos = fh['clusters'].keys()
                return np.array([fh['clusters'][next_halo][native_quantity].value for next_halo in halos])
            yield native_quantity_getter


    def _generate_native_quantity_list(self):
        with h5py.File(self._file, 'r') as fh:
            # quantities are the same for all halos - we use the first here
            hgroup = fh['clusters']['halo_1']
            hobjects = []
            #get all the names of objects in this tree
            hgroup.visit(hobjects.append)
            #filter out the group objects and keep the dataset objects
            hdatasets = [hobject for hobject in hobjects if type(hgroup[hobject]) != h5py.Group]
            native_quantities = set(hdatasets)
        return native_quantities


    def _get_native_quantity_info_dict(self, quantity, default=None):
        with h5py.File(self._file,'r') as fh:
            
            fhGroup = fh['clusters']['halo_1']
            if quantity not in fhGroup:
                return default
            else:
                info_dict = dict()
                for key in fhGroup[quantity].attrs:
                    info_dict[key] = fhGroup[quantity].attrs[key]
                return info_dict


    def _get_quantities_iter(self, quantities, filters, native_filters):
        for native_quantity_getter in self._iter_native_dataset(native_filters):
            data = self._load_quantities(quantities.union(set(filters.variable_names)),
                                         native_quantity_getter)
            try:
                data = filters.filter(data)
            except ValueError:
                filteredData = dict()
                for quantity in data.keys():
                    if(data[quantity].dtype != 'O'): filteredData[quantity] = filters.filter(data[quantity])
                    else:
                        haloLists = data[quantity]
                        filteredData[quantity] = np.array([filters.filter({quantity: haloVals})[quantity] for haloVals in haloLists])
                data = filteredData

            for q in set(data).difference(quantities):
                del data[q]
            yield data
            del data
   
