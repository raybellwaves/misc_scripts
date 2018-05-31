"""
Created on Apr 25 2017

seasonalskillmetrics.py

Seasonal forecasting skill metrics

See also: https://github.com/WFRT/verif
          https://github.com/PeterRochford/SkillMetrics

@author: RayBell
"""
import xarray as xr
import numpy as np
import scipy.stats as st
import math
import sys

def errorchecks_da_or_ds(a):
    """Error checks if DataArray or DataSet

       Args:

       * a:

       Returns:

           Raises an error is array is not compatible

       Example:
       
           errorchecks_da_or_ds(a)
    """

    # Check variable is a xarray.DataArray
    tmp = isinstance(a, xr.core.dataarray.DataArray)
    if tmp == True:
        errorchecks1var(a)
    else:
        tmp2 = isinstance(a, xr.core.dataset.Dataset)
        if tmp2 == True:
            errorchecksvars(a)
        else:
            raise TypeError('a is: ' + str(type(a)) + '. Should be xr.DataArray or xr.Dataset')

def errorchecks1var(a):
    """Error checks to ensure the subsequent functions work correctly

       Args:

       * a:
                          data - xarray.DataArray of size (time) or
                                                (ensemble, time)
                                                    (time, data)
                                     (time, latitude, longitude)
			   (ensemble, time, latitude, longitude)
                    (model, ensemble, time, latitude, longitude)

       Returns:

           Raises an error is variable is not compatible

       Example:
       
           errorchecks1var(a)
    """

    # Check variable is a xarray.DataArray
    tmp = isinstance( a, xr.core.dataarray.DataArray )
    if tmp == False:
       raise TypeError('a is: ' + str(type(a)) + '. Should be xr.DataArray')

    # Check the sizes of the dataset is expected
    tmp = len(a.dims)
    if tmp == 1:
        strcheck = 'time'
        if strcheck not in a.dims: 
            raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be time')
    if tmp == 2:
        strcheck = ['ensemble', 'time', 'data']
        good = 0
        if strcheck in a.dims:
            good = 1
        if good == 0:
            raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be ensemble, time or data')

    if tmp >= 3: 
        strcheck = ['time', 'latitude', 'longitude']
        for i in strcheck:
            if i not in a.dims: 
                raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be time, latitude and longitude')

def errorchecksvars(a):
    """Error checks to ensure the subsequent functions work correctly

       Args:

       * a:
                          data - xarray.Dataset of size (time) or
                                              (ensemble, time) or
                                                  (time, data) or
                                       (model, ensemble, time) or
                                   (time, latitude, longitude) or
			 (ensemble, time, latitude, longitude)
                  (model, ensemble, time, latitude, longitude)

       Returns:

           Raises an error is variable is not compatible

       Example:
       
           errorchecksvars(a)
    """
    # Check the sizes of the dataset is expected
    tmp = len(a.dims)
    if tmp == 1:
        strcheck = 'time'
        if strcheck not in a.dims: 
            raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be time')
    if tmp == 2:
        #strcheck = ['ensemble', 'time', 'data']
        #good = 0
        #if strcheck in a.dims:
        #    good = 1
        #if good == 0:
        #    raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be ensemble, time or data')
        strcheck = ['ensemble', 'time']
        for i in strcheck:
            if i not in a.dims:
                raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be ensemble and time')
    if tmp == 3: 
        strcheck = ['time', 'latitude', 'longitude']
        strcheck2 = ['model','ensemble','time']
        for i in strcheck:
            if i not in a.dims:
                for j in strcheck2:
                    if j not in a.dims:  
                        raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be time, latitude and longitude or model, ensemble and time')
    if tmp == 4:
        strcheck = ['model', 'time', 'latitude', 'longitude']
        strcheck2 = ['ensemble', 'time', 'latitude', 'longitude']
        for i in strcheck:
            if i not in a.dims: 
                for j in strcheck2:
                    if j not in a.dims:
                        raise NameError('The dimensions in a are: ' + str(a.dims) + '; should be model, time, latitude and longitude or ensemble, time, latitude and longitude')

def errorchecks(mod, obs):
    """Error checks to ensure the subsequent functions work correctly

       Args:

       * mod:
           Model data - xarray.DataArray/Dataset of size (time) or
                          (ensemble, time, latitude, longitude)
                   (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray/Dataset of size (time) or
                                            (time, latitude, longitude)

       Returns:

           Raises an error is variables are not compatible

       Example:
       
           errorchecks(mod,obs)
    """

    # Check both variables are xarray.DataArray
    a = isinstance( mod, xr.core.dataarray.DataArray )
    b = isinstance( obs, xr.core.dataarray.DataArray )
    if a == False or b == False:
        aa = isinstance( mod, xr.core.dataset.Dataset )
        bb = isinstance( mod, xr.core.dataset.Dataset )
        if aa == False or bb == False:
            raise TypeError('mod is: ' + str(type(mod)) + '; obs is: ' + str(type(obs)) + '. Both should be xr.DataArray/xr.Dataset')

    # Check the sizes of the datasets are compatible
    a = len(mod.dims)
    b = len(obs.dims)
    if a == 1 or b == 1:
        # Check that the names of the dimensions for 'time' is labelled correctly
        strcheck = 'time'
        if strcheck not in obs.dims: 
            raise NameError('The dimensions in obs are: ' + str(obs.dims) + '; should be time')
        if strcheck not in mod.dims: 
            raise NameError('The dimensions in mod are: ' + str(mod.dims) + '; should be time')

        # Check that the values for 'time' are the same
        b = mod.coords['time'].values
        c = obs.coords['time'].values
        if np.array_equal(b,c) is False:
           raise ValueError('The value of ' + a + ' in mod are: ' + str(b) + '; in obs are: ' + str(c))

    else:
        # Check that the names of the dimensions for 'time', 'latitude' and 'longitude' are labelled correctly
        strcheck = ['time', 'latitude', 'longitude']
        for a in strcheck:
            if a not in obs.dims: 
                raise NameError('The dimensions in obs are: ' + str(obs.dims) + '; should be time, latitude and longitude')
        for a in strcheck:
            if a not in mod.dims: 
                raise NameError('The dimensions in mod are: ' + str(mod.dims) + '; should be time, latitude and longitude')

        # Check that the values for 'time', 'latitude' and 'longitude' are the same
        for a in strcheck:
            b = mod.coords[a].values
            c = obs.coords[a].values
            if np.array_equal(b,c) is False:
                raise ValueError('The value of ' + a + ' in mod are: ' + str(b) + '; in obs are: ' + str(c))

def climatology(a):
    """Climatology

       Calculate the climatology. Average over different dimensions depending on how my dimensions the DataArray is.

       Args:

       * a:
           xarray.DataArray/Dataset of size (time, latitude, longitude) or
                                  (ensemble, time, latitude, longitude) or
                            (model, ensemble, time, latitude, longitude) 

       Returns:

       * b:
           Climatology
           xarray.DataArray/Dataset of size (latitude, longitude) or
                                     (model, latitude, longitude)

       Example:

       modclim = climatology(mod)

    """

    errorchecks_da_or_ds(a)
    # Number of dimensions
    ndims = len(a.dims)
    if ndims <= 3:
        b = a.mean(dim=('time'))
    else:
        b = a.mean(dim=('time', 'ensemble'))
    return b

def standarddeviation(a):
    """Standard deviation

       Calculate the standard deviation. Average over different dimensions depending on how my dimensions the DataArray is.

       Args:

       * a:
           xarray.DataArray/Dataset of size (time, latitude, longitude) or
                                  (ensemble, time, latitude, longitude) or
                            (model, ensemble, time, latitude, longitude) 

       Returns:

       * b:
           Climatology
           xarray.DataArray/Dataset of size (latitude, longitude) or
                                     (model, latitude, longitude)

       Example:

       modstd = standarddeviation(std)

    """

    errorchecks_da_or_ds(a)
    # Number of dimensions
    ndims = len(a.dims)
    if ndims == 3:
        b = a.std(dim=('time'))
    else:
        b = a.std(dim=('time', 'ensemble'))
    return b

def meanstatebias(modclim, obsclim):
    """Mean-state Bias.

       Calculate the mean-state bias

       Args:

       * modclim:
           xarray.DataArray/Dataset of size (latitude, longitude) or
                                     (model, latitude, longitude)

       * obsclim:
           xarray.DataArray/Dataset of size (latitude, longitude)

       Returns:

       * b:
           Mean-state Bias
           xarray.DataArray/Dataset of size (latitude, longitude)

       Example:

           modbias = meanstatebias(modclim, obsclim)

    """

    b = modclim - obsclim
    return b

def absolute_error(mod, obs):
    """Absolute error.

       Calualte the absolute error (modulus of mod minus obs)

       Args:

       * mod:
           xarray.DataArray of size (time)

       * obs:
           xarray.DataArray of size (time)

       Returns:

       * b:
           absolute error
           xarray.DataArray of size (time)

       Example:

           ae = absolute_error(mod, obs)

    """

    errorchecks(mod, obs)

    b = np.fabs(mod - obs)
    return b    

def observationanomalies(obs, obsclim):
    """Observational anomalies.

       Calculate observational anomalies

       Args:

       * obs:
           xarray.DataArray of size (time, latitude, longitude)

       * obsclim:
           xarray.DataArray of size (latitude, longitude)

       Returns:

       * b:
           Observational anomalies
           xarray.DataArray of size (time, latitude, longitude)

       Example:

           obsanom = observationanomalies(obs, obsclim)

    """

    b = obs - obsclim
    return b

def modelanomalies(mod, obsclim, modbias, **kwargs):
    """Model anomalies.

       Calculate model anomalies. 
       Option to remove mean-state bias

       Args:

       * mod:
           xarray.DataArray of size (time, latitude, longitude) or 
                          (ensemble, time, latitude, longitude) or
                   (model, ensemble, time, latitude, longitude)

       * obsclim:
           xarray.DataArray of size (latitude, longitude)

       * modbias:
           xarray.DataArray of size (latitude, longitude)

       * kwargs:
           removebias = 1 - remove the model bias

       Returns:

       * b:
           Model anomalies
           xarray.DataArray of size (time, latitude, longitude) or
                          (ensemble, time, latitude, longitude) or
                   (model, ensemble, time, latitude, longitude)

       Example:

           modanom = modelanomalies(mod, obsclim, modclim, removebias = 1)

    """

    if 'removebias' in kwargs:
        b = mod - obsclim - modbias
    return b

def standardizeobsanom(obsanom):
    """Standardize Observational anomalies.

       Divide observational anomalies by the standard devations

       Args:

       * obsanom:
           Observational anomalies
           xarray.DataArray of size (time, latitude, longitude)

       Returns:

       * b:
           Standardized Observational anomalies
           xarray.DataArray of size (time, latitude, longitude)

       Example:

           stdzobsanom = standardizeobsanom(obsanom)

    """

    b = obsanom / obsanom.std(dim=('time'))
    return b

def standardizemodanom(modanom):
    """Standardize Model anomalies.

       Divide mode anomalies by the standard devations

       Args:

       * modanom:
           Model anomalies
           xarray.DataArray of size (time, latitude, longitude) or
                          (ensemble, time, latitude, longitude)

       Returns:

       * b:
           Standardized model anomalies
           xarray.DataArray of size (time, latitude, longitude) or
                          (ensemble, time, latitude, longitude)

       Example:

           stdzmodanom = standardizemodanom(modanom)

    """

    #b = xr.full_like(modanom, np.nan)

    if len(modanom.dims) == 2:
        b = modanom / modanom.std(dim=('ensemble','time'))
    if len(modanom.dims) == 3:
        b = modanom / modanom.std(dim=('time'))
    if len(modanom.dims) == 4:
        b = modanom / modanom.std(dim=('ensemble','time'))

    return b

def combinemodels(a, nens):
    """Combine models.

       Combine all model ensembles. 

       Args:

       * a:
           xarray.DataArray/Dataset of size (model, ensemble, time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           Model anomalies
           xarray.DataArray of size (ensemble, time, latitude, longitude)

       Example:

           modanom4d = combinemodels(modanom)

    """

    len(a.coords['time'].values)
    if len(a.dims) == 5:    
        # model, ensemble, time, lat, lon)
        #b = a.stack(allensembles=['model','ensemble'])
        if isinstance( a, xr.core.dataarray.DataArray ) == True:
            # Build data array back up removing empty ensembles as models have different ensemble lengths
            datatmp = np.full([nens, len(a.coords['time'].values), len(a.coords['latitude'].values), len(a.coords['longitude'].values)], np.nan)
            b = xr.DataArray(datatmp, coords=[range(1,nens + 1), a.coords['time'], a.coords['latitude'], a.coords['longitude']], dims=['ensemble', 'time', 'latitude', 'longitude'])  
     
            counter = 0
            # Loop over model
            for i in range(0,len(a.coords['model'].values)):
                # Loop over ensemble
                for j in range(0,len(a.coords['ensemble'].values)):
                    # Check that data exists for that ensemble
                    if not np.all(a[i,j,:,:,:].isnull()):
                        b[counter,:,:,:] = a[i,j,:,:,:]
                        counter += 1
        else:  
            datavars = a.data_vars
            datatmp = np.full([len(datavars), nens, len(a.coords['time'].values), len(a.coords['latitude'].values), len(a.coords['longitude'].values)], np.nan)

            # Loop over data variables
            for i in range(0, len(datavars)):
                counter = 0
                # Loop over model
                for j in range(0,len(a.coords['model'].values)):
                    # Loop over ensembles
                    for k in range(0,len(a.coords['ensemble'].values)):
                        if not np.all(a[list(datavars.keys())[i]][j,k,:,:,:].isnull()):
                            datatmp[i,counter,:,:,:] = a[list(datavars.keys())[i]][j,k,:,:,:].values   
                            counter += 1
            bda = xr.DataArray(datatmp, coords=[datavars, range(1,nens + 1), a.coords['time'], a.coords['latitude'], a.coords['longitude']], dims=['datavar', 'ensemble', 'time', 'latitude', 'longitude'])
            b = bda.to_dataset(dim='datavar')
    if len(a.dims) == 3:
    # model, ensemble, time
        if isinstance( a, xr.core.dataarray.DataArray ) == True:
            pass
            print('code this')
        else:
            datavars = a.data_vars
            datatmp = np.full([len(datavars), nens, len(a.coords['time'].values)], np.nan)

            # Loop over data variables
            for i in range(0, len(datavars)):
                counter = 0

                for j in range(0,len(a.coords['model'].values)):
                    # Loop over ensembles
                    for k in range(0,len(a.coords['ensemble'].values)):
                        if not np.all(a[list(datavars.keys())[i]][j,k,:].isnull()):
                            datatmp[i,counter,:] = a[list(datavars.keys())[i]][j,k,:].values   
                            counter += 1
            bda = xr.DataArray(datatmp, coords=[datavars, range(1,nens + 1), a.coords['time']], dims=['datavar', 'ensemble', 'time'])
            b = bda.to_dataset(dim='datavar')
    return b

def combinemodelsNAO(a, nens):
    """Combine models.

       Combine all model ensembles. 

       Args:

       * a:
           xarray.DataArray/Dataset of size (latitude, longitude, model, ensemble, time)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           Model anomalies
           xarray.DataArray/Dataset of size (latitude, longitude, ensemble, time)

       Example:

           stdanom = ssm.combinemodelsNAO(stdanom, nens))

    """
 
    datatmp = np.full([len(a.coords['latitude'].values), len(a.coords['longitude'].values), nens, len(a.coords['time'].values)], np.nan)
    b = xr.DataArray(datatmp, coords=[a.coords['latitude'], a.coords['longitude'], range(1,nens + 1), a.coords['time']], dims=['latitude', 'longitude','ensemble', 'time'])
     
    counter = 0
    # Loop over model
    for i in range(0,len(a.coords['model'].values)):
        # Loop over ensemble
        for j in range(0,len(a.coords['ensemble'].values)):
            # Check that data exists for that ensemble
            if not np.all(a[:,:,i,j,:].isnull()):
            #if not np.isnan(np.sum(a[i,j,:,:,:].values)):
               b[:,:,counter,:] = a[:,:,i,j,:]
               counter += 1      
    return b


def rootmeansquareerror(mod, obs, nens):
    """Root Mean Square Error.

       Calculates the Root Mean Square Error of model dataset compared to an observed dataset.
       See line 242 in Bell et al (2017).
       Values are 0 -> +.

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           Root Mean Square Error 
           xarray.DataArray of size (latitude, longitude)

       Example:

       Calculate the Root Mean Square Error of a model dataset compared to an observed dataset

       rmse = rootmeansquareerror(mod, obs, nens)
    """
    
    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
   
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)

    tmp = (modanom - obsanom)**2
    if isinstance( mod, xr.core.dataarray.DataArray ) == True: 
        b = np.sqrt( (1. / (nens * len(modanom.coords['time'].values))) * tmp.sum(dim=('ensemble','time')) )
        b = b.rename('rmse') # add name to b
    else:
        tmp2 = (1. / (nens * len(modanom.coords['time'].values))) * tmp.sum(dim=('ensemble','time'))
        b = tmp2.apply(np.sqrt)
    return b 

def saturationrmse(mod, obs, nens):
    """Saturation Root Mean Square Error.

       Root Mean Square Error normalized by the sun of the variances of the model and the observed.
       It shows how important error growth is relative to the system's own natural variability.
       A forecast is said to be saturated with eror if Saturation RMSE > 1.
       This occurs when the anomaly correlation is <= 0. 
       Calculates the Root Mean Square Error of model dataset compared to an observed dataset.
       See line 250 in Bell et al (2017).
       Values are 0.5 -> 1.5.

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           Saturation Root Mean Square Error 
           xarray.DataArray of size (latitude, longitude)

       Example:

       Calculate the Saturation Root Mean Square Error of a model dataset compared to an observed dataset

       satrmse = saturationrmse(mod, obs, nens)
    """

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
   
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)

    rmse = rootmeansquareerror(mod, obs, nens)

    tmp1 = (modanom.std(dim=('ensemble','time')))**2
    tmp2 = (obsanom.std(dim=('time')))**2
    b = rmse / np.sqrt(tmp1 + tmp2)
    b = b.rename('saturationrmse') # add name to b
    return b 

def anomalycorrelation(mod, obs, nens, siglevel=None):
    """Anomaly Correlation.

       See line 236 in Bell et al (2017).
       Values are -1 -> 1.

       Args:

       * mod:
           Model data - xarray.DataArray/Dataset of size (time, latitude, longitude) or 
                                               (ensemble, time, latitude, longitude) or
                                        (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       * siglevel:
           Significant level for which to copute t-test to. 0-1. float.

       Returns:

       * b:
           Anomaly Correlation
           xarray.DataArray/Dataset of size (latitude, longitude)

       * b1:
           Significance map
           boolean xarray.DataArray/Dataset of size (latitude, longitude) - 1 or nan

       Example:

       Calculate the Anomaly correlation of a model dataset compared to an observed dataset

       anomcorr = anomalycorrelation(mod, obs, nens)
    """

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
    
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)   
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)
    modanommean = modanom.mean(dim=('ensemble'))

    # see https://hrishichandanpurkar.blogspot.com/2017/09/vectorized-functions-for-correlation.html
    x,y = xr.align(obsanom,modanommean)
    n = len(x.coords['time'])
    xmean = x.mean(dim=('time'))
    ymean = y.mean(dim=('time'))
    xstd = x.std(dim=('time'))
    ystd = y.std(dim=('time'))
    cov = np.sum((x - xmean)*(y - ymean), axis=0)/(n)
    b = cov/(xstd*ystd)

    if isinstance( mod, xr.core.dataarray.DataArray ) == True:
        b = b.rename('anomalycorrelation') # add name to b

    if siglevel is not None:
        dof = len(modanom.coords['time'].values) - 2
        #dof = (len(modanom.coords['time'].values) * nens) - 2
        if isinstance( mod, xr.core.dataarray.DataArray ) == True:
            tmp = np.full([len(b.coords['latitude'].values), len(b.coords['longitude'].values)], 1) 
            b1 = xr.DataArray(tmp, coords=[b.coords['latitude'].values, b.coords['longitude'].values], dims=['latitude', 'longitude'])
            # Make the distribution normal by putting through a fisher r to z transform
            # https://en.wikipedia.org/wiki/Fisher_transformation
            zarr = np.arctanh(b.values)
            # If zarr is greater or equal to 1.0 set t 0.99 at it messes up t_arr calc
            zarr[zarr >= 1] = 0.99
            # Calculate the t statisic 
            # http://vassarstats.net/rsig.html
            tarr = zarr / (np.sqrt((1.0 - np.square(zarr))/(dof)))
            # Calculate associated probability
            pvalues = st.norm.cdf(tarr)
            b1 = b1.where(pvalues >= siglevel)
        else:
            zarr = b.apply(np.arctanh)
            zarr = zarr.where(zarr <= 1).fillna(0.99) # Where will keep values
            tmp = zarr.apply(np.square)
            tmp2 = (1.0 - tmp)/dof
            tmp3 = tmp2.apply(np.sqrt)
            tarr = zarr / tmp3
            pvalues = tarr.apply(st.norm.cdf)
            b1 = pvalues.where(pvalues <= siglevel).fillna(1)
            b1 = b1.where(pvalues > siglevel)
        # Add sig level as an attriube
        b1.attrs['sig_level'] = siglevel
        b1.attrs['dof'] = dof
        return b, b1
    else:
        return b 

def spread(mod, obs, nens):
    """Spread.

       How widely the ensembles vary.
       See line 244 in Bell et al (2017).

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           Spread 
           xarray.DataArray of size (latitude, longitude)

       Example:

       Calculate the Spread of a model

       spread = spread(mod, obs, nens)
    """

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
    
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)    
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)

    p = math.factorial(nens) / (2 * math.factorial(nens - 2))

    # Setup index array for spread calculation
    refarr = np.zeros(int(p))
    start = 0
    for i in range(0,nens-1):
        end = start + nens - i - 2
        refarr[start:end+1] = i
        start = end + 1
    refarr2 = np.zeros(int(p))
    start = 0
    for i in range(0,nens-1):
        end = start + nens - i - 2
        refarr2[start:end+1] = np.arange(i,nens-1) + 1
        start = end + 1

    tmp = np.zeros((int(p), len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
    # Loop over ensemble combinations
    for i in range(0,int(p)):
        ref1 = int(refarr[i])
        ref2 = int(refarr2[i])
        tmp[i,:,:,:] = (modanom[ref1,:,:,:].values - modanom[ref2,:,:,:].values)**2 
    tmp2 = np.sqrt( (1./(p*len(mod.coords['time'].values))) * np.sum(tmp, axis=(0,1)) )
    b = xr.DataArray(tmp2, coords=[mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['latitude', 'longitude']) 
    return b   


def rpss(mod, obs, nens):
    """Rank Probability Skill score.

       Probabilistic skill relative to a climatological forecast.
       A value greater than 0 is better than climatology (Weigel et al. 2007)
       See line 258 in Bell et al (2017).

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * b:
           rpss 
           xarray.DataArray of size (latitude, longitude)

       Example:

       Calculate the RPSS of a forecast compared to observed

       rpss = rpss(mod, obs, nens)
    """

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
    
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)

    # Standardize anomalies
    stdzobsanom = standardizeobsanom(obsanom)
    stdzmodanom = standardizemodanom(modanom)
    
    event = 0.43 # Needed as normalized by standard deviation so mean = 1 and std = 0.
    # For the normal distribution terciles can be defined y the stand z score like in t-statistics.
    # So the upper tercile is 0.43 * std and lower -0.43 * std. 
    # 0.43 is the PDF into 3
    lowerevent = event * -1

    climo = np.array((1./3,2./3,1.0))

    if isinstance( mod, xr.core.dataarray.DataArray ) == True:
        # Create new DataArray for cumulative probability
        datatmp = np.full([3, len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)], 0, dtype='float')
        pcumobs = xr.DataArray(datatmp, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])
        pcumfct = xr.zeros_like(pcumobs)

        a = pcumobs.values
        # Observations 

        b = stdzobsanom.values
        # Lower tercile
        a[0,:,:,:][b <= lowerevent] = 1.0
        a[1,:,:,:][b <= lowerevent] = 1.0
        a[2,:,:,:][b <= lowerevent] = 1.0
        # Middle tercile
        a[0,:,:,:][(b > lowerevent) & (b < event) ] = 0
        a[1,:,:,:][(b > lowerevent) & (b < event) ] = 1.0
        a[2,:,:,:][(b > lowerevent) & (b < event) ] = 1.0
        # Upper tercile
        a[0,:,:,:][b >= event] = 0
        a[1,:,:,:][b >= event] = 0
        a[2,:,:,:][b >= event] = 1.0

        pcumobs = xr.DataArray(a, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])
    
        # Model
        a = pcumfct.values
        b = stdzmodanom.values
        c = np.zeros((3, len(stdzmodanom.coords['ensemble'].values), len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
        # Count how many ensembles in each tercile  
        for i in range(0,nens):
            # Lower tercile
            c[0,i,:,:,:][b[i,:,:,:] <= lowerevent] = (1.0/nens)
            # Middle tercile
            c[1,i,:,:,:][(b[i,:,:,:] > lowerevent) & (b[i,:,:,:] < event)] = (1.0/nens)
            # Upper tercile
            c[2,i,:,:,:][b[i,:,:,:] >= event] = (1.0/nens)

        # Sum along ensemble dimension in c to create a
        a[0,:,:,:] = c[0,:,:,:,:].sum(axis=0)
        a[1,:,:,:] = c[1,:,:,:,:].sum(axis=0) + a[0,:,:,:]
        a[2,:,:,:] = c[2,:,:,:,:].sum(axis=0) + a[1,:,:,:]
        pcumfct = xr.DataArray(a, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])

        # Setup RPS(C) array
        rpsc = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
        rps = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))

        # Loop over time
        for i in range(0, len(mod.coords['time'].values)):
            # Loop over terciles
            for j in range(0, 3):
                rpsc = rpsc + (climo[j] - pcumobs[j,i,:,:].values)**2
                rps = rps + (pcumfct[j,i,:,:].values - pcumobs[j,i,:,:].values)**2

        rpss = 1.0 - (rps/rpsc)
        b = xr.DataArray(rpss, coords=[mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['latitude', 'longitude'])
    else:
        # Create new Dataset for cumulative probability
        datavars = obs.data_vars

        datatmp = np.full([len(datavars), 3, len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)], 0, dtype='float')
        # Loop over variables
        for i in range(0, len(datavars)):
            # Lower tercile
            datatmp[i,0,:,:,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 1.0
            datatmp[i,1,:,:,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 1.0
            datatmp[i,2,:,:,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 1.0
            # Middle tercile
            datatmp[i,0,:,:,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 0
            datatmp[i,1,:,:,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 1.0
            datatmp[i,2,:,:,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 1.0
            # Upper tercile
            datatmp[i,0,:,:,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 0
            datatmp[i,1,:,:,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 0
            datatmp[i,2,:,:,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 1.0

        pcumobsda = xr.DataArray(datatmp, coords=[datavars, np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['datavar', 'tercile', 'time', 'latitude', 'longitude'])
        pcumobs = pcumobsda.to_dataset(dim='datavar')

        pcumfct = xr.zeros_like(pcumobs)

        datatmp2 = np.full([len(datavars), 3, len(stdzmodanom.coords['ensemble'].values), len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)], 0, dtype='float')
        # Loop over variables
        for i in range(0, len(datavars)):
            # Loop over ensembles
            for j in range(0,nens):
                tmp = stdzmodanom[list(datavars.keys())[i]].sel(ensemble=j+1).expand_dims('ensemble')
                # Lower tercile
                datatmp2[i,0,j,:,:,:][tmp <= lowerevent] = (1.0/nens)
                # Middle tercile
                datatmp2[i,1,j,:,:,:][(tmp > lowerevent) & (tmp < event)] = (1.0/nens)
                # Upper tercile
                datatmp2[i,2,j,:,:,:][tmp > event] = (1.0/nens)
            # Sum along ensemble dimension
            pcumfct[list(datavars.keys())[i]][dict(tercile=0)] = datatmp2[i,0,:,:,:,:].sum(axis=0)
            pcumfct[list(datavars.keys())[i]][dict(tercile=1)] = datatmp2[i,1,:,:,:,:].sum(axis=0) + pcumfct[list(datavars.keys())[i]][dict(tercile=0)]
            pcumfct[list(datavars.keys())[i]][dict(tercile=2)] = datatmp2[i,2,:,:,:,:].sum(axis=0) + pcumfct[list(datavars.keys())[i]][dict(tercile=1)]


        # Setup RPS(C) array
        datatmp3 = np.zeros((len(datavars), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
        rpscda = xr.DataArray(datatmp3, coords=[datavars, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['datavar', 'latitude', 'longitude'])
        rpsc = rpscda.to_dataset(dim='datavar')
        rps = xr.zeros_like(rpsc)
        # Loop over time
        for i, time in enumerate(mod.coords['time'].values):
            # Loop over terciles
            for j in range(0, 3):
                rpsc = rpsc + (climo[j] - pcumobs.sel(tercile=j, time=time))**2
                #sys.exit()
                #tmp1 = pcumfct.sel(tercile=j, time=tile)
                #tmp2 = pcumobs.sel(tercile=j, time=time)
                #tmp = tmp1 - tmp2
                rps = rps + (pcumfct.sel(tercile=j, time=time) - pcumobs.sel(tercile=j, time=time))**2
        b = 1.0 - (rps/rpsc)
    return b


def mc_rpss_time(mod, obs, nens, mask=None):
    """Monte Carlo Rank Probability Skill Score for time.

       Probabilistic skill relative to a climatological forecast.
       A value greater than 0 is better than climatology (Weigel et al. 2007)
       See line 258 in Bell et al (2017).

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       * mask:
           Mask - xarray.DataArray of size (latitude, longitude)

       Returns:

       * da:
           mc_rpss 
           xarray.DataArray of size (time, nmodels, maxcombinations)

       Example:

       Calculate the RPSS of a forecast compared to observed

       rpss = rpss_time(mod, obs, nens)
    """
    # Only do for 2 model (CFSV2-2011 and GloSea5) and multi-model mean
    # Setup dataarray
    nmodels = len(mod.coords['model'])
    nmodelscoords = np.linspace(1,nmodels,nmodels, dtype=np.int)
    maxcomb = int(math.factorial(nmodels) / (2 * math.factorial(nmodels - 2)))
    maxcombcoords = np.linspace(1,maxcomb,maxcomb, dtype=np.int)
    ntime = len(mod.coords['time'])
    datatmp = np.full([ntime, nmodels, maxcomb], np.nan, dtype=np.float)
    sys.exit()
    #da = xr.DataArray(datatmp, coords=[nmodelscoords, maxcombcoords], dims=['nmodels', 'ncombinations'])


def mc_rpss(mod, obs, nens, mask=None):
    """Monte Carlo Rank Probability Skill Score.

       Calculate all possible model conbinations Rank Probability Skill Score.

       For example if mod is (model: 7, ensemble: 10, time: 17, latitude: 70, longitiude: 130)
       Setup a DataArray of n model combinations (7) x max conbinations 7!/2(7-2)! (21).
       Requires a mininum of 3 models.
       Loop through n model combinations:
       e.g. for 4 models
       1 model combinations = 4!/1!(4-1)! = 4
       2 model combinations = 4!/2!(4-2)! = 6
       3 model combinatinos = 4!/3!(4-3)! = 4
       4 model combinations = 4!/4!(4-4)! = 1
       e.g. for 7 models
       1 model combinations = 7!/1!(7-1)! = 7
       2 model combinations = 7!/2!(7-2)! = 21
       3 model combinatinos = 7!/3!(7-3)! = 35
       4 model combinations = 7!/4!(7-4)! = 35
       5 model combinations = 7!/5!(7-5)! = 21
       6 model combinations = 7!/6!(7-6)! = 7
       7 model combinatinos = 7!/7!(7-7)! = 1

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Number of ensembles (integer) - has to be same for all models

       * mask:
           Mask - xarray.DataArray of size (latitude, longitude)

       Returns:

       * da:
           mc_rpss 
           xarray.DataArray of size (nmodels, maxcombinations)

       Example:

       mc_rpss = ssm.mc_rpss(sub_mod, sub_obs, 10, mask=sub_reg_mask)

    """

    # Setup dataarray
    nmodels = len(mod.coords['model'])
    nmodelscoords = np.linspace(1,nmodels,nmodels, dtype=np.int)
    # Create combination arrary
    modcombarr = np.linspace(1,nmodels,nmodels, dtype=np.int)
    for i in range(1, nmodels+1):
        modcombarr[i-1] = int(math.factorial(nmodels) / (math.factorial(i) * math.factorial(nmodels - i)))
    maxcomb = np.max(modcombarr)
    maxcombcoords = np.linspace(1,maxcomb,maxcomb, dtype=np.int)
    datatmp = np.full([nmodels, maxcomb], np.nan, dtype=np.float)
    da = xr.DataArray(datatmp, coords=[nmodelscoords, maxcombcoords], dims=['nmodels', 'ncombinations'])
    
    # Loop over nmodelcombinations
    for i in range(1, nmodels+1):
        nmodelcombs = modcombarr[i-1]
        modelcomblist = [None] * nmodelcombs
        # 1 model combinations
        if i == 1:
            for j in range(0, nmodelcombs):
                _mod = mod[j,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
        # 2 model combinations
        elif i == 2:
            # Hard code this later...
            ix0 = 0
            ix1 = 1
            for j in range(0, nmodelcombs):
                # Need to loop through all possible combinations
                ix = np.array([ix0, ix1])
                _mod = mod[ix,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens*i)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
                # Increase the indices
                ix1 = ix1 + 1
                if ix1 == nmodels:
                    ix0 = ix0 + 1
                    ix1 = ix0 + 1
        # 3 model combinations
        elif i == 3:
            ix0 = 0
            ix1 = 1
            ix2 = 2
            for j in range(0, nmodelcombs):
                # Need to loop through all possible combinations
                ix = np.array([ix0, ix1, ix2])
                _mod = mod[ix,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens*i)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
                # Increase the indices
                ix2 = ix2 + 1
                if ix2 == nmodels:
                    ix1 = ix1 + 1
                    ix2 = ix1 + 1
                if ix1 == nmodels - 1:
                    ix0 = ix0 + 1
                    ix1 = ix0 + 1
                    ix2 = ix1 + 1
        # 4 model combinations
        elif i == 4:
            ix0 = 0
            ix1 = 1
            ix2 = 2
            ix3 = 3
            for j in range(0, nmodelcombs):
                # Need to loop through all possible combinations
                ix = np.array([ix0, ix1, ix2, ix3])
                _mod = mod[ix,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens*i)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
                # Increase the indices
                ix3 = ix3 + 1
                if ix3 == nmodels:
                    ix2 = ix2 + 1
                    ix3 = ix2 + 1
                if ix2 == nmodels - 1:
                    ix1 = ix1 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
                if ix1 == nmodels - 2:
                    ix0 = ix0 + 1
                    ix1 = ix0 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
        # 5 model combinations
        elif i == 5:
            ix0 = 0
            ix1 = 1
            ix2 = 2
            ix3 = 3
            ix4 = 4
            for j in range(0, nmodelcombs):
                # Need to loop through all possible combinations
                ix = np.array([ix0, ix1, ix2, ix3, ix4])
                _mod = mod[ix,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens*i)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
                # Increase the indices
                ix4 = ix4 + 1
                if ix4 == nmodels:
                    ix3 = ix3 + 1
                    ix4 = ix3 + 1
                if ix3 == nmodels - 1:
                    ix2 = ix2 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
                if ix2 == nmodels - 2:
                    ix1 = ix1 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
                if ix1 == nmodels - 3:
                    ix0 = ix0 + 1
                    ix1 = ix0 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
        # 6 model combination
        elif i == 6:
            ix0 = 0; ix1 = 1; ix2 = 2; ix3 = 3; ix4 = 4; ix5 = 5
            for j in range(0, nmodelcombs):
                # Need to loop through all possible combinations
                ix = np.array([ix0, ix1, ix2, ix3, ix4, ix5])            
                _mod = mod[ix,:,:,:,:]
                modelcomblist[j] = str(_mod.coords['model'].values)
                # Put through rpss
                _rpss = rpss(_mod, obs, nens*i)
                # Get spatial mean
                rpsswmask = _rpss.where(mask > 0)
                _rpssmean = rpsswmask.mean()
                da[i-1,j] = _rpssmean
                # Increase the indices
                ix5 = ix5 + 1
                if ix5 == nmodels:
                    ix4 = ix4 + 1
                    ix5 = ix4 + 1
                if ix4 == nmodels - 1:
                    ix3 = ix3 + 1
                    ix4 = ix3 + 1
                    ix5 = ix4 + 1
                if ix3 == nmodels - 2:
                    ix2 = ix2 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
                    ix5 = ix4 + 1
                if ix2 == nmodels - 3:
                    ix1 = ix1 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
                    ix5 = ix4 + 1
                if ix1 == nmodels - 4:
                    ix0 = ix0 + 1
                    ix1 = ix0 + 1
                    ix2 = ix1 + 1
                    ix3 = ix2 + 1
                    ix4 = ix3 + 1
                    ix5 = ix4 + 1
                    # 7 model combination          
        else:
            _mod = mod
            modelcomblist = str(_mod.coords['model'].values)
            # Put through rpss
            _rpss = rpss(_mod, obs, nens*nmodels)
            # Get spatial mean
            rpsswmask = _rpss.where(mask > 0)
            _rpssmean = rpsswmask.mean()
            da[i-1,0] = _rpssmean
        da.attrs[str(i)+'modelcombs'] = modelcomblist
    print(da)
    return da


def rocscore(mod, obs, nens, tercile):
    """Relative Operating Characteristics Score.

       Work out hit rate % and false alarm % for the bins starting from highest to lowest.
       two-by-two contingency tables (events (E), non-event (E'), warning (W) and no warning
       (W'). 
       See table 1 in Mason & Graham (1999).
       The following outcomes are possible: hit (h) of the event occured and a warning was provided, 
       a false alarm (f) if an event did not occur but a warning was provided; a miss (m) if an 
       event occurred but a warning was not provided; a correct rejection, if an event did not
       occur and warning was not provided (c).
       This set of hit rates is plotted against false-alarm rates to generate the ROC curve.
       The ROC score is the area under the curve.
       For skill, the ROC curve will lay above the 45o line and the area will be > 0.5.       
       See line 275 in Bell et al (2017).

       Args:

       * mod:
           Model data - xarray.DataArray of size (time, latitude, longitude) or 
                                       (ensemble, time, latitude, longitude) or
                                (model, ensemble, time, latitude, longitude)

       * obs:
           Observational data - xarray.DataArray of size (time, latitude, longitude)

       * nens:
           Combined number of ensembles (integer)

       * tercile:
           Tercile of interest (0,1 or 2; lower, middle or upper) (integer)

       Returns:

       * b:
           rocscore 
           xarray.DataArray of size (latitude, longitude)

       Example:

       Calculate the ROC score of a forecast compared to observed

       rocscore = rocscore(mod, obs, nens)
    """     

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
    
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 5:
        modanom = combinemodels(modanom, nens)

    stdzobsanom = standardizeobsanom(obsanom)
    stdzmodanom = standardizemodanom(modanom)

    event = 0.43 # Needed as normalized by standard deviation so mean = 1 and std = 0.
    # For the normal distribution terciles can be defined y the stand z score like in t-statistics.
    # So the upper tercile is 0.43 * std and lower -0.43 * std. 
    # 0.43 is the PDF into 3
    lowerevent = event * -1

    # Create new DataArray for cumulative probability
    datatmp = np.full([3, len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)], 0, dtype='float')
    pobs = xr.DataArray(datatmp, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])
    pfct = xr.zeros_like(pobs)

    # Observations 
    a = pobs.values
    b = stdzobsanom.values
    # Lower tercile
    a[0,:,:,:][b <= lowerevent] = 1.0
    a[1,:,:,:][b <= lowerevent] = 0.0
    a[2,:,:,:][b <= lowerevent] = 0.0
    # Middle tercile
    a[0,:,:,:][(b > lowerevent) & (b < event) ] = 0.0
    a[1,:,:,:][(b > lowerevent) & (b < event) ] = 1.0
    a[2,:,:,:][(b > lowerevent) & (b < event) ] = 0.0
    # Upper tercile
    a[0,:,:,:][b >= event] = 0.0
    a[1,:,:,:][b >= event] = 0.0
    a[2,:,:,:][b >= event] = 1.0

    pobs = xr.DataArray(a, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])

    # Model
    a = pfct.values
    b = stdzmodanom.values
    c = np.zeros((3, len(stdzmodanom.coords['ensemble'].values), len(mod.coords['time'].values), len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
    # Count how many ensembles in each tercile  
    for i in range(0,nens):
        # Lower tercile
        c[0,i,:,:,:][b[i,:,:,:] <= lowerevent] = (1.0/nens)
        # Middle tercile
        c[1,i,:,:,:][(b[i,:,:,:] > lowerevent) & (b[i,:,:,:] < event)] = (1.0/nens)
        # Upper tercile
        c[2,i,:,:,:][b[i,:,:,:] >= event] = (1.0/nens)

    # Sum along ensemble dimension in c to create a
    a[0,:,:,:] = c[0,:,:,:,:].sum(axis=0)
    a[1,:,:,:] = c[1,:,:,:,:].sum(axis=0)
    a[2,:,:,:] = c[2,:,:,:,:].sum(axis=0)
    pfct = xr.DataArray(a, coords=[np.arange(0,3), mod.coords['time'].values, mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['tercile', 'time', 'latitude', 'longitude'])

    # Setup bins.
    # Want mid-point somewhere around 0.33
    # There are a couple of methods for ROC score: can either choose 10 bins or use the number of ensembles as bins
    # Choosing 10 bin method here
    nrocbins = 10
    rocbins = np.linspace(0.15, 0.45, num = nrocbins - 1) # Percentage of ensembles

    a = pobs.values
    b = pfct.values

    # Setup arrays
    hitrate = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values), len(rocbins)+2))
    falsealarm = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values), len(rocbins)+2))
    # Last value is always 1
    hitrate[:,:,-1] = 1.0
    falsealarm[:,:,-1] = 1.0
    auc = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))
    normauc = np.zeros((len(mod.coords['latitude'].values), len(mod.coords['longitude'].values)))

    #---- loop testing ----#
    for i in range(0, len(mod.coords['latitude'].values)):
        for j in range(0, len(mod.coords['longitude'].values)):  
           # How many events were observed in the tercile?  
           nevents = np.sum(a[tercile,:,i,j])
           nnonevents = len(mod.coords['time'].values) - nevents
           # Extract the forecasts for these events
           pfct_tercile_ref = np.where(a[tercile,:,i,j] == 1)
           pfct_tercile = b[tercile,:,i,j][pfct_tercile_ref]
           # Extract data for the other years
           pfct_nontercile_ref = np.where(a[tercile,:,i,j] != 1)
           pfct_nontercile = b[tercile,:,i,j][pfct_nontercile_ref]

           for k in range(1, len(rocbins)+1):
               if k == 1:    
                   # How many years have forecast the events with > max bin (%) of the ensembles?
                   hitrate[i,j,k] = sum(l > rocbins[-k] for l in pfct_tercile) / nevents
                   falsealarm[i,j,k] = sum(l > rocbins[-k] for l in pfct_nontercile) / nnonevents
               elif k == len(rocbins)+1:
                   # How many years have forecast the events with < min bin % of the ensembles?
                   hitrate[i,j,k] = (sum(l < rocbins[0] for l in pfct_tercile) / nevents) + hitrate[i,j,k-1]
                   falsealarm[i,j,k] = (sum(l < rocbins[0] for l in pfct_nontercile) / nnonevents) + falsealarm[i,j,k-1]
               else: 
                   hitrate[i,j,k] = (sum(l < rocbins[(k-1)*-1] and l > rocbins[-k] for l in pfct_tercile) / nevents) + hitrate[i,j,k-1]
                   falsealarm[i,j,k] = (sum(l < rocbins[(k-1)*-1] and l > rocbins[-k] for l in pfct_nontercile) / nnonevents) + falsealarm[i,j,k-1]

           # Area under the curve
           auc[i,j] = np.trapz(hitrate[i,j,:], x=falsealarm[i,j,:]) 
           normauc[i,j] = 2 * (auc[i,j] - 0.5) 

    b = xr.DataArray(normauc, coords=[mod.coords['latitude'].values, mod.coords['longitude'].values], dims=['latitude', 'longitude']) 
    return b

def rocscore_1d(mod, obs, nens):
    """Relative Operating Characteristics Score.

       Work out hit rate % and false alarm % for the bins starting from highest to lowest.
       two-by-two contingency tables (events (E), non-event (E'), warning (W) and no warning
       (W'). 
       See table 1 in Mason & Graham (1999).
       The following outcomes are possible: hit (h) of the event occured and a warning was provided, 
       a false alarm (f) if an event did not occur but a warning was provided; a miss (m) if an 
       event occurred but a warning was not provided; a correct rejection, if an event did not
       occur and warning was not provided (c).
       This set of hit rates is plotted against false-alarm rates to generate the ROC curve.
       The ROC score is the area under the curve.
       For skill, the ROC curve will lay above the 45o line and the area will be > 0.5.       
       See line 275 in Bell et al (2017).

       Args:

       * mod:
           Model data - xarray.Dataset of size (time) or 
                                     (ensemble, time) or
                              (model, ensemble, time)

       * obs:
           Observational data - xarray.Dataset of size (time)

       * nens:
           Combined number of ensembles (integer)

       Returns:

       * rocscore:
           hit rate
           flase alarm 
           xarray.Dataset of size (tercile, 10+1)

       Example:

       Calculate the area under the curve, hit rate and false alarm rate of a forecast compared to observed
       rocscore = rocscore_1d(mod, obs, nens)
    """

    errorchecks(mod, obs)

    modclim = climatology(mod)
    obsclim = climatology(obs)

    modbias = meanstatebias(modclim, obsclim)
    
    obsanom = observationanomalies(obs, obsclim)

    modanom = modelanomalies(mod, obsclim, modbias, removebias = 1)
    # If mod is 5 dimensions combine all ensembles to make 4 dimensions
    if len(mod.dims) == 3:
        modanom = combinemodels(modanom, nens)
    stdzobsanom = standardizeobsanom(obsanom)
    stdzmodanom = standardizemodanom(modanom)

    event = 0.43 # Needed as normalized by standard deviation so mean = 1 and std = 0.
    # For the normal distribution terciles can be defined y the stand z score like in t-statistics.
    # So the upper tercile is 0.43 * std and lower -0.43 * std. 
    # 0.43 is the PDF into 3
    lowerevent = event * -1

    # Create new Dataset for cumulative probability
    datavars = obs.data_vars

    datatmp = np.full([len(datavars), 3, len(mod.coords['time'].values)], 0, dtype='float')
    # Loop over variables
    for i in range(0, len(datavars)):
        # Lower tercile
        datatmp[i,0,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 1.0
        datatmp[i,1,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 0.0
        datatmp[i,2,:][stdzobsanom[list(datavars.keys())[i]].values <= lowerevent] = 0.0
        # Middle tercile
        datatmp[i,0,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 0.0
        datatmp[i,1,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 1.0
        datatmp[i,2,:][(stdzobsanom[list(datavars.keys())[i]].values > lowerevent) & (stdzobsanom[list(datavars.keys())[i]].values < event) ] = 0.0
        # Upper tercile
        datatmp[i,0,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 0
        datatmp[i,1,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 0
        datatmp[i,2,:][stdzobsanom[list(datavars.keys())[i]].values >= event] = 1.0

    pobsda = xr.DataArray(datatmp, coords=[datavars, np.arange(0,3), mod.coords['time'].values], dims=['datavar', 'tercile', 'time'])
    pobs = pobsda.to_dataset(dim='datavar')

    pfct = xr.zeros_like(pobs)

    datatmp2 = np.full([len(datavars), 3, len(stdzmodanom.coords['ensemble'].values), len(mod.coords['time'].values)], 0, dtype='float')
    # Loop over variables
    for i in range(0, len(datavars)):
        # Loop over ensembles
        for j in range(0,nens):
            tmp = stdzmodanom[list(datavars.keys())[i]].sel(ensemble=j+1).expand_dims('ensemble')
            # Lower tercile
            datatmp2[i,0,j,:][tmp <= lowerevent] = (1.0/nens)
            # Middle tercile
            datatmp2[i,1,j,:][(tmp > lowerevent) & (tmp < event)] = (1.0/nens)
            # Upper tercile
            datatmp2[i,2,j,:][tmp > event] = (1.0/nens)
        # Sum along ensemble dimension
        pfct[list(datavars.keys())[i]][dict(tercile=0)] = datatmp2[i,0,:,:].sum(axis=0)
        pfct[list(datavars.keys())[i]][dict(tercile=1)] = datatmp2[i,1,:,:].sum(axis=0)
        pfct[list(datavars.keys())[i]][dict(tercile=2)] = datatmp2[i,2,:,:].sum(axis=0)

    # Setup bins.
    # Want mid-point somewhere around 0.33
    # There are a couple of methods for ROC score: can either choose 10 bins or use the number of ensembles as bins
    # Choosing 10 bin method here
    nrocbins = 10
    rocbins = np.linspace(0.15, 0.45, num = nrocbins - 1) # Percentage of ensembles

    # Setup arrays
    datatmp3 = np.zeros((len(datavars), 3, len(rocbins)+2))
    hitrateda = xr.DataArray(datatmp3, coords=[datavars, np.arange(0,3), np.arange(0, len(rocbins)+2)], dims=['datavar', 'tercile', 'rocbin'])
    # Last value is always 1
    hitrateda[:, :, -1] = 1.0
    falsealarmda = hitrateda.copy(deep=True)

    # Loop over variables
    for i in range(0, len(datavars)):
        # Loop over terciles
        for j in range(0, len(hitrateda.coords['tercile'].values)):
            # How many events were in observed in the tercile?
            nevents = pobs[list(datavars.keys())[i]][dict(tercile=j)].sum(dim='time')
            nnonevents = len(mod.coords['time'].values) - nevents
            # Extract the forecasts for these events (time indexs)
            ind = np.where(pobs[list(datavars.keys())[i]][dict(tercile=j)].values == 1)
            pfct_tercile = pfct[list(datavars.keys())[i]][dict(tercile=j)][ind]
            # Extract data for the other years
            ind = np.where(pobs[list(datavars.keys())[i]][dict(tercile=j)].values != 1)
            pfct_nontercile = pfct[list(datavars.keys())[i]][dict(tercile=j)][ind]
            # Loop over rocbins
            for k in range(1, len(rocbins)+1):
                if k == 1:    
                    # How many years have forecast the events with > max bin (%) of the ensembles?
                    hitrateda[i,j,k] = sum(l > rocbins[-k] for l in pfct_tercile) / nevents
                    falsealarmda[i,j,k] = sum(l > rocbins[-k] for l in pfct_nontercile) / nnonevents
                elif k == len(rocbins)+1:
                    # How many years have forecast the events with < min bin % of the ensembles?
                    hitrateda[i,j,k] = (sum(l < rocbins[0] for l in pfct_tercile) / nevents) + hitrateda[i,j,k-1]
                    falsealarmda[i,j,k] = (sum(l < rocbins[0] for l in pfct_nontercile) / nnonevents) + falsealarmda[i,j,k-1]
                else: 
                    hitrateda[i,j,k] = (sum(l < rocbins[(k-1)*-1] and l > rocbins[-k] for l in pfct_tercile) / nevents) + hitrateda[i,j,k-1]
                    falsealarmda[i,j,k] = (sum(l < rocbins[(k-1)*-1] and l > rocbins[-k] for l in pfct_nontercile) / nnonevents) + falsealarmda[i,j,k-1]

    b = hitrateda.to_dataset(dim='datavar')
    b1 = falsealarmda.to_dataset(dim='datavar')
    return b, b1
